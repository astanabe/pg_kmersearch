/*
 * kmersearch_gin.c - GIN index support functions for pg_kmersearch
 *
 * This module contains all GIN-related functionality including:
 * - extract_value functions for DNA2 and DNA4 types
 * - extract_query function for query processing  
 * - consistent function for index consistency checking
 * - compare_partial function for partial key comparison
 * - Supporting utility functions for k-mer extraction and processing
 */

#include "kmersearch.h"

/* PostgreSQL function info declarations for GIN functions */
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna2);
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna4);
PG_FUNCTION_INFO_V1(kmersearch_extract_query);
PG_FUNCTION_INFO_V1(kmersearch_consistent);
PG_FUNCTION_INFO_V1(kmersearch_compare_partial);

/*
 * Forward declarations for static functions
 */
static bool kmersearch_is_highfreq_kmer_parallel(VarBit *ngram_key);

/*
 * Get index information from index OID
 */
bool
kmersearch_get_index_info(Oid index_oid, Oid *table_oid, char **column_name, int *k_size)
{
    int ret;
    bool found = false;
    StringInfoData query;
    
    /* Connect to SPI */
    if (SPI_connect() != SPI_OK_CONNECT)
        return false;
    
    /* Build query to get index information */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT table_oid, column_name, kmer_size FROM kmersearch_index_info "
        "WHERE index_oid = %u",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    if (ret != SPI_OK_SELECT)
        ereport(ERROR, (errmsg("SPI_execute failed with code %d", ret)));
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        Datum table_oid_datum, column_name_datum, k_size_datum;
        bool isnull;
        
        table_oid_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull && table_oid)
            *table_oid = DatumGetObjectId(table_oid_datum);
        
        column_name_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull);
        (void) column_name_datum;  /* Suppress unused variable warning */
        if (!isnull && column_name)
        {
            char *col_name = SPI_getvalue(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2);
            *column_name = pstrdup(col_name);
        }
        
        k_size_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull);
        if (!isnull && k_size)
            *k_size = DatumGetInt32(k_size_datum);
        
        found = true;
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return found;
}


/*
 * GIN extract_value function for DNA2
 * Note: Exclusion filtering is applied separately after index creation
 * via kmersearch_analyze_table() and related functions
 */
Datum
kmersearch_extract_value_dna2(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *dna = (kmersearch_dna2 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    void *kmer_uint_array;
    int kmer_uint_count;
    KmerOccurrence *occurrences;
    int max_occurrences;
    int final_count = 0;
    int i;
    
    if (kmersearch_kmer_size < 4 || kmersearch_kmer_size > 32)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 32")));
    
    if (kmersearch_preclude_highfreq_kmer) {
        /* New optimized flow: extract uint k-mers and filter before VarBit creation */
        kmersearch_extract_dna2_kmer2_as_uint_direct((VarBit *)dna, kmersearch_kmer_size, &kmer_uint_array, &kmer_uint_count);
        
        if (kmer_uint_count == 0) {
            *nkeys = 0;
            PG_RETURN_POINTER(NULL);
        }
        
        /* Initialize occurrence tracking */
        max_occurrences = kmer_uint_count;
        occurrences = (KmerOccurrence *) palloc(max_occurrences * sizeof(KmerOccurrence));
        
        /* Process each uint k-mer */
        for (i = 0; i < kmer_uint_count; i++) {
            uint64 kmer_uint;
            bool is_high_frequency = false;
            
            /* Extract k-mer based on size */
            if (kmersearch_kmer_size <= 8) {
                kmer_uint = ((uint16 *)kmer_uint_array)[i];
            } else if (kmersearch_kmer_size <= 16) {
                kmer_uint = ((uint32 *)kmer_uint_array)[i];
            } else {
                kmer_uint = ((uint64 *)kmer_uint_array)[i];
            }
            
            /* Check cache */
            if (kmersearch_force_use_parallel_highfreq_kmer_cache || IsParallelWorker()) {
                if (parallel_highfreq_cache && parallel_highfreq_cache->is_initialized) {
                    is_high_frequency = kmersearch_lookup_kmer2_as_uint_in_parallel_cache(kmer_uint, NULL, NULL);
                } else {
                    ereport(ERROR,
                            (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                             errmsg("parallel high-frequency k-mer cache is not initialized"),
                             errhint("Use kmersearch_parallel_highfreq_kmers_cache_load() to create the cache first.")));
                }
            } else {
                if (global_highfreq_cache.is_valid) {
                    is_high_frequency = kmersearch_lookup_kmer2_as_uint_in_global_cache(kmer_uint, NULL, NULL);
                } else {
                    ereport(ERROR,
                            (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                             errmsg("global high-frequency k-mer cache is not initialized"),
                             errhint("Use kmersearch_highfreq_kmers_cache_load() to create the cache first.")));
                }
            }
            
            /* Skip high-frequency k-mers */
            if (is_high_frequency)
                continue;
            
            /* Track occurrence count */
            kmersearch_find_or_add_kmer_occurrence(occurrences, &final_count, kmer_uint, max_occurrences);
        }
        
        /* Create final ngram_key2 array */
        if (final_count > 0) {
            keys = (Datum *) palloc(final_count * sizeof(Datum));
            for (i = 0; i < final_count; i++) {
                VarBit *ngram_key = kmersearch_create_ngram_key2_from_kmer2_as_uint(
                    occurrences[i].kmer_value, kmersearch_kmer_size, occurrences[i].count);
                keys[i] = PointerGetDatum(ngram_key);
            }
        } else {
            keys = NULL;
        }
        
        /* Cleanup */
        pfree(kmer_uint_array);
        pfree(occurrences);
        
        *nkeys = final_count;
    } else {
        /* Original flow: extract ngram_key2 directly without filtering */
        keys = kmersearch_extract_dna2_ngram_key2_direct((VarBit *)dna, kmersearch_kmer_size, nkeys);
    }
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
}

/*
 * GIN extract_value function for DNA4
 * Note: Exclusion filtering is applied separately after index creation
 * via kmersearch_analyze_table() and related functions
 */
Datum
kmersearch_extract_value_dna4(PG_FUNCTION_ARGS)
{
    kmersearch_dna4 *dna = (kmersearch_dna4 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    void *kmer_uint_array;
    int kmer_uint_count;
    KmerOccurrence *occurrences;
    int max_occurrences;
    int final_count = 0;
    int i;
    
    if (kmersearch_kmer_size < 4 || kmersearch_kmer_size > 32)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 32")));
    
    if (kmersearch_preclude_highfreq_kmer) {
        /* New optimized flow: extract uint k-mers and filter before VarBit creation */
        kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct((VarBit *)dna, kmersearch_kmer_size, &kmer_uint_array, &kmer_uint_count);
        
        if (kmer_uint_count == 0) {
            *nkeys = 0;
            PG_RETURN_POINTER(NULL);
        }
        
        /* Initialize occurrence tracking */
        max_occurrences = kmer_uint_count;
        occurrences = (KmerOccurrence *) palloc(max_occurrences * sizeof(KmerOccurrence));
        
        /* Process each expanded uint k-mer */
        for (i = 0; i < kmer_uint_count; i++) {
            uint64 kmer_uint;
            bool is_high_frequency = false;
            
            /* Extract k-mer based on size */
            if (kmersearch_kmer_size <= 8) {
                kmer_uint = ((uint16 *)kmer_uint_array)[i];
            } else if (kmersearch_kmer_size <= 16) {
                kmer_uint = ((uint32 *)kmer_uint_array)[i];
            } else {
                kmer_uint = ((uint64 *)kmer_uint_array)[i];
            }
            
            /* Check cache */
            if (kmersearch_force_use_parallel_highfreq_kmer_cache || IsParallelWorker()) {
                if (parallel_highfreq_cache && parallel_highfreq_cache->is_initialized) {
                    is_high_frequency = kmersearch_lookup_kmer2_as_uint_in_parallel_cache(kmer_uint, NULL, NULL);
                } else {
                    ereport(ERROR,
                            (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                             errmsg("parallel high-frequency k-mer cache is not initialized"),
                             errhint("Use kmersearch_parallel_highfreq_kmers_cache_load() to create the cache first.")));
                }
            } else {
                if (global_highfreq_cache.is_valid) {
                    is_high_frequency = kmersearch_lookup_kmer2_as_uint_in_global_cache(kmer_uint, NULL, NULL);
                } else {
                    ereport(ERROR,
                            (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                             errmsg("global high-frequency k-mer cache is not initialized"),
                             errhint("Use kmersearch_highfreq_kmers_cache_load() to create the cache first.")));
                }
            }
            
            /* Skip high-frequency k-mers */
            if (is_high_frequency)
                continue;
            
            /* Track occurrence count */
            kmersearch_find_or_add_kmer_occurrence(occurrences, &final_count, kmer_uint, max_occurrences);
        }
        
        /* Create final ngram_key2 array */
        if (final_count > 0) {
            keys = (Datum *) palloc(final_count * sizeof(Datum));
            for (i = 0; i < final_count; i++) {
                VarBit *ngram_key = kmersearch_create_ngram_key2_from_kmer2_as_uint(
                    occurrences[i].kmer_value, kmersearch_kmer_size, occurrences[i].count);
                keys[i] = PointerGetDatum(ngram_key);
            }
        } else {
            keys = NULL;
        }
        
        /* Cleanup */
        pfree(kmer_uint_array);
        pfree(occurrences);
        
        *nkeys = final_count;
    } else {
        /* Original flow: extract ngram_key2 directly without filtering */
        keys = kmersearch_extract_dna4_ngram_key2_with_expansion_direct((VarBit *)dna, kmersearch_kmer_size, nkeys);
    }
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
}

/*
 * Filter high-frequency k-mers from query keys and set actual_min_score
 * 
 * Input: query_keys containing high-frequency k-mers
 * Output: filtered query_keys with high-frequency k-mers removed
 * Side effect: actual_min_score is cached with filtered keys as cache key
 */
VarBit **
filter_ngram_key2_and_set_actual_min_score(VarBit **query_keys, int *nkeys, 
                                           const char *query_string)
{
    VarBit **filtered_keys;
    int filtered_count = 0;
    int original_nkeys = *nkeys;
    int actual_min_score;
    int i;
    
    if (!query_keys || *nkeys <= 0)
        return query_keys;
    
    /* Step 1: Calculate actual_min_score with original keys and cache it */
    /* Note: get_cached_actual_min_score will internally create cache key from filtered keys */
    actual_min_score = get_cached_actual_min_score(query_keys, *nkeys);
    
    /* Step 2: Filter out high-frequency k-mers if enabled */
    if (!kmersearch_preclude_highfreq_kmer) {
        /* High-frequency k-mer exclusion disabled - return original keys */
        return query_keys;
    }
    
    filtered_keys = (VarBit **) palloc(*nkeys * sizeof(VarBit *));
    
    for (i = 0; i < *nkeys; i++) {
        if (!kmersearch_is_kmer_highfreq(query_keys[i])) {
            filtered_keys[filtered_count++] = query_keys[i];
        }
    }
    
    *nkeys = filtered_count;
    
    elog(DEBUG1, "filter_ngram_key2_and_set_actual_min_score: filtered %d high-freq k-mers from %d total, "
                 "cached actual_min_score=%d", 
         original_nkeys - filtered_count, original_nkeys, actual_min_score);
    
    if (filtered_count == 0) {
        pfree(filtered_keys);
        return NULL;
    }
    
    return filtered_keys;
}

/*
 * GIN extract_query function
 */
Datum
kmersearch_extract_query(PG_FUNCTION_ARGS)
{
    Datum query = PG_GETARG_DATUM(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    StrategyNumber strategy = PG_GETARG_UINT16(2);
    bool **pmatch = (bool **) PG_GETARG_POINTER(3);
    Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4);
    bool **nullFlags = (bool **) PG_GETARG_POINTER(5);
    int32 *searchMode = (int32 *) PG_GETARG_POINTER(6);
    
    text *query_text = DatumGetTextP(query);
    char *query_string = text_to_cstring(query_text);
    int query_len = strlen(query_string);
    VarBit **varbit_keys;
    Datum *keys;
    int i;
    
    
    if (query_len < kmersearch_kmer_size)
        ereport(ERROR, (errmsg("Query sequence must be at least %d bases long", kmersearch_kmer_size)));
    
    if (kmersearch_kmer_size < 4 || kmersearch_kmer_size > 32)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 32")));
    
    /* Use cached query pattern extraction */
    varbit_keys = get_cached_query_kmer(query_string, kmersearch_kmer_size, nkeys);
    
    /* Filter high-frequency k-mers and cache actual_min_score */
    if (varbit_keys != NULL && *nkeys > 0) {
        varbit_keys = filter_ngram_key2_and_set_actual_min_score(varbit_keys, nkeys, query_string);
    }
    
    if (varbit_keys == NULL || *nkeys == 0) {
        keys = NULL;
    } else {
        keys = (Datum *) palloc(*nkeys * sizeof(Datum));
        for (i = 0; i < *nkeys; i++) {
            keys[i] = PointerGetDatum(varbit_keys[i]);
        }
        /* NOTE: varbit_keys may be newly allocated by filter function, don't free it */
    }
    
    *pmatch = NULL;
    *extra_data = NULL;
    *nullFlags = NULL;
    *searchMode = GIN_SEARCH_MODE_DEFAULT;
    
    pfree(query_string);
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
}

/*
 * GIN consistent function
 */
Datum
kmersearch_consistent(PG_FUNCTION_ARGS)
{
    bool *check = (bool *) PG_GETARG_POINTER(0);
    StrategyNumber strategy = PG_GETARG_UINT16(1);
    Datum query = PG_GETARG_DATUM(2);
    int32 nkeys = PG_GETARG_INT32(3);
    Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4);
    bool *recheck = (bool *) PG_GETARG_POINTER(5);
    Datum *queryKeys = (Datum *) PG_GETARG_POINTER(6);
    bool *nullFlags = (bool *) PG_GETARG_POINTER(7);
    
    int match_count = 0;
    int actual_min_score;
    int i;
    VarBit **query_key_array;
    
    /* High-frequency k-mer cache checking is not needed during search operations */
    
    *recheck = false;  /* No recheck needed - actual_min_score accounts for high-frequency k-mer exclusion */
    
    /* Count matching keys */
    for (i = 0; i < nkeys; i++)
    {
        if (check[i])
            match_count++;
    }
    
    /* Convert queryKeys to VarBit array for actual min score calculation */
    query_key_array = (VarBit **) palloc(nkeys * sizeof(VarBit *));
    for (i = 0; i < nkeys; i++)
    {
        query_key_array[i] = DatumGetVarBitP(queryKeys[i]);
        /* Debug: Log query key information */
        elog(DEBUG1, "kmersearch_gin_consistent: query_key[%d] bit length = %d", 
             i, VARBITLEN(query_key_array[i]));
    }
    
    /* Debug: Log before calling get_cached_actual_min_score_or_error */
    elog(DEBUG1, "kmersearch_gin_consistent: calling get_cached_actual_min_score_or_error with nkeys = %d", nkeys);
    
    /* Query keys are already filtered - use error version to ensure cache hit */
    actual_min_score = get_cached_actual_min_score_or_error(query_key_array, nkeys);
    
    pfree(query_key_array);
    
    /* Return true if match count meets actual minimum score */
    PG_RETURN_BOOL(match_count >= actual_min_score);
}

/*
 * GIN compare_partial function - simple byte comparison for varbit
 */
Datum
kmersearch_compare_partial(PG_FUNCTION_ARGS)
{
    VarBit *a = DatumGetVarBitP(PG_GETARG_DATUM(0));
    VarBit *b = DatumGetVarBitP(PG_GETARG_DATUM(1));
    int result;
    
    int32 len_a = VARBITLEN(a);
    int32 len_b = VARBITLEN(b);
    
    if (len_a < len_b)
        result = -1;
    else if (len_a > len_b)
        result = 1;
    else
    {
        result = memcmp(VARBITS(a), VARBITS(b), VARBITBYTES(a));
    }
    
    PG_RETURN_INT32(result);
}

/*
 * High-frequency k-mer filtering functions implementation
 */
Datum *
kmersearch_filter_highfreq_ngram_key2(Datum *original_keys, int *nkeys, HTAB *highfreq_hash, int k)
{
    Datum *filtered_keys;
    int original_count;
    int filtered_count;
    int i;
    
    if (!original_keys || !nkeys || *nkeys <= 0)
        return NULL;
    
    /* If no high-frequency hash, return original keys unchanged */
    if (!highfreq_hash)
        return original_keys;
    
    original_count = *nkeys;
    filtered_keys = (Datum *) palloc(original_count * sizeof(Datum));
    filtered_count = 0;
    
    /* Filter out high-frequency k-mers */
    for (i = 0; i < original_count; i++)
    {
        VarBit *ngram_key;
        bool found;
        
        ngram_key = DatumGetVarBitP(original_keys[i]);
        if (!ngram_key)
            continue;
        
        /* Use ngram_key2 directly for high-frequency lookup - no occurrence bits removal needed */
        {
            /* Calculate hash using complete ngram_key2 (kmer2 + occurrence bits) */
            int bit_length = VARBITLEN(ngram_key);
            int byte_count = (bit_length + 7) / 8;  /* Round up to next byte */
            uint64 hash_value = DatumGetUInt64(hash_any((unsigned char *) VARBITS(ngram_key), byte_count));
            
            hash_search(highfreq_hash,
                       (void *) &hash_value,
                       HASH_FIND,
                       &found);
        }
        
        if (!found)
        {
            /* Not a high-frequency k-mer, keep it */
            filtered_keys[filtered_count++] = original_keys[i];
        }
        
        /* No need to free ngram_key since it's managed by caller */
    }
    
    /* Update the count */
    *nkeys = filtered_count;
    
    /* If no keys left, return NULL */
    if (filtered_count == 0)
    {
        pfree(filtered_keys);
        return NULL;
    }
    
    /* Resize the array if significantly smaller */
    if (filtered_count < original_count / 2)
    {
        filtered_keys = (Datum *) repalloc(filtered_keys, filtered_count * sizeof(Datum));
    }
    
    return filtered_keys;
}

/*
 * Filter high-frequency k-mers from keys using parallel cache
 */
Datum *
kmersearch_filter_highfreq_ngram_key2_parallel(Datum *original_keys, int *nkeys, int k)
{
    Datum *filtered_keys;
    int filtered_count = 0;
    int i;
    
    if (!original_keys || *nkeys <= 0)
        return original_keys;
    
    /* Allocate space for filtered keys */
    filtered_keys = (Datum *) palloc(*nkeys * sizeof(Datum));
    
    /* Filter out high-frequency k-mers using direct ngram_key2 comparison */
    for (i = 0; i < *nkeys; i++) {
        VarBit *ngram_key = (VarBit *) DatumGetPointer(original_keys[i]);
        
        /* Use ngram_key2 directly for high-frequency check - no occurrence bits removal needed */
        if (!kmersearch_is_highfreq_kmer_parallel(ngram_key)) {
            /* Keep this k-mer */
            filtered_keys[filtered_count++] = original_keys[i];
        } else {
            /* Free the filtered k-mer */
            pfree(ngram_key);
        }
        
        /* No need to clean up temporary key since we use ngram_key directly */
    }
    
    /* Free original keys array */
    pfree(original_keys);
    
    /* Update the key count */
    *nkeys = filtered_count;
    
    /* Return filtered keys (or NULL if no keys remain) */
    if (filtered_count == 0) {
        pfree(filtered_keys);
        return NULL;
    }
    
    return filtered_keys;
}

/*
 * Check if k-mer is highly frequent using parallel cache
 */
static bool
kmersearch_is_highfreq_kmer_parallel(VarBit *ngram_key)
{
    uint64 ngram_hash;
    
    /* If parallel cache is not available, return false */
    if (parallel_cache_hash == NULL)
        return false;
    
    /* Calculate hash for the ngram_key2 (kmer2 + occurrence bits) */
    ngram_hash = hash_any((unsigned char *) VARDATA(ngram_key), 
                         VARSIZE(ngram_key) - VARHDRSZ);
    
    /* Look up in parallel cache using complete ngram_key2 */
    return kmersearch_parallel_cache_lookup(ngram_hash);
}

/*
 * Evaluate optimized match condition
 */
bool
evaluate_optimized_match_condition(VarBit **query_keys, int nkeys, int shared_count, const char *query_string, int query_total_kmers)
{
    int actual_min_score;
    
    
    /* Get cached actual min score (with TopMemoryContext caching for performance) */
    actual_min_score = get_cached_actual_min_score(query_keys, nkeys);
    elog(LOG, "evaluate_optimized_match_condition: get_cached_actual_min_score returned %d", actual_min_score);
    
    /* Use optimized condition check with cached actual_min_score */
    return (shared_count >= actual_min_score);
}

