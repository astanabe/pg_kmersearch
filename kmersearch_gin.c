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
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    /* Extract ngram_key2 (kmer2 + occurrence bits) directly from DNA2 */
    keys = kmersearch_extract_dna2_kmer2_direct((VarBit *)dna, k, nkeys);
    
    /* Apply high-frequency k-mer filtering using ngram_key2 direct comparison */
    if (keys && *nkeys > 0 && kmersearch_preclude_highfreq_kmer) {
        if (kmersearch_force_use_parallel_highfreq_kmer_cache || IsParallelWorker()) {
            /* Use parallel cache for worker processes or when forcing dshash */
            if (parallel_highfreq_cache && parallel_highfreq_cache->is_initialized) {
                keys = kmersearch_filter_highfreq_ngram_key2_parallel(keys, nkeys, k);
            } else {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("parallel high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_parallel_highfreq_kmers_cache_load() to create the cache first.")));
            }
        } else {
            /* Use global cache for main process */
            if (global_highfreq_cache.is_valid) {
                keys = kmersearch_filter_highfreq_ngram_key2(keys, nkeys, global_highfreq_cache.highfreq_hash, k);
            } else {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("global high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_highfreq_kmers_cache_load() to create the cache first.")));
            }
        }
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
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    /* Extract ngram_key2 (kmer2 + occurrence bits) from DNA4 with degenerate expansion */
    keys = kmersearch_extract_dna4_kmer2_with_expansion_direct((VarBit *)dna, k, nkeys);
    
    /* Apply high-frequency k-mer filtering using ngram_key2 direct comparison */
    if (keys && *nkeys > 0 && kmersearch_preclude_highfreq_kmer) {
        if (kmersearch_force_use_parallel_highfreq_kmer_cache || IsParallelWorker()) {
            /* Use parallel cache for worker processes or when forcing dshash */
            if (parallel_highfreq_cache && parallel_highfreq_cache->is_initialized) {
                keys = kmersearch_filter_highfreq_ngram_key2_parallel(keys, nkeys, k);
            } else {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("parallel high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_parallel_highfreq_kmers_cache_load() to create the cache first.")));
            }
        } else {
            /* Use global cache for main process */
            if (global_highfreq_cache.is_valid) {
                keys = kmersearch_filter_highfreq_ngram_key2(keys, nkeys, global_highfreq_cache.highfreq_hash, k);
            } else {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("global high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_highfreq_kmers_cache_load() to create the cache first.")));
            }
        }
    }
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
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
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    text *query_text = DatumGetTextP(query);
    char *query_string = text_to_cstring(query_text);
    int query_len = strlen(query_string);
    Datum *keys;
    
    if (query_len < k)
        ereport(ERROR, (errmsg("Query sequence must be at least %d bases long", k)));
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    keys = kmersearch_extract_kmers(query_string, query_len, k, nkeys);
    
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
    
    *recheck = true;  /* Always recheck for scoring */
    
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
    }
    
    /* Use cached actual minimum score for better performance */
    actual_min_score = get_cached_actual_min_score(query_key_array, nkeys);
    
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