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

/*
 * Extract k-mers from DNA sequence and create n-gram keys
 */
Datum *
kmersearch_extract_kmers(const char *sequence, int seq_len, int k, int *nkeys)
{
    Datum *keys;
    int max_keys = seq_len - k + 1;
    int key_count = 0;
    bool has_degenerate = false;
    int i, j;
    
    if (max_keys <= 0)
    {
        *nkeys = 0;
        return NULL;
    }
    
    keys = (Datum *) palloc(sizeof(Datum) * max_keys * 10);  /* Extra space for degenerate expansion */
    
    for (i = 0; i <= seq_len - k; i++)
    {
        char kmer[65];  /* Max k=64 + null terminator */
        strncpy(kmer, sequence + i, k);
        kmer[k] = '\0';
        
        /* Check if this k-mer has degenerate codes */
        has_degenerate = false;
        for (j = 0; j < k; j++)
        {
            char c = toupper(kmer[j]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            {
                has_degenerate = true;
                break;
            }
        }
        
        if (has_degenerate)
        {
            /* Expand degenerate codes */
            char *expanded[10];
            int expand_count;
            
            kmersearch_expand_degenerate_sequence(kmer, k, expanded, &expand_count);
            
            for (j = 0; j < expand_count; j++)
            {
                /* Count occurrences of this expanded k-mer */
                int occurrence = 1;
                VarBit *ngram_key;
                int prev;
                
                for (prev = 0; prev < key_count; prev++)
                {
                    /* This is simplified - in reality we'd need to compare the actual k-mer */
                    /* For now, assume each k-mer appears once per position */
                }
                
                ngram_key = kmersearch_create_ngram_key2(expanded[j], k, occurrence);
                keys[key_count++] = PointerGetDatum(ngram_key);
                
                pfree(expanded[j]);
            }
        }
        else
        {
            /* Simple case - no degenerate codes */
            int occurrence = 1;
            VarBit *ngram_key = kmersearch_create_ngram_key2(kmer, k, occurrence);
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
    }
    
    *nkeys = key_count;
    return keys;
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
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    /* Extract ngram_key2 (kmer2 + occurrence bits) directly from DNA2 */
    keys = kmersearch_extract_dna2_kmer2_direct((VarBit *)dna, k, nkeys);
    
    /* Apply high-frequency k-mer filtering using ngram_key2 direct comparison */
    if (keys && *nkeys > 0 && kmersearch_preclude_highfreq_kmer) {
        if (kmersearch_force_use_dshash || IsParallelWorker()) {
            /* Use parallel cache for worker processes or when forcing dshash */
            if (parallel_highfreq_cache && parallel_highfreq_cache->is_initialized) {
                keys = kmersearch_filter_highfreq_kmers_from_keys_parallel(keys, nkeys, k);
            } else {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("parallel high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_parallel_highfreq_kmers_cache_load() to create the cache first.")));
            }
        } else {
            /* Use global cache for main process */
            if (global_highfreq_cache.is_valid) {
                keys = kmersearch_filter_highfreq_kmers_from_keys(keys, nkeys, global_highfreq_cache.highfreq_hash, k);
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
        if (kmersearch_force_use_dshash || IsParallelWorker()) {
            /* Use parallel cache for worker processes or when forcing dshash */
            if (parallel_highfreq_cache && parallel_highfreq_cache->is_initialized) {
                keys = kmersearch_filter_highfreq_kmers_from_keys_parallel(keys, nkeys, k);
            } else {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("parallel high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_parallel_highfreq_kmers_cache_load() to create the cache first.")));
            }
        } else {
            /* Use global cache for main process */
            if (global_highfreq_cache.is_valid) {
                keys = kmersearch_filter_highfreq_kmers_from_keys(keys, nkeys, global_highfreq_cache.highfreq_hash, k);
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
    
    if (query_len < 8)
        ereport(ERROR, (errmsg("Query sequence must be at least 8 bases long")));
    
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
    
    /* Handle high-frequency k-mer cache only if enabled */
    if (kmersearch_preclude_highfreq_kmer) {
        if (kmersearch_force_use_dshash || IsParallelWorker()) {
            /* Use parallel cache for worker processes or when forcing dshash */
            if (!(parallel_highfreq_cache && parallel_highfreq_cache->is_initialized)) {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("parallel high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_parallel_highfreq_kmers_cache_load() to create the cache first.")));
            }
        } else {
            /* Use global cache for main process */
            if (!global_highfreq_cache.is_valid) {
                ereport(ERROR,
                        (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                         errmsg("global high-frequency k-mer cache is not initialized"),
                         errhint("Use kmersearch_highfreq_kmers_cache_load() to create the cache first.")));
            }
        }
    }
    
    *recheck = true;  /* Always recheck for scoring */
    
    /* Count matching keys */
    for (i = 0; i < nkeys; i++)
    {
        if (check[i])
            match_count++;
    }
    
    /* Convert queryKeys to VarBit array for excluded k-mer checking */
    query_key_array = (VarBit **) palloc(nkeys * sizeof(VarBit *));
    for (i = 0; i < nkeys; i++)
    {
        query_key_array[i] = DatumGetVarBitP(queryKeys[i]);
    }
    
    /* Calculate actual minimum score using comprehensive scoring logic */
    actual_min_score = calculate_actual_min_score(query_key_array, nkeys, nkeys);
    
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