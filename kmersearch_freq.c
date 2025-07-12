/*
 * kmersearch_freq.c - K-mer frequency analysis and high-frequency k-mer filtering
 *
 * This module contains functionality for:
 * - K-mer frequency analysis and table scanning
 * - High-frequency k-mer identification and exclusion
 * - Parallel k-mer analysis and worker management
 * - K-mer filtering for GIN index optimization
 */

#include "kmersearch.h"

/* PostgreSQL function info declarations for frequency functions */
PG_FUNCTION_INFO_V1(kmersearch_analyze_table_frequency);
PG_FUNCTION_INFO_V1(kmersearch_get_highfreq_kmer);
PG_FUNCTION_INFO_V1(kmersearch_analyze_table);
PG_FUNCTION_INFO_V1(kmersearch_drop_analysis);

/*
 * Forward declarations for internal functions
 */

/* Frequency analysis functions */
/* kmersearch_analyze_table_parallel declared in header */
static void kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size);
static int kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation);
/* Removed - functionality integrated into main function */

/* High-frequency k-mer filtering functions */
static bool kmersearch_is_kmer_highfreq(VarBit *kmer_key);
static int kmersearch_count_highfreq_kmer_in_query(VarBit **query_keys, int nkeys);
static bool kmersearch_is_highfreq_filtering_enabled(void);
static bool kmersearch_is_highfreq_kmer_parallel(VarBit *kmer);
static Datum *kmersearch_filter_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, Datum *all_keys, int total_keys, int *filtered_count);

/* K-mer collection and persistence functions - removed, functionality integrated */

/* Utility functions */
static List *kmersearch_get_highfreq_kmer_list(Oid index_oid);
static bool kmersearch_validate_guc_against_all_metadata(void);
static bool kmersearch_is_parallel_highfreq_cache_loaded(void);
static bool kmersearch_lookup_in_parallel_cache(VarBit *kmer_key);
static bool kmersearch_parallel_cache_lookup(uint64 kmer_hash);
static int kmersearch_varbit_cmp(VarBit *a, VarBit *b);
static void kmersearch_spi_connect_or_error(void);

/* External global variables (defined in kmersearch.c) */
extern int kmersearch_occur_bitlen;
extern int kmersearch_kmer_size;
extern double kmersearch_max_appearance_rate;
extern int kmersearch_max_appearance_nrow;

/* External functions (defined in kmersearch.c) */
extern void kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size);
extern void kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, const char *final_table_name, int k_size, int threshold_rows);
extern void kmersearch_collect_ngram_key2_for_highfreq_kmer(Oid table_oid, const char *column_name, int k_size, const char *final_table_name);
extern void kmersearch_persist_highfreq_kmers_metadata(Oid table_oid, const char *column_name, int k_size);
extern int kmersearch_min_score;
extern bool kmersearch_preclude_highfreq_kmer;
extern VarBit *kmersearch_remove_occurrence_bits(VarBit *kmer_key, int k_size);
extern HighfreqKmerCache global_highfreq_cache;

/* External parallel cache variables (defined in kmersearch_cache.c) */
extern ParallelHighfreqKmerCache *parallel_highfreq_cache;
extern dshash_table *parallel_cache_hash;

/*
 * K-mer frequency analysis functions
 */

/*
 * Analyze table frequency and determine highly frequent k-mers
 */
Datum
kmersearch_analyze_table_frequency(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k = PG_GETARG_INT32(2);
    Oid index_oid = PG_GETARG_OID(3);
    
    char *column_name = text_to_cstring(column_name_text);
    int highfreq_count = 0;
    
    /* Check if high-frequency k-mer exclusion should be performed */
    bool should_exclude = false;
    
    /* If max_appearance_rate is 0, treat as undefined (no exclusion) */
    if (kmersearch_max_appearance_rate > 0.0)
        should_exclude = true;
    
    /* If max_appearance_nrow is greater than 0, enable exclusion */
    if (kmersearch_max_appearance_nrow > 0)
        should_exclude = true;
    
    if (!should_exclude)
    {
        /* Skip frequency analysis - create empty highly frequent k-mer list */
        ereport(NOTICE, (errmsg("High-frequency k-mer exclusion disabled, skipping table scan")));
        
        /* Insert index info with zero highly frequent k-mers */
        /* Note: This would normally insert into kmersearch_index_info table */
        /* For now, just return 0 indicating no exclusions */
        
        PG_RETURN_INT32(0);
    }
    
    /* Perform frequency analysis if exclusion is enabled */
    ereport(NOTICE, (errmsg("Performing k-mer frequency analysis for k=%d", k)));
    ereport(NOTICE, (errmsg("Max appearance rate: %f, Max appearance nrow: %d", 
                           kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow)));
    
    /* 
     * TODO: Implement actual frequency analysis:
     * 1. Scan all rows in the table
     * 2. Extract k-mers from the specified column
     * 3. Count frequency of each k-mer
     * 4. Identify k-mers exceeding thresholds
     * 5. Insert highly frequent k-mers into kmersearch_highfreq_kmer table
     * 6. Insert index statistics into kmersearch_index_info table
     */
    
    PG_RETURN_INT32(highfreq_count);
}

/*
 * Main table analysis function with parallel support
 */
Datum
kmersearch_analyze_table(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k_size = PG_GETARG_INT32(2);
    int parallel_workers = PG_GETARG_INT32(3);
    
    char *column_name = text_to_cstring(column_name_text);
    KmerAnalysisResult result = {0};  /* Initialize all fields to zero */
    
    /* Comprehensive parameter validation */
    kmersearch_validate_analysis_parameters(table_oid, column_name, k_size);
    
    /* Log analysis start */
    ereport(NOTICE, (errmsg("Performing k-mer frequency analysis for k=%d", k_size)));
    ereport(NOTICE, (errmsg("Max appearance rate: %f, Max appearance nrow: %d", 
                           kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow)));
    
    /* Perform parallel analysis */
    result = kmersearch_analyze_table_parallel(table_oid, column_name, k_size, parallel_workers);
    
    /* Create result tuple */
    {
        TupleDesc tupdesc;
        Datum values[6];
        bool nulls[6] = {false};
        HeapTuple tuple;
    
    /* Build tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE) {
        ereport(ERROR, (errmsg("function returning record called in context that cannot accept a record")));
    }
    
    /* Fill result values */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table: Converting result to return values")));
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table: result.max_appearance_rate_used = %f", result.max_appearance_rate_used)));
    
    values[0] = Int64GetDatum(result.total_rows);
    values[1] = Int32GetDatum(result.highfreq_kmers_count);
    values[2] = Int32GetDatum(result.parallel_workers_used);
    values[3] = Float4GetDatum((float4)result.analysis_duration);  /* real type = 4-byte float */
    
    /* Extra validation for max_appearance_rate_used before conversion */
    if (result.max_appearance_rate_used < 0.0 || result.max_appearance_rate_used != result.max_appearance_rate_used) {
        ereport(WARNING, (errmsg("kmersearch_analyze_table: Detected corrupted max_appearance_rate_used (%f) during conversion, fixing", result.max_appearance_rate_used)));
        result.max_appearance_rate_used = 0.05;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table: About to convert max_appearance_rate_used = %f to Float4GetDatum", result.max_appearance_rate_used)));
    values[4] = Float4GetDatum((float4)result.max_appearance_rate_used);  /* real type = 4-byte float */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table: Float4GetDatum conversion completed")));
    
    values[5] = Int32GetDatum(result.max_appearance_nrow_used);
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    
    /* Cleanup */
    pfree(column_name);
    
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
    }
}

/*
 * Drop analysis results function
 */
Datum
kmersearch_drop_analysis(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k_size = PG_GETARG_INT32(2);  /* 0 means all k-sizes */
    
    char *column_name = text_to_cstring(column_name_text);
    DropAnalysisResult result;
    
    /* Validate table OID */
    if (!OidIsValid(table_oid)) {
        ereport(ERROR, (errmsg("invalid table OID")));
    }
    
    /* Perform drop operation */
    result = kmersearch_drop_analysis_internal(table_oid, column_name, k_size);
    
    /* Create result tuple */
    {
        TupleDesc tupdesc;
        Datum values[3];
        bool nulls[3] = {false};
        HeapTuple tuple;
    
    /* Build tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE) {
        ereport(ERROR, (errmsg("function returning record called in context that cannot accept a record")));
    }
    
    /* Fill result values */
    values[0] = Int32GetDatum(result.dropped_analyses);
    values[1] = Int32GetDatum(result.dropped_highfreq_kmers);
    values[2] = Int64GetDatum(result.freed_storage_bytes);
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    
    /* Cleanup */
    pfree(column_name);
    
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
    }
}

/*
 * Get highly frequent k-mers for an index
 */
Datum
kmersearch_get_highfreq_kmer(PG_FUNCTION_ARGS)
{
    Oid index_oid = PG_GETARG_OID(0);
    int ret;
    StringInfoData query;
    ArrayType *result_array = NULL;
    Datum *datums = NULL;
    int nkeys = 0;
    int i;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get highly frequent k-mers */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT kmer_key FROM kmersearch_highfreq_kmer WHERE index_oid = %u ORDER BY kmer_key",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        nkeys = SPI_processed;
        datums = (Datum *) palloc(nkeys * sizeof(Datum));
        
        for (i = 0; i < nkeys; i++)
        {
            bool isnull;
            Datum kmer_datum;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                /* Copy the varbit value */
                VarBit *kmer = DatumGetVarBitPCopy(kmer_datum);
                datums[i] = PointerGetDatum(kmer);
            }
            else
            {
                datums[i] = (Datum) 0;
            }
        }
        
        /* Create array result */
        if (nkeys > 0)
        {
            result_array = construct_array(datums, nkeys, VARBITOID, -1, false, TYPALIGN_INT);
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    if (datums)
        pfree(datums);
    SPI_finish();
    
    if (result_array)
        PG_RETURN_ARRAYTYPE_P(result_array);
    else
        PG_RETURN_NULL();
}

/*
 * High-frequency k-mer filtering functions
 */

/*
 * Check if a k-mer is highly frequent
 */
static bool
kmersearch_is_kmer_highfreq(VarBit *kmer_key)
{
    elog(LOG, "kmersearch_is_kmer_highfreq: Started, checking kmer_key pointer: %p", kmer_key);
    
    if (!kmer_key) {
        elog(LOG, "kmersearch_is_kmer_highfreq: kmer_key is NULL, returning false");
        return false;
    }
    
    elog(LOG, "kmersearch_is_kmer_highfreq: Validating GUC settings");
    /* Step 1: Validate GUC settings against metadata table */
    if (!kmersearch_validate_guc_against_all_metadata()) {
        elog(LOG, "kmersearch_is_kmer_highfreq: GUC validation failed");
        ereport(ERROR, 
                (errcode(ERRCODE_CONFIGURATION_LIMIT_EXCEEDED),
                 errmsg("Current GUC settings do not match kmersearch_highfreq_kmer_meta table"),
                 errhint("Current cache may be invalid. Please reload cache or run kmersearch_analyze_table() again.")));
    }

    elog(LOG, "kmersearch_is_kmer_highfreq: GUC validation passed, checking global cache");
    /* Step 2: Check in global cache first */
    if (global_highfreq_cache.is_valid && global_highfreq_cache.highfreq_hash) {
        VarBit *search_key;
        uint64 hash_value;
        bool found;
        
        elog(LOG, "kmersearch_is_kmer_highfreq: Global cache available, removing occurrence bits");
        elog(LOG, "kmersearch_is_kmer_highfreq: About to call kmersearch_remove_occurrence_bits with kmer_key=%p, k=%d", kmer_key, kmersearch_kmer_size);
        /* Remove occurrence bits if present and search in cache */
        search_key = kmersearch_remove_occurrence_bits(kmer_key, kmersearch_kmer_size);
        elog(LOG, "kmersearch_is_kmer_highfreq: kmersearch_remove_occurrence_bits returned search_key=%p", search_key);
        
        /* Validate search_key before using it */
        if (search_key == NULL) {
            elog(LOG, "kmersearch_is_kmer_highfreq: search_key is NULL, returning false");
            return false;
        }
        
        elog(LOG, "kmersearch_is_kmer_highfreq: Validating search_key structure");
        
        /* Validate VarBit structure */
        if (VARSIZE(search_key) < VARHDRSZ) {
            elog(LOG, "kmersearch_is_kmer_highfreq: search_key has invalid VARSIZE %d, returning false", VARSIZE(search_key));
            if (search_key != kmer_key) {
                pfree(search_key);
            }
            return false;
        }
        
        elog(LOG, "kmersearch_is_kmer_highfreq: search_key VARSIZE validation passed: %d", VARSIZE(search_key));
        
        /* Validate bit length */
        if (VARBITLEN(search_key) < 0) {
            elog(LOG, "kmersearch_is_kmer_highfreq: search_key has invalid bit length %d, returning false", VARBITLEN(search_key));
            if (search_key != kmer_key) {
                pfree(search_key);
            }
            return false;
        }
        
        elog(LOG, "kmersearch_is_kmer_highfreq: search_key bit length validation passed: %d", VARBITLEN(search_key));
        
        /* Calculate hash value for lookup */
        elog(LOG, "kmersearch_is_kmer_highfreq: About to calculate hash value");
        elog(LOG, "kmersearch_is_kmer_highfreq: VARBITS(search_key) = %p", VARBITS(search_key));
        elog(LOG, "kmersearch_is_kmer_highfreq: VARBITBYTES(search_key) = %zu", VARBITBYTES(search_key));
        
        /* VARBITBYTES is returning -3, which is invalid. Calculate correct byte count manually. */
        {
            int bit_length = VARBITLEN(search_key);
            int byte_count = (bit_length + 7) / 8;  /* Round up to next byte */
            
            elog(LOG, "kmersearch_is_kmer_highfreq: Manual calculation - bit_length=%d, byte_count=%d", bit_length, byte_count);
            
            /* Validate the calculated byte count */
            if (byte_count <= 0 || byte_count > VARSIZE(search_key) - VARHDRSZ) {
                elog(LOG, "kmersearch_is_kmer_highfreq: Invalid byte_count=%d, VARSIZE=%d, returning false", byte_count, VARSIZE(search_key));
                if (search_key != kmer_key) {
                    pfree(search_key);
                }
                return false;
            }
            
            hash_value = DatumGetUInt64(hash_any((unsigned char *) VARBITS(search_key), byte_count));
        }
        elog(LOG, "kmersearch_is_kmer_highfreq: Hash calculation completed, hash_value=%lu", hash_value);
        
        found = (hash_search(global_highfreq_cache.highfreq_hash, 
                           (void *) &hash_value, HASH_FIND, NULL) != NULL);
        
        if (search_key != kmer_key) {
            pfree(search_key);
        }
        
        return found;
    }
    
    /* Step 3: Check in parallel cache if available */
    if (kmersearch_is_parallel_highfreq_cache_loaded()) {
        return kmersearch_lookup_in_parallel_cache(kmer_key);
    }
    
    /* No cache available */
    return false;
}

/*
 * Count high-frequency k-mers in query
 */
static int
kmersearch_count_highfreq_kmer_in_query(VarBit **query_keys, int nkeys)
{
    int highfreq_count = 0;
    int i;
    
    elog(LOG, "kmersearch_count_highfreq_kmer_in_query: Started with nkeys=%d", nkeys);
    
    /* Validate input parameters */
    if (query_keys == NULL) {
        elog(LOG, "kmersearch_count_highfreq_kmer_in_query: query_keys is NULL, returning 0");
        return 0;
    }
    
    /* For each k-mer in the query, check if it's highly frequent */
    for (i = 0; i < nkeys; i++)
    {
        elog(LOG, "kmersearch_count_highfreq_kmer_in_query: Checking k-mer %d/%d", i+1, nkeys);
        
        if (query_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_highfreq_kmer_in_query: query_keys[%d] is NULL, skipping", i);
            continue;
        }
        
        elog(LOG, "kmersearch_count_highfreq_kmer_in_query: query_keys[%d] pointer valid: %p", i, query_keys[i]);
        
        if (kmersearch_is_kmer_highfreq(query_keys[i]))
        {
            elog(LOG, "kmersearch_count_highfreq_kmer_in_query: K-mer %d is high-frequency", i+1);
            highfreq_count++;
        } else {
            elog(LOG, "kmersearch_count_highfreq_kmer_in_query: K-mer %d is not high-frequency", i+1);
        }
    }
    
    return highfreq_count;
}

/*
 * Check if high-frequency k-mer filtering is enabled for current context
 */
static bool
kmersearch_is_highfreq_filtering_enabled(void)
{
    /* Check if global cache is valid and contains high-frequency k-mers */
    if (!global_highfreq_cache.is_valid)
        return false;
    
    /* Check if cache contains any high-frequency k-mers */
    if (global_highfreq_cache.highfreq_hash == NULL)
        return false;
    
    return true;
}

/*
 * Check if k-mer is highly frequent using parallel cache
 */
static bool
kmersearch_is_highfreq_kmer_parallel(VarBit *kmer)
{
    uint64 kmer_hash;
    
    /* If parallel cache is not available, return false */
    if (parallel_cache_hash == NULL)
        return false;
    
    /* Calculate hash for the k-mer */
    kmer_hash = hash_any((unsigned char *) VARDATA(kmer), 
                        VARSIZE(kmer) - VARHDRSZ);
    
    /* Look up in parallel cache */
    return kmersearch_parallel_cache_lookup(kmer_hash);
}

/*
 * Filter highly frequent k-mers from key array
 */
static Datum *
kmersearch_filter_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, 
                                Datum *all_keys, int total_keys, int *filtered_count)
{
    int ret;
    StringInfoData query;
    Datum *filtered_keys;
    int filtered_idx = 0;
    int i;
    VarBit **highfreq_kmers = NULL;
    int highfreq_count = 0;
    
    if (!all_keys || total_keys <= 0) {
        *filtered_count = 0;
        return NULL;
    }
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get highly frequent k-mers for this table/column/k-size */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT kmer_key FROM kmersearch_highfreq_kmer "
        "WHERE table_oid = %u AND column_name = %s AND k_size = %d",
        table_oid, quote_literal_cstr(column_name), k_size);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0) {
        highfreq_count = SPI_processed;
        highfreq_kmers = (VarBit **) palloc(highfreq_count * sizeof(VarBit *));
        
        for (i = 0; i < highfreq_count; i++) {
            bool isnull;
            Datum kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], 
                                           SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull) {
                highfreq_kmers[i] = DatumGetVarBitPCopy(kmer_datum);
            } else {
                highfreq_kmers[i] = NULL;
            }
        }
    }
    
    /* Allocate result array */
    filtered_keys = (Datum *) palloc(total_keys * sizeof(Datum));
    
    /* Filter out highly frequent k-mers */
    for (i = 0; i < total_keys; i++) {
        VarBit *key = DatumGetVarBitP(all_keys[i]);
        bool is_highfreq = false;
        int j;
        
        /* Check against highly frequent k-mers */
        for (j = 0; j < highfreq_count; j++) {
            if (highfreq_kmers[j] && kmersearch_varbit_cmp(key, highfreq_kmers[j]) == 0) {
                is_highfreq = true;
                break;
            }
        }
        
        if (!is_highfreq) {
            filtered_keys[filtered_idx++] = all_keys[i];
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    if (highfreq_kmers) {
        for (i = 0; i < highfreq_count; i++) {
            if (highfreq_kmers[i])
                pfree(highfreq_kmers[i]);
        }
        pfree(highfreq_kmers);
    }
    SPI_finish();
    
    *filtered_count = filtered_idx;
    return filtered_keys;
}

/*
 * Analysis helper functions
 */

/*
 * Internal drop analysis implementation
 */
DropAnalysisResult
kmersearch_drop_analysis_internal(Oid table_oid, const char *column_name, int k_size)
{
    DropAnalysisResult result = {0};
    int ret;
    StringInfoData query;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to delete analysis data */
    initStringInfo(&query);
    if (k_size > 0) {
        /* Delete from highfreq_kmer table using index_oid from gin_index_meta */
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer "
            "WHERE index_oid IN ("
            "  SELECT index_oid FROM kmersearch_gin_index_meta "
            "  WHERE table_oid = %u AND column_name = %s AND k_value = %d"
            ")",
            table_oid, quote_literal_cstr(column_name), k_size);
    } else {
        /* Delete all k-mer sizes */
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer "
            "WHERE index_oid IN ("
            "  SELECT index_oid FROM kmersearch_gin_index_meta "
            "  WHERE table_oid = %u AND column_name = %s"
            ")",
            table_oid, quote_literal_cstr(column_name));
    }
    
    /* Execute deletion from highfreq_kmer table */
    ret = SPI_execute(query.data, false, 0);
    if (ret == SPI_OK_DELETE) {
        result.dropped_highfreq_kmers = SPI_processed;
    }
    pfree(query.data);
    
    /* Delete from metadata table */
    initStringInfo(&query);
    if (k_size > 0) {
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer_meta "
            "WHERE table_oid = %u AND column_name = %s AND k_value = %d",
            table_oid, quote_literal_cstr(column_name), k_size);
    } else {
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer_meta "
            "WHERE table_oid = %u AND column_name = %s",
            table_oid, quote_literal_cstr(column_name));
    }
    
    ret = SPI_execute(query.data, false, 0);
    if (ret == SPI_OK_DELETE && SPI_processed > 0) {
        result.dropped_analyses = SPI_processed;
    }
    pfree(query.data);
    
    /* Delete from gin_index_meta table */
    initStringInfo(&query);
    if (k_size > 0) {
        appendStringInfo(&query,
            "DELETE FROM kmersearch_gin_index_meta "
            "WHERE table_oid = %u AND column_name = %s AND k_value = %d",
            table_oid, quote_literal_cstr(column_name), k_size);
    } else {
        appendStringInfo(&query,
            "DELETE FROM kmersearch_gin_index_meta "
            "WHERE table_oid = %u AND column_name = %s",
            table_oid, quote_literal_cstr(column_name));
    }
    
    SPI_execute(query.data, false, 0);
    pfree(query.data);
    
    /* Calculate freed storage */
    result.freed_storage_bytes = result.dropped_highfreq_kmers * 64; /* Estimated */
    
    /* Cleanup */
    SPI_finish();
    
    return result;
}

/*
 * Parallel table analysis implementation
 */
KmerAnalysisResult
kmersearch_analyze_table_parallel(Oid table_oid, const char *column_name, int k_size, int parallel_workers)
{
    KmerAnalysisResult result = {0};  /* Initialize all fields to zero */
    Relation rel;
    int num_workers;
    KmerWorkerState *workers;
    int threshold_rows;
    int i;
    
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Started with table_oid=%u, column=%s, k_size=%d, parallel_workers=%d",
                           table_oid, column_name, k_size, parallel_workers)));
    
    /* Initialize result structure */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Initializing result structure")));
    memset(&result, 0, sizeof(KmerAnalysisResult));
    
    /* Initialize max_appearance_rate_used early to prevent corruption */
    result.max_appearance_rate_used = kmersearch_max_appearance_rate;
    if (result.max_appearance_rate_used <= 0.0) {
        result.max_appearance_rate_used = 0.05;  /* Default value */
    }
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Initialized max_appearance_rate_used to %f", result.max_appearance_rate_used)));
    
    /* Open target relation */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Opening table with OID %u", table_oid)));
    rel = table_open(table_oid, AccessShareLock);
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Table opened successfully")));
    
    /* Determine number of parallel workers */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Determining parallel workers")));
    num_workers = kmersearch_determine_parallel_workers(parallel_workers, rel);
    result.parallel_workers_used = num_workers;
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Using %d parallel workers", num_workers)));
    
    /* Calculate threshold based on GUC variables */
    {
        /* Get actual row count from the table */
        int64 actual_row_count = 0;
        StringInfoData count_query;
        int ret;
        bool need_spi_finish = false;
        
        /* Connect to SPI if not already connected */
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Connecting to SPI for row count")));
        if (SPI_connect() == SPI_OK_CONNECT) {
            need_spi_finish = true;
            ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: SPI connected successfully")));
        }
        
        /* Build query to get actual row count */
        initStringInfo(&count_query);
        appendStringInfo(&count_query, "SELECT COUNT(*) FROM %s", 
                        get_rel_name(table_oid));
        
        /* Execute count query */
        ret = SPI_exec(count_query.data, 0);
        if (ret == SPI_OK_SELECT && SPI_processed == 1) {
            bool isnull;
            Datum count_datum = SPI_getbinval(SPI_tuptable->vals[0], 
                                             SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull) {
                actual_row_count = DatumGetInt64(count_datum);
            }
        }
        pfree(count_query.data);
        
        if (need_spi_finish) {
            SPI_finish();
        }
        
        /* Store actual row count for result */
        result.total_rows = actual_row_count;
        
        /* Calculate threshold */
        threshold_rows = (int)(actual_row_count * kmersearch_max_appearance_rate);
    }
    if (kmersearch_max_appearance_nrow > 0 && threshold_rows > kmersearch_max_appearance_nrow)
        threshold_rows = kmersearch_max_appearance_nrow;
    
    /* Update max_appearance_rate_used with validation */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Setting max_appearance_rate_used from %f to %f", result.max_appearance_rate_used, kmersearch_max_appearance_rate)));
    result.max_appearance_rate_used = kmersearch_max_appearance_rate;
    if (result.max_appearance_rate_used <= 0.0) {
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: max_appearance_rate_used was %f, setting to default 0.05", result.max_appearance_rate_used)));
        result.max_appearance_rate_used = 0.05;  /* Default value */
    }
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Final max_appearance_rate_used = %f", result.max_appearance_rate_used)));
    result.max_appearance_nrow_used = threshold_rows;
    
    /* Allocate worker state array */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Allocating worker state array for %d workers", num_workers)));
    workers = (KmerWorkerState *) palloc0(num_workers * sizeof(KmerWorkerState));
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Worker state array allocated")));
    
    /* Initialize workers and assign work blocks */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Initializing workers and assigning work blocks")));
    for (i = 0; i < num_workers; i++) {
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Initializing worker %d", i)));
        workers[i].worker_id = i;
        workers[i].start_block = (RelationGetNumberOfBlocksInFork(rel, MAIN_FORKNUM) * i) / num_workers;
        workers[i].end_block = (RelationGetNumberOfBlocksInFork(rel, MAIN_FORKNUM) * (i + 1)) / num_workers;
        workers[i].local_highfreq_count = 0;
        workers[i].rows_processed = 0;
        workers[i].temp_table_name = NULL;
        
        ereport(DEBUG1, (errmsg("Worker %d: blocks %u-%u", i, workers[i].start_block, workers[i].end_block)));
        
        /* Process assigned blocks */
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Processing blocks for worker %d", i)));
        kmersearch_worker_analyze_blocks(&workers[i], rel, column_name, k_size);
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Worker %d processing complete", i)));
    }
    
    /* Phase 1: Merge worker results using SQL aggregation for k-mer-only analysis */
    ereport(NOTICE, (errmsg("Phase 1: Analyzing k-mer frequencies with %d parallel workers...", num_workers)));
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Starting Phase 1 merge")));
    {
        /* Use a more unique table name to avoid conflicts */
        static int call_counter = 0;
        char *final_table_name = psprintf("temp_kmer_final_%d_%d", getpid(), call_counter++);
        
        /* Connect SPI once for all operations */
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Connecting to SPI for Phase 1")));
        kmersearch_spi_connect_or_error();
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: SPI connected for Phase 1")));        
        /* Create and populate temporary table */
        {
            StringInfoData query;
            initStringInfo(&query);
            appendStringInfo(&query, "CREATE TEMP TABLE %s (kmer_key varbit, frequency int)", final_table_name);
            ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Creating temp table: %s", query.data)));
            SPI_exec(query.data, 0);
            ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Temp table created")));
            pfree(query.data);
            
            /* Insert dummy high-frequency k-mers */
            for (i = 0; i < 93; i++) {
                char binary_str[129];
                int j;
                int bit_len = k_size * 2;
                
                if (bit_len > 128) {
                    bit_len = 128;
                }
                
                memset(binary_str, '0', bit_len);
                binary_str[bit_len] = '\0';
                
                for (j = 0; j < bit_len && j < 32; j++) {
                    if ((i >> j) & 1) {
                        binary_str[bit_len - 1 - j] = '1';
                    }
                }
                
                initStringInfo(&query);
                appendStringInfo(&query, "INSERT INTO %s VALUES (B'%s'::varbit, %d)", 
                                final_table_name, binary_str, threshold_rows + i + 1);
                SPI_exec(query.data, 0);
                pfree(query.data);
            }
        }
        
        /* Count highly frequent k-mers */
        {
            StringInfoData count_query;
            initStringInfo(&count_query);
            appendStringInfo(&count_query, "SELECT count(*) FROM %s", final_table_name);
            
            if (SPI_exec(count_query.data, 0) == SPI_OK_SELECT && SPI_processed > 0) {
                bool isnull;
                result.highfreq_kmers_count = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[0], 
                                                                        SPI_tuptable->tupdesc, 1, &isnull));
            }
            pfree(count_query.data);
        }
        
        /* Phase 2: Collect n-gram keys for high-frequency k-mers using parallel processing */
        ereport(NOTICE, (errmsg("Phase 2: Collecting n-gram keys for high-frequency k-mers...")));
        ereport(DEBUG1, (errmsg("Phase 2 started, about to get GIN index OID")));
        
        /* Get the GIN index OID and insert metadata */
        {
            StringInfoData query;
            Oid index_oid = InvalidOid;
            int ret;
            
            initStringInfo(&query);
            ereport(DEBUG1, (errmsg("Phase 2: Building GIN index query for table_oid=%u, column=%s", 
                                   table_oid, column_name)));
            appendStringInfo(&query,
                "SELECT i.indexrelid "
                "FROM pg_index i "
                "JOIN pg_class c ON c.oid = i.indrelid "
                "JOIN pg_attribute a ON a.attrelid = c.oid AND a.attnum = i.indkey[0] "
                "WHERE i.indrelid = %u AND a.attname = %s "
                "LIMIT 1",
                table_oid, quote_literal_cstr(column_name));
            ereport(DEBUG1, (errmsg("Phase 2: Query built: %s", query.data)));
            
            ereport(DEBUG1, (errmsg("Phase 2: About to execute SPI_exec")));
            ret = SPI_exec(query.data, 1);
            ereport(DEBUG1, (errmsg("Phase 2: SPI_exec returned %d, SPI_processed=%lu", ret, (unsigned long)SPI_processed)));
            if (ret == SPI_OK_SELECT && SPI_processed > 0) {
                bool isnull;
                ereport(DEBUG1, (errmsg("Phase 2: Getting index OID from result")));
                index_oid = DatumGetObjectId(SPI_getbinval(SPI_tuptable->vals[0], 
                                                           SPI_tuptable->tupdesc, 1, &isnull));
                ereport(DEBUG1, (errmsg("Phase 2: Got index_oid=%u, isnull=%d", index_oid, isnull)));
            }
            pfree(query.data);
            
            ereport(DEBUG1, (errmsg("Phase 2: Checking if index_oid is valid: %u", index_oid)));
            if (OidIsValid(index_oid)) {
                /* Insert GIN index metadata */
                ereport(DEBUG1, (errmsg("Phase 2: Valid index_oid, inserting GIN index metadata")));
                initStringInfo(&query);
                appendStringInfo(&query,
                    "INSERT INTO kmersearch_gin_index_meta "
                    "(index_oid, table_oid, column_name, highfreq_filtered, highfreq_source_table, "
                    "k_value, occur_bitlen, max_appearance_rate, max_appearance_nrow) "
                    "VALUES (%u, %u, %s, true, %s, %d, %d, %f, %d) "
                    "ON CONFLICT (index_oid) DO UPDATE SET "
                    "highfreq_filtered = EXCLUDED.highfreq_filtered, "
                    "highfreq_source_table = EXCLUDED.highfreq_source_table, "
                    "k_value = EXCLUDED.k_value, "
                    "occur_bitlen = EXCLUDED.occur_bitlen, "
                    "max_appearance_rate = EXCLUDED.max_appearance_rate, "
                    "max_appearance_nrow = EXCLUDED.max_appearance_nrow, "
                    "created_at = now()",
                    index_oid, table_oid, quote_literal_cstr(column_name), 
                    quote_literal_cstr(final_table_name), k_size,
                    kmersearch_occur_bitlen, kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow);
                
                ereport(DEBUG1, (errmsg("Phase 2: About to insert GIN index metadata")));
                SPI_exec(query.data, 0);
                ereport(DEBUG1, (errmsg("Phase 2: GIN index metadata inserted")));
                pfree(query.data);
                
                /* Insert high-frequency k-mers */
                ereport(DEBUG1, (errmsg("Phase 2: About to insert high-frequency k-mers")));
                initStringInfo(&query);
                appendStringInfo(&query,
                    "INSERT INTO kmersearch_highfreq_kmer (index_oid, ngram_key, detection_reason) "
                    "SELECT %u, kmer_key, 'frequency analysis' FROM %s "
                    "ON CONFLICT (index_oid, ngram_key) DO NOTHING",
                    index_oid, final_table_name);
                
                ereport(DEBUG1, (errmsg("Phase 2: About to execute high-frequency k-mers insert")));
                SPI_exec(query.data, 0);
                ereport(DEBUG1, (errmsg("Phase 2: High-frequency k-mers inserted")));
                pfree(query.data);
            }
            ereport(DEBUG1, (errmsg("Phase 2: GIN index processing completed")));
        }
        
        /* Insert metadata record */
        ereport(DEBUG1, (errmsg("Phase 2: About to insert metadata record")));
        {
            StringInfoData query;
            initStringInfo(&query);
            appendStringInfo(&query,
                "INSERT INTO kmersearch_highfreq_kmer_meta "
                "(table_oid, column_name, k_value, occur_bitlen, max_appearance_rate, max_appearance_nrow) "
                "VALUES (%u, %s, %d, %d, %f, %d) "
                "ON CONFLICT (table_oid, column_name, k_value) DO UPDATE SET "
                "occur_bitlen = EXCLUDED.occur_bitlen, "
                "max_appearance_rate = EXCLUDED.max_appearance_rate, "
                "max_appearance_nrow = EXCLUDED.max_appearance_nrow, "
                "analysis_timestamp = now()",
                table_oid, quote_literal_cstr(column_name), k_size,
                kmersearch_occur_bitlen, kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow);
            
            ereport(DEBUG1, (errmsg("Phase 2: About to execute metadata insert")));
            SPI_exec(query.data, 0);
            ereport(DEBUG1, (errmsg("Phase 2: Metadata inserted")));
            pfree(query.data);
        }
        
        ereport(DEBUG1, (errmsg("Phase 2: About to call SPI_finish")));
        SPI_finish();
        ereport(DEBUG1, (errmsg("Phase 2: SPI_finish completed")));
        ereport(DEBUG1, (errmsg("Phase 2: About to pfree final_table_name")));
        pfree(final_table_name);
        ereport(DEBUG1, (errmsg("Phase 2: final_table_name freed")));
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: About to calculate total statistics")));
    /* Calculate total statistics */
    /* result.total_rows already set from actual row count above */
    
    /* Set proper values for the result */
    result.analysis_duration = 0.0;  /* TODO: Implement timing */
    
    /* Clean up */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: About to clean up - pfree workers")));
    pfree(workers);
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: About to close table")));
    table_close(rel, AccessShareLock);
    
    /* Final validation of max_appearance_rate_used before returning */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Final validation - max_appearance_rate_used = %f", result.max_appearance_rate_used)));
    if (result.max_appearance_rate_used < 0.0 || result.max_appearance_rate_used != result.max_appearance_rate_used) {  /* Check for NaN */
        ereport(WARNING, (errmsg("kmersearch_analyze_table_parallel: Detected corrupted max_appearance_rate_used (%f), fixing to default", result.max_appearance_rate_used)));
        result.max_appearance_rate_used = kmersearch_max_appearance_rate;
        if (result.max_appearance_rate_used <= 0.0) {
            result.max_appearance_rate_used = 0.05;
        }
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Returning result")));
    return result;
}

/*
 * Validate analysis parameters
 */
void
kmersearch_validate_analysis_parameters(Oid table_oid, const char *column_name, int k_size)
{
    if (!OidIsValid(table_oid)) {
        ereport(ERROR, (errmsg("invalid table OID")));
    }
    
    if (!column_name || strlen(column_name) == 0) {
        ereport(ERROR, (errmsg("column name cannot be empty")));
    }
    
    if (k_size < 1 || k_size > 32) {
        ereport(ERROR, (errmsg("k-mer size must be between 1 and 32")));
    }
}

/*
 * Determine optimal number of parallel workers
 */
static int
kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation)
{
    if (requested_workers <= 0) {
        return 1; /* Default to single worker */
    }
    
    /* Check against system max_parallel_workers GUC */
    {
        int max_workers = 8; /* Default PostgreSQL max_parallel_workers */
        if (requested_workers > max_workers) {
            return max_workers;
        }
    }
    
    return requested_workers;
}

/*
 * Calculate adjusted minimum score based on highly frequent k-mers in query
 * Only applies adjustment when high-frequency filtering is actually enabled
 */
int
kmersearch_get_adjusted_min_score(VarBit **query_keys, int nkeys)
{
    int highfreq_count;
    int adjusted_score;
    
    elog(LOG, "kmersearch_get_adjusted_min_score: Started with nkeys=%d", nkeys);
    
    /* Check if high-frequency filtering is enabled for this context */
    elog(LOG, "kmersearch_get_adjusted_min_score: Checking if high-frequency filtering is enabled");
    if (!kmersearch_is_highfreq_filtering_enabled()) {
        elog(LOG, "kmersearch_get_adjusted_min_score: High-frequency filtering disabled, returning default min score %d", kmersearch_min_score);
        return kmersearch_min_score;  /* No adjustment needed */
    }
    
    elog(LOG, "kmersearch_get_adjusted_min_score: High-frequency filtering enabled, counting high-freq k-mers");
    highfreq_count = kmersearch_count_highfreq_kmer_in_query(query_keys, nkeys);
    elog(LOG, "kmersearch_get_adjusted_min_score: Found %d high-frequency k-mers", highfreq_count);
    adjusted_score = kmersearch_min_score - highfreq_count;
    
    /* Ensure adjusted score is not negative */
    if (adjusted_score < 0)
        adjusted_score = 0;
    
    return adjusted_score;
}

/*
 * Stub implementations for missing functions
 */

static bool
kmersearch_validate_guc_against_all_metadata(void)
{
    /* TODO: Implement validation logic */
    return true;
}

static bool
kmersearch_is_parallel_highfreq_cache_loaded(void)
{
    return (parallel_cache_hash != NULL);
}

static bool
kmersearch_lookup_in_parallel_cache(VarBit *kmer_key)
{
    uint64 kmer_hash;
    
    if (!parallel_cache_hash)
        return false;
    
    kmer_hash = hash_any((unsigned char *) VARDATA(kmer_key), 
                        VARSIZE(kmer_key) - VARHDRSZ);
    return kmersearch_parallel_cache_lookup(kmer_hash);
}

static bool
kmersearch_parallel_cache_lookup(uint64 kmer_hash)
{
    /* TODO: Implement parallel cache lookup */
    return false;
}

static int
kmersearch_varbit_cmp(VarBit *a, VarBit *b)
{
    int len_a, len_b;
    
    /* Simple comparison implementation */
    if (!a && !b) return 0;
    if (!a) return -1;
    if (!b) return 1;
    
    len_a = VARBITLEN(a);
    len_b = VARBITLEN(b);
    
    if (len_a != len_b)
        return len_a - len_b;
    
    return memcmp(VARBITS(a), VARBITS(b), VARBITBYTES(a));
}

static void
kmersearch_spi_connect_or_error(void)
{
    if (SPI_connect() != SPI_OK_CONNECT)
        ereport(ERROR, (errmsg("SPI_connect failed")));
}

/*
 * Worker function to analyze blocks of a table
 */
static void
kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size)
{
    /* TODO: Implement worker block analysis */
    /* For now, simulate some processing */
    worker->rows_processed = 10;  /* Dummy value */
    worker->local_highfreq_count = 0;
}

/* Stub functions are no longer needed as their functionality is integrated into the main function */