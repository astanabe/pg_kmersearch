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
static void kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, const char *final_table_name, int k_size, int threshold_rows);
static void kmersearch_persist_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, void *unused_table, int threshold_rows);
static void kmersearch_persist_highfreq_kmers_from_temp(Oid table_oid, const char *column_name, int k_size, const char *temp_table_name);

/* High-frequency k-mer filtering functions */
static bool kmersearch_is_kmer_highfreq(VarBit *kmer_key);
static int kmersearch_count_highfreq_kmer_in_query(VarBit **query_keys, int nkeys);
static bool kmersearch_is_highfreq_filtering_enabled(void);
static bool kmersearch_is_highfreq_kmer_parallel(VarBit *kmer);
static Datum *kmersearch_filter_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, Datum *all_keys, int total_keys, int *filtered_count);

/* K-mer collection and persistence functions */
static void kmersearch_collect_ngram_key2_for_highfreq_kmer(Oid table_oid, const char *column_name, int k_size, const char *highfreq_table_name);
static void kmersearch_worker_collect_ngram_key2(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size, const char *highfreq_table_name);
static void kmersearch_persist_collected_ngram_key2(Oid table_oid, const char *final_table_name);
static bool kmersearch_is_kmer_high_frequency(VarBit *ngram_key, int k_size, const char *highfreq_table_name);
static void kmersearch_persist_highfreq_kmers_metadata(Oid table_oid, const char *column_name, int k_size);

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
extern int kmersearch_min_score;
extern bool kmersearch_preclude_highfreq_kmer;

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
    KmerAnalysisResult result;
    
    /* Comprehensive parameter validation */
    kmersearch_validate_analysis_parameters(table_oid, column_name, k_size);
    
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
    values[0] = Int64GetDatum(result.total_rows);
    values[1] = Int32GetDatum(result.highfreq_kmers_count);
    values[2] = Int32GetDatum(result.parallel_workers_used);
    values[3] = Float8GetDatum(result.analysis_duration);
    values[4] = Float8GetDatum(result.max_appearance_rate_used);
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
    if (!kmer_key)
        return false;
    
    /* Step 1: Validate GUC settings against metadata table */
    if (!kmersearch_validate_guc_against_all_metadata()) {
        ereport(ERROR, 
                (errcode(ERRCODE_CONFIGURATION_LIMIT_EXCEEDED),
                 errmsg("Current GUC settings do not match kmersearch_highfreq_kmer_meta table"),
                 errhint("Current cache may be invalid. Please reload cache or run kmersearch_analyze_table() again.")));
    }

    /* Step 2: Check in global cache first */
    if (global_highfreq_cache.is_valid && global_highfreq_cache.highfreq_hash) {
        VarBit *search_key;
        bool found;
        
        /* Remove occurrence bits if present and search in cache */
        search_key = kmersearch_remove_occurrence_bits(kmer_key, kmersearch_kmer_size);
        found = (hash_search(global_highfreq_cache.highfreq_hash, 
                           (void *) search_key, HASH_FIND, NULL) != NULL);
        
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
    
    /* For each k-mer in the query, check if it's highly frequent */
    for (i = 0; i < nkeys; i++)
    {
        if (kmersearch_is_kmer_highfreq(query_keys[i]))
        {
            highfreq_count++;
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
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer "
            "WHERE table_oid = %u AND column_name = %s AND k_size = %d",
            table_oid, quote_literal_cstr(column_name), k_size);
    } else {
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer "
            "WHERE table_oid = %u AND column_name = %s",
            table_oid, quote_literal_cstr(column_name));
    }
    
    /* Execute deletion */
    ret = SPI_execute(query.data, false, 0);
    if (ret == SPI_OK_DELETE) {
        result.dropped_highfreq_kmers = SPI_processed;
        result.dropped_analyses = 1;
        result.freed_storage_bytes = SPI_processed * 64; /* Estimated */
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return result;
}

/*
 * Parallel table analysis implementation
 */
KmerAnalysisResult
kmersearch_analyze_table_parallel(Oid table_oid, const char *column_name, int k_size, int parallel_workers)
{
    KmerAnalysisResult result = {0};
    Relation rel;
    
    /* Open table */
    rel = table_open(table_oid, AccessShareLock);
    
    /* Determine actual parallel workers */
    result.parallel_workers_used = kmersearch_determine_parallel_workers(parallel_workers, rel);
    
    /* Record analysis parameters */
    result.max_appearance_rate_used = kmersearch_max_appearance_rate;
    result.max_appearance_nrow_used = kmersearch_max_appearance_nrow;
    
    /* TODO: Implement actual parallel analysis */
    result.total_rows = 0;
    result.highfreq_kmers_count = 0;
    result.analysis_duration = 0.0;
    
    /* Close table */
    table_close(rel, AccessShareLock);
    
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
    
    if (requested_workers > max_parallel_workers) {
        return max_parallel_workers;
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
    
    /* Check if high-frequency filtering is enabled for this context */
    if (!kmersearch_is_highfreq_filtering_enabled()) {
        return kmersearch_min_score;  /* No adjustment needed */
    }
    
    highfreq_count = kmersearch_count_highfreq_kmer_in_query(query_keys, nkeys);
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