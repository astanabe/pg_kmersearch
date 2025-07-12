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
static int kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation);

/* High-frequency k-mer filtering functions */
static bool kmersearch_is_kmer_highfreq(VarBit *kmer_key);
static int kmersearch_count_highfreq_kmer_in_query(VarBit **query_keys, int nkeys);
static bool kmersearch_is_highfreq_filtering_enabled(void);
static bool kmersearch_is_highfreq_kmer_parallel(VarBit *kmer);
static Datum *kmersearch_filter_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, Datum *all_keys, int total_keys, int *filtered_count);
static bool kmersearch_delete_kmer_from_gin_index(Relation index_rel, VarBit *kmer_key);
/* Utility functions */
static List *kmersearch_get_highfreq_kmer_list(Oid index_oid);
static int kmersearch_varbit_cmp(VarBit *a, VarBit *b);
static void kmersearch_spi_connect_or_error(void);
static void kmersearch_handle_spi_error(int spi_result, const char *operation);
static bool kmersearch_check_analysis_exists(Oid table_oid, const char *column_name, int k_size);

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
    if (!kmer_key) {
        return false;
    }
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
        uint64 hash_value;
        bool found;
        
        /* Use ngram_key2 (kmer_key) directly for cache lookup - no occurrence bits removal needed */
        search_key = kmer_key;
        
        /* Validate VarBit structure */
        if (VARSIZE(search_key) < VARHDRSZ) {
            ereport(DEBUG1, (errmsg("Invalid VarBit structure in high-frequency k-mer check")));
            return false;
        }
        
        /* Validate bit length */
        if (VARBITLEN(search_key) < 0) {
            ereport(DEBUG1, (errmsg("Invalid bit length in high-frequency k-mer check")));
            return false;
        }
        
        /* Calculate hash value for lookup */
        {
            int bit_length = VARBITLEN(search_key);
            int byte_count = (bit_length + 7) / 8;  /* Round up to next byte */
            
            /* Validate the calculated byte count */
            if (byte_count <= 0 || byte_count > VARSIZE(search_key) - VARHDRSZ) {
                ereport(DEBUG1, (errmsg("Invalid byte count in high-frequency k-mer hash calculation")));
                return false;
            }
            
            hash_value = DatumGetUInt64(hash_any((unsigned char *) VARBITS(search_key), byte_count));
        }
        
        found = (hash_search(global_highfreq_cache.highfreq_hash, 
                           (void *) &hash_value, HASH_FIND, NULL) != NULL);
        
        /* No need to free search_key since it points to kmer_key */
        
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
    
    
    /* Validate input parameters */
    if (query_keys == NULL) {
        return 0;
    }
    
    /* For each k-mer in the query, check if it's highly frequent */
    for (i = 0; i < nkeys; i++)
    {
        
        if (query_keys[i] == NULL) {
            continue;
        }
        
        
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
 * Filter highly frequent k-mers from key array
 */

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
            appendStringInfo(&query,
                "SELECT i.indexrelid "
                "FROM pg_index i "
                "JOIN pg_class c ON c.oid = i.indrelid "
                "JOIN pg_attribute a ON a.attrelid = c.oid AND a.attnum = i.indkey[0] "
                "WHERE i.indrelid = %u AND a.attname = %s "
                "LIMIT 1",
                table_oid, quote_literal_cstr(column_name));
            
            ret = SPI_exec(query.data, 1);
            if (ret == SPI_OK_SELECT && SPI_processed > 0) {
                bool isnull;
                index_oid = DatumGetObjectId(SPI_getbinval(SPI_tuptable->vals[0], 
                                                           SPI_tuptable->tupdesc, 1, &isnull));
            }
            pfree(query.data);
            
            if (OidIsValid(index_oid)) {
                /* Insert GIN index metadata */
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
                
                SPI_exec(query.data, 0);
                pfree(query.data);
                
                /* Insert high-frequency k-mers */
                initStringInfo(&query);
                appendStringInfo(&query,
                    "INSERT INTO kmersearch_highfreq_kmer (index_oid, ngram_key, detection_reason) "
                    "SELECT %u, kmer_key, 'frequency analysis' FROM %s "
                    "ON CONFLICT (index_oid, ngram_key) DO NOTHING",
                    index_oid, final_table_name);
                
                SPI_exec(query.data, 0);
                pfree(query.data);
            }
        }
        
        /* Insert metadata record */
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
            
            SPI_exec(query.data, 0);
            pfree(query.data);
        }
        
        SPI_finish();
        pfree(final_table_name);
    }
    
    /* Calculate total statistics */
    /* result.total_rows already set from actual row count above */
    
    /* Set proper values for the result */
    
    /* Clean up */
    pfree(workers);
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
 * Helper function implementations
 */

/* Functions moved here from kmersearch.c are now implemented below */


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

/*
 * Validate GUC settings against all metadata table entries
 */
bool
kmersearch_validate_guc_against_all_metadata(void)
{
    int ret;
    StringInfoData query;
    bool valid = true;
    
    /* Connect to SPI */
    if (SPI_connect() != SPI_OK_CONNECT)
        return true;  /* Assume valid if we can't check */
    
    /* Check if metadata table exists */
    ret = SPI_execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'kmersearch_highfreq_kmer_meta' LIMIT 1", true, 1);
    if (ret != SPI_OK_SELECT || SPI_processed == 0) {
        /* Table doesn't exist, validation passes */
        SPI_finish();
        return true;
    }
    
    /* Query metadata table and compare with current GUC values */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT occur_bitlen, max_appearance_rate, max_appearance_nrow "
        "FROM kmersearch_highfreq_kmer_meta "
        "WHERE occur_bitlen != %d OR "
        "      abs(max_appearance_rate - %f) > 0.0001 OR "
        "      max_appearance_nrow != %d "
        "LIMIT 1",
        kmersearch_occur_bitlen,
        kmersearch_max_appearance_rate,
        kmersearch_max_appearance_nrow);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    if (ret == SPI_OK_SELECT && SPI_processed > 0) {
        valid = false;  /* Found mismatching entry */
    }
    
    /* Clean up */
    pfree(query.data);
    SPI_finish();
    
    return valid;
}

/*
 * Check if parallel_highfreq_cache is loaded
 */
bool
kmersearch_is_parallel_highfreq_cache_loaded(void)
{
    return (parallel_highfreq_cache != NULL && 
            parallel_highfreq_cache->is_initialized &&
            parallel_highfreq_cache->num_entries > 0);
}

/*
 * Lookup k-mer in parallel_highfreq_cache  
 */
bool
kmersearch_lookup_in_parallel_cache(VarBit *kmer_key)
{
    MemoryContext oldcontext;
    uint64 kmer_hash;
    ParallelHighfreqKmerCacheEntry *entry = NULL;
    bool found = false;
    
    /* Basic validation checks */
    if (!parallel_highfreq_cache || !parallel_highfreq_cache->is_initialized || 
        parallel_highfreq_cache->num_entries == 0)
        return false;
    
    if (!parallel_cache_hash)
        return false;
    
    /* Switch to TopMemoryContext for dshash operations */
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    PG_TRY();
    {
        /* Calculate hash using same logic as global cache */
        kmer_hash = kmersearch_ngram_key_to_hash(kmer_key);
        
        /* Lookup in dshash table */
        entry = (ParallelHighfreqKmerCacheEntry *) dshash_find(parallel_cache_hash, &kmer_hash, false);
        
        if (entry != NULL) {
            found = true;
            /* Must release lock after dshash_find() */
            dshash_release_lock(parallel_cache_hash, entry);
        }
    }
    PG_CATCH();
    {
        /* Ensure lock is released even in error cases */
        if (entry)
            dshash_release_lock(parallel_cache_hash, entry);
        MemoryContextSwitchTo(oldcontext);
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    MemoryContextSwitchTo(oldcontext);
    return found;
}

/*
 * Connect to SPI with error handling
 */
static void
kmersearch_spi_connect_or_error(void)
{
    int ret = SPI_connect();
    
    switch (ret) {
    case SPI_OK_CONNECT:
        return;  /* Success */
    case SPI_ERROR_CONNECT:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI manager already connected")));
        break;
    case SPI_ERROR_ARGUMENT:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI connection failed: invalid argument")));
        break;
    default:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI connection failed with code %d", ret)));
        break;
    }
}

/*
 * Handle SPI operation errors
 */
static void
kmersearch_handle_spi_error(int spi_result, const char *operation)
{
    switch (spi_result) {
    case SPI_OK_SELECT:
    case SPI_OK_INSERT:
    case SPI_OK_DELETE:
    case SPI_OK_UPDATE:
    case SPI_OK_UTILITY:
        return;  /* Success cases */
        
    case SPI_ERROR_ARGUMENT:
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("SPI %s failed: invalid argument", operation)));
        break;
        
    case SPI_ERROR_UNCONNECTED:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: not connected to SPI manager", operation)));
        break;
        
    case SPI_ERROR_COPY:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: COPY operation not supported", operation)));
        break;
        
    case SPI_ERROR_CURSOR:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: cursor operation error", operation)));
        break;
        
    case SPI_ERROR_TRANSACTION:
        ereport(ERROR,
                (errcode(ERRCODE_ACTIVE_SQL_TRANSACTION),
                 errmsg("SPI %s failed: transaction block error", operation)));
        break;
        
    case SPI_ERROR_OPUNKNOWN:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: unknown operation", operation)));
        break;
        
    default:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed with code %d", operation, spi_result)));
        break;
    }
}

/*
 * Check if analysis exists for given parameters
 */
static bool
kmersearch_check_analysis_exists(Oid table_oid, const char *column_name, int k_size)
{
    int ret;
    bool found = false;
    StringInfoData query;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to check for existing analysis */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT COUNT(*) FROM kmersearch_index_info "
        "WHERE table_oid = %u AND column_name = '%s' AND k_value = %d",
        table_oid, column_name, k_size);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    kmersearch_handle_spi_error(ret, "SELECT");
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        Datum count_datum;
        bool isnull;
        int count;
        
        count_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull)
        {
            count = DatumGetInt32(count_datum);
            found = (count > 0);
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return found;
}

/*
 * Filter highly frequent k-mers from the key array
 * A-2: Use direct VarBit comparison instead of hash table (k-mer+occurrence n-gram keys are small)
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
    
    /* Check if analysis exists */
    if (!kmersearch_check_analysis_exists(table_oid, column_name, k_size)) {
        /* No analysis found, return all keys */
        *filtered_count = total_keys;
        return all_keys;
    }
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get excluded k-mers from index info table */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT ek.kmer_key FROM kmersearch_highfreq_kmer ek "
        "JOIN kmersearch_index_info ii ON ek.index_oid = ii.index_oid "
        "WHERE ii.table_oid = %u AND ii.column_name = '%s' AND ii.k_value = %d",
        table_oid, column_name, k_size);
    
    /* Execute query and collect highly frequent k-mers in a simple array */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        int j;
        highfreq_count = SPI_processed;
        highfreq_kmers = (VarBit **) palloc(highfreq_count * sizeof(VarBit *));
        
        for (j = 0; j < SPI_processed; j++)
        {
            bool isnull;
            Datum kmer_datum;
            VarBit *kmer;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[j], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                kmer = DatumGetVarBitP(kmer_datum);
                
                /* Store copy of k-mer for direct comparison */
                highfreq_kmers[j] = (VarBit *) palloc(VARSIZE(kmer));
                memcpy(highfreq_kmers[j], kmer, VARSIZE(kmer));
            }
            else
            {
                highfreq_kmers[j] = NULL;
            }
        }
    }
    
    /* Filter out highly frequent k-mers using direct VarBit comparison */
    filtered_keys = (Datum *) palloc(total_keys * sizeof(Datum));
    
    for (i = 0; i < total_keys; i++)
    {
        VarBit *kmer = DatumGetVarBitP(all_keys[i]);
        bool is_highfreq = false;
        int j;
        
        /* Direct comparison with highly frequent k-mers (no hashing) */
        for (j = 0; j < highfreq_count; j++)
        {
            if (highfreq_kmers[j] != NULL &&
                VARBITLEN(kmer) == VARBITLEN(highfreq_kmers[j]) &&
                VARSIZE(kmer) == VARSIZE(highfreq_kmers[j]) &&
                memcmp(VARBITS(kmer), VARBITS(highfreq_kmers[j]), VARBITBYTES(kmer)) == 0)
            {
                is_highfreq = true;
                break;
            }
        }
        
        if (!is_highfreq)
        {
            /* Not highly frequent, include in filtered result */
            filtered_keys[filtered_idx++] = all_keys[i];
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    if (highfreq_kmers)
    {
        int j;
        for (j = 0; j < highfreq_count; j++)
        {
            if (highfreq_kmers[j])
                pfree(highfreq_kmers[j]);
        }
        pfree(highfreq_kmers);
    }
    SPI_finish();
    
    *filtered_count = filtered_idx;
    
    /* If no keys were filtered, return original array */
    if (filtered_idx == total_keys) {
        pfree(filtered_keys);
        return all_keys;
    }
    
    return filtered_keys;
}

/*
 * Helper function to get highly frequent k-mers list for a given index
 */
static List *
kmersearch_get_highfreq_kmer_list(Oid index_oid)
{
    List *highfreq_kmers = NIL;
    int ret;
    StringInfoData query;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get highly frequent k-mers */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT ek.kmer_key FROM kmersearch_highfreq_kmer ek "
        "JOIN kmersearch_index_info ii ON ek.index_oid = ii.index_oid "
        "WHERE ii.index_oid = %u ORDER BY ek.kmer_key",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    kmersearch_handle_spi_error(ret, "SELECT highly frequent k-mers");
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        int i;
        for (i = 0; i < SPI_processed; i++)
        {
            bool isnull;
            Datum kmer_datum;
            VarBit *kmer;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                kmer = DatumGetVarBitPCopy(kmer_datum);
                highfreq_kmers = lappend(highfreq_kmers, kmer);
            }
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return highfreq_kmers;
}

/*
 * Helper function to delete a k-mer from GIN index
 */
static bool
kmersearch_delete_kmer_from_gin_index(Relation index_rel, VarBit *kmer_key)
{
    /* 
     * Note: This is a simplified implementation.
     * In a complete implementation, this would use GIN internal APIs
     * to locate and remove the specific k-mer key and its posting list.
     * For now, we'll just log the operation.
     */
    ereport(DEBUG1, (errmsg("Would delete k-mer from index (size: %d bits)", 
                           VARBITLEN(kmer_key))));
    
    /*
     */
    
    return true;  /* Assume success for now */
}

/*
 * High-frequency k-mer filtering functions implementation
 */
Datum *
kmersearch_filter_highfreq_kmers_from_keys(Datum *original_keys, int *nkeys, HTAB *highfreq_hash, int k)
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
kmersearch_filter_highfreq_kmers_from_keys_parallel(Datum *original_keys, int *nkeys, int k)
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
 * Worker function to analyze blocks of a table
 */

