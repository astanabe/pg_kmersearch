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
#include <sys/time.h>

/* PostgreSQL function info declarations for frequency functions */
PG_FUNCTION_INFO_V1(kmersearch_perform_highfreq_analysis);
PG_FUNCTION_INFO_V1(kmersearch_undo_highfreq_analysis);

/*
 * Forward declarations for internal functions
 */

/* Frequency analysis functions */
/* kmersearch_perform_highfreq_analysis_parallel declared in header */
static int kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation);

/* High-frequency k-mer filtering functions */
static bool kmersearch_is_kmer_highfreq(VarBit *kmer_key);
static int kmersearch_count_highfreq_kmer_in_query(VarBit **query_keys, int nkeys);
static bool kmersearch_is_highfreq_filtering_enabled(void);
/* Utility functions */
static int kmersearch_varbit_cmp(VarBit *a, VarBit *b);
static void kmersearch_spi_connect_or_error(void);

/* External functions are now declared in kmersearch.h */
extern HighfreqKmerCache global_highfreq_cache;

/* External parallel cache variables (defined in kmersearch_cache.c) */
extern ParallelHighfreqKmerCache *parallel_highfreq_cache;
extern dshash_table *parallel_cache_hash;

/*
 * K-mer frequency analysis functions
 */


/*
 * Main table analysis function with parallel support
 */
Datum
kmersearch_perform_highfreq_analysis(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_P(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    char *table_name = text_to_cstring(table_name_text);
    char *column_name = text_to_cstring(column_name_text);
    KmerAnalysisResult result = {0};  /* Initialize all fields to zero */
    Oid table_oid;
    int k_size;
    int parallel_workers;
    
    check_guc_initialization();
    
    /* Get table OID from table name */
    table_oid = RelnameGetRelid(table_name);
    if (!OidIsValid(table_oid))
    {
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("relation \"%s\" does not exist", table_name)));
    }
    
    /* Get configuration from GUC variables */
    k_size = kmersearch_kmer_size;
    parallel_workers = max_parallel_maintenance_workers;
    
    /* Comprehensive parameter validation */
    kmersearch_validate_analysis_parameters(table_oid, column_name, k_size);
    
    /* Log analysis start */
    ereport(NOTICE, (errmsg("Performing k-mer frequency analysis for k=%d", k_size)));
    ereport(NOTICE, (errmsg("Max appearance rate: %f, Max appearance nrow: %d", 
                           kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow)));
    
    /* Perform parallel analysis */
    result = kmersearch_perform_highfreq_analysis_parallel(table_oid, column_name, k_size, parallel_workers);
    
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
    ereport(DEBUG1, (errmsg("kmersearch_perform_highfreq_analysis: Converting result to return values")));
    ereport(DEBUG1, (errmsg("kmersearch_perform_highfreq_analysis: result.max_appearance_rate_used = %f", result.max_appearance_rate_used)));
    
    values[0] = Int64GetDatum(result.total_rows);
    values[1] = Int32GetDatum(result.highfreq_kmers_count);
    values[2] = Int32GetDatum(result.parallel_workers_used);
    values[3] = Float4GetDatum((float4)result.analysis_duration);  /* real type = 4-byte float */
    
    /* Extra validation for max_appearance_rate_used before conversion */
    if (result.max_appearance_rate_used < 0.0 || result.max_appearance_rate_used != result.max_appearance_rate_used) {
        ereport(WARNING, (errmsg("kmersearch_perform_highfreq_analysis: Detected corrupted max_appearance_rate_used (%f) during conversion, fixing", result.max_appearance_rate_used)));
        result.max_appearance_rate_used = 0.5;
    }
    
    values[4] = Float4GetDatum((float4)result.max_appearance_rate_used);  /* real type = 4-byte float */
    ereport(DEBUG1, (errmsg("kmersearch_perform_highfreq_analysis: Float4GetDatum conversion completed")));
    
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
kmersearch_undo_highfreq_analysis(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_P(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    char *table_name = text_to_cstring(table_name_text);
    char *column_name = text_to_cstring(column_name_text);
    DropAnalysisResult result;
    Oid table_oid;
    
    check_guc_initialization();
    
    /* Get table OID from table name */
    table_oid = RelnameGetRelid(table_name);
    if (!OidIsValid(table_oid))
    {
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("relation \"%s\" does not exist", table_name)));
    }
    
    /* Perform drop operation (delete all k-mer sizes for this table/column) */
    result = kmersearch_undo_highfreq_analysis_internal(table_oid, column_name, 0);
    
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
 * High-frequency k-mer filtering functions
 */


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
 * Filter highly frequent k-mers from key array
 */

/*
 * Analysis helper functions
 */

/*
 * Internal drop analysis implementation
 */
DropAnalysisResult
kmersearch_undo_highfreq_analysis_internal(Oid table_oid, const char *column_name, int k_size)
{
    DropAnalysisResult result = {0};
    int ret;
    StringInfoData query;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to delete high-frequency k-mer data */
    initStringInfo(&query);
    if (k_size > 0) {
        /* Delete specific k-mer size from highfreq_kmer table */
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer "
            "WHERE table_oid = %u AND column_name = %s "
            "AND EXISTS ("
            "  SELECT 1 FROM kmersearch_highfreq_kmer_meta "
            "  WHERE table_oid = %u AND column_name = %s AND kmer_size = %d"
            ")",
            table_oid, quote_literal_cstr(column_name), 
            table_oid, quote_literal_cstr(column_name), k_size);
    } else {
        /* Delete all k-mer sizes for this table/column */
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer "
            "WHERE table_oid = %u AND column_name = %s",
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
        /* Delete specific k-mer size from metadata table */
        appendStringInfo(&query,
            "DELETE FROM kmersearch_highfreq_kmer_meta "
            "WHERE table_oid = %u AND column_name = %s AND kmer_size = %d",
            table_oid, quote_literal_cstr(column_name), k_size);
    } else {
        /* Delete all k-mer sizes for this table/column from metadata table */
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
            "WHERE table_oid = %u AND column_name = %s AND kmer_size = %d",
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
kmersearch_perform_highfreq_analysis_parallel(Oid table_oid, const char *column_name, int k_size, int parallel_workers)
{
    KmerAnalysisResult result = {0};  /* Initialize all fields to zero */
    Relation rel;
    int num_workers;
    KmerWorkerState *workers;
    int threshold_rows;
    int rate_threshold;
    int i;
    
    /* Initialize result structure */
    ereport(DEBUG1, (errmsg("kmersearch_perform_highfreq_analysis_parallel: Initializing result structure")));
    memset(&result, 0, sizeof(KmerAnalysisResult));
    
    /* Initialize max_appearance_rate_used early to prevent corruption */
    result.max_appearance_rate_used = kmersearch_max_appearance_rate;
    if (result.max_appearance_rate_used <= 0.0) {
        result.max_appearance_rate_used = 0.5;  /* Default value */
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
        {
            char *table_name = get_rel_name(table_oid);
            const char *escaped_table_name = quote_identifier(table_name);
            appendStringInfo(&count_query, "SELECT COUNT(*) FROM %s", escaped_table_name);
        }
        
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
        
        /* Check for empty table */
        if (actual_row_count == 0) {
            ereport(WARNING,
                    (errcode(ERRCODE_DATA_EXCEPTION),
                     errmsg("Cannot perform frequency analysis on empty table"),
                     errhint("Insert data into the table before analysis")));
            result.highfreq_kmers_count = 0;
            return result;
        }
        
        /* Calculate threshold using both rate and nrow limits */
        rate_threshold = (int)(actual_row_count * kmersearch_max_appearance_rate);
        
        /* Validate rate-based threshold calculation */
        if (kmersearch_max_appearance_rate > 0.0 && rate_threshold == 0) {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("Rate-based threshold calculation resulted in 0 with non-zero max_appearance_rate"),
                     errdetail("total_rows=%ld, max_appearance_rate=%f, calculated_threshold=%d",
                               actual_row_count, kmersearch_max_appearance_rate, rate_threshold),
                     errhint("Increase table size or adjust kmersearch.max_appearance_rate to a larger value. "
                             "For %ld rows, minimum rate should be approximately %f",
                             actual_row_count, 1.0 / actual_row_count)));
        }
        
        if (kmersearch_max_appearance_nrow > 0) {
            /* Use the more restrictive (smaller) threshold between rate and nrow limits */
            threshold_rows = (rate_threshold < kmersearch_max_appearance_nrow) ? 
                             rate_threshold : kmersearch_max_appearance_nrow;
        } else {
            /* No nrow limit (0 means unlimited), use rate-based threshold only */
            threshold_rows = rate_threshold;
        }
    }
    
    /* Update max_appearance_rate_used with validation */
    ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Setting max_appearance_rate_used from %f to %f", result.max_appearance_rate_used, kmersearch_max_appearance_rate)));
    result.max_appearance_rate_used = kmersearch_max_appearance_rate;
    if (result.max_appearance_rate_used <= 0.0) {
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: max_appearance_rate_used was %f, setting to default 0.5", result.max_appearance_rate_used)));
        result.max_appearance_rate_used = 0.5;  /* Default value */
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
        /* Use a globally unique table name to avoid conflicts */
        char *final_table_name = kmersearch_generate_unique_temp_table_name("temp_kmer_final", -1);
        
        /* Connect SPI once for all operations with transaction management */
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Connecting to SPI for Phase 1")));
        kmersearch_spi_connect_or_error();
        
        /* Begin transaction for consistent k-mer analysis data */
        PG_TRY();
        {
            ereport(DEBUG1, (errmsg("Starting transaction for high-frequency k-mer analysis")));
            SPI_execute("BEGIN", false, 0);
        }
        PG_CATCH();
        {
            SPI_finish();
            PG_RE_THROW();
        }
        PG_END_TRY();
        
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: SPI connected for Phase 1")));        
        /* Create final aggregation table with structure matching worker tables */
        {
            StringInfoData query;
            const char *data_type;
            int ret;
            
            /* Always use bigint for k-mer hash values for consistency with merge function */
            data_type = "bigint";  /* Always store hash values as bigint for consistency */
            
            initStringInfo(&query);
            /* Drop table if it exists to avoid conflicts */
            appendStringInfo(&query, "DROP TABLE IF EXISTS %s CASCADE", final_table_name);
            ret = SPI_exec(query.data, 0);
            ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Dropped existing table (if any): %s, result: %d", final_table_name, ret)));
            pfree(query.data);
            
            /* Now create the table */
            initStringInfo(&query);
            appendStringInfo(&query, "CREATE TABLE %s (kmer_data %s PRIMARY KEY, frequency_count integer)", final_table_name, data_type);
            ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Creating temp table: %s", query.data)));
            
            ret = SPI_exec(query.data, 0);
            if (ret != SPI_OK_UTILITY) {
                ereport(ERROR, (errmsg("Failed to create temporary table %s, SPI result: %d", final_table_name, ret)));
            }
            ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Temp table created")));
            pfree(query.data);
        }
        
        /* Phase 1.5: Merge worker results using SQL aggregation */
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Starting worker results merge")));
        kmersearch_merge_worker_results_sql(workers, num_workers, final_table_name, k_size, threshold_rows);
        ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Worker results merged")));
        
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
            
            /* Insert GIN index metadata only if a GIN index exists */
            if (OidIsValid(index_oid)) {
                int gin_ret;
                
                initStringInfo(&query);
                appendStringInfo(&query,
                    "INSERT INTO kmersearch_gin_index_meta "
                    "(index_oid, table_oid, column_name, highfreq_filtered, highfreq_source_table, "
                    "kmer_size, occur_bitlen, max_appearance_rate, max_appearance_nrow) "
                    "VALUES (%u, %u, %s, true, %s, %d, %d, %f, %d) "
                    "ON CONFLICT (index_oid) DO UPDATE SET "
                    "highfreq_filtered = EXCLUDED.highfreq_filtered, "
                    "highfreq_source_table = EXCLUDED.highfreq_source_table, "
                    "kmer_size = EXCLUDED.kmer_size, "
                    "occur_bitlen = EXCLUDED.occur_bitlen, "
                    "max_appearance_rate = EXCLUDED.max_appearance_rate, "
                    "max_appearance_nrow = EXCLUDED.max_appearance_nrow, "
                    "created_at = now()",
                    index_oid, table_oid, quote_literal_cstr(column_name), 
                    quote_literal_cstr(final_table_name), k_size,
                    kmersearch_occur_bitlen, kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow);
                
                gin_ret = SPI_exec(query.data, 0);
                if (gin_ret != SPI_OK_INSERT && gin_ret != SPI_OK_UPDATE && gin_ret != SPI_OK_INSERT_RETURNING) {
                    ereport(ERROR, 
                            (errcode(ERRCODE_INTERNAL_ERROR),
                             errmsg("Failed to insert/update GIN index metadata"),
                             errdetail("SPI_exec returned %d for GIN metadata query", gin_ret),
                             errhint("Check kmersearch_gin_index_meta table structure and permissions")));
                } else {
                    ereport(DEBUG1, (errmsg("Successfully inserted/updated %lu GIN index metadata records", SPI_processed)));
                }
                pfree(query.data);
            }
            
            /* Always collect n-gram keys for high-frequency k-mers, regardless of GIN index existence */
            ereport(DEBUG1, (errmsg("kmersearch_analyze_table_parallel: Collecting n-gram keys for high-frequency k-mers from final table: %s", final_table_name)));
            kmersearch_collect_ngram_key2_for_highfreq_kmer(table_oid, column_name, k_size, final_table_name);
        }
        
        /* Insert metadata record */
        {
            StringInfoData query;
            int ret;
            
            initStringInfo(&query);
            appendStringInfo(&query,
                "INSERT INTO kmersearch_highfreq_kmer_meta "
                "(table_oid, column_name, kmer_size, occur_bitlen, max_appearance_rate, max_appearance_nrow) "
                "VALUES (%u, %s, %d, %d, %f, %d) "
                "ON CONFLICT (table_oid, column_name, kmer_size) DO UPDATE SET "
                "occur_bitlen = EXCLUDED.occur_bitlen, "
                "max_appearance_rate = EXCLUDED.max_appearance_rate, "
                "max_appearance_nrow = EXCLUDED.max_appearance_nrow, "
                "analysis_timestamp = now()",
                table_oid, quote_literal_cstr(column_name), k_size,
                kmersearch_occur_bitlen, kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow);
            
            ret = SPI_exec(query.data, 0);
            if (ret != SPI_OK_INSERT && ret != SPI_OK_UPDATE && ret != SPI_OK_INSERT_RETURNING) {
                ereport(ERROR, 
                        (errcode(ERRCODE_INTERNAL_ERROR),
                         errmsg("Failed to insert/update high-frequency k-mer metadata"),
                         errdetail("SPI_exec returned %d for metadata INSERT/UPDATE query", ret),
                         errhint("Check kmersearch_highfreq_kmer_meta table structure and permissions")));
            } else {
                ereport(DEBUG1, (errmsg("Successfully inserted/updated %lu high-frequency k-mer metadata records", SPI_processed)));
            }
            pfree(query.data);
        }
        
        /* Clean up the temporary table */
        {
            StringInfoData cleanup_query;
            initStringInfo(&cleanup_query);
            appendStringInfo(&cleanup_query, "DROP TABLE IF EXISTS %s", final_table_name);
            SPI_exec(cleanup_query.data, 0);
            pfree(cleanup_query.data);
        }
        
        /* Commit transaction for successful analysis completion */
        PG_TRY();
        {
            ereport(DEBUG1, (errmsg("Committing high-frequency k-mer analysis transaction")));
            SPI_execute("COMMIT", false, 0);
            ereport(NOTICE, (errmsg("High-frequency k-mer analysis transaction committed successfully")));
        }
        PG_CATCH();
        {
            ereport(WARNING, (errmsg("Failed to commit transaction, attempting rollback")));
            SPI_execute("ROLLBACK", false, 0);
            SPI_finish();
            pfree(final_table_name);
            PG_RE_THROW();
        }
        PG_END_TRY();
        
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
            result.max_appearance_rate_used = 0.5;
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
                 errhint("Current cache may be invalid. Please reload cache or run kmersearch_perform_highfreq_analysis() again.")));
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
 * Worker function to analyze blocks of a table
 */

