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
PG_FUNCTION_INFO_V1(kmersearch_perform_highfreq_analysis);
PG_FUNCTION_INFO_V1(kmersearch_undo_highfreq_analysis);

/*
 * Forward declarations for internal functions
 */

/* Frequency analysis functions */
/* kmersearch_perform_highfreq_analysis_parallel declared in header */

/* Utility functions */

/* Analysis dshash functions */
static int analysis_kmer_hash_compare(const void *a, const void *b, size_t size, void *arg);
static uint32 analysis_kmer_hash_hash(const void *key, size_t size, void *arg);
static uint32 kmersearch_uint16_identity_hash(const void *key, size_t keysize, void *arg);
static uint32 kmersearch_uint32_identity_hash(const void *key, size_t keysize, void *arg);
static bool kmersearch_create_analysis_dshash(int estimated_entries, int kmer_size);
static void kmersearch_cleanup_analysis_dshash(void);
/* kmersearch_is_kmer_hash_in_analysis_dshash is now extern and declared in kmersearch.h */
static void kmersearch_populate_analysis_dshash_from_workers(KmerWorkerState *workers, int num_workers, uint64 threshold_rows);
static int kmersearch_get_analysis_dshash_count(int k_size, uint64 threshold_rows);
static void kmersearch_insert_kmer2_as_uint_from_dshash(Oid table_oid, const char *column_name, int k_size, uint64 threshold_rows);

/* External functions are now declared in kmersearch.h */

/*
 * Create worker temporary table for n-gram keys
 */
void
create_worker_ngram_temp_table(const char *table_name)
{
    StringInfoData query;
    int ret;
    
    initStringInfo(&query);
    
    /* First try to drop table if it exists to avoid conflicts */
    appendStringInfo(&query, "DROP TABLE IF EXISTS %s", table_name);
    
    ret = SPI_exec(query.data, 0);
    if (ret < 0)
        ereport(ERROR, (errmsg("Failed to drop existing temp table %s", table_name)));
    
    /* Create the temp table */
    resetStringInfo(&query);
    appendStringInfo(&query,
        "CREATE TEMP TABLE %s ("
        "ngram_key varbit NOT NULL"
        ")", table_name);
    
    ret = SPI_exec(query.data, 0);
    if (ret < 0)
        ereport(ERROR, (errmsg("Failed to create temp table %s", table_name)));
    
    pfree(query.data);
}

/*
 * Persist high-frequency k-mers from temporary table to permanent table
 */
void
kmersearch_persist_highfreq_kmers_from_temp(Oid table_oid, const char *column_name, int k_size,
                                           const char *temp_table_name)
{
    StringInfoData query;
    int ret;
    
    initStringInfo(&query);
    
    /* Insert highly frequent k-mers into permanent table */
    /* Note: Converting integer values to bit strings for ngram_key */
    ereport(DEBUG1, (errmsg("Building INSERT query for %s with occur_bitlen=%d", temp_table_name, kmersearch_occur_bitlen)));
    
    if (k_size > 32) {
        /* k > 32 not supported */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k_size)));
    }
    
    /* For k <= 32, use single column approach */
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, ngram_key, detection_reason) "
        "SELECT %u, '%s', "
        "  CASE "
        "    WHEN %d <= 8 THEN (kmer_data::integer)::bit(32) || (frequency_count::integer)::bit(%d) "
        "    WHEN %d <= 16 THEN (kmer_data::bigint)::bit(64) || (frequency_count::integer)::bit(%d) "
        "    ELSE (kmer_data::bigint)::bit(64) || (frequency_count::integer)::bit(%d) "
        "  END AS ngram_key, "
        "  'high_frequency' "
        "FROM %s "
        "WHERE kmer_data IS NOT NULL AND frequency_count > 0",
        table_oid, column_name, 
        k_size, kmersearch_occur_bitlen,
        k_size, kmersearch_occur_bitlen,
        kmersearch_occur_bitlen,
        temp_table_name);
    
    ereport(DEBUG1, (errmsg("Generated INSERT query: %s", query.data)));
    
    SPI_connect();
    
    /* Validate source table has data before insertion */
    {
        StringInfoData count_query;
        int count_ret;
        
        initStringInfo(&count_query);
        appendStringInfo(&count_query, "SELECT COUNT(*) FROM %s", temp_table_name);
        
        count_ret = SPI_exec(count_query.data, 0);
        if (count_ret == SPI_OK_SELECT && SPI_processed == 1) {
            bool isnull;
            int64 source_count = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], 
                                                             SPI_tuptable->tupdesc, 1, &isnull));
            ereport(NOTICE, (errmsg("Source table %s contains %ld high-frequency k-mer records", 
                                    temp_table_name, source_count)));
            
            if (source_count == 0) {
                ereport(WARNING, (errmsg("No high-frequency k-mers found in source table %s", temp_table_name)));
            }
        }
        pfree(count_query.data);
    }
    
    /* Execute high-frequency k-mer insertion with error checking */
    ereport(DEBUG1, (errmsg("Executing high-frequency k-mer insertion query")));
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_INSERT && ret != SPI_OK_INSERT_RETURNING) {
        ereport(ERROR, 
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to insert high-frequency k-mers into kmersearch_highfreq_kmer"),
                 errdetail("SPI_exec returned %d for INSERT query", ret),
                 errhint("Check table structure and permissions")));
    }
    
    /* Validate insertion results */
    if (SPI_processed == 0) {
        ereport(WARNING, (errmsg("No records were inserted into kmersearch_highfreq_kmer. Check data compatibility and SQL query syntax.")));
    } else {
        ereport(NOTICE, (errmsg("Successfully inserted %lu high-frequency k-mer records", SPI_processed)));
    }
    
    /* Insert metadata record */
    pfree(query.data);
    initStringInfo(&query);
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer_meta "
        "(table_oid, column_name, kmer_size, occur_bitlen, max_appearance_rate, max_appearance_nrow) "
        "VALUES (%u, '%s', %d, %d, %f, %d) "
        "ON CONFLICT (table_oid, column_name, kmer_size) DO UPDATE SET "
        "occur_bitlen = EXCLUDED.occur_bitlen, "
        "max_appearance_rate = EXCLUDED.max_appearance_rate, "
        "max_appearance_nrow = EXCLUDED.max_appearance_nrow, "
        "analysis_timestamp = now()",
        table_oid, column_name, k_size, kmersearch_occur_bitlen, 
        kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow);
    
    /* Execute metadata insertion with error checking */
    ereport(DEBUG1, (errmsg("Executing metadata insertion query")));
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_INSERT && ret != SPI_OK_UPDATE && ret != SPI_OK_INSERT_RETURNING) {
        ereport(ERROR, 
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to insert/update metadata in kmersearch_highfreq_kmer_meta"),
                 errdetail("SPI_exec returned %d for metadata INSERT/UPDATE query", ret),
                 errhint("Check table structure and permissions")));
    }
    ereport(DEBUG1, (errmsg("Successfully inserted/updated %lu metadata records", SPI_processed)));
    
    SPI_finish();
    
    pfree(query.data);
}

/*
 * Add k-mer to buffer, flush to temp table if full
 */
void
kmersearch_add_to_buffer(KmerBuffer *buffer, KmerData kmer_data, const char *temp_table_name)
{
    CompactKmerFreq *entry;
    
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_buffer_to_table(buffer, temp_table_name);
    }
    
    /* Add new entry */
    entry = &buffer->entries[buffer->count];
    entry->kmer_data = kmer_data;
    entry->frequency_count = 1;
    entry->is_highfreq = false;
    buffer->count++;
}

/*
 * Flush buffer contents to temporary table
 */
void
kmersearch_flush_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    
    if (buffer->count == 0) return;
    
    /* Aggregate entries with same kmer_data before insertion */
    kmersearch_aggregate_buffer_entries(buffer);
    
    initStringInfo(&query);
    
    /* Build bulk INSERT statement */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        if (buffer->kmer_size <= 8) {
            appendStringInfo(&query, "(%u, %d)", 
                           buffer->entries[i].kmer_data.k8_data,
                           buffer->entries[i].frequency_count);
        } else if (buffer->kmer_size <= 16) {
            appendStringInfo(&query, "(%u, %d)", 
                           buffer->entries[i].kmer_data.k16_data,
                           buffer->entries[i].frequency_count);
        } else if (buffer->kmer_size <= 32) {
            appendStringInfo(&query, "(%lu, %d)", 
                           buffer->entries[i].kmer_data.k32_data,
                           buffer->entries[i].frequency_count);
        } else {
            /* k > 32 not supported */
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("k-mer length must be between 4 and 32"),
                     errdetail("Provided k-mer length: %d", buffer->kmer_size)));
        }
    }
    
    /* Handle conflict resolution */
    appendStringInfo(&query, " ON CONFLICT (kmer_data) DO UPDATE SET frequency_count = %s.frequency_count + EXCLUDED.frequency_count", temp_table_name);
    
    /* Execute the query */
    SPI_connect();
    SPI_exec(query.data, 0);
    SPI_finish();
    
    /* Reset buffer */
    buffer->count = 0;
    
    pfree(query.data);
}

/*
 * Aggregate duplicate k-mer entries in buffer
 */
void
kmersearch_aggregate_buffer_entries(KmerBuffer *buffer)
{
    int i, j;
    int write_pos = 0;
    bool merged;
    
    if (buffer->count <= 1) return;
    
    for (i = 0; i < buffer->count; i++) {
        merged = false;
        
        /* Check if this entry can be merged with any previous entry */
        for (j = 0; j < write_pos; j++) {
            bool same_kmer = false;
            
            if (buffer->kmer_size <= 8) {
                same_kmer = (buffer->entries[i].kmer_data.k8_data == buffer->entries[j].kmer_data.k8_data);
            } else if (buffer->kmer_size <= 16) {
                same_kmer = (buffer->entries[i].kmer_data.k16_data == buffer->entries[j].kmer_data.k16_data);
            } else if (buffer->kmer_size <= 32) {
                same_kmer = (buffer->entries[i].kmer_data.k32_data == buffer->entries[j].kmer_data.k32_data);
            } else {
                /* k > 32 not supported */
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                         errmsg("k-mer length must be between 4 and 32"),
                         errdetail("Provided k-mer length: %d", buffer->kmer_size)));
            }
            
            if (same_kmer) {
                /* Merge frequency counts */
                buffer->entries[j].frequency_count += buffer->entries[i].frequency_count;
                merged = true;
                break;
            }
        }
        
        /* If not merged, add to write position */
        if (!merged) {
            if (write_pos != i) {
                buffer->entries[write_pos] = buffer->entries[i];
            }
            write_pos++;
        }
    }
    
    /* Update count to reflect aggregated entries */
    buffer->count = write_pos;
}

/* Analysis context for managing dshash resources */
typedef struct KmerAnalysisContext
{
    dsm_segment *dsm_seg;
    dsa_area *dsa;
    dshash_table *hash;
    dsm_handle dsm_handle;
    dshash_table_handle hash_handle;
} KmerAnalysisContext;

/* Handles structure for passing to workers via shm_toc */
typedef struct KmerAnalysisHandles
{
    dsm_handle dsm_handle;
    dshash_table_handle hash_handle;
} KmerAnalysisHandles;

/* Forward declaration for block-based k-mer extraction */
static void kmersearch_extract_kmers_from_block(Oid table_oid, AttrNumber column_attnum, Oid column_type_oid,
                                               BlockNumber block, int kmer_size, dshash_table *hash);
static void kmersearch_update_kmer_counts_in_dshash(Datum sequence_datum, int kmer_size, dshash_table *hash, Oid column_type_oid);
static void create_analysis_dshash_resources(KmerAnalysisContext *ctx, int estimated_entries, int kmer_size);

/* Analysis-specific dshash resources for temp_kmer_final replacement */
/* IMPORTANT: These are only valid in the main process, NOT in parallel workers */
static dsm_segment *analysis_dsm_segment = NULL;
static dsa_area *analysis_dsa = NULL;
static dshash_table *analysis_highfreq_hash = NULL;

/* High-frequency k-mer dshash entry structure */
typedef struct AnalysisHighfreqKmerEntry {
    uint64 kmer_hash;     /* Hash key generated by kmersearch_get_kmer_hash */
    int frequency_count;  /* K-mer frequency count */
} AnalysisHighfreqKmerEntry;

/*
 * K-mer frequency analysis functions
 */


/*
 * Main table analysis function with parallel support
 */
Datum
kmersearch_perform_highfreq_analysis(PG_FUNCTION_ARGS)
{
    text *table_name_or_oid_text = PG_GETARG_TEXT_P(0);
    text *column_name_or_attnum_text = PG_GETARG_TEXT_P(1);
    KmerAnalysisResult result = {0};  /* Initialize all fields to zero */
    Oid table_oid;
    char *column_name;
    AttrNumber column_attnum;
    int parallel_workers;
    char *table_str;
    char *column_str;
    char *endptr;
    unsigned long oid_val;
    long attnum_val;
    
    /* This function must not be called from parallel workers */
    if (IsParallelWorker())
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
                 errmsg("cannot execute kmersearch_perform_highfreq_analysis() in a parallel worker")));
    
    check_guc_initialization();
    
    /* Parse table identifier (name or OID) */
    table_str = text_to_cstring(table_name_or_oid_text);
    oid_val = strtoul(table_str, &endptr, 10);
    
    if (*endptr == '\0' && oid_val != 0 && OidIsValid(oid_val))
    {
        /* Input is OID */
        table_oid = (Oid)oid_val;
    }
    else
    {
        /* Input is table name */
        table_oid = RelnameGetRelid(table_str);
        if (!OidIsValid(table_oid))
        {
            ereport(ERROR,
                    (errcode(ERRCODE_UNDEFINED_TABLE),
                     errmsg("relation \"%s\" does not exist", table_str)));
        }
    }
    
    /* Parse column identifier (name or attnum) */
    column_str = text_to_cstring(column_name_or_attnum_text);
    attnum_val = strtol(column_str, &endptr, 10);
    
    if (*endptr == '\0' && attnum_val > 0)
    {
        /* Input is attnum */
        Relation rel;
        TupleDesc tupdesc;
        
        column_attnum = (AttrNumber)attnum_val;
        rel = table_open(table_oid, AccessShareLock);
        tupdesc = RelationGetDescr(rel);
        
        if (column_attnum <= 0 || column_attnum > tupdesc->natts)
        {
            table_close(rel, AccessShareLock);
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_COLUMN_REFERENCE),
                     errmsg("invalid column number %ld", attnum_val)));
        }
        
        column_name = pstrdup(NameStr(TupleDescAttr(tupdesc, column_attnum - 1)->attname));
        table_close(rel, AccessShareLock);
    }
    else
    {
        /* Input is column name */
        column_name = column_str;
        /* Verify column exists and get attnum */
        column_attnum = get_attnum(table_oid, column_name);
        if (column_attnum == InvalidAttrNumber)
        {
            ereport(ERROR,
                    (errcode(ERRCODE_UNDEFINED_COLUMN),
                     errmsg("column \"%s\" does not exist", column_name)));
        }
    }
    
    /* Get configuration from GUC variables */
    /* PostgreSQL will automatically limit based on max_parallel_workers and max_parallel_maintenance_workers */
    parallel_workers = max_parallel_maintenance_workers;
    
    /* Comprehensive parameter validation */
    kmersearch_validate_analysis_parameters(table_oid, column_name, kmersearch_kmer_size);
    
    /* Log analysis start */
    ereport(NOTICE, (errmsg("Performing k-mer frequency analysis for k=%d", kmersearch_kmer_size)));
    ereport(NOTICE, (errmsg("Max appearance rate: %f, Max appearance nrow: %d", 
                           kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow)));
    
    /* Perform parallel analysis */
    result = kmersearch_perform_highfreq_analysis_parallel(table_oid, column_name, kmersearch_kmer_size, parallel_workers);
    
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
    
    /* Extra validation for max_appearance_rate_used before conversion */
    if (result.max_appearance_rate_used < 0.0 || result.max_appearance_rate_used != result.max_appearance_rate_used) {
        ereport(WARNING, (errmsg("kmersearch_perform_highfreq_analysis: Detected corrupted max_appearance_rate_used (%f) during conversion, fixing", result.max_appearance_rate_used)));
        result.max_appearance_rate_used = 0.5;
    }
    
    values[3] = Float4GetDatum((float4)result.max_appearance_rate_used);  /* real type = 4-byte float */
    ereport(DEBUG1, (errmsg("kmersearch_perform_highfreq_analysis: Float4GetDatum conversion completed")));
    
    values[4] = Int32GetDatum(result.max_appearance_nrow_used);
    
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
    
    /* This function must not be called from parallel workers */
    if (IsParallelWorker())
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
                 errmsg("cannot execute kmersearch_undo_highfreq_analysis() in a parallel worker")));
    
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
int
kmersearch_count_highfreq_kmer_in_query(VarBit **query_keys, int nkeys)
{
    int highfreq_count = 0;
    int i;
    
    /* Debug: Log function entry */
    elog(DEBUG1, "kmersearch_count_highfreq_kmer_in_query: called with nkeys = %d", nkeys);
    
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
        
        /* Debug: Log each key */
        elog(DEBUG1, "kmersearch_count_highfreq_kmer_in_query: processing key[%d], bit length = %d", 
             i, VARBITLEN(query_keys[i]));
        
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
bool
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
kmersearch_perform_highfreq_analysis_parallel(Oid table_oid, const char *column_name, int k_size, int requested_workers)
{
    KmerAnalysisResult result = {0};  /* Initialize all fields to zero */
    ParallelContext *pcxt = NULL;
    KmerAnalysisSharedState *shared_state = NULL;
    KmerAnalysisContext analysis_ctx = {0};
    KmerAnalysisHandles *handles = NULL;
    shm_toc *toc;
    char *shm_pointer;
    bool table_locked = false;
    int64 total_rows;
    uint64 threshold_rows;
    uint64 rate_based_threshold;
    
    elog(DEBUG1, "kmersearch_perform_highfreq_analysis_parallel: Starting parallel analysis");
    
    
    PG_TRY();
    {
        /* Pre-validation */
        Relation rel;
        AttrNumber column_attnum;
        Form_pg_attribute attr;
        Oid column_type_oid;
        BlockNumber total_blocks;
        StringInfoData count_query;
        int ret;
        bool isnull;
        Datum count_datum;
        
        total_rows = 0; /* Initialize */
        
        rel = table_open(table_oid, AccessShareLock);
        column_attnum = get_attnum(table_oid, column_name);
        if (column_attnum == InvalidAttrNumber)
            elog(ERROR, "Column \"%s\" does not exist", column_name);
        
        attr = TupleDescAttr(RelationGetDescr(rel), column_attnum - 1);
        column_type_oid = attr->atttypid;
        if (column_type_oid != TypenameGetTypid("dna2") && 
            column_type_oid != TypenameGetTypid("dna4"))
            elog(ERROR, "Column must be DNA2 or DNA4 type");
        
        total_blocks = RelationGetNumberOfBlocks(rel);
        table_close(rel, AccessShareLock);
        
        /* Get row count */
        initStringInfo(&count_query);
        {
            char *table_name_str = get_rel_name(table_oid);
            const char *escaped_table_name = quote_identifier(table_name_str);
            appendStringInfo(&count_query, "SELECT COUNT(*) FROM %s", escaped_table_name);
        }
        
        if (SPI_connect() == SPI_OK_CONNECT) {
            ret = SPI_exec(count_query.data, 0);
            if (ret == SPI_OK_SELECT && SPI_processed == 1) {
                count_datum = SPI_getbinval(SPI_tuptable->vals[0], 
                                                 SPI_tuptable->tupdesc, 1, &isnull);
                if (!isnull) {
                    total_rows = DatumGetInt64(count_datum);
                }
            }
            SPI_finish();
        }
        pfree(count_query.data);
        
        result.total_rows = total_rows;
        result.max_appearance_rate_used = kmersearch_max_appearance_rate;
        result.max_appearance_nrow_used = kmersearch_max_appearance_nrow;
        
        /* Table lock acquisition */
        LockRelationOid(table_oid, ExclusiveLock);
        table_locked = true;
        
        /* Create dshash resources BEFORE entering parallel mode */
        /* Ensure global variables are NULL before entering parallel mode */
        analysis_dsm_segment = NULL;
        analysis_dsa = NULL;
        analysis_highfreq_hash = NULL;
        
        create_analysis_dshash_resources(&analysis_ctx, 100000, k_size);
        
        /* Enter parallel mode */
        EnterParallelMode();
        
        /* Now set global pointers AFTER entering parallel mode so workers don't inherit them */
        analysis_dsm_segment = analysis_ctx.dsm_seg;
        analysis_dsa = analysis_ctx.dsa;
        analysis_highfreq_hash = analysis_ctx.hash;
        
        /* Create parallel context */
        pcxt = CreateParallelContext("pg_kmersearch", "kmersearch_analysis_worker", requested_workers);
        
        /* DSM size estimation - estimate each chunk separately as per PostgreSQL best practice */
        shm_toc_estimate_chunk(&pcxt->estimator, MAXALIGN(sizeof(KmerAnalysisSharedState)));
        shm_toc_estimate_chunk(&pcxt->estimator, MAXALIGN(sizeof(KmerAnalysisHandles)));
        shm_toc_estimate_keys(&pcxt->estimator, 2); /* SHARED_STATE, HANDLES */
        
        elog(DEBUG1, "kmersearch_perform_highfreq_analysis_parallel: DSM estimation details:");
        elog(DEBUG1, "  - sizeof(KmerAnalysisSharedState) = %zu, MAXALIGN = %zu", 
             sizeof(KmerAnalysisSharedState), MAXALIGN(sizeof(KmerAnalysisSharedState)));
        elog(DEBUG1, "  - sizeof(dsm_handle) = %zu, MAXALIGN = %zu", 
             sizeof(dsm_handle), MAXALIGN(sizeof(dsm_handle)));
        elog(DEBUG1, "  - sizeof(dshash_table_handle) = %zu, MAXALIGN = %zu", 
             sizeof(dshash_table_handle), MAXALIGN(sizeof(dshash_table_handle)));
        
        /* Initialize DSM */
        elog(DEBUG1, "kmersearch_perform_highfreq_analysis_parallel: Initializing parallel DSM");
        InitializeParallelDSM(pcxt);
        toc = pcxt->toc;
        
        /* Set up shared state */
        shm_pointer = shm_toc_allocate(toc, sizeof(KmerAnalysisSharedState));
        if (!shm_pointer)
        {
            elog(ERROR, "Failed to allocate memory for shared state in shm_toc");
        }
        shared_state = (KmerAnalysisSharedState *)shm_pointer;
        memset(shared_state, 0, sizeof(KmerAnalysisSharedState));  /* Clear all fields first */
        LWLockInitialize(&shared_state->mutex.lock, LWTRANCHE_KMERSEARCH_ANALYSIS);
        shared_state->num_workers = requested_workers;
        shared_state->table_oid = table_oid;
        shared_state->column_attnum = column_attnum;
        shared_state->column_type_oid = column_type_oid;
        shared_state->kmer_size = k_size;
        shared_state->batch_size = kmersearch_highfreq_analysis_batch_size;
        shared_state->all_processed = false;
        shared_state->next_block = 0;
        shared_state->total_blocks = total_blocks;
        shared_state->worker_error_occurred = false;
        elog(DEBUG1, "kmersearch_perform_highfreq_analysis_parallel: Initialized shared state - next_block=%u, total_blocks=%u",
             shared_state->next_block, shared_state->total_blocks);
        shm_toc_insert(toc, KMERSEARCH_KEY_SHARED_STATE, shared_state);
        
        /* Store handles for workers */
        handles = (KmerAnalysisHandles *)shm_toc_allocate(toc, sizeof(KmerAnalysisHandles));
        if (!handles)
        {
            elog(ERROR, "Failed to allocate memory for analysis handles in shm_toc");
        }
        handles->dsm_handle = analysis_ctx.dsm_handle;
        handles->hash_handle = analysis_ctx.hash_handle;
        shm_toc_insert(toc, KMERSEARCH_KEY_HANDLES, handles);
        
        elog(DEBUG1, "kmersearch_perform_highfreq_analysis_parallel: Handles stored in shm_toc");
        
        /* Launch parallel workers */
        LaunchParallelWorkers(pcxt);
        result.parallel_workers_used = pcxt->nworkers_launched;
        elog(NOTICE, "kmersearch_perform_highfreq_analysis: Launched %d parallel workers", 
             pcxt->nworkers_launched);
        
        /* Wait for workers to complete */
        WaitForParallelWorkersToFinish(pcxt);
        
        /* Check for worker errors */
        if (shared_state->worker_error_occurred)
            elog(ERROR, "Parallel worker error: %s", shared_state->error_message);
        
        /* Calculate threshold based on GUC variables */
        rate_based_threshold = (uint64)(result.total_rows * kmersearch_max_appearance_rate);
        
        /* Determine final threshold: min of rate-based and nrow-based (excluding 0) */
        if (kmersearch_max_appearance_nrow > 0) {
            /* Both are set, use the smaller one */
            threshold_rows = (rate_based_threshold < kmersearch_max_appearance_nrow) ? 
                            rate_based_threshold : kmersearch_max_appearance_nrow;
        } else {
            /* Only rate-based threshold */
            threshold_rows = rate_based_threshold;
        }
        elog(DEBUG1, "Threshold calculation: total_rows=%ld, rate=%.2f, rate_based=%lu, nrow=%d, final=%lu",
             result.total_rows, kmersearch_max_appearance_rate, (unsigned long)rate_based_threshold, 
             kmersearch_max_appearance_nrow, (unsigned long)threshold_rows);
        
        /* Update max_appearance_nrow_used to the actual threshold used */
        result.max_appearance_nrow_used = threshold_rows;
        
        /* Get high-frequency k-mer count before exiting parallel mode */
        result.highfreq_kmers_count = kmersearch_get_analysis_dshash_count(k_size, threshold_rows);
        
        /* Clean up parallel context and exit parallel mode BEFORE any SQL operations */
        DestroyParallelContext(pcxt);
        pcxt = NULL; /* Mark as destroyed to prevent double-free */
        ExitParallelMode();
        
        /* Now we can safely execute SQL operations */
        /* Save results - only in main process (though we should always be main process here) */
        kmersearch_spi_connect_or_error();
        kmersearch_insert_kmer2_as_uint_from_dshash(table_oid, column_name, k_size, threshold_rows);
        
        /* Insert metadata */
        {
            StringInfoData query;
            int meta_ret;
            
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
            
            meta_ret = SPI_exec(query.data, 0);
            if (meta_ret != SPI_OK_INSERT && meta_ret != SPI_OK_UPDATE)
                elog(ERROR, "Failed to insert/update metadata");
            
            pfree(query.data);
        }
        
        SPI_finish();
    }
    PG_CATCH();
    {
        /* Error cleanup */
        kmersearch_cleanup_analysis_dshash();
        if (pcxt) {
            DestroyParallelContext(pcxt);
        }
        if (table_locked) {
            UnlockRelationOid(table_oid, ExclusiveLock);
        }
        ExitParallelMode();
        
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    /* Normal cleanup - already done above, just cleanup dshash and unlock */
    kmersearch_cleanup_analysis_dshash();
    UnlockRelationOid(table_oid, ExclusiveLock);
    
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
 * Helper function implementations
 */

/* Functions moved here from kmersearch.c are now implemented below */



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
void
kmersearch_spi_connect_or_error(void)
{
    int ret;
    
    /* Never connect SPI in parallel workers */
    if (IsParallelWorker())
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
                 errmsg("cannot connect to SPI in a parallel worker")));
    
    ret = SPI_connect();
    
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
bool
kmersearch_is_kmer_highfreq(VarBit *kmer_key)
{
    uint64 kmer_uint = 0;
    bool is_highfreq = false;
    int ret;
    int kmer_bit_len;
    
    if (!kmer_key) {
        return false;
    }
    
    /* Debug: Log input key information */
    kmer_bit_len = VARBITLEN(kmer_key);
    elog(DEBUG2, "kmersearch_is_kmer_highfreq: input key bit length = %d, kmersearch_occur_bitlen = %d, kmersearch_kmer_size = %d",
         kmer_bit_len, kmersearch_occur_bitlen, kmersearch_kmer_size);
    
    /* Determine if input is ngram_key2 (with occurrence count) or kmer2 (without occurrence count)
     * by checking if the bit length is exactly k-mer size * 2 (kmer2) or larger (ngram_key2) */
    if (kmer_bit_len == kmersearch_kmer_size * 2) {
        /* Input is kmer2 format (without occurrence count) */
        elog(DEBUG2, "kmersearch_is_kmer_highfreq: input is kmer2 format (no occurrence count)");
        
        /* For kmer2, check that k-mer bit length is at least 8 */
        if (kmer_bit_len < 8) {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("k-mer key is too short"),
                     errdetail("k-mer bit length is %d, but must be at least 8 bits.", kmer_bit_len)));
        }
    } else if (kmer_bit_len == kmersearch_kmer_size * 2 + kmersearch_occur_bitlen) {
        /* Input is ngram_key2 format (with occurrence count) */
        elog(DEBUG2, "kmersearch_is_kmer_highfreq: input is ngram_key2 format (with occurrence count)");
        
        /* For ngram_key2, check that k-mer portion is at least 8 bits */
        if (kmersearch_kmer_size * 2 < 8) {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("k-mer portion of key is too short"),
                     errdetail("k-mer size is %d, resulting in %d-bit k-mer, but must be at least 8 bits.", 
                              kmersearch_kmer_size, kmersearch_kmer_size * 2)));
        }
        
        /* Extract kmer2 portion by creating a new VarBit with only k-mer bits */
        {
            int kmer2_bits = kmersearch_kmer_size * 2;
            int kmer2_bytes = (kmer2_bits + 7) / 8;
            VarBit *kmer2_only;
            bits8 *src_data, *dst_data;
            int i;
            
            /* Allocate new VarBit for kmer2 only */
            kmer2_only = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer2_bytes);
            SET_VARSIZE(kmer2_only, VARHDRSZ + sizeof(int32) + kmer2_bytes);
            VARBITLEN(kmer2_only) = kmer2_bits;
            
            src_data = VARBITS(kmer_key);
            dst_data = VARBITS(kmer2_only);
            
            /* Copy only the k-mer bits (first kmer2_bits) */
            for (i = 0; i < kmer2_bits; i++) {
                int byte_pos = i / 8;
                int bit_offset = i % 8;
                
                if (src_data[byte_pos] & (0x80 >> bit_offset)) {
                    dst_data[byte_pos] |= (0x80 >> bit_offset);
                }
            }
            
            /* Recursively call with kmer2 format */
            is_highfreq = kmersearch_is_kmer_highfreq(kmer2_only);
            
            /* Clean up */
            pfree(kmer2_only);
            
            return is_highfreq;
        }
    } else {
        /* Unexpected bit length */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("invalid k-mer key format"),
                 errdetail("Key bit length is %d, expected either %d (kmer2) or %d (ngram_key2).", 
                          kmer_bit_len, kmersearch_kmer_size * 2, 
                          kmersearch_kmer_size * 2 + kmersearch_occur_bitlen)));
    }
    
    /* Step 1: Validate GUC settings against metadata table */
    if (!kmersearch_validate_guc_against_all_metadata()) {
        ereport(ERROR, 
                (errcode(ERRCODE_CONFIGURATION_LIMIT_EXCEEDED),
                 errmsg("Current GUC settings do not match kmersearch_highfreq_kmer_meta table"),
                 errhint("Current cache may be invalid. Please reload cache or run kmersearch_perform_highfreq_analysis() again.")));
    }
    
    /* kmer_key is kmer2 (without occurrence count), convert directly to uint */
    if (kmersearch_kmer_size <= 8) {
        uint16 value = 0;
        kmersearch_convert_kmer2_to_uint16(kmer_key, &value);
        kmer_uint = value;
    } else if (kmersearch_kmer_size <= 16) {
        uint32 value = 0;
        kmersearch_convert_kmer2_to_uint32(kmer_key, &value);
        kmer_uint = value;
    } else {
        kmersearch_convert_kmer2_to_uint64(kmer_key, &kmer_uint);
    }
    
    /* Priority 1: Check in global cache */
    if (global_highfreq_cache.is_valid && global_highfreq_cache.highfreq_hash) {
        is_highfreq = kmersearch_lookup_kmer2_as_uint_in_global_cache(kmer_uint, NULL, NULL);
        return is_highfreq;
    }
    
    /* Priority 2: Check in parallel cache */
    if (kmersearch_is_parallel_highfreq_cache_loaded()) {
        is_highfreq = kmersearch_lookup_kmer2_as_uint_in_parallel_cache(kmer_uint, NULL, NULL);
        return is_highfreq;
    }
    
    /* Priority 3: Check kmersearch_highfreq_kmer table directly */
    ret = SPI_connect();
    if (ret == SPI_OK_CONNECT) {
        StringInfoData query;
        
        initStringInfo(&query);
        
        /* Build query based on k-mer size */
        if (kmersearch_kmer_size <= 8) {
            appendStringInfo(&query,
                "SELECT 1 FROM kmersearch_highfreq_kmer "
                "WHERE kmer2_as_uint = %u "
                "LIMIT 1",
                (unsigned int)kmer_uint);
        } else if (kmersearch_kmer_size <= 16) {
            appendStringInfo(&query,
                "SELECT 1 FROM kmersearch_highfreq_kmer "
                "WHERE kmer2_as_uint = %u "
                "LIMIT 1", 
                (unsigned int)kmer_uint);
        } else {
            appendStringInfo(&query,
                "SELECT 1 FROM kmersearch_highfreq_kmer "
                "WHERE kmer2_as_uint = %lu "
                "LIMIT 1",
                kmer_uint);
        }
        
        ret = SPI_execute(query.data, true, 1);
        if (ret == SPI_OK_SELECT && SPI_processed > 0) {
            is_highfreq = true;
        }
        
        pfree(query.data);
        SPI_finish();
    }
    
    return is_highfreq;
}

/*
 * Analysis dshash functions implementation
 */

/*
 * Hash function for uint64 k-mer hash values
 */
static uint32
analysis_kmer_hash_hash(const void *key, size_t size, void *arg)
{
    const uint64 *hash_key = (const uint64 *) key;
    return (uint32) (*hash_key ^ (*hash_key >> 32));
}

/*
 * Compare function for uint64 k-mer hash values
 */
static int
analysis_kmer_hash_compare(const void *a, const void *b, size_t size, void *arg)
{
    const uint64 *hash_a = (const uint64 *) a;
    const uint64 *hash_b = (const uint64 *) b;
    
    if (*hash_a < *hash_b)
        return -1;
    else if (*hash_a > *hash_b)
        return 1;
    else
        return 0;
}

/*
 * Identity hash function for uint16 (kmer2_as_uint is already a hash value)
 */
static uint32
kmersearch_uint16_identity_hash(const void *key, size_t keysize, void *arg)
{
    return (uint32)(*(const uint16 *)key);
}

/*
 * Identity hash function for uint32 (kmer2_as_uint is already a hash value)
 */
static uint32
kmersearch_uint32_identity_hash(const void *key, size_t keysize, void *arg)
{
    return *(const uint32 *)key;
}

/*
 * Create DSM segment, DSA area, and dshash table for high-frequency k-mer analysis
 * This must be called before parallel workers are launched
 */
static bool
kmersearch_create_analysis_dshash(int estimated_entries, int kmer_size)
{
    dshash_parameters params;
    Size segment_size;
    
    /* Clean up any existing resources first */
    if (analysis_dsm_segment != NULL || analysis_dsa != NULL || analysis_highfreq_hash != NULL) {
        kmersearch_cleanup_analysis_dshash();
    }
    
    /* Calculate required DSM segment size */
    segment_size = 1024 * 1024;  /* Start with 1MB */
    if (estimated_entries > 10000) {
        segment_size = estimated_entries * 64;  /* Estimate 64 bytes per entry */
    }
    
    /* Limit maximum size to avoid out of shared memory */
    /* With 128MB shared_buffers, we need to be very conservative */
    if (segment_size > 4 * 1024 * 1024) {
        segment_size = 4 * 1024 * 1024;  /* Maximum 4MB for 128MB shared_buffers */
    }
    /* Ensure minimum size for DSA overhead */
    if (segment_size < 1 * 1024 * 1024) {
        segment_size = 1 * 1024 * 1024;  /* Minimum 1MB */
    }
    
    /* Step 1: Create DSM segment */
    PG_TRY();
    {
        analysis_dsm_segment = dsm_create(segment_size, 0);
        if (analysis_dsm_segment == NULL) {
            ereport(ERROR, (errmsg("Failed to create DSM segment for analysis")));
        }
        
        /* Pin the DSM segment to prevent automatic cleanup */
        dsm_pin_mapping(analysis_dsm_segment);
    }
    PG_CATCH();
    {
        ereport(ERROR, (errmsg("Failed to create or pin DSM segment for analysis")));
        return false;
    }
    PG_END_TRY();
    
    /* Step 2: Create DSA area */
    PG_TRY();
    {
        analysis_dsa = dsa_create_in_place(dsm_segment_address(analysis_dsm_segment), 
                                          segment_size, 
                                          LWTRANCHE_KMERSEARCH_ANALYSIS, 
                                          analysis_dsm_segment);
        if (analysis_dsa == NULL) {
            ereport(ERROR, (errmsg("Failed to create DSA area for analysis")));
        }
        
        /* Pin the DSA area */
        dsa_pin(analysis_dsa);
    }
    PG_CATCH();
    {
        if (analysis_dsm_segment != NULL) {
            dsm_unpin_mapping(analysis_dsm_segment);
            dsm_detach(analysis_dsm_segment);
            analysis_dsm_segment = NULL;
        }
        ereport(ERROR, (errmsg("Failed to create or pin DSA area for analysis")));
        return false;
    }
    PG_END_TRY();
    
    /* Step 3: Set up dshash parameters based on k-mer size */
    memset(&params, 0, sizeof(params));
    if (kmer_size <= 8) {
        params.key_size = sizeof(uint16);
        params.entry_size = sizeof(KmerEntry16);
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint16_identity_hash;
    } else if (kmer_size <= 16) {
        params.key_size = sizeof(uint32);
        params.entry_size = sizeof(KmerEntry32);
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint32_identity_hash;
    } else {
        params.key_size = sizeof(uint64);
        params.entry_size = sizeof(KmerEntry64);
        params.compare_function = dshash_memcmp;
        params.hash_function = dshash_memhash;
    }
    params.tranche_id = LWTRANCHE_KMERSEARCH_ANALYSIS;
    
    /* Step 4: Create dshash table */
    PG_TRY();
    {
        analysis_highfreq_hash = dshash_create(analysis_dsa, &params, NULL);
        if (analysis_highfreq_hash == NULL) {
            ereport(ERROR, (errmsg("Failed to create dshash table for analysis")));
        }
    }
    PG_CATCH();
    {
        if (analysis_dsa != NULL) {
            dsa_unpin(analysis_dsa);
            dsa_detach(analysis_dsa);
            analysis_dsa = NULL;
        }
        if (analysis_dsm_segment != NULL) {
            dsm_unpin_mapping(analysis_dsm_segment);
            dsm_detach(analysis_dsm_segment);
            analysis_dsm_segment = NULL;
        }
        ereport(ERROR, (errmsg("Failed to create dshash table for analysis")));
        return false;
    }
    PG_END_TRY();
    
    ereport(NOTICE, (errmsg("Analysis dshash table created successfully with estimated capacity for %d entries", estimated_entries)));
    return true;
}

/*
 * Clean up analysis DSM/DSA/dshash resources
 */
static void
kmersearch_cleanup_analysis_dshash(void)
{
    /* Never perform cleanup in parallel workers */
    if (IsParallelWorker()) {
        return;
    }
    
    /* Step 1: Clean up dshash table */
    if (analysis_highfreq_hash != NULL) {
        dshash_destroy(analysis_highfreq_hash);
        analysis_highfreq_hash = NULL;
    }
    
    /* Step 2: Clean up DSA area */
    if (analysis_dsa != NULL) {
        dsa_unpin(analysis_dsa);
        dsa_detach(analysis_dsa);
        analysis_dsa = NULL;
    }
    
    /* Step 3: Clean up DSM segment */
    if (analysis_dsm_segment != NULL) {
        dsm_detach(analysis_dsm_segment);
        analysis_dsm_segment = NULL;
    }
}

/*
 * Check if a k-mer hash exists in the analysis dshash table
 * This function is called from kmersearch.c, so it's not static
 */
bool
kmersearch_is_kmer_hash_in_analysis_dshash(uint64 kmer_hash)
{
    AnalysisHighfreqKmerEntry *entry;
    bool found = false;
    
    /* This function should not be called from parallel workers */
    if (IsParallelWorker() || analysis_highfreq_hash == NULL) {
        return false;
    }
    
    /* Look up the k-mer hash in the dshash table */
    entry = (AnalysisHighfreqKmerEntry *) dshash_find(analysis_highfreq_hash, 
                                                     &kmer_hash, 
                                                     false);
    
    if (entry != NULL) {
        found = true;
        dshash_release_lock(analysis_highfreq_hash, entry);
    }
    
    return found;
}

/*
 * Populate analysis dshash using hybrid SQL+dshash approach
 * First creates temp_kmer_final table with SQL aggregation, then loads into dshash
 */
static void
kmersearch_populate_analysis_dshash_from_workers(KmerWorkerState *workers, int num_workers, uint64 threshold_rows)
{
    StringInfoData query;
    char *temp_kmer_final_name;
    int worker_id;
    int ret;
    int processed_entries = 0;
    StringInfoData debug_query;
    bool first_table;
    
    if (analysis_highfreq_hash == NULL) {
        ereport(ERROR, (errmsg("Analysis dshash table not initialized")));
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_populate_analysis_dshash_from_workers: Starting hybrid SQL+dshash merge with threshold %lu", (unsigned long)threshold_rows)));
    
    /* Generate unique name for temp_kmer_final table */
    temp_kmer_final_name = kmersearch_generate_unique_temp_table_name("temp_kmer_final", -1);
    
    /* Step 1: Create temp_kmer_final table with SQL aggregation (like original implementation) */
    initStringInfo(&query);
    appendStringInfo(&query, 
        "CREATE TEMP TABLE %s ("
        "kmer_data bigint PRIMARY KEY, "
        "frequency_count integer"
        ")", temp_kmer_final_name);
    
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_UTILITY) {
        ereport(ERROR, (errmsg("Failed to create temp_kmer_final table")));
    }
    pfree(query.data);
    
    /* Debug: Check worker table contents before merging */
    for (worker_id = 0; worker_id < num_workers; worker_id++) {
        if (workers[worker_id].temp_table_name == NULL) {
            continue;
        }
        
        initStringInfo(&debug_query);
        appendStringInfo(&debug_query, "SELECT count(*) FROM %s", workers[worker_id].temp_table_name);
        
        ret = SPI_exec(debug_query.data, 0);
        if (ret == SPI_OK_SELECT && SPI_processed > 0) {
            bool isnull;
            int row_count = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull));
            ereport(DEBUG1, (errmsg("Worker %d table %s contains %d k-mer entries", 
                                   worker_id, workers[worker_id].temp_table_name, row_count)));
        }
        pfree(debug_query.data);
    }
    
    /* Step 2: Build UNION ALL query to combine all worker tables with SQL aggregation */
    initStringInfo(&query);
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) ", temp_kmer_final_name);
    appendStringInfoString(&query, "SELECT kmer_data, sum(frequency_count) FROM (");
    
    first_table = true;
    for (worker_id = 0; worker_id < num_workers; worker_id++) {
        if (workers[worker_id].temp_table_name == NULL) {
            continue;
        }
        
        if (!first_table) {
            appendStringInfoString(&query, " UNION ALL ");
        }
        appendStringInfo(&query, "SELECT kmer_data, frequency_count FROM %s", 
                        workers[worker_id].temp_table_name);
        first_table = false;
    }
    
    appendStringInfo(&query, ") AS combined GROUP BY kmer_data HAVING sum(frequency_count) > %lu", 
                    (unsigned long)threshold_rows);
    
    /* Execute SQL aggregation query */
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_INSERT) {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to execute k-mer aggregation query"),
                 errdetail("SPI_exec returned %d", ret),
                 errhint("Query was: %s", query.data)));
    }
    
    ereport(DEBUG1, (errmsg("SQL aggregation completed, %lu high-frequency k-mers found", SPI_processed)));
    pfree(query.data);
    
    /* Step 3: Load aggregated results from temp_kmer_final into dshash */
    initStringInfo(&query);
    appendStringInfo(&query, "SELECT kmer_data, frequency_count FROM %s", temp_kmer_final_name);
    
    ret = SPI_exec(query.data, 0);
    if (ret == SPI_OK_SELECT) {
        int i;
        
        ereport(DEBUG1, (errmsg("Loading %lu k-mers from temp_kmer_final into dshash", SPI_processed)));
        
        /* Insert each high-frequency k-mer into the dshash table */
        for (i = 0; i < SPI_processed; i++) {
            bool isnull_kmer, isnull_freq;
            uint64 kmer_hash = DatumGetUInt64(SPI_getbinval(SPI_tuptable->vals[i], 
                                                            SPI_tuptable->tupdesc, 1, &isnull_kmer));
            int frequency = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[i], 
                                                       SPI_tuptable->tupdesc, 2, &isnull_freq));
            
            if (!isnull_kmer && !isnull_freq) {
                AnalysisHighfreqKmerEntry *entry;
                bool found;
                
                /* Insert entry in dshash table (should be new since we're loading from aggregated table) */
                entry = (AnalysisHighfreqKmerEntry *) dshash_find_or_insert(analysis_highfreq_hash,
                                                                           &kmer_hash,
                                                                           &found);
                if (entry != NULL) {
                    entry->kmer_hash = kmer_hash;
                    entry->frequency_count = frequency;
                    processed_entries++;
                    dshash_release_lock(analysis_highfreq_hash, entry);
                    
                    if (found) {
                        ereport(WARNING, (errmsg("Duplicate k-mer hash found during dshash loading: %lu", kmer_hash)));
                    }
                }
            }
        }
    }
    pfree(query.data);
    
    /* Step 4: Clean up worker temp tables and temp_kmer_final */
    for (worker_id = 0; worker_id < num_workers; worker_id++) {
        if (workers[worker_id].temp_table_name == NULL) {
            continue;
        }
        
        initStringInfo(&query);
        appendStringInfo(&query, "DROP TABLE IF EXISTS %s", workers[worker_id].temp_table_name);
        SPI_exec(query.data, 0);
        pfree(query.data);
    }
    
    /* Clean up temp_kmer_final table */
    initStringInfo(&query);
    appendStringInfo(&query, "DROP TABLE IF EXISTS %s", temp_kmer_final_name);
    SPI_exec(query.data, 0);
    pfree(query.data);
    
    pfree(temp_kmer_final_name);
    
    ereport(NOTICE, (errmsg("Hybrid SQL+dshash approach completed: loaded %d high-frequency k-mer entries into dshash", processed_entries)));
}

/*
 * Get count of high-frequency k-mers (those exceeding threshold)
 */
static int
kmersearch_get_analysis_dshash_count(int k_size, uint64 threshold_rows)
{
    int count = 0;
    dshash_seq_status status;
    void *entry;
    
    if (analysis_highfreq_hash == NULL) {
        return 0;
    }
    
    /* Iterate through all entries and count those exceeding threshold */
    dshash_seq_init(&status, analysis_highfreq_hash, false);
    while ((entry = dshash_seq_next(&status)) != NULL) {
        int kmer_count;
        
        /* Extract count based on k-mer size */
        if (k_size <= 8) {
            KmerEntry16 *entry16 = (KmerEntry16 *)entry;
            kmer_count = entry16->count;
        } else if (k_size <= 16) {
            KmerEntry32 *entry32 = (KmerEntry32 *)entry;
            kmer_count = entry32->count;
        } else {
            KmerEntry64 *entry64 = (KmerEntry64 *)entry;
            kmer_count = entry64->count;
        }
        
        /* Count only if exceeds threshold */
        if (kmer_count > threshold_rows) {
            count++;
        }
    }
    dshash_seq_term(&status);
    
    return count;
}

/*
 * Insert kmer2_as_uint values directly from dshash into kmersearch_highfreq_kmer table
 */
static void
kmersearch_insert_kmer2_as_uint_from_dshash(Oid table_oid, const char *column_name, int k_size, uint64 threshold_rows)
{
    dshash_seq_status status;
    void *entry;
    StringInfoData query;
    int processed_entries = 0;
    bool first_entry = true;
    int total_entries = 0;
    int high_freq_entries = 0;
    
    /* Ensure this function is never called from parallel workers */
    if (IsParallelWorker())
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
                 errmsg("cannot execute kmersearch_insert_kmer2_as_uint_from_dshash() in a parallel worker")));
    
    if (analysis_highfreq_hash == NULL) {
        ereport(ERROR, (errmsg("Analysis dshash table not initialized")));
    }
    
    /* Double-check we're not in a parallel worker */
    if (IsParallelWorker()) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
                 errmsg("analysis_highfreq_hash access detected in parallel worker"),
                 errdetail("This should never happen after initial check")));
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_insert_kmer2_as_uint_from_dshash: Starting direct insertion")));
    
    /* Threshold already calculated in the caller - just use it */
    elog(DEBUG1, "Using threshold_rows = %lu for filtering high-frequency k-mers", (unsigned long)threshold_rows);
    
    /* Build base insert query */
    initStringInfo(&query);
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, kmer2_as_uint, detection_reason) VALUES ");
    
    /* Iterate through all entries in the dshash table */
    dshash_seq_init(&status, analysis_highfreq_hash, false);
    
    while ((entry = dshash_seq_next(&status)) != NULL) {
        uint64 kmer2_as_uint;
        int count;
        
        total_entries++;
        
        /* Extract kmer and count based on k-mer size */
        if (k_size <= 8) {
            KmerEntry16 *entry16 = (KmerEntry16 *)entry;
            kmer2_as_uint = (uint64)entry16->kmer;
            count = entry16->count;
        } else if (k_size <= 16) {
            KmerEntry32 *entry32 = (KmerEntry32 *)entry;
            kmer2_as_uint = (uint64)entry32->kmer;
            count = entry32->count;
        } else {
            KmerEntry64 *entry64 = (KmerEntry64 *)entry;
            kmer2_as_uint = entry64->kmer;
            count = entry64->count;
        }
        
        /* Debug: log k-mer counts */
        if (count > 0) {
            elog(DEBUG1, "K-mer %lu has count=%d, threshold=%lu", (unsigned long)kmer2_as_uint, count, (unsigned long)threshold_rows);
        }
        
        /* Only insert if count exceeds threshold */
        if (count > threshold_rows) {
            high_freq_entries++;
            if (!first_entry) {
                appendStringInfoString(&query, ", ");
            }
            
            /* Add entry to batch insert */
            appendStringInfo(&query, "(%u, %s, %ld, 'frequency_analysis')",
                table_oid, 
                quote_literal_cstr(column_name),
                (int64)kmer2_as_uint);
            
            processed_entries++;
            first_entry = false;
            
            /* Execute batch insert every 1000 entries to avoid query length limits */
            if (processed_entries % 1000 == 0) {
                int ret;
                
                /* Safety check before SQL execution */
                if (IsParallelWorker()) {
                    ereport(ERROR,
                            (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
                             errmsg("About to execute SQL in parallel worker - this should never happen"),
                             errdetail("Query: %s", query.data)));
                }
                
                ret = SPI_exec(query.data, 0);
                if (ret != SPI_OK_INSERT) {
                    ereport(ERROR, (errmsg("Failed to insert batch of kmer2_as_uint values")));
                }
                
                /* Reset query for next batch */
                pfree(query.data);
                initStringInfo(&query);
                appendStringInfo(&query,
                    "INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, kmer2_as_uint, detection_reason) VALUES ");
                first_entry = true;
            }
        }
    }
    
    /* Execute final batch if there are remaining entries */
    if (!first_entry) {
        int ret;
        
        /* Final safety check before SQL execution */
        if (IsParallelWorker()) {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
                     errmsg("About to execute SQL in parallel worker - this should never happen"),
                     errdetail("Query: %s", query.data)));
        }
        
        ret = SPI_exec(query.data, 0);
        if (ret != SPI_OK_INSERT) {
            ereport(ERROR, (errmsg("Failed to insert final batch of kmer2_as_uint values")));
        }
    }
    
    dshash_seq_term(&status);
    pfree(query.data);
    
    ereport(NOTICE, (errmsg("K-mer analysis summary: total_entries=%d, threshold=%lu, high_freq_entries=%d, processed_entries=%d", 
                           total_entries, (unsigned long)threshold_rows, high_freq_entries, processed_entries)));
    ereport(NOTICE, (errmsg("Direct storage completed: inserted %d kmer2_as_uint values into kmersearch_highfreq_kmer table", processed_entries)));
}

/*
 * Calculate optimal buffer size based on memory constraints
 */
static int
kmersearch_calculate_buffer_size(int k_size)
{
    const int TARGET_MEMORY_MB = 50;  /* Target 50MB per worker */
    const int MIN_BUFFER_SIZE = 1000;
    const int MAX_BUFFER_SIZE = 100000;
    
    size_t entry_size = sizeof(CompactKmerFreq);
    int max_entries = (TARGET_MEMORY_MB * 1024 * 1024) / entry_size;
    
    if (max_entries < MIN_BUFFER_SIZE) return MIN_BUFFER_SIZE;
    if (max_entries > MAX_BUFFER_SIZE) return MAX_BUFFER_SIZE;
    return max_entries;
}

/*
 * Initialize k-mer buffer
 */
static void
kmersearch_init_buffer(KmerBuffer *buffer, int k_size)
{
    buffer->capacity = kmersearch_calculate_buffer_size(k_size);
    buffer->entries = (CompactKmerFreq *) palloc0(buffer->capacity * sizeof(CompactKmerFreq));
    buffer->count = 0;
    buffer->kmer_size = k_size;
}

/*
 * Flush hash buffer to temporary table (simplified for Phase 1)
 */
static void
kmersearch_flush_hash_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    CompactKmerFreq *entry;
    int ret;
    
    if (buffer->count == 0) return;
    
    initStringInfo(&query);
    
    /* Build INSERT statement with multiple VALUES */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        entry = &buffer->entries[i];
        
        if (i > 0) appendStringInfoString(&query, ", ");
        
        /* Use k32_data field to access the stored hash value */
        appendStringInfo(&query, "(%lu, %d)",
                        (unsigned long)entry->kmer_data.k32_data,  /* Hash stored as k32_data */
                        entry->frequency_count);
    }
    
    appendStringInfoString(&query, " ON CONFLICT (kmer_data) DO UPDATE SET frequency_count = ");
    appendStringInfo(&query, "%s.frequency_count + EXCLUDED.frequency_count", temp_table_name);
    
    /* Execute the INSERT */
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_INSERT) {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("failed to insert k-mer data into temporary table %s", temp_table_name)));
    }
    
    pfree(query.data);
    
    /* Reset buffer */
    buffer->count = 0;
}

/*
 * Add hash value to buffer for Phase 1 processing
 */
static void
kmersearch_add_hash_to_buffer(KmerBuffer *buffer, uint64_t kmer_hash, const char *temp_table_name)
{
    CompactKmerFreq *entry;
    int i;
    
    /* No buffer-level deduplication - each row contribution should be counted */
    /* Row-level deduplication is already handled in the calling function */
    
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_hash_buffer_to_table(buffer, temp_table_name);
        
        /* After flush, continue to add the new entry */
    }
    
    /* Add new entry using k32_data field to store uint64_t hash */
    entry = &buffer->entries[buffer->count];
    entry->kmer_data.k32_data = kmer_hash;  /* Store hash in k32_data field */
    entry->frequency_count = 1;
    entry->is_highfreq = false;
    buffer->count++;
}


/*
 * Worker function to analyze blocks and extract k-mer frequencies
 */
void
kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, 
                                const char *column_name, int k_size, int target_attno, bool is_dna4_type)
{
    TableScanDesc scan;
    HeapTuple tuple;
    TupleDesc tupdesc;
    int i;
    Datum *kmer_datums;
    BlockNumber current_block;
    bool isnull;
    Datum value;
    VarBit *sequence;
    VarBit **kmers;
    int nkeys;
    int j;
    
    /* Use passed parameters instead of determining them again */
    tupdesc = RelationGetDescr(rel);
    
    /* Initialize buffer */
    kmersearch_init_buffer(&worker->buffer, k_size);
    
    /* Create temporary table for this worker with unique name */
    worker->temp_table_name = kmersearch_generate_unique_temp_table_name("temp_kmer_worker", worker->worker_id);
    
    /* Scan assigned blocks only */
    
    for (current_block = worker->start_block; current_block < worker->end_block; current_block++) {
        Buffer buffer;
        Page page;
        OffsetNumber max_offset;
        OffsetNumber offset;
        
        /* Read the block */
        buffer = ReadBufferExtended(rel, MAIN_FORKNUM, current_block, RBM_NORMAL, NULL);
        LockBuffer(buffer, BUFFER_LOCK_SHARE);
        page = BufferGetPage(buffer);
        max_offset = PageGetMaxOffsetNumber(page);
        
        /* Process all tuples in this block */
        for (offset = FirstOffsetNumber; offset <= max_offset; offset = OffsetNumberNext(offset)) {
            ItemId item_id;
            HeapTupleData tuple_data;
            HeapTuple block_tuple = &tuple_data;
            
            /* Get the tuple */
            item_id = PageGetItemId(page, offset);
            if (!ItemIdIsNormal(item_id))
                continue;
                
            block_tuple->t_len = ItemIdGetLength(item_id);
            block_tuple->t_data = (HeapTupleHeader) PageGetItem(page, item_id);
            block_tuple->t_tableOid = RelationGetRelid(rel);
            block_tuple->t_self.ip_blkid.bi_hi = current_block >> 16;
            block_tuple->t_self.ip_blkid.bi_lo = current_block & 0xFFFF;
            block_tuple->t_self.ip_posid = offset;
            
            worker->rows_processed++;
            
            /* Extract the DNA sequence value */
            value = heap_getattr(block_tuple, target_attno, tupdesc, &isnull);
        if (isnull) {
            continue;  /* Skip NULL values */
        }
        
        /* Convert DNA data to VarBit representation */
        sequence = DatumGetVarBitP(value);
        
        /* Extract k-mers from the sequence using SIMD-optimized function based on DNA type */
        if (is_dna4_type) {
            kmer_datums = kmersearch_extract_dna4_kmer2_with_expansion_direct(sequence, k_size, &nkeys);
        } else {
            kmer_datums = kmersearch_extract_dna2_kmer2_direct(sequence, k_size, &nkeys);
        }
        
        if (kmer_datums == NULL || nkeys == 0) {
            continue;
        }
        
        /* Use hash set to track unique k-mers in this row to avoid counting duplicates */
        {
            HTAB *row_kmer_set = NULL;
            HASHCTL hash_ctl;
            
            /* Create hash table for unique k-mers in this row */
            memset(&hash_ctl, 0, sizeof(hash_ctl));
            hash_ctl.keysize = sizeof(uint64_t);
            hash_ctl.entrysize = sizeof(uint64_t);
            hash_ctl.hcxt = CurrentMemoryContext;
            
            row_kmer_set = hash_create("row_kmer_set", nkeys, &hash_ctl, 
                                       HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
            
            /* Process each k-mer in this row - convert to hash and deduplicate */
            for (j = 0; j < nkeys; j++) {
                uint64_t kmer_hash;
                VarBit *kmer;
                bool found;
                
                /* Convert Datum to VarBit */
                kmer = DatumGetVarBitP(kmer_datums[j]);
                if (kmer == NULL) {
                    continue;
                }
                
                /* Convert k-mer to consistent hash value */
                kmer_hash = kmersearch_get_kmer_hash(kmer, 0, k_size);
                
                /* Only add to buffer if not already seen in this row */
                hash_search(row_kmer_set, (void *) &kmer_hash, HASH_ENTER, &found);
                if (!found) {
                    /* Add to buffer (will flush to temp table if full) */
                    kmersearch_add_hash_to_buffer(&worker->buffer, kmer_hash, worker->temp_table_name);
                }
            }
            
            /* Clean up row hash table */
            hash_destroy(row_kmer_set);
        }
        
        /* Cleanup k-mer array */
        if (kmer_datums) {
            pfree(kmer_datums);
        }
        }
        
        /* Release buffer and lock */
        UnlockReleaseBuffer(buffer);
    }
    
    /* Flush any remaining buffer contents */
    kmersearch_flush_hash_buffer_to_table(&worker->buffer, worker->temp_table_name);
    
    /* Cleanup buffer */
    if (worker->buffer.entries) {
        pfree(worker->buffer.entries);
    }
}

/*
 * Parallel worker function for k-mer analysis
 */
PGDLLEXPORT void
kmersearch_analysis_worker(dsm_segment *seg, shm_toc *toc)
{
    KmerAnalysisSharedState *shared_state = NULL;
    KmerAnalysisHandles *handles = NULL;
    dsa_area *dsa = NULL;
    dshash_table *hash = NULL;
    dshash_parameters params;
    dsm_segment *analysis_seg = NULL;
    BlockNumber current_block;
    bool has_work = true;
    
    /* Ensure global variables are NULL in worker processes */
    analysis_dsm_segment = NULL;
    analysis_dsa = NULL;
    analysis_highfreq_hash = NULL;
    
    /* Double-check we're in a parallel worker */
    if (!IsParallelWorker())
    {
        elog(ERROR, "kmersearch_analysis_worker called from non-parallel context");
    }
    
    elog(DEBUG1, "kmersearch_analysis_worker: Worker started - PID %d, IsParallelWorker=%d", 
         MyProcPid, IsParallelWorker());
    
    /* Get shared state from shm_toc */
    shared_state = (KmerAnalysisSharedState *)shm_toc_lookup(toc, KMERSEARCH_KEY_SHARED_STATE, false);
    if (!shared_state) {
        elog(ERROR, "kmersearch_analysis_worker: Failed to get shared state");
    }
    elog(DEBUG1, "kmersearch_analysis_worker: Got shared state, table_oid=%u, kmer_size=%d", 
         shared_state->table_oid, shared_state->kmer_size);
    
    /* Get handles from shm_toc */
    handles = (KmerAnalysisHandles *)shm_toc_lookup(toc, KMERSEARCH_KEY_HANDLES, false);
    if (!handles) {
        elog(ERROR, "kmersearch_analysis_worker: Failed to get handles");
    }
    
    /* Attach to DSM segment containing the DSA */
    analysis_seg = dsm_attach(handles->dsm_handle);
    if (!analysis_seg) {
        elog(ERROR, "kmersearch_analysis_worker: Failed to attach to DSM segment");
    }
    
    /* Get DSA from the DSM segment (DSA was created in place) */
    dsa = dsa_attach_in_place(dsm_segment_address(analysis_seg), analysis_seg);
    if (!dsa) {
        elog(ERROR, "kmersearch_analysis_worker: Failed to attach to DSA in place");
    }
    
    /* Set up dshash parameters based on k-mer size */
    memset(&params, 0, sizeof(params));
    if (shared_state->kmer_size <= 8) {
        params.key_size = sizeof(uint16);
        params.entry_size = sizeof(KmerEntry16);
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint16_identity_hash;
    } else if (shared_state->kmer_size <= 16) {
        params.key_size = sizeof(uint32);
        params.entry_size = sizeof(KmerEntry32);
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint32_identity_hash;
    } else {
        params.key_size = sizeof(uint64);
        params.entry_size = sizeof(KmerEntry64);
        params.compare_function = dshash_memcmp;
        params.hash_function = dshash_memhash;
    }
    params.tranche_id = LWTRANCHE_KMERSEARCH_ANALYSIS;
    
    /* Attach to dshash table */
    hash = dshash_attach(dsa, &params, handles->hash_handle, NULL);
    if (!hash) {
        elog(ERROR, "kmersearch_analysis_worker: Failed to attach to dshash");
    }
    elog(DEBUG1, "kmersearch_analysis_worker: Successfully attached to dshash");
    
    /* Dynamic work acquisition loop */
    while (has_work)
    {
        /* Acquire block number with short exclusive lock */
        LWLockAcquire(&shared_state->mutex.lock, LW_EXCLUSIVE);
        
        if (shared_state->next_block < shared_state->total_blocks)
        {
            /* Get block to process */
            current_block = shared_state->next_block;
            shared_state->next_block++;
            LWLockRelease(&shared_state->mutex.lock);
            
            /* Process block (parallel execution outside lock) */
            kmersearch_extract_kmers_from_block(shared_state->table_oid,
                                               shared_state->column_attnum,
                                               shared_state->column_type_oid,
                                               current_block,
                                               shared_state->kmer_size,
                                               hash);
            elog(DEBUG2, "kmersearch_analysis_worker: Finished processing block %u", current_block);
        }
        else
        {
            /* All blocks processed */
            shared_state->all_processed = true;
            LWLockRelease(&shared_state->mutex.lock);
            has_work = false;
        }
    }
    
    /* Normal cleanup */
    if (hash)
        dshash_detach(hash);
    
    /* Detach from DSA and DSM - important to prevent leaks */
    if (dsa)
        dsa_detach(dsa);
    
    if (analysis_seg)
        dsm_detach(analysis_seg);
}

/*
 * Extract k-mers from a specific block
 */
static void
kmersearch_extract_kmers_from_block(Oid table_oid, AttrNumber column_attnum, Oid column_type_oid,
                                   BlockNumber block, int kmer_size, dshash_table *hash)
{
    Relation rel = NULL;
    Buffer buffer;
    Page page;
    OffsetNumber offset;
    int ntuples;
    bool isnull;
    Datum datum;
    
    elog(DEBUG2, "kmersearch_extract_kmers_from_block: Starting block %u", block);
    
    /* Open table */
    rel = table_open(table_oid, AccessShareLock);
    
    /* Read block */
    buffer = ReadBuffer(rel, block);
    LockBuffer(buffer, BUFFER_LOCK_SHARE);
    page = BufferGetPage(buffer);
    ntuples = PageGetMaxOffsetNumber(page);
    
    /* Process each tuple in block */
    for (offset = FirstOffsetNumber; offset <= ntuples; offset++)
    {
        ItemId itemid = PageGetItemId(page, offset);
        HeapTupleData tuple;
        
        if (!ItemIdIsNormal(itemid))
            continue;
            
        tuple.t_data = (HeapTupleHeader) PageGetItem(page, itemid);
        tuple.t_len = ItemIdGetLength(itemid);
        ItemPointerSet(&(tuple.t_self), block, offset);
        
        /* Get column value */
        datum = heap_getattr(&tuple, column_attnum, RelationGetDescr(rel), &isnull);
        if (isnull)
            continue;
            
        /* Extract k-mers and update dshash */
        kmersearch_update_kmer_counts_in_dshash(datum, kmer_size, hash, column_type_oid);
    }
    
    UnlockReleaseBuffer(buffer);
    table_close(rel, AccessShareLock);
}

/*
 * Update k-mer counts in dshash
 * Each unique k-mer within a row is counted only once
 */
static void
kmersearch_update_kmer_counts_in_dshash(Datum sequence_datum, int kmer_size, dshash_table *hash, Oid column_type_oid)
{
    void *kmer_array = NULL;
    int kmer_count;
    bool found;
    Oid dna2_oid;
    Oid dna4_oid;
    int i;
    HTAB *local_kmer_set = NULL;
    HASHCTL hashctl;
    
    /* Extract k-mers based on data type */
    dna2_oid = TypenameGetTypid("dna2");
    dna4_oid = TypenameGetTypid("dna4");
    if (column_type_oid == dna2_oid)
    {
        VarBit *seq;
        
        /* Get properly detoasted sequence */
        seq = DatumGetVarBitP(sequence_datum);
        
        /* Direct extraction from DNA2 */
        kmersearch_extract_dna2_kmer2_as_uint_direct(seq, kmer_size, &kmer_array, &kmer_count);
        
        /* Free detoasted copy if needed */
        if ((void *)seq != DatumGetPointer(sequence_datum))
            pfree(seq);
    }
    else if (column_type_oid == dna4_oid)
    {
        VarBit *seq;
        
        /* Get properly detoasted sequence */
        seq = DatumGetVarBitP(sequence_datum);
        
        /* Direct extraction with expansion from DNA4 */
        kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(seq, kmer_size, &kmer_array, &kmer_count);
        
        /* Free detoasted copy if needed */  
        if ((void *)seq != DatumGetPointer(sequence_datum))
            pfree(seq);
    }
    else
    {
        elog(ERROR, "Column must be DNA2 or DNA4 type");
    }
    
    if (kmer_array == NULL || kmer_count <= 0)
    {
        if (kmer_array)
            pfree(kmer_array);
        return; /* No valid k-mers */
    }
    
    /* Create local hash table to track unique k-mers within this row */
    memset(&hashctl, 0, sizeof(hashctl));
    if (kmer_size <= 8) {
        hashctl.keysize = sizeof(uint16);
    } else if (kmer_size <= 16) {
        hashctl.keysize = sizeof(uint32);
    } else {
        hashctl.keysize = sizeof(uint64);
    }
    hashctl.entrysize = hashctl.keysize;  /* We only need the key, no value */
    hashctl.hcxt = CurrentMemoryContext;
    local_kmer_set = hash_create("row k-mer set", 
                                 (kmer_count < 1000) ? kmer_count : 1000,  /* Initial size - most rows won't have too many unique k-mers */
                                 &hashctl,
                                 HASH_ELEM | HASH_CONTEXT);
    
    /* Update counts for each unique k-mer in this row */
    elog(DEBUG2, "Processing row with %d k-mers", kmer_count);
    for (i = 0; i < kmer_count; i++)
    {
        bool local_found;
        
        /* Process based on k-mer size */
        if (kmer_size <= 8)
        {
            uint16 *kmer16_array = (uint16 *)kmer_array;
            uint16 kmer16 = kmer16_array[i];
            
            /* Check if we've already seen this k-mer in this row */
            hash_search(local_kmer_set, &kmer16, HASH_ENTER, &local_found);
            if (!local_found)  /* First time seeing this k-mer in this row */
            {
                KmerEntry16 *entry16;
                
                elog(DEBUG2, "First occurrence of k-mer %u in this row", kmer16);
                entry16 = (KmerEntry16 *)dshash_find_or_insert(hash, &kmer16, &found);
                if (!found)
                {
                    entry16->kmer = kmer16;
                    entry16->count = 1;
                    elog(DEBUG2, "New k-mer %u, count=1", kmer16);
                }
                else
                {
                    entry16->count++;
                    elog(DEBUG2, "Existing k-mer %u, count=%d", kmer16, entry16->count);
                }
                dshash_release_lock(hash, entry16);
            }
            else
            {
                elog(DEBUG2, "Duplicate k-mer %u in this row, skipping", kmer16);
            }
        }
        else if (kmer_size <= 16)
        {
            uint32 *kmer32_array = (uint32 *)kmer_array;
            uint32 kmer32 = kmer32_array[i];
            
            /* Check if we've already seen this k-mer in this row */
            hash_search(local_kmer_set, &kmer32, HASH_ENTER, &local_found);
            if (!local_found)  /* First time seeing this k-mer in this row */
            {
                KmerEntry32 *entry32;
                
                entry32 = (KmerEntry32 *)dshash_find_or_insert(hash, &kmer32, &found);
                if (!found)
                {
                    entry32->kmer = kmer32;
                    entry32->count = 1;
                }
                else
                {
                    if (entry32->count < INT_MAX)
                        entry32->count++;
                }
                dshash_release_lock(hash, entry32);
            }
        }
        else /* kmer_size <= 32 */
        {
            uint64 *kmer64_array = (uint64 *)kmer_array;
            uint64 kmer64 = kmer64_array[i];
            
            /* Check if we've already seen this k-mer in this row */
            hash_search(local_kmer_set, &kmer64, HASH_ENTER, &local_found);
            if (!local_found)  /* First time seeing this k-mer in this row */
            {
                KmerEntry64 *entry64;
                
                entry64 = (KmerEntry64 *)dshash_find_or_insert(hash, &kmer64, &found);
                if (!found)
                {
                    entry64->kmer = kmer64;
                    entry64->count = 1;
                }
                else
                {
                    if (entry64->count < INT_MAX)
                        entry64->count++;
                }
                dshash_release_lock(hash, entry64);
            }
        }
    }
    
    /* Normal cleanup */
    if (local_kmer_set)
        hash_destroy(local_kmer_set);
    if (kmer_array)
        pfree(kmer_array);
}

/*
 * Create dshash resources for analysis (must be called before EnterParallelMode)
 */
static void
create_analysis_dshash_resources(KmerAnalysisContext *ctx, int estimated_entries, int kmer_size)
{
    dshash_parameters params;
    Size segment_size;
    
    /* Calculate DSM segment size */
    segment_size = 1024 * 1024;  /* Start with 1MB */
    if (estimated_entries > 10000) {
        segment_size = estimated_entries * 64;  /* Estimate 64 bytes per entry */
    }
    
    /* Limit maximum size to avoid out of shared memory */
    if (segment_size > 4 * 1024 * 1024) {
        segment_size = 4 * 1024 * 1024;  /* Maximum 4MB for 128MB shared_buffers */
    }
    /* Ensure minimum size for DSA overhead */
    if (segment_size < 1 * 1024 * 1024) {
        segment_size = 1 * 1024 * 1024;  /* Minimum 1MB */
    }
    
    elog(DEBUG1, "create_analysis_dshash_resources: Creating DSM segment with size %zu for %d estimated entries", 
         segment_size, estimated_entries);
    
    /* Create DSM segment */
    ctx->dsm_seg = dsm_create(segment_size, 0);
    if (ctx->dsm_seg == NULL) {
        ereport(ERROR, (errmsg("Failed to create DSM segment for analysis")));
    }
    /* Pin both the segment and mapping to prevent cleanup issues */
    dsm_pin_segment(ctx->dsm_seg);
    dsm_pin_mapping(ctx->dsm_seg);
    ctx->dsm_handle = dsm_segment_handle(ctx->dsm_seg);
    
    /* Create DSA area */
    ctx->dsa = dsa_create_in_place(dsm_segment_address(ctx->dsm_seg), 
                                   segment_size, 
                                   LWTRANCHE_KMERSEARCH_ANALYSIS, 
                                   ctx->dsm_seg);
    if (ctx->dsa == NULL) {
        dsm_detach(ctx->dsm_seg);
        ereport(ERROR, (errmsg("Failed to create DSA area for analysis")));
    }
    dsa_pin(ctx->dsa);
    
    /* Set up dshash parameters based on k-mer size */
    memset(&params, 0, sizeof(params));
    if (kmer_size <= 8) {
        params.key_size = sizeof(uint16);
        params.entry_size = sizeof(KmerEntry16);
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint16_identity_hash;
    } else if (kmer_size <= 16) {
        params.key_size = sizeof(uint32);
        params.entry_size = sizeof(KmerEntry32);
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint32_identity_hash;
    } else {
        params.key_size = sizeof(uint64);
        params.entry_size = sizeof(KmerEntry64);
        params.compare_function = dshash_memcmp;
        params.hash_function = dshash_memhash;
    }
    params.tranche_id = LWTRANCHE_KMERSEARCH_ANALYSIS;
    
    /* Create dshash table */
    ctx->hash = dshash_create(ctx->dsa, &params, NULL);
    if (ctx->hash == NULL) {
        dsa_unpin(ctx->dsa);
        dsa_detach(ctx->dsa);
        dsm_unpin_mapping(ctx->dsm_seg);
        dsm_detach(ctx->dsm_seg);
        ereport(ERROR, (errmsg("Failed to create dshash table for analysis")));
    }
    ctx->hash_handle = dshash_get_hash_table_handle(ctx->hash);
    
    /* IMPORTANT: Do NOT set global pointers here - they should only be set after
     * EnterParallelMode() to prevent parallel workers from inheriting them */
    
    elog(DEBUG1, "create_analysis_dshash_resources: Resources created successfully");
}

/*
 * Merge worker results using SQL aggregation
 */
void
kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, 
                                   const char *final_table_name, int k_size, int threshold_rows)
{
    StringInfoData query;
    StringInfoData union_query;
    const char *data_type;
    int i;
    int ret;
    
    initStringInfo(&query);
    initStringInfo(&union_query);
    
    /* Create final aggregation table */
    if (k_size > 32) {
        /* k > 32 not supported */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k_size)));
    }
    
    appendStringInfo(&query, 
        "CREATE TEMP TABLE %s ("
        "kmer_data bigint PRIMARY KEY, "
        "frequency_count integer"
        ")", final_table_name);
    
    /* Use existing SPI connection from main function - don't call SPI_connect() again */
    SPI_exec(query.data, 0);
    
    /* Debug: Check worker table contents before merging */
    for (i = 0; i < num_workers; i++) {
        StringInfoData debug_query;
        initStringInfo(&debug_query);
        appendStringInfo(&debug_query, "SELECT count(*) FROM %s", workers[i].temp_table_name);
        
        ret = SPI_exec(debug_query.data, 0);
        if (ret == SPI_OK_SELECT && SPI_processed > 0) {
            bool isnull;
            int row_count = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull));
        }
        pfree(debug_query.data);
    }

    /* Build UNION ALL query to combine all worker tables */
    resetStringInfo(&query);
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) ", final_table_name);
    appendStringInfoString(&query, "SELECT kmer_data, sum(frequency_count) FROM (");
    
    for (i = 0; i < num_workers; i++) {
        if (i > 0) appendStringInfoString(&query, " UNION ALL ");
        appendStringInfo(&query, "SELECT kmer_data, frequency_count FROM %s", 
                        workers[i].temp_table_name);
    }
    
    appendStringInfo(&query, ") AS combined GROUP BY kmer_data HAVING sum(frequency_count) > %d", 
                    threshold_rows);
    
    /* Execute aggregation query with detailed error handling */
    ret = SPI_exec(query.data, 0);
    
    if (ret != SPI_OK_INSERT) {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to execute k-mer aggregation query"),
                 errdetail("SPI_exec returned %d", ret),
                 errhint("Query was: %s", query.data)));
    }
    
    /* Don't call SPI_finish() - leave connection open for main function */
    
    pfree(query.data);
    pfree(union_query.data);
}

/*
 * Check if kmer2 value exists in high-frequency temporary table
 */
static bool
is_kmer2_in_highfreq_table(uint64_t kmer2_value, const char *highfreq_table_name)
{
    StringInfoData query;
    bool result = false;
    int ret;
    
    initStringInfo(&query);
    
    /* Check if kmer2_value exists in temp_kmer_final table */
    if (strstr(highfreq_table_name, "temp_kmer_final_") != NULL) {
        /* Final aggregation table uses kmer_data column (bigint) */
        appendStringInfo(&query,
            "SELECT 1 FROM %s WHERE kmer_data = %lu LIMIT 1",
            highfreq_table_name, kmer2_value);
    } else {
        /* 
         * Permanent table uses different structure - not supported yet
         * NOTE: This path is never reached in current implementation since
         * highfreq_table_name is always generated with "temp_kmer_final_" prefix
         * via kmersearch_generate_unique_temp_table_name()
         */
        ereport(ERROR, (errmsg("Permanent table lookup not supported in is_kmer2_in_highfreq_table")));
    }
    
    SPI_connect();
    ret = SPI_exec(query.data, 0);
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
        result = true;
        
    SPI_finish();
    pfree(query.data);
    
    return result;
}

/*
 * Check if kmer2 value corresponds to a high-frequency k-mer using analysis dshash
 * This function uses kmersearch_get_kmer_hash for consistent hash calculation
 */
static bool
is_kmer2_in_analysis_dshash(uint64_t kmer2_value)
{
    /* The kmer2_value is already a hash generated by kmersearch_get_kmer_hash(),
     * so we can directly use it for dshash lookup */
    return kmersearch_is_kmer_hash_in_analysis_dshash(kmer2_value);
}


