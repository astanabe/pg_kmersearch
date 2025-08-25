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
#include "partitioning/partdesc.h"
#include <sqlite3.h>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include "storage/fd.h"

/* PostgreSQL function info declarations for frequency functions */
PG_FUNCTION_INFO_V1(kmersearch_perform_highfreq_analysis);
PG_FUNCTION_INFO_V1(kmersearch_undo_highfreq_analysis);
PG_FUNCTION_INFO_V1(kmersearch_delete_tempfiles);

/* Temporary k-mer frequency entry structures for batch processing */
typedef struct TempKmerFreqEntry16
{
    uint16      uintkey;
    uint64      nrow;      /* Number of rows containing this k-mer */
} TempKmerFreqEntry16;

typedef struct TempKmerFreqEntry32
{
    uint32      uintkey;
    uint64      nrow;      /* Number of rows containing this k-mer */
} TempKmerFreqEntry32;

typedef struct TempKmerFreqEntry64
{
    uint64      uintkey;
    uint64      nrow;      /* Number of rows containing this k-mer */
} TempKmerFreqEntry64;

/* SQLite3-based worker context */
typedef struct SQLiteWorkerContext
{
    sqlite3     *db;                    /* SQLite3 database handle */
    sqlite3_stmt *insert_stmt;          /* Prepared INSERT statement */
    char        db_path[MAXPGPATH];     /* Database file path */
    HTAB        *batch_hash;            /* Batch hash table for aggregation */
    int         batch_count;            /* Current batch count */
    int         total_bits;             /* Total bits for k-mer + occurrence */
    Oid         dna2_oid;               /* Cached OID for dna2 type */
    Oid         dna4_oid;               /* Cached OID for dna4 type */
    Oid         column_type_oid;        /* Column data type OID */
} SQLiteWorkerContext;

static int analysis_kmer_hash_compare(const void *a, const void *b, size_t size, void *arg);
static uint32 analysis_kmer_hash_hash(const void *key, size_t size, void *arg);
static uint32 kmersearch_uint16_identity_hash(const void *key, size_t keysize, void *arg);
static uint32 kmersearch_uint32_identity_hash(const void *key, size_t keysize, void *arg);

/* SQLite3 helper functions */
static void kmersearch_register_worker_temp_file(KmerAnalysisSharedState *shared_state, 
                                                const char *file_path, int worker_id);
static void kmersearch_flush_batch_to_sqlite(HTAB *batch_hash, sqlite3 *db, 
                                            sqlite3_stmt *stmt, int total_bits);
static void kmersearch_process_block_with_batch(BlockNumber block,
                                               Oid table_oid,
                                               KmerAnalysisSharedState *shared_state,
                                               SQLiteWorkerContext *ctx);

/*
 * Create worker temporary table for k-mer uintkeys
 */
void
kmersearch_create_worker_kmer_temp_table(const char *table_name)
{
    StringInfoData query;
    int ret;
    int total_bits;
    const char *kmer_type;
    const char *freq_type;
    
    initStringInfo(&query);
    
    /* Calculate total bits needed */
    total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;
    
    /* Determine appropriate data types based on total bits */
    if (total_bits <= 16) {
        kmer_type = "smallint";  /* int2 */
    } else if (total_bits <= 32) {
        kmer_type = "integer";   /* int4 */
    } else {
        kmer_type = "bigint";    /* int8 */
    }
    
    /* Frequency count always uses integer (int4) as it represents occurrence count */
    freq_type = "integer";
    
    /* First try to drop table if it exists to avoid conflicts */
    appendStringInfo(&query, "DROP TABLE IF EXISTS %s", table_name);
    
    ret = SPI_exec(query.data, 0);
    if (ret < 0)
        ereport(ERROR, (errmsg("Failed to drop existing temp table %s", table_name)));
    
    /* Create the temp table with optimized data types */
    resetStringInfo(&query);
    appendStringInfo(&query,
        "CREATE TEMP TABLE %s ("
        "kmer_data %s, "
        "frequency_count %s"
        ")", table_name, kmer_type, freq_type);
    
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
    /* Note: Combining kmer_data and frequency_count into single uintkey */
    
    if (k_size > 32) {
        /* k > 32 not supported */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k_size)));
    }
    
    /* Combine kmer_data and frequency_count into single bigint uintkey */
    /* Cast to bigint to ensure compatibility with target table */
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, uintkey, detection_reason) "
        "SELECT %u, '%s', "
        "  (kmer_data::bigint << %d) | frequency_count::bigint AS uintkey, "
        "  'high_frequency' "
        "FROM %s "
        "WHERE kmer_data IS NOT NULL AND frequency_count > 0",
        table_oid, column_name, 
        kmersearch_occur_bitlen,
        temp_table_name);
    
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
            
            if (source_count == 0) {
                ereport(WARNING, (errmsg("No high-frequency k-mers found in source table %s", temp_table_name)));
            }
        }
        pfree(count_query.data);
    }
    
    /* Execute high-frequency k-mer insertion with error checking */
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
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_INSERT && ret != SPI_OK_UPDATE && ret != SPI_OK_INSERT_RETURNING) {
        ereport(ERROR, 
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to insert/update metadata in kmersearch_highfreq_kmer_meta"),
                 errdetail("SPI_exec returned %d for metadata INSERT/UPDATE query", ret),
                 errhint("Check table structure and permissions")));
    }
    
    SPI_finish();
    
    pfree(query.data);
}

/*
 * Buffer management functions for 16-bit uintkeys
 */
void
kmersearch_add_to_buffer16(UintkeyBuffer16 *buffer, uint16 uintkey, const char *temp_table_name)
{
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_buffer16_to_table(buffer, temp_table_name);
    }
    
    /* Add new uintkey */
    buffer->uintkeys[buffer->count] = uintkey;
    buffer->count++;
}

/*
 * Buffer management functions for 32-bit uintkeys
 */
void
kmersearch_add_to_buffer32(UintkeyBuffer32 *buffer, uint32 uintkey, const char *temp_table_name)
{
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_buffer32_to_table(buffer, temp_table_name);
    }
    
    /* Add new uintkey */
    buffer->uintkeys[buffer->count] = uintkey;
    buffer->count++;
}

/*
 * Buffer management functions for 64-bit uintkeys
 */
void
kmersearch_add_to_buffer64(UintkeyBuffer64 *buffer, uint64 uintkey, const char *temp_table_name)
{
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_buffer64_to_table(buffer, temp_table_name);
    }
    
    /* Add new uintkey */
    buffer->uintkeys[buffer->count] = uintkey;
    buffer->count++;
}

/*
 * Flush 16-bit buffer contents to temporary table
 */
void
kmersearch_flush_buffer16_to_table(UintkeyBuffer16 *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    
    if (buffer->count == 0) return;
    
    initStringInfo(&query);
    
    /* Build bulk INSERT statement */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        /* Each uintkey has frequency_count of 1 (one row occurrence) */
        appendStringInfo(&query, "(%lu, 1)", 
                        (unsigned long)buffer->uintkeys[i]);
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
 * Flush 32-bit buffer contents to temporary table
 */
void
kmersearch_flush_buffer32_to_table(UintkeyBuffer32 *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    
    if (buffer->count == 0) return;
    
    initStringInfo(&query);
    
    /* Build bulk INSERT statement */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        /* Each uintkey has frequency_count of 1 (one row occurrence) */
        appendStringInfo(&query, "(%u, 1)", buffer->uintkeys[i]);
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
 * Flush 64-bit buffer contents to temporary table
 */
void
kmersearch_flush_buffer64_to_table(UintkeyBuffer64 *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i, j;
    int write_pos = 0;
    bool merged;
    
    if (buffer->count == 0) return;
    
    initStringInfo(&query);
    
    /* Build bulk INSERT statement */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        /* Each uintkey has frequency_count of 1 (one row occurrence) */
        appendStringInfo(&query, "(%lu, 1)", (unsigned long)buffer->uintkeys[i]);
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


/* Forward declaration for partition block mapping */
static PartitionBlockMapping kmersearch_map_global_to_partition_block(BlockNumber global_block, 
                                        KmerAnalysisSharedState *state);

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
    
    kmersearch_check_guc_initialization();
    
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
    /* Use max_parallel_maintenance_workers for analysis/maintenance operations like ANALYZE */
    parallel_workers = max_parallel_maintenance_workers;
    
    /* Comprehensive parameter validation */
    kmersearch_validate_analysis_parameters(table_oid, column_name, kmersearch_kmer_size);
    
    /* Log analysis start */
    
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
    
    values[0] = Int64GetDatum(result.total_rows);
    values[1] = Int32GetDatum(result.highfreq_kmers_count);
    values[2] = Int32GetDatum(result.parallel_workers_used);
    
    /* Extra validation for max_appearance_rate_used before conversion */
    if (result.max_appearance_rate_used < 0.0 || result.max_appearance_rate_used != result.max_appearance_rate_used) {
        ereport(WARNING, (errmsg("kmersearch_perform_highfreq_analysis: Detected corrupted max_appearance_rate_used (%f) during conversion, fixing", result.max_appearance_rate_used)));
        result.max_appearance_rate_used = 0.5;
    }
    
    values[3] = Float4GetDatum((float4)result.max_appearance_rate_used);  /* real type = 4-byte float */
    
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
    
    kmersearch_check_guc_initialization();
    
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
    shm_toc *toc;
    char *shm_pointer;
    bool table_locked = false;
    int64 total_rows;
    uint64 threshold_rows;
    uint64 rate_based_threshold;
    char temp_dir_to_delete[MAXPGPATH] = {0};  /* Initialize to empty string */
    KmerSearchTableType table_type;
    List *partition_oids = NIL;
    ListCell *lc;
    int partition_idx = 0;
    BlockNumber total_blocks_all_partitions = 0;
    PartitionBlockInfo *partition_blocks = NULL;
    int num_partitions = 0;
    
    
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
        
        /* Check table type */
        table_type = kmersearch_get_table_type(table_oid);
        
        if (table_type == KMERSEARCH_TABLE_PARTITIONED)
        {
            /* Get partition list */
            partition_oids = kmersearch_get_partition_oids(table_oid);
            num_partitions = list_length(partition_oids);
            
            if (num_partitions == 0)
            {
                elog(ERROR, "Partitioned table has no partitions");
            }
            
            /* Allocate partition block info array */
            partition_blocks = palloc(sizeof(PartitionBlockInfo) * num_partitions);
            
            /* Calculate total blocks from all partitions */
            partition_idx = 0;
            foreach(lc, partition_oids)
            {
                Oid part_oid = lfirst_oid(lc);
                Relation part_rel = table_open(part_oid, AccessShareLock);
                BlockNumber part_blocks = RelationGetNumberOfBlocks(part_rel);
                
                partition_blocks[partition_idx].partition_oid = part_oid;
                partition_blocks[partition_idx].start_block = total_blocks_all_partitions;
                partition_blocks[partition_idx].end_block = total_blocks_all_partitions + part_blocks - 1;
                
                
                total_blocks_all_partitions += part_blocks;
                partition_idx++;
                
                table_close(part_rel, AccessShareLock);
            }
            
            total_blocks = total_blocks_all_partitions;
            
            /* Validate column exists in first partition */
            {
                Oid first_part_oid = linitial_oid(partition_oids);
                column_attnum = get_attnum(first_part_oid, column_name);
                if (column_attnum == InvalidAttrNumber)
                    elog(ERROR, "Column \"%s\" does not exist in partitions", column_name);
                
                rel = table_open(first_part_oid, AccessShareLock);
                attr = TupleDescAttr(RelationGetDescr(rel), column_attnum - 1);
                column_type_oid = attr->atttypid;
                if (column_type_oid != TypenameGetTypid("dna2") && 
                    column_type_oid != TypenameGetTypid("dna4"))
                    elog(ERROR, "Column must be DNA2 or DNA4 type");
                table_close(rel, AccessShareLock);
            }
        }
        else
        {
            /* Regular table processing */
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
        }
        
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
        
        
        /* Enter parallel mode */
        EnterParallelMode();
        
        /* Create parallel context */
        pcxt = CreateParallelContext("pg_kmersearch", "kmersearch_analysis_worker", requested_workers);
        
        /* DSM size estimation - estimate each chunk separately as per PostgreSQL best practice */
        shm_toc_estimate_chunk(&pcxt->estimator, MAXALIGN(sizeof(KmerAnalysisSharedState)));
        if (table_type == KMERSEARCH_TABLE_PARTITIONED)
        {
            /* Additional space for partition block info array */
            shm_toc_estimate_chunk(&pcxt->estimator, MAXALIGN(sizeof(PartitionBlockInfo) * num_partitions));
            shm_toc_estimate_keys(&pcxt->estimator, 2); /* SHARED_STATE, PARTITION_BLOCKS */
        }
        else
        {
            shm_toc_estimate_keys(&pcxt->estimator, 1); /* SHARED_STATE */
        }
        elog(DEBUG1, "  - sizeof(KmerAnalysisSharedState) = %zu, MAXALIGN = %zu", 
             sizeof(KmerAnalysisSharedState), MAXALIGN(sizeof(KmerAnalysisSharedState)));
        
        /* Initialize DSM */
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
        shared_state->all_processed = false;
        shared_state->next_block = 0;
        shared_state->total_blocks = total_blocks;
        
        /* Initialize progress tracking atomics */
        pg_atomic_init_u64(&shared_state->total_rows_processed, 0);
        pg_atomic_init_u64(&shared_state->total_batches_committed, 0);
        
        /* Initialize temporary directory path for SQLite3 files */
        {
            Oid tablespace_oid;
            char base_temp_dir[MAXPGPATH];
            
            tablespace_oid = GetNextTempTableSpace();
            
            /* If no temp_tablespaces configured, use database's default tablespace */
            if (tablespace_oid == InvalidOid && OidIsValid(MyDatabaseTableSpace))
            {
                tablespace_oid = MyDatabaseTableSpace;
                elog(DEBUG1, "Using database default tablespace OID: %u", tablespace_oid);
            }
            
            /* Get the base temporary directory path (pgsql_tmp) */
            TempTablespacePath(base_temp_dir, tablespace_oid);
            
            /* Create a dedicated temporary directory for this analysis session */
            snprintf(shared_state->temp_dir_path, MAXPGPATH, "%s/pg_kmersearch_%d_%ld",
                     base_temp_dir, MyProcPid, GetCurrentTimestamp());
            
            /* Create the temporary directory (this will also create pgsql_tmp if it doesn't exist) */
            PathNameCreateTemporaryDir(base_temp_dir, shared_state->temp_dir_path);
            
            /* Save the temp directory path for cleanup */
            strlcpy(temp_dir_to_delete, shared_state->temp_dir_path, MAXPGPATH);
            
            elog(DEBUG1, "Created dedicated temp directory: %s (tablespace OID: %u)",
                 shared_state->temp_dir_path, tablespace_oid);
        }
        
        /* Initialize worker file lock */
        SpinLockInit(&shared_state->worker_file_lock);
        shared_state->worker_error_occurred = false;
        shared_state->total_rows = total_rows;
        
        /* Set partition-specific fields */
        shared_state->is_partitioned = (table_type == KMERSEARCH_TABLE_PARTITIONED);
        shared_state->num_partitions = num_partitions;
        
        if (shared_state->is_partitioned)
        {
            /* Allocate and copy partition block info to shared memory */
            PartitionBlockInfo *shm_partition_blocks;
            shm_partition_blocks = (PartitionBlockInfo *)shm_toc_allocate(toc, 
                sizeof(PartitionBlockInfo) * num_partitions);
            if (!shm_partition_blocks)
            {
                elog(ERROR, "Failed to allocate memory for partition blocks in shm_toc");
            }
            memcpy(shm_partition_blocks, partition_blocks, sizeof(PartitionBlockInfo) * num_partitions);
            shared_state->partition_blocks = shm_partition_blocks;
            shared_state->total_blocks_all_partitions = total_blocks_all_partitions;
            pg_atomic_init_u32(&shared_state->next_global_block, 0);
            shm_toc_insert(toc, KMERSEARCH_KEY_PARTITION_BLOCKS, shm_partition_blocks);
        }
        else
        {
            shared_state->partition_blocks = NULL;
            shared_state->total_blocks_all_partitions = 0;
        }
        
        elog(DEBUG1, "kmersearch_perform_highfreq_analysis_parallel: Initialized shared state - next_block=%u, total_blocks=%u, is_partitioned=%d",
             shared_state->next_block, shared_state->total_blocks, shared_state->is_partitioned);
        shm_toc_insert(toc, KMERSEARCH_KEY_SHARED_STATE, shared_state);
        
        /* Launch parallel workers */
        LaunchParallelWorkers(pcxt);
        result.parallel_workers_used = pcxt->nworkers_launched;
        
        /* Show initial progress */
        {
            BlockNumber total_blocks_to_process;
            
            if (shared_state->is_partitioned) {
                total_blocks_to_process = shared_state->total_blocks_all_partitions;
            } else {
                total_blocks_to_process = shared_state->total_blocks;
            }
            
            ereport(INFO, 
                    (errmsg("Starting high-frequency k-mer analysis: %lu rows in %u blocks with %d parallel workers", 
                            result.total_rows, total_blocks_to_process, result.parallel_workers_used)));
        }
        
        /* Wait for workers to complete */
        WaitForParallelWorkersToFinish(pcxt);
        
        /* Check for worker errors */
        if (shared_state->worker_error_occurred)
            elog(ERROR, "Parallel worker error: %s", shared_state->error_message);
        
        /* Announce start of parent aggregation */
        ereport(INFO,
                (errmsg("Parallel scan completed. Starting aggregation of results.")));
        
        /* Aggregate results from worker SQLite3 files */
        {
            sqlite3 *parent_db = NULL;
            sqlite3_stmt *select_stmt = NULL;
            char parent_db_path[MAXPGPATH];
            int rc;
            int fd;
            int worker_files_processed = 0;
            int total_worker_files = 0;
            
            /* Create parent SQLite3 database */
            snprintf(parent_db_path, MAXPGPATH, "%s/pg_kmersearch_%d_XXXXXX",
                     shared_state->temp_dir_path, MyProcPid);
            fd = mkstemp(parent_db_path);
            if (fd < 0)
            {
                ereport(ERROR,
                        (errcode_for_file_access(),
                         errmsg("could not create temporary file \"%s\": %m", parent_db_path)));
            }
            close(fd);
            
            rc = sqlite3_open(parent_db_path, &parent_db);
            if (rc != SQLITE_OK)
            {
                ereport(ERROR,
                        (errmsg("could not open parent SQLite3 database: %s", sqlite3_errmsg(parent_db))));
            }
            
            /* Configure parent database for better performance - MUST be done before creating any tables */
            sqlite3_exec(parent_db, "PRAGMA page_size = 8192", NULL, NULL, NULL);
            sqlite3_exec(parent_db, "PRAGMA cache_size = 50000", NULL, NULL, NULL);  /* 50000 * 8KB = ~400MB cache */
            sqlite3_exec(parent_db, "PRAGMA temp_store = MEMORY", NULL, NULL, NULL);
            sqlite3_exec(parent_db, "PRAGMA synchronous = OFF", NULL, NULL, NULL);
            sqlite3_exec(parent_db, "PRAGMA journal_mode = OFF", NULL, NULL, NULL);
            
            /* Create table in parent database */
            rc = sqlite3_exec(parent_db,
                             "CREATE TABLE kmersearch_highfreq_kmer ("
                             "uintkey INTEGER PRIMARY KEY, "
                             "nrow INTEGER)",
                             NULL, NULL, NULL);
            if (rc != SQLITE_OK)
            {
                ereport(ERROR,
                        (errmsg("could not create parent SQLite3 table: %s", sqlite3_errmsg(parent_db))));
            }
            
            /* Aggregate data from each worker */
            /* Count total worker files first */
            for (int i = 0; i < MAX_PARALLEL_WORKERS; i++)
            {
                if (strlen(shared_state->worker_temp_files[i]) > 0)
                    total_worker_files++;
            }
            
            for (int i = 0; i < MAX_PARALLEL_WORKERS; i++)
            {
                if (strlen(shared_state->worker_temp_files[i]) > 0)
                {
                    char attach_sql[MAXPGPATH + 100];
                    
                    worker_files_processed++;
                    
                    /* Report progress for each temp file */
                    ereport(INFO,
                            (errmsg("Processing temporary file %d of %d",
                                    worker_files_processed, total_worker_files)));
                    
                    /* Attach worker database */
                    snprintf(attach_sql, sizeof(attach_sql),
                            "ATTACH DATABASE '%s' AS worker_db",
                            shared_state->worker_temp_files[i]);
                    rc = sqlite3_exec(parent_db, attach_sql, NULL, NULL, NULL);
                    if (rc != SQLITE_OK)
                    {
                        elog(WARNING, "Failed to attach worker database: %s", sqlite3_errmsg(parent_db));
                        continue;
                    }
                    
                    /* Merge worker data into parent - use transaction for batch insert */
                    {
                        sqlite3_stmt *merge_stmt;
                        int merge_count = 0;
                        
                        /* Begin transaction for batch insert */
                        sqlite3_exec(parent_db, "BEGIN TRANSACTION", NULL, NULL, NULL);
                        
                        rc = sqlite3_prepare_v2(parent_db,
                                              "INSERT INTO main.kmersearch_highfreq_kmer (uintkey, nrow) "
                                              "VALUES (?, ?) "
                                              "ON CONFLICT(uintkey) DO UPDATE SET nrow = nrow + excluded.nrow",
                                              -1, &merge_stmt, NULL);
                        if (rc == SQLITE_OK)
                        {
                            rc = sqlite3_prepare_v2(parent_db,
                                                  "SELECT uintkey, nrow FROM worker_db.kmersearch_highfreq_kmer",
                                                  -1, &select_stmt, NULL);
                            if (rc == SQLITE_OK)
                            {
                                while (sqlite3_step(select_stmt) == SQLITE_ROW)
                                {
                                    int64 uintkey = sqlite3_column_int64(select_stmt, 0);
                                    int64 nrow = sqlite3_column_int64(select_stmt, 1);
                                    
                                    sqlite3_bind_int64(merge_stmt, 1, uintkey);
                                    sqlite3_bind_int64(merge_stmt, 2, nrow);
                                    sqlite3_step(merge_stmt);
                                    sqlite3_reset(merge_stmt);
                                    merge_count++;
                                }
                                sqlite3_finalize(select_stmt);
                            }
                            sqlite3_finalize(merge_stmt);
                        }
                        else
                        {
                            elog(WARNING, "Failed to prepare merge statement: %s", sqlite3_errmsg(parent_db));
                        }
                        
                        /* Commit transaction */
                        sqlite3_exec(parent_db, "COMMIT", NULL, NULL, NULL);
                    }
                    
                    /* Detach and delete worker database */
                    sqlite3_exec(parent_db, "DETACH DATABASE worker_db", NULL, NULL, NULL);
                    unlink(shared_state->worker_temp_files[i]);
                }
            }
            
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
            
            /* Count high-frequency k-mers from SQLite3 */
            {
                sqlite3_stmt *count_stmt;
                
                rc = sqlite3_prepare_v2(parent_db,
                                       "SELECT COUNT(*) FROM kmersearch_highfreq_kmer WHERE nrow > ?",
                                       -1, &count_stmt, NULL);
                if (rc == SQLITE_OK)
                {
                    sqlite3_bind_int64(count_stmt, 1, (int64)threshold_rows);
                    if (sqlite3_step(count_stmt) == SQLITE_ROW)
                    {
                        result.highfreq_kmers_count = sqlite3_column_int(count_stmt, 0);
                    }
                    sqlite3_finalize(count_stmt);
                }
            }
            
            /* temp_dir_to_delete already contains the path from initialization */
            
            /* Clean up parallel context and exit parallel mode BEFORE any SQL operations */
            DestroyParallelContext(pcxt);
            pcxt = NULL; /* Mark as destroyed to prevent double-free */
            ExitParallelMode();
            
            /* Now we can safely execute SQL operations */
            /* Save results to PostgreSQL */
            kmersearch_spi_connect_or_error();
            
            /* Prepare to extract high-frequency k-mers from SQLite3 */
            rc = sqlite3_prepare_v2(parent_db,
                                   "SELECT uintkey FROM kmersearch_highfreq_kmer WHERE nrow > ?",
                                   -1, &select_stmt, NULL);
            if (rc != SQLITE_OK)
            {
                sqlite3_close(parent_db);
                unlink(parent_db_path);
                ereport(ERROR,
                        (errmsg("could not prepare SQLite3 select statement: %s", sqlite3_errmsg(parent_db))));
            }
            
            sqlite3_bind_int64(select_stmt, 1, (int64)threshold_rows);
            
            /* Announce start of writing to PostgreSQL table */
            ereport(INFO,
                    (errmsg("Writing %d high-frequency k-mers to kmersearch_highfreq_kmer table...",
                            result.highfreq_kmers_count)));
            
            /* Insert high-frequency k-mers into PostgreSQL */
            {
                int kmers_written = 0;
                
                while (sqlite3_step(select_stmt) == SQLITE_ROW)
                {
                    uint64 uintkey;
                    StringInfoData insert_query;
                    int insert_ret;
                    int total_bits;
                    
                    /* Get uintkey based on total_bits */
                    total_bits = k_size * 2 + kmersearch_occur_bitlen;
                    if (total_bits <= 32)
                    {
                        uintkey = (uint64)sqlite3_column_int(select_stmt, 0);
                    }
                    else
                    {
                        uintkey = (uint64)sqlite3_column_int64(select_stmt, 0);
                    }
                    
                    /* Insert into PostgreSQL table - cast uint64 to int64 for bigint column */
                    initStringInfo(&insert_query);
                    appendStringInfo(&insert_query,
                        "INSERT INTO kmersearch_highfreq_kmer "
                        "(table_oid, column_name, uintkey, detection_reason) "
                        "VALUES (%u, %s, %ld, 'threshold') "
                        "ON CONFLICT (table_oid, column_name, uintkey) DO NOTHING",
                        table_oid, quote_literal_cstr(column_name), (int64)uintkey);
                    
                    insert_ret = SPI_exec(insert_query.data, 0);
                    if (insert_ret != SPI_OK_INSERT && insert_ret != SPI_OK_UPDATE)
                    {
                        elog(WARNING, "Failed to insert high-frequency k-mer");
                    }
                    pfree(insert_query.data);
                    
                    kmers_written++;
                    /* Report progress every 1000 k-mers */
                    if (kmers_written % 1000 == 0)
                    {
                        elog(INFO, "Written %d/%d high-frequency k-mers...",
                             kmers_written, result.highfreq_kmers_count);
                    }
                }
                
                /* Clean up SQLite3 resources */
                sqlite3_finalize(select_stmt);
                sqlite3_close(parent_db);
                unlink(parent_db_path);
                
                /* Report completion of writing */
                ereport(INFO,
                        (errmsg("Successfully wrote %d high-frequency k-mers to database.",
                                kmers_written)));
            }
        }
        
        /* Note: result.highfreq_kmers_count is already set above */
        
        /* Note: DestroyParallelContext and ExitParallelMode already called above */
        
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
    
    /* Normal cleanup - already done above, just unlock */
    UnlockRelationOid(table_oid, ExclusiveLock);
    
    /* Free partition blocks if allocated */
    if (partition_blocks)
    {
        pfree(partition_blocks);
    }
    
    /* Free partition OID list */
    if (partition_oids)
    {
        list_free(partition_oids);
    }
    
    /* Clean up the temporary directory if it was created */
    if (temp_dir_to_delete[0] != '\0')
    {
        PathNameDeleteTemporaryDir(temp_dir_to_delete);
        elog(DEBUG1, "Deleted temporary directory: %s", temp_dir_to_delete);
    }
    
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
 * Identity hash function for uint16 (uintkey is already a hash value)
 */
static uint32
kmersearch_uint16_identity_hash(const void *key, size_t keysize, void *arg)
{
    return (uint32)(*(const uint16 *)key);
}

/*
 * Identity hash function for uint32 (uintkey is already a hash value)
 */
static uint32
kmersearch_uint32_identity_hash(const void *key, size_t keysize, void *arg)
{
    return *(const uint32 *)key;
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
    
    /* Entry size depends on k-mer size */
    size_t entry_size;
    int total_bits = k_size * 2 + kmersearch_occur_bitlen;
    int max_entries;
    
    if (total_bits <= 16)
        entry_size = sizeof(uint16);
    else if (total_bits <= 32)
        entry_size = sizeof(uint32);
    else
        entry_size = sizeof(uint64);
    
    max_entries = (TARGET_MEMORY_MB * 1024 * 1024) / entry_size;
    
    if (max_entries < MIN_BUFFER_SIZE) return MIN_BUFFER_SIZE;
    if (max_entries > MAX_BUFFER_SIZE) return MAX_BUFFER_SIZE;
    return max_entries;
}

/*
 * Initialize k-mer buffer
 */
/*
 * Initialize 16-bit buffer
 */
void
kmersearch_init_buffer16(UintkeyBuffer16 *buffer, int k_size)
{
    buffer->capacity = kmersearch_calculate_buffer_size(k_size);
    buffer->uintkeys = (uint16 *) palloc0(buffer->capacity * sizeof(uint16));
    buffer->count = 0;
    buffer->kmer_size = k_size;
}

/*
 * Initialize 32-bit buffer
 */
void
kmersearch_init_buffer32(UintkeyBuffer32 *buffer, int k_size)
{
    buffer->capacity = kmersearch_calculate_buffer_size(k_size);
    buffer->uintkeys = (uint32 *) palloc0(buffer->capacity * sizeof(uint32));
    buffer->count = 0;
    buffer->kmer_size = k_size;
}

/*
 * Initialize 64-bit buffer
 */
void
kmersearch_init_buffer64(UintkeyBuffer64 *buffer, int k_size)
{
    buffer->capacity = kmersearch_calculate_buffer_size(k_size);
    buffer->uintkeys = (uint64 *) palloc0(buffer->capacity * sizeof(uint64));
    buffer->count = 0;
    buffer->kmer_size = k_size;
}

/*
 * Flush 16-bit hash buffer to temporary table
 */
static void
kmersearch_flush_hash_buffer16_to_table(UintkeyBuffer16 *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    int ret;
    
    if (buffer->count == 0) return;
    
    initStringInfo(&query);
    
    /* Build INSERT statement with multiple VALUES */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        /* Each uintkey has frequency_count of 1 (one row occurrence) */
        appendStringInfo(&query, "(%u, 1)", buffer->uintkeys[i]);
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
 * Flush 32-bit hash buffer to temporary table
 */
static void
kmersearch_flush_hash_buffer32_to_table(UintkeyBuffer32 *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    int ret;
    
    if (buffer->count == 0) return;
    
    initStringInfo(&query);
    
    /* Build INSERT statement with multiple VALUES */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        /* Each uintkey has frequency_count of 1 (one row occurrence) */
        appendStringInfo(&query, "(%u, 1)", buffer->uintkeys[i]);
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
 * Flush 64-bit hash buffer to temporary table
 */
static void
kmersearch_flush_hash_buffer64_to_table(UintkeyBuffer64 *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    int ret;
    
    if (buffer->count == 0) return;
    
    initStringInfo(&query);
    
    /* Build INSERT statement with multiple VALUES */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        /* Each uintkey has frequency_count of 1 (one row occurrence) */
        appendStringInfo(&query, "(%lu, 1)", (unsigned long)buffer->uintkeys[i]);
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
 * Add 16-bit hash value to buffer
 */
static void
kmersearch_add_hash_to_buffer16(UintkeyBuffer16 *buffer, uint16 uintkey_value, const char *temp_table_name)
{
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_hash_buffer16_to_table(buffer, temp_table_name);
    }
    
    /* Add new uintkey */
    buffer->uintkeys[buffer->count] = uintkey_value;
    buffer->count++;
}

/*
 * Add 32-bit hash value to buffer
 */
static void
kmersearch_add_hash_to_buffer32(UintkeyBuffer32 *buffer, uint32 uintkey_value, const char *temp_table_name)
{
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_hash_buffer32_to_table(buffer, temp_table_name);
    }
    
    /* Add new uintkey */
    buffer->uintkeys[buffer->count] = uintkey_value;
    buffer->count++;
}

/*
 * Add 64-bit hash value to buffer
 */
static void
kmersearch_add_hash_to_buffer64(UintkeyBuffer64 *buffer, uint64 uintkey_value, const char *temp_table_name)
{
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_hash_buffer64_to_table(buffer, temp_table_name);
    }
    
    /* Add new uintkey */
    buffer->uintkeys[buffer->count] = uintkey_value;
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
    void *uintkeys = NULL;
    BlockNumber current_block;
    bool isnull;
    Datum value;
    VarBit *sequence;
    VarBit **kmers;
    int nkeys;
    int j;
    int total_bits;
    
    /* Use passed parameters instead of determining them again */
    tupdesc = RelationGetDescr(rel);
    
    /* Initialize buffer based on k-mer size */
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    if (total_bits <= 16) {
        worker->buffer_type = 0;  /* 16-bit */
        worker->buffer = palloc0(sizeof(UintkeyBuffer16));
        kmersearch_init_buffer16((UintkeyBuffer16 *)worker->buffer, k_size);
    } else if (total_bits <= 32) {
        worker->buffer_type = 1;  /* 32-bit */
        worker->buffer = palloc0(sizeof(UintkeyBuffer32));
        kmersearch_init_buffer32((UintkeyBuffer32 *)worker->buffer, k_size);
    } else {
        worker->buffer_type = 2;  /* 64-bit */
        worker->buffer = palloc0(sizeof(UintkeyBuffer64));
        kmersearch_init_buffer64((UintkeyBuffer64 *)worker->buffer, k_size);
    }
    
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
        
        /* Extract k-mers from the sequence as uintkey format */
        {
            int ui;
            
            if (is_dna4_type) {
                kmersearch_extract_uintkey_from_dna4(sequence, &uintkeys, &nkeys);
            } else {
                kmersearch_extract_uintkey_from_dna2(sequence, &uintkeys, &nkeys);
            }
            
            /* No need to convert to Datum array - use uintkeys directly */
        }
        
        if (uintkeys == NULL || nkeys == 0) {
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
            
            /* Process each uintkey in this row - deduplicate within row */
            for (j = 0; j < nkeys; j++) {
                uint64_t uintkey_value;
                bool found;
                
                /* Get uintkey value directly from array */
                if (total_bits <= 16) {
                    uintkey_value = ((uint16 *)uintkeys)[j];
                } else if (total_bits <= 32) {
                    uintkey_value = ((uint32 *)uintkeys)[j];
                } else {
                    uintkey_value = ((uint64 *)uintkeys)[j];
                }
                
                /* Only add to buffer if not already seen in this row */
                hash_search(row_kmer_set, (void *) &uintkey_value, HASH_ENTER, &found);
                if (!found) {
                    /* Add to buffer (will flush to temp table if full) */
                    if (worker->buffer_type == 0) {
                        kmersearch_add_hash_to_buffer16((UintkeyBuffer16 *)worker->buffer, 
                                                       (uint16)uintkey_value, worker->temp_table_name);
                    } else if (worker->buffer_type == 1) {
                        kmersearch_add_hash_to_buffer32((UintkeyBuffer32 *)worker->buffer, 
                                                       (uint32)uintkey_value, worker->temp_table_name);
                    } else {
                        kmersearch_add_hash_to_buffer64((UintkeyBuffer64 *)worker->buffer, 
                                                       uintkey_value, worker->temp_table_name);
                    }
                }
            }
            
            /* Clean up row hash table */
            hash_destroy(row_kmer_set);
        }
        
        /* Cleanup uintkey array */
        if (uintkeys) {
            pfree(uintkeys);
        }
        }
        
        /* Release buffer and lock */
        UnlockReleaseBuffer(buffer);
    }
    
    /* Flush any remaining buffer contents */
    if (worker->buffer_type == 0) {
        kmersearch_flush_hash_buffer16_to_table((UintkeyBuffer16 *)worker->buffer, worker->temp_table_name);
        if (((UintkeyBuffer16 *)worker->buffer)->uintkeys) {
            pfree(((UintkeyBuffer16 *)worker->buffer)->uintkeys);
        }
    } else if (worker->buffer_type == 1) {
        kmersearch_flush_hash_buffer32_to_table((UintkeyBuffer32 *)worker->buffer, worker->temp_table_name);
        if (((UintkeyBuffer32 *)worker->buffer)->uintkeys) {
            pfree(((UintkeyBuffer32 *)worker->buffer)->uintkeys);
        }
    } else {
        kmersearch_flush_hash_buffer64_to_table((UintkeyBuffer64 *)worker->buffer, worker->temp_table_name);
        if (((UintkeyBuffer64 *)worker->buffer)->uintkeys) {
            pfree(((UintkeyBuffer64 *)worker->buffer)->uintkeys);
        }
    }
    
    /* Cleanup buffer structure itself */
    if (worker->buffer) {
        pfree(worker->buffer);
    }
}

/*
 * Parallel worker function for k-mer analysis (SQLite3-based)
 */
PGDLLEXPORT void
kmersearch_analysis_worker(dsm_segment *seg, shm_toc *toc)
{
    KmerAnalysisSharedState *shared_state = NULL;
    SQLiteWorkerContext ctx;
    bool has_work = true;
    int worker_id;
    int rc;
    int fd;
    
    /* Initialize context */
    memset(&ctx, 0, sizeof(SQLiteWorkerContext));
    
    /* Double-check we're in a parallel worker */
    if (!IsParallelWorker())
    {
        elog(ERROR, "kmersearch_analysis_worker called from non-parallel context");
    }
    
    /* Get shared state from shm_toc */
    shared_state = (KmerAnalysisSharedState *)shm_toc_lookup(toc, KMERSEARCH_KEY_SHARED_STATE, false);
    if (!shared_state) {
        elog(ERROR, "kmersearch_analysis_worker: Failed to get shared state");
    }
    
    /* Get worker ID (simplified - use PID modulo MAX_PARALLEL_WORKERS) */
    worker_id = MyProcPid % MAX_PARALLEL_WORKERS;
    
    /* If partitioned table, get partition block info from shared memory */
    if (shared_state->is_partitioned)
    {
        shared_state->partition_blocks = (PartitionBlockInfo *)shm_toc_lookup(toc, 
            KMERSEARCH_KEY_PARTITION_BLOCKS, false);
        if (!shared_state->partition_blocks)
        {
            elog(ERROR, "kmersearch_analysis_worker: Failed to get partition blocks");
        }
    }
    
    /* Calculate total bits for data type selection */
    ctx.total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;
    
    /* Cache type OIDs at worker start (avoid repeated lookups) */
    ctx.dna2_oid = TypenameGetTypid("dna2");
    ctx.dna4_oid = TypenameGetTypid("dna4");
    
    /* Get column type from table metadata */
    {
        Relation rel;
        TupleDesc tupdesc;
        Form_pg_attribute attr;
        Oid column_type_oid;
        
        rel = table_open(shared_state->table_oid, AccessShareLock);
        tupdesc = RelationGetDescr(rel);
        attr = TupleDescAttr(tupdesc, shared_state->column_attnum - 1);
        column_type_oid = attr->atttypid;
        table_close(rel, AccessShareLock);
        
        /* Validate column type */
        if (column_type_oid != ctx.dna2_oid && column_type_oid != ctx.dna4_oid)
        {
            elog(ERROR, "Column must be DNA2 or DNA4 type");
        }
        
        /* Store column type for later use */
        ctx.column_type_oid = column_type_oid;
    }
    
    /* Create temporary SQLite3 database */
    snprintf(ctx.db_path, MAXPGPATH, "%s/pg_kmersearch_%d_XXXXXX",
             shared_state->temp_dir_path, MyProcPid);
    
    fd = mkstemp(ctx.db_path);
    if (fd < 0)
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not create temporary file \"%s\": %m", ctx.db_path)));
    }
    close(fd);  /* SQLite3 will open it */
    
    /* Open SQLite3 database */
    rc = sqlite3_open(ctx.db_path, &ctx.db);
    if (rc != SQLITE_OK)
    {
        ereport(ERROR,
                (errmsg("could not open SQLite3 database: %s", sqlite3_errmsg(ctx.db))));
    }
    
    /* Configure SQLite3 for maximum performance - page_size MUST be set before creating tables */
    sqlite3_exec(ctx.db, "PRAGMA page_size = 8192", NULL, NULL, NULL);     /* Larger page size for bulk operations */
    sqlite3_exec(ctx.db, "PRAGMA cache_size = 20000", NULL, NULL, NULL);   /* 20000 * 8KB = ~160MB cache per worker */
    sqlite3_exec(ctx.db, "PRAGMA temp_store = MEMORY", NULL, NULL, NULL);
    sqlite3_exec(ctx.db, "PRAGMA synchronous = OFF", NULL, NULL, NULL);
    sqlite3_exec(ctx.db, "PRAGMA journal_mode = OFF", NULL, NULL, NULL);  /* No journaling for temp data */
    sqlite3_exec(ctx.db, "PRAGMA locking_mode = EXCLUSIVE", NULL, NULL, NULL);  /* No locking overhead */
    
    /* Create table */
    rc = sqlite3_exec(ctx.db, 
                     "CREATE TABLE kmersearch_highfreq_kmer ("
                     "uintkey INTEGER PRIMARY KEY, "
                     "nrow INTEGER)",
                     NULL, NULL, NULL);
    if (rc != SQLITE_OK)
    {
        ereport(ERROR,
                (errmsg("could not create SQLite3 table: %s", sqlite3_errmsg(ctx.db))));
    }
    
    /* Prepare INSERT statement */
    rc = sqlite3_prepare_v2(ctx.db,
                           "INSERT INTO kmersearch_highfreq_kmer (uintkey, nrow) "
                           "VALUES (?, ?) "
                           "ON CONFLICT(uintkey) DO UPDATE SET nrow = nrow + excluded.nrow",
                           -1, &ctx.insert_stmt, NULL);
    if (rc != SQLITE_OK)
    {
        ereport(ERROR,
                (errmsg("could not prepare SQLite3 statement: %s", sqlite3_errmsg(ctx.db))));
    }
    
    /* Dynamic work acquisition loop */
    {
        int blocks_processed = 0;
        
        while (has_work)
    {
        BlockNumber current_block;
        bool got_block = false;
        
        if (shared_state->is_partitioned)
        {
            /* Partitioned table: use global block counter */
            BlockNumber global_block;
            
            /* Atomic increment to get next global block */
            global_block = pg_atomic_fetch_add_u32(&shared_state->next_global_block, 1);
            
            if (global_block < shared_state->total_blocks_all_partitions)
            {
                /* Map global block to partition and local block */
                PartitionBlockMapping mapping = kmersearch_map_global_to_partition_block(
                    global_block, shared_state);
                
                
                /* Process block with batch - pass partition OID for partitioned tables */
                kmersearch_process_block_with_batch(mapping.local_block_number,
                                                   mapping.partition_oid,
                                                   shared_state, &ctx);
                blocks_processed++;
                
                got_block = true;
            }
        }
        else
        {
            /* Regular table: use existing logic */
            LWLockAcquire(&shared_state->mutex.lock, LW_EXCLUSIVE);
            
            if (shared_state->next_block < shared_state->total_blocks)
            {
                current_block = shared_state->next_block;
                shared_state->next_block++;
                LWLockRelease(&shared_state->mutex.lock);
                
                
                /* Process block with batch - use table_oid for regular tables */
                kmersearch_process_block_with_batch(current_block, 
                                                   shared_state->table_oid,
                                                   shared_state, &ctx);
                blocks_processed++;
                
                got_block = true;
            }
            else
            {
                shared_state->all_processed = true;
                LWLockRelease(&shared_state->mutex.lock);
            }
        }
        
            if (!got_block)
                has_work = false;
        }
    }
    
    /* Flush any remaining batch data */
    if (ctx.batch_hash && ctx.batch_count > 0)
    {
        uint64 total_rows;
        uint64 batch_num;
        
        kmersearch_flush_batch_to_sqlite(ctx.batch_hash, ctx.db, ctx.insert_stmt, ctx.total_bits);
        
        /* Update shared progress counters for final batch */
        total_rows = pg_atomic_add_fetch_u64(&shared_state->total_rows_processed, ctx.batch_count);
        batch_num = pg_atomic_add_fetch_u64(&shared_state->total_batches_committed, 1);
        
        /* Report final cumulative progress (deterministic) */
        ereport(INFO,
                (errmsg("Batch %lu completed: %lu total rows processed (final)",
                        (unsigned long)batch_num, (unsigned long)total_rows)));
    }
    
    /* Destroy hash table only once at the end */
    if (ctx.batch_hash)
    {
        hash_destroy(ctx.batch_hash);
        ctx.batch_hash = NULL;
    }
    
    /* Get final statistics before closing DB */
    if (ctx.db)
    {
        sqlite3_stmt *final_stats;
        rc = sqlite3_prepare_v2(ctx.db,
                                    "SELECT COUNT(*), SUM(nrow) FROM kmersearch_highfreq_kmer",
                                    -1, &final_stats, NULL);
        if (rc == SQLITE_OK && sqlite3_step(final_stats) == SQLITE_ROW)
        {
            int final_count = sqlite3_column_int(final_stats, 0);
            int64 final_nrow = sqlite3_column_int64(final_stats, 1);
            sqlite3_finalize(final_stats);
        }
    }
    
    /* Cleanup SQLite3 resources */
    if (ctx.insert_stmt)
        sqlite3_finalize(ctx.insert_stmt);
    if (ctx.db)
        sqlite3_close(ctx.db);
    
    /* Register temporary file path for parent process */
    kmersearch_register_worker_temp_file(shared_state, ctx.db_path, worker_id);
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
    int total_bits;
    const char *kmer_type;
    const char *freq_type;
    
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
    
    /* Calculate total bits to determine appropriate data types */
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    
    /* Determine appropriate data types based on total bits */
    if (total_bits <= 16) {
        kmer_type = "smallint";  /* int2 */
    } else if (total_bits <= 32) {
        kmer_type = "integer";   /* int4 */
    } else {
        kmer_type = "bigint";    /* int8 */
    }
    
    /* Frequency count always uses integer (int4) as it represents occurrence count */
    freq_type = "integer";
    
    appendStringInfo(&query, 
        "CREATE TEMP TABLE %s ("
        "kmer_data %s PRIMARY KEY, "
        "frequency_count %s"
        ")", final_table_name, kmer_type, freq_type);
    
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
 * Map global block number to partition and local block
 */
static PartitionBlockMapping
kmersearch_map_global_to_partition_block(BlockNumber global_block, 
                                        KmerAnalysisSharedState *state)
{
    PartitionBlockMapping result;
    int i;
    
    /* Progress reporting at specific block milestones (every 5%) */
    if (state->total_blocks_all_partitions > 20 && state->total_rows > 0) {
        BlockNumber milestone_interval = state->total_blocks_all_partitions / 20; /* 5% intervals */
        
        if (global_block > 0 && global_block % milestone_interval == 0) {
            int current_percentage = (int)((global_block * 100) / state->total_blocks_all_partitions);
            uint64 estimated_rows = ((uint64)state->total_rows * global_block) / state->total_blocks_all_partitions;
            
            if (current_percentage < 100) {
                ereport(INFO,
                        (errmsg("Progress: %d%% (estimated %lu / %lu rows processed)",
                                current_percentage, estimated_rows, state->total_rows)));
            }
        }
    }
    
    for (i = 0; i < state->num_partitions; i++)
    {
        if (global_block >= state->partition_blocks[i].start_block &&
            global_block <= state->partition_blocks[i].end_block)
        {
            result.partition_oid = state->partition_blocks[i].partition_oid;
            result.local_block_number = global_block - state->partition_blocks[i].start_block;
            return result;
        }
    }
    
    /* Should not happen if everything is correct */
    elog(ERROR, "Invalid global block number %u", global_block);
}

/*
 * Partition detection functions
 */

/*
 * Determine if a table is a regular table, partitioned table, or partition child
 */
KmerSearchTableType
kmersearch_get_table_type(Oid table_oid)
{
    Relation rel;
    KmerSearchTableType result;
    
    rel = table_open(table_oid, AccessShareLock);
    
    if (rel->rd_rel->relkind == RELKIND_PARTITIONED_TABLE)
    {
        result = KMERSEARCH_TABLE_PARTITIONED;
    }
    else if (rel->rd_rel->relispartition)
    {
        result = KMERSEARCH_TABLE_PARTITION_CHILD;
    }
    else
    {
        result = KMERSEARCH_TABLE_REGULAR;
    }
    
    table_close(rel, AccessShareLock);
    
    return result;
}

/*
 * Get list of partition OIDs for a partitioned table
 * Returns empty list if table is not partitioned
 */
List *
kmersearch_get_partition_oids(Oid parent_oid)
{
    Relation parent_rel;
    PartitionDesc partdesc;
    List *partition_oids = NIL;
    int i;
    
    /* Open parent table */
    parent_rel = table_open(parent_oid, AccessShareLock);
    
    /* Check if it's actually a partitioned table */
    if (parent_rel->rd_rel->relkind != RELKIND_PARTITIONED_TABLE)
    {
        table_close(parent_rel, AccessShareLock);
        return NIL;
    }
    
    /* Get partition descriptor - include detached partitions */
    partdesc = RelationGetPartitionDesc(parent_rel, false);
    
    /* Build list of partition OIDs */
    if (partdesc != NULL)
    {
        for (i = 0; i < partdesc->nparts; i++)
        {
            partition_oids = lappend_oid(partition_oids, partdesc->oids[i]);
        }
    }
    
    table_close(parent_rel, AccessShareLock);
    
    return partition_oids;
}

/*
 * Delete temporary SQLite3 files created by high-frequency k-mer analysis
 */
Datum
kmersearch_delete_tempfiles(PG_FUNCTION_ARGS)
{
    int deleted_count = 0;
    int64 deleted_size = 0;
    int error_count = 0;
    List *tablespace_oids = NIL;
    ListCell *lc;
    TupleDesc tupdesc;
    Datum values[3];
    bool nulls[3] = {false};
    HeapTuple tuple;
    int num_tablespaces;
    Oid *tablespace_array;
    time_t current_time;
    int64 file_size;
    
    /* Build result tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
    {
        ereport(ERROR,
                (errmsg("function returning record called in context "
                        "that cannot accept a record")));
    }
    
    /* Get temp_tablespaces */
    tablespace_array = palloc(sizeof(Oid) * 256);  /* Max 256 tablespaces */
    num_tablespaces = GetTempTablespaces(tablespace_array, 256);
    
    if (num_tablespaces == 0)
    {
        /* No temp_tablespaces configured, use MyDatabaseTableSpace */
        /* This will automatically fall back to system default if MyDatabaseTableSpace is InvalidOid */
        tablespace_oids = lappend_oid(tablespace_oids, MyDatabaseTableSpace);
    }
    else
    {
        /* Add temp_tablespaces to list */
        for (int i = 0; i < num_tablespaces; i++)
        {
            tablespace_oids = lappend_oid(tablespace_oids, tablespace_array[i]);
        }
    }
    
    /* Scan each tablespace's pgsql_tmp directory */
    foreach(lc, tablespace_oids)
    {
        Oid tablespace_oid = lfirst_oid(lc);
        char temp_path[MAXPGPATH];
        DIR *dir;
        struct dirent *de;
        
        /* Get temp directory path for this tablespace */
        TempTablespacePath(temp_path, tablespace_oid);
        
        /* Open directory */
        dir = opendir(temp_path);
        if (dir == NULL)
        {
            /* Directory doesn't exist, skip with warning */
            ereport(WARNING,
                    (errmsg("could not open temp directory \"%s\": %m", temp_path)));
            continue;
        }
        
        /* Scan directory for pgsql_tmp.kmersearch_* directories and pg_kmersearch_* files */
        while ((de = readdir(dir)) != NULL)
        {
            char full_path[MAXPGPATH];
            struct stat st;
            
            /* Check if it's a kmersearch temp directory */
            if (strncmp(de->d_name, "pgsql_tmp.kmersearch_", 21) == 0)
            {
                /* This is a kmersearch temporary directory */
                char subdir_path[MAXPGPATH];
                DIR *subdir;
                struct dirent *subde;
                
                snprintf(subdir_path, MAXPGPATH, "%s/%s", temp_path, de->d_name);
                
                /* Check directory age - skip directories created within last 60 seconds */
                if (stat(subdir_path, &st) < 0)
                {
                    ereport(WARNING,
                            (errmsg("could not stat directory \"%s\": %m", subdir_path)));
                    error_count++;
                    continue;
                }
                
                current_time = time(NULL);
                if (current_time - st.st_mtime < 60)
                {
                    ereport(INFO,
                            (errmsg("skipping recently created directory \"%s\"", subdir_path)));
                    continue;
                }
                
                /* Scan and delete files inside the directory */
                subdir = opendir(subdir_path);
                if (subdir == NULL)
                {
                    ereport(WARNING,
                            (errmsg("could not open directory \"%s\": %m", subdir_path)));
                    error_count++;
                    continue;
                }
                
                while ((subde = readdir(subdir)) != NULL)
                {
                    char file_path[MAXPGPATH];
                    
                    if (strcmp(subde->d_name, ".") == 0 || strcmp(subde->d_name, "..") == 0)
                        continue;
                    
                    snprintf(file_path, MAXPGPATH, "%s/%s", subdir_path, subde->d_name);
                    
                    if (stat(file_path, &st) < 0)
                    {
                        ereport(WARNING,
                                (errmsg("could not stat file \"%s\": %m", file_path)));
                        error_count++;
                        continue;
                    }
                    
                    if (S_ISREG(st.st_mode))
                    {
                        file_size = st.st_size;
                        if (unlink(file_path) < 0)
                        {
                            ereport(WARNING,
                                    (errmsg("could not delete file \"%s\": %m", file_path)));
                            error_count++;
                        }
                        else
                        {
                            deleted_count++;
                            deleted_size += file_size;
                            ereport(INFO,
                                    (errmsg("deleted temporary file \"%s\" (size: %ld bytes)",
                                            file_path, file_size)));
                        }
                    }
                }
                
                closedir(subdir);
                
                /* Note: We do NOT use PathNameDeleteTemporaryDir here */
                /* We only delete the files inside, leaving directory removal to PostgreSQL */
            }
            /* Also check for old-style pg_kmersearch_* files directly in pgsql_tmp */
            else if (strncmp(de->d_name, "pg_kmersearch_", 14) == 0)
            {
                /* Build full path */
                snprintf(full_path, MAXPGPATH, "%s/%s", temp_path, de->d_name);
                
                /* Get file info */
                if (stat(full_path, &st) < 0)
                {
                    ereport(WARNING,
                            (errmsg("could not stat file \"%s\": %m", full_path)));
                    error_count++;
                    continue;
                }
                
                /* Check if it's a regular file */
                if (!S_ISREG(st.st_mode))
                    continue;
                
                /* Check file age - skip files created within last 60 seconds */
                current_time = time(NULL);
                if (current_time - st.st_mtime < 60)
                {
                    ereport(INFO,
                            (errmsg("skipping recently created file \"%s\"", full_path)));
                    continue;
                }
                
                /* Record file size */
                file_size = st.st_size;
                
                /* Delete the file */
                if (unlink(full_path) < 0)
                {
                    ereport(WARNING,
                            (errmsg("could not delete file \"%s\": %m", full_path)));
                    error_count++;
                }
                else
                {
                    deleted_count++;
                    deleted_size += file_size;
                    ereport(INFO,
                            (errmsg("deleted temporary file \"%s\" (size: %ld bytes)",
                                    full_path, file_size)));
                }
            }
        }
        
        closedir(dir);
    }
    
    /* Free allocated memory */
    pfree(tablespace_array);
    list_free(tablespace_oids);
    
    /* Build result tuple */
    values[0] = Int32GetDatum(deleted_count);
    values[1] = Int64GetDatum(deleted_size);
    values[2] = Int32GetDatum(error_count);
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}

/*
 * Register worker temporary file path in shared memory
 */
static void
kmersearch_register_worker_temp_file(KmerAnalysisSharedState *shared_state,
                                    const char *file_path, int worker_id)
{
    SpinLockAcquire(&shared_state->worker_file_lock);
    strlcpy(shared_state->worker_temp_files[worker_id], file_path, MAXPGPATH);
    SpinLockRelease(&shared_state->worker_file_lock);
}

/*
 * Clear hash table for reuse
 */
static void
kmersearch_clear_batch_hash(HTAB *batch_hash, int total_bits)
{
    HASH_SEQ_STATUS status;
    void *entry;
    List *keys_to_remove = NIL;
    ListCell *lc;
    
    /* First pass: collect keys to remove */
    hash_seq_init(&status, batch_hash);
    while ((entry = hash_seq_search(&status)) != NULL)
    {
        if (total_bits <= 16)
        {
            TempKmerFreqEntry16 *e = (TempKmerFreqEntry16 *)entry;
            uint16 *key_copy = (uint16 *)palloc(sizeof(uint16));
            *key_copy = e->uintkey;
            keys_to_remove = lappend(keys_to_remove, key_copy);
        }
        else if (total_bits <= 32)
        {
            TempKmerFreqEntry32 *e = (TempKmerFreqEntry32 *)entry;
            uint32 *key_copy = (uint32 *)palloc(sizeof(uint32));
            *key_copy = e->uintkey;
            keys_to_remove = lappend(keys_to_remove, key_copy);
        }
        else
        {
            TempKmerFreqEntry64 *e = (TempKmerFreqEntry64 *)entry;
            uint64 *key_copy = (uint64 *)palloc(sizeof(uint64));
            *key_copy = e->uintkey;
            keys_to_remove = lappend(keys_to_remove, key_copy);
        }
    }
    
    /* Second pass: remove entries using the collected keys */
    foreach(lc, keys_to_remove)
    {
        void *key = lfirst(lc);
        hash_search(batch_hash, key, HASH_REMOVE, NULL);
        pfree(key);
    }
    
    /* Clean up the list */
    list_free(keys_to_remove);
}

/*
 * Flush batch hash table to SQLite3
 */
static void
kmersearch_flush_batch_to_sqlite(HTAB *batch_hash, sqlite3 *db,
                                sqlite3_stmt *stmt, int total_bits)
{
    HASH_SEQ_STATUS status;
    void *entry;
    int rc;
    
    /* Begin transaction */
    rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    if (rc != SQLITE_OK)
    {
        ereport(ERROR,
                (errmsg("SQLite3 BEGIN TRANSACTION failed: %s", sqlite3_errmsg(db))));
    }
    
    /* Scan hash table and write to SQLite3 */
    hash_seq_init(&status, batch_hash);
    
    while ((entry = hash_seq_search(&status)) != NULL)
    {
        if (total_bits <= 16)
        {
            TempKmerFreqEntry16 *e = (TempKmerFreqEntry16 *)entry;
            sqlite3_bind_int(stmt, 1, (int16)e->uintkey);
            sqlite3_bind_int64(stmt, 2, (int64)e->nrow);
        }
        else if (total_bits <= 32)
        {
            TempKmerFreqEntry32 *e = (TempKmerFreqEntry32 *)entry;
            sqlite3_bind_int(stmt, 1, (int32)e->uintkey);
            sqlite3_bind_int64(stmt, 2, (int64)e->nrow);
        }
        else
        {
            TempKmerFreqEntry64 *e = (TempKmerFreqEntry64 *)entry;
            sqlite3_bind_int64(stmt, 1, (int64)e->uintkey);
            sqlite3_bind_int64(stmt, 2, (int64)e->nrow);
        }
        
        rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE)
        {
            ereport(ERROR,
                    (errmsg("SQLite3 INSERT failed: %s", sqlite3_errmsg(db))));
        }
        sqlite3_reset(stmt);
    }
    
    /* Commit transaction */
    rc = sqlite3_exec(db, "COMMIT", NULL, NULL, NULL);
    if (rc != SQLITE_OK)
    {
        ereport(ERROR,
                (errmsg("SQLite3 COMMIT failed: %s", sqlite3_errmsg(db))));
    }
}

/*
 * Process a block with batch aggregation
 */
static void
kmersearch_process_block_with_batch(BlockNumber block,
                                   Oid table_oid,
                                   KmerAnalysisSharedState *shared_state,
                                   SQLiteWorkerContext *ctx)
{
    Relation rel;
    Buffer buffer;
    Page page;
    OffsetNumber maxoff;
    
    /* Initialize batch hash if needed */
    if (ctx->batch_hash == NULL)
    {
        HASHCTL hashctl;
        memset(&hashctl, 0, sizeof(hashctl));
        
        if (ctx->total_bits <= 16)
        {
            hashctl.keysize = sizeof(uint16);
            hashctl.entrysize = sizeof(TempKmerFreqEntry16);
        }
        else if (ctx->total_bits <= 32)
        {
            hashctl.keysize = sizeof(uint32);
            hashctl.entrysize = sizeof(TempKmerFreqEntry32);
        }
        else
        {
            hashctl.keysize = sizeof(uint64);
            hashctl.entrysize = sizeof(TempKmerFreqEntry64);
        }
        
        ctx->batch_hash = hash_create("KmerBatchHash",
                                     kmersearch_highfreq_analysis_hashtable_size,
                                     &hashctl, HASH_ELEM | HASH_BLOBS);
    }
    
    /* Open table and read block */
    rel = table_open(table_oid, AccessShareLock);
    buffer = ReadBuffer(rel, block);
    LockBuffer(buffer, BUFFER_LOCK_SHARE);
    page = BufferGetPage(buffer);
    maxoff = PageGetMaxOffsetNumber(page);
    
    /* Process each tuple in the page */
    for (OffsetNumber offnum = FirstOffsetNumber;
         offnum <= maxoff;
         offnum = OffsetNumberNext(offnum))
    {
        ItemId itemid = PageGetItemId(page, offnum);
        HeapTupleData tuple;
        bool isnull;
        Datum datum;
        void *kmer_array = NULL;
        int kmer_count;
        
        if (!ItemIdIsNormal(itemid))
            continue;
        
        tuple.t_data = (HeapTupleHeader) PageGetItem(page, itemid);
        tuple.t_len = ItemIdGetLength(itemid);
        
        /* Get sequence data */
        datum = heap_getattr(&tuple, shared_state->column_attnum,
                           RelationGetDescr(rel), &isnull);
        if (isnull)
            continue;
        
        /* Extract k-mers based on data type (type OIDs are cached at worker start) */
        if (ctx->column_type_oid == ctx->dna2_oid)
        {
            VarBit *seq = DatumGetVarBitP(datum);
            kmersearch_extract_uintkey_from_dna2(seq, &kmer_array, &kmer_count);
        }
        else if (ctx->column_type_oid == ctx->dna4_oid)
        {
            VarBit *seq = DatumGetVarBitP(datum);
            kmersearch_extract_uintkey_from_dna4(seq, &kmer_array, &kmer_count);
        }
        else
        {
            elog(ERROR, "Unsupported column type");
            continue;
        }
        
        
        /* Add k-mers to batch hash table */
        for (int i = 0; i < kmer_count; i++)
        {
            bool found;
            void *entry;
            
            if (ctx->total_bits <= 16)
            {
                uint16 uintkey = ((uint16 *)kmer_array)[i];
                TempKmerFreqEntry16 *freq_entry;
                
                entry = hash_search(ctx->batch_hash, &uintkey, HASH_ENTER, &found);
                freq_entry = (TempKmerFreqEntry16 *)entry;
                if (!found)
                {
                    freq_entry->uintkey = uintkey;
                    freq_entry->nrow = 0;
                }
                freq_entry->nrow++;
            }
            else if (ctx->total_bits <= 32)
            {
                uint32 uintkey = ((uint32 *)kmer_array)[i];
                TempKmerFreqEntry32 *freq_entry;
                
                entry = hash_search(ctx->batch_hash, &uintkey, HASH_ENTER, &found);
                freq_entry = (TempKmerFreqEntry32 *)entry;
                if (!found)
                {
                    freq_entry->uintkey = uintkey;
                    freq_entry->nrow = 0;
                }
                freq_entry->nrow++;
            }
            else
            {
                uint64 uintkey = ((uint64 *)kmer_array)[i];
                TempKmerFreqEntry64 *freq_entry;
                
                entry = hash_search(ctx->batch_hash, &uintkey, HASH_ENTER, &found);
                freq_entry = (TempKmerFreqEntry64 *)entry;
                if (!found)
                {
                    freq_entry->uintkey = uintkey;
                    freq_entry->nrow = 0;
                }
                freq_entry->nrow++;
            }
        }
        
        if (kmer_array)
            pfree(kmer_array);
        
        ctx->batch_count++;
        
        /* Flush batch if size limit reached */
        if (ctx->batch_count >= kmersearch_highfreq_analysis_batch_size)
        {
            uint64 total_rows;
            uint64 batch_num;
            
            kmersearch_flush_batch_to_sqlite(ctx->batch_hash, ctx->db,
                                           ctx->insert_stmt, ctx->total_bits);
            /* Clear hash table for reuse instead of destroying */
            kmersearch_clear_batch_hash(ctx->batch_hash, ctx->total_bits);
            
            /* Update shared progress counters */
            total_rows = pg_atomic_add_fetch_u64(&shared_state->total_rows_processed, ctx->batch_count);
            batch_num = pg_atomic_add_fetch_u64(&shared_state->total_batches_committed, 1);
            
            /* Report cumulative progress (deterministic) */
            ereport(INFO,
                    (errmsg("Batch %lu completed: %lu total rows processed",
                            (unsigned long)batch_num, (unsigned long)total_rows)));
            
            ctx->batch_count = 0;
        }
    }
    
    UnlockReleaseBuffer(buffer);
    table_close(rel, AccessShareLock);
}

