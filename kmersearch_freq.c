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
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include "storage/fd.h"
#include "storage/bufmgr.h"
#include "access/genam.h"
#include "utils/lsyscache.h"

/* PostgreSQL function info declarations for frequency functions */
PG_FUNCTION_INFO_V1(kmersearch_perform_highfreq_analysis);
PG_FUNCTION_INFO_V1(kmersearch_undo_highfreq_analysis);
PG_FUNCTION_INFO_V1(kmersearch_delete_tempfiles);

/* Temporary k-mer frequency entry structures for batch processing */
typedef struct TempKmerFreqEntry16
{
    uint16      uintkey;
    uint64      appearance_nrow;      /* Number of rows containing this k-mer */
} TempKmerFreqEntry16;

typedef struct TempKmerFreqEntry32
{
    uint32      uintkey;
    uint64      appearance_nrow;      /* Number of rows containing this k-mer */
} TempKmerFreqEntry32;

typedef struct TempKmerFreqEntry64
{
    uint64      uintkey;
    uint64      appearance_nrow;      /* Number of rows containing this k-mer */
} TempKmerFreqEntry64;

/* File hash table helper functions */
static void kmersearch_register_worker_temp_file(KmerAnalysisSharedState *shared_state,
                                                const char *file_path, int worker_id);
static void kmersearch_flush_batch_to_fht(FileHashWorkerContext *ctx);
static void kmersearch_process_block_with_batch(BlockNumber block,
                                               Oid table_oid,
                                               KmerAnalysisSharedState *shared_state,
                                               FileHashWorkerContext *ctx);

/* Parallel merge functions */
PGDLLEXPORT void kmersearch_parallel_merge_worker(dsm_segment *seg, shm_toc *toc);
static void kmersearch_aggregate_temp_files_parallel(char file_paths[][MAXPGPATH],
                                                    int num_files,
                                                    const char *temp_dir_path,
                                                    int total_bits);

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
    /* Combining kmer_data and frequency_count into single uintkey */
    
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
}/* Forward declaration for partition block mapping */
static PartitionBlockMapping kmersearch_map_global_to_partition_block(BlockNumber global_block, 
                                        KmerAnalysisSharedState *state);/*
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
        
        /* Store table name and column name for progress reporting */
        {
            char *table_name_str = get_rel_name(table_oid);
            if (table_name_str)
            {
                strlcpy(shared_state->table_name, table_name_str, NAMEDATALEN);
                pfree(table_name_str);
            }
            else
            {
                snprintf(shared_state->table_name, NAMEDATALEN, "oid_%u", table_oid);
            }
        }
        strlcpy(shared_state->column_name, column_name, NAMEDATALEN);
        
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
            snprintf(shared_state->temp_dir_path, MAXPGPATH, "%s/pg_kmersearch_db%u_%d_%ld",
                     base_temp_dir, MyDatabaseId, MyProcPid, GetCurrentTimestamp());
            
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

        /* Aggregate results from worker file hash tables */
        {
            char worker_files[MAX_PARALLEL_WORKERS][MAXPGPATH];
            char *aggregated_file_path;
            int num_worker_files = 0;
            int total_bits;

            total_bits = k_size * 2 + kmersearch_occur_bitlen;

            /* Collect worker temp files for parallel aggregation */
            for (int i = 0; i < MAX_PARALLEL_WORKERS; i++)
            {
                if (strlen(shared_state->worker_temp_files[i]) > 0)
                {
                    strlcpy(worker_files[num_worker_files], shared_state->worker_temp_files[i], MAXPGPATH);
                    num_worker_files++;
                }
            }

            if (num_worker_files == 0)
            {
                ereport(ERROR,
                        (errmsg("No worker temporary files found")));
            }

            /* If we have multiple worker files, use parallel aggregation */
            if (num_worker_files > 1)
            {
                ereport(INFO,
                        (errmsg("Starting parallel aggregation of %d temporary files", num_worker_files)));

                kmersearch_aggregate_temp_files_parallel(worker_files, num_worker_files,
                                                        shared_state->temp_dir_path, total_bits);
            }

            /* Use the first worker file (aggregated if multiple, original if single) */
            aggregated_file_path = worker_files[0];

            /* Calculate threshold based on GUC variables */
            rate_based_threshold = (uint64)(result.total_rows * kmersearch_max_appearance_rate);

            /* Determine final threshold: min of rate-based and nrow-based (excluding 0) */
            if (kmersearch_max_appearance_nrow > 0) {
                threshold_rows = (rate_based_threshold < (uint64)kmersearch_max_appearance_nrow) ?
                                rate_based_threshold : (uint64)kmersearch_max_appearance_nrow;
            } else {
                threshold_rows = rate_based_threshold;
            }
            elog(DEBUG1, "Threshold calculation: total_rows=%ld, rate=%.2f, rate_based=%lu, nrow=%d, final=%lu",
                 result.total_rows, kmersearch_max_appearance_rate, (unsigned long)rate_based_threshold,
                 kmersearch_max_appearance_nrow, (unsigned long)threshold_rows);

            result.max_appearance_nrow_used = threshold_rows;

            /* Count and extract high-frequency k-mers from file hash table */
            {
                int highfreq_count = 0;

                if (total_bits <= 16)
                {
                    FileHashTable16Context *fht_ctx = kmersearch_fht16_open(aggregated_file_path);
                    FileHashTableIterator16 iter;
                    uint16 uintkey;
                    uint64 appearance_nrow;

                    kmersearch_fht16_iterator_init(&iter, fht_ctx);
                    while (kmersearch_fht16_iterate(&iter, &uintkey, &appearance_nrow))
                    {
                        if (appearance_nrow > threshold_rows)
                            highfreq_count++;
                    }
                    kmersearch_fht16_close(fht_ctx);
                }
                else if (total_bits <= 32)
                {
                    FileHashTable32Context *fht_ctx = kmersearch_fht32_open(aggregated_file_path);
                    FileHashTableIterator32 iter;
                    uint32 uintkey;
                    uint64 appearance_nrow;

                    kmersearch_fht32_iterator_init(&iter, fht_ctx);
                    while (kmersearch_fht32_iterate(&iter, &uintkey, &appearance_nrow))
                    {
                        if (appearance_nrow > threshold_rows)
                            highfreq_count++;
                    }
                    kmersearch_fht32_close(fht_ctx);
                }
                else
                {
                    FileHashTable64Context *fht_ctx = kmersearch_fht64_open(aggregated_file_path);
                    FileHashTableIterator64 iter;
                    uint64 uintkey;
                    uint64 appearance_nrow;

                    kmersearch_fht64_iterator_init(&iter, fht_ctx);
                    while (kmersearch_fht64_iterate(&iter, &uintkey, &appearance_nrow))
                    {
                        if (appearance_nrow > threshold_rows)
                            highfreq_count++;
                    }
                    kmersearch_fht64_close(fht_ctx);
                }

                result.highfreq_kmers_count = highfreq_count;
            }

            /* Clean up parallel context and exit parallel mode BEFORE any SQL operations */
            DestroyParallelContext(pcxt);
            pcxt = NULL;
            ExitParallelMode();

            /* Save results to PostgreSQL */
            kmersearch_spi_connect_or_error();

            ereport(INFO,
                    (errmsg("Writing %d high-frequency k-mers to kmersearch_highfreq_kmer table...",
                            result.highfreq_kmers_count)));

            /* Insert high-frequency k-mers into PostgreSQL */
            {
                int kmers_written = 0;

                if (total_bits <= 16)
                {
                    FileHashTable16Context *fht_ctx = kmersearch_fht16_open(aggregated_file_path);
                    FileHashTableIterator16 iter;
                    uint16 uintkey;
                    uint64 appearance_nrow;

                    kmersearch_fht16_iterator_init(&iter, fht_ctx);
                    while (kmersearch_fht16_iterate(&iter, &uintkey, &appearance_nrow))
                    {
                        if (appearance_nrow > threshold_rows)
                        {
                            StringInfoData insert_query;
                            int insert_ret;

                            initStringInfo(&insert_query);
                            appendStringInfo(&insert_query,
                                "INSERT INTO kmersearch_highfreq_kmer "
                                "(table_oid, column_name, uintkey, appearance_nrow, detection_reason) "
                                "VALUES (%u, %s, %ld, %lu, 'threshold') "
                                "ON CONFLICT (table_oid, column_name, uintkey) DO UPDATE SET "
                                "appearance_nrow = EXCLUDED.appearance_nrow",
                                table_oid, quote_literal_cstr(column_name),
                                (int64)uintkey, (unsigned long)appearance_nrow);

                            insert_ret = SPI_exec(insert_query.data, 0);
                            if (insert_ret != SPI_OK_INSERT && insert_ret != SPI_OK_UPDATE)
                                elog(WARNING, "Failed to insert high-frequency k-mer");
                            pfree(insert_query.data);

                            kmers_written++;
                            if (kmers_written % 1000 == 0)
                                elog(INFO, "Written %d/%d high-frequency k-mers...",
                                     kmers_written, result.highfreq_kmers_count);
                        }
                    }
                    kmersearch_fht16_close(fht_ctx);
                }
                else if (total_bits <= 32)
                {
                    FileHashTable32Context *fht_ctx = kmersearch_fht32_open(aggregated_file_path);
                    FileHashTableIterator32 iter;
                    uint32 uintkey;
                    uint64 appearance_nrow;

                    kmersearch_fht32_iterator_init(&iter, fht_ctx);
                    while (kmersearch_fht32_iterate(&iter, &uintkey, &appearance_nrow))
                    {
                        if (appearance_nrow > threshold_rows)
                        {
                            StringInfoData insert_query;
                            int insert_ret;

                            initStringInfo(&insert_query);
                            appendStringInfo(&insert_query,
                                "INSERT INTO kmersearch_highfreq_kmer "
                                "(table_oid, column_name, uintkey, appearance_nrow, detection_reason) "
                                "VALUES (%u, %s, %ld, %lu, 'threshold') "
                                "ON CONFLICT (table_oid, column_name, uintkey) DO UPDATE SET "
                                "appearance_nrow = EXCLUDED.appearance_nrow",
                                table_oid, quote_literal_cstr(column_name),
                                (int64)uintkey, (unsigned long)appearance_nrow);

                            insert_ret = SPI_exec(insert_query.data, 0);
                            if (insert_ret != SPI_OK_INSERT && insert_ret != SPI_OK_UPDATE)
                                elog(WARNING, "Failed to insert high-frequency k-mer");
                            pfree(insert_query.data);

                            kmers_written++;
                            if (kmers_written % 1000 == 0)
                                elog(INFO, "Written %d/%d high-frequency k-mers...",
                                     kmers_written, result.highfreq_kmers_count);
                        }
                    }
                    kmersearch_fht32_close(fht_ctx);
                }
                else
                {
                    FileHashTable64Context *fht_ctx = kmersearch_fht64_open(aggregated_file_path);
                    FileHashTableIterator64 iter;
                    uint64 uintkey;
                    uint64 appearance_nrow;

                    kmersearch_fht64_iterator_init(&iter, fht_ctx);
                    while (kmersearch_fht64_iterate(&iter, &uintkey, &appearance_nrow))
                    {
                        if (appearance_nrow > threshold_rows)
                        {
                            StringInfoData insert_query;
                            int insert_ret;

                            initStringInfo(&insert_query);
                            appendStringInfo(&insert_query,
                                "INSERT INTO kmersearch_highfreq_kmer "
                                "(table_oid, column_name, uintkey, appearance_nrow, detection_reason) "
                                "VALUES (%u, %s, %ld, %lu, 'threshold') "
                                "ON CONFLICT (table_oid, column_name, uintkey) DO UPDATE SET "
                                "appearance_nrow = EXCLUDED.appearance_nrow",
                                table_oid, quote_literal_cstr(column_name),
                                (int64)uintkey, (unsigned long)appearance_nrow);

                            insert_ret = SPI_exec(insert_query.data, 0);
                            if (insert_ret != SPI_OK_INSERT && insert_ret != SPI_OK_UPDATE)
                                elog(WARNING, "Failed to insert high-frequency k-mer");
                            pfree(insert_query.data);

                            kmers_written++;
                            if (kmers_written % 1000 == 0)
                                elog(INFO, "Written %d/%d high-frequency k-mers...",
                                     kmers_written, result.highfreq_kmers_count);
                        }
                    }
                    kmersearch_fht64_close(fht_ctx);
                }

                /* Delete aggregated file */
                unlink(aggregated_file_path);

                ereport(INFO,
                        (errmsg("Successfully wrote %d high-frequency k-mers to database.",
                                kmers_written)));
            }
        }
        
        /* result.highfreq_kmers_count is already set above */
        
        /* DestroyParallelContext and ExitParallelMode already called above */
        
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
 * Parallel worker function for k-mer analysis (file hash table based)
 */
PGDLLEXPORT void
kmersearch_analysis_worker(dsm_segment *seg, shm_toc *toc)
{
    KmerAnalysisSharedState *shared_state = NULL;
    FileHashWorkerContext ctx;
    bool has_work = true;
    int worker_id;
    int fd;

    memset(&ctx, 0, sizeof(FileHashWorkerContext));

    if (!IsParallelWorker())
    {
        elog(ERROR, "kmersearch_analysis_worker called from non-parallel context");
    }

    shared_state = (KmerAnalysisSharedState *)shm_toc_lookup(toc, KMERSEARCH_KEY_SHARED_STATE, false);
    if (!shared_state) {
        elog(ERROR, "kmersearch_analysis_worker: Failed to get shared state");
    }

    worker_id = MyProcPid % MAX_PARALLEL_WORKERS;

    if (shared_state->is_partitioned)
    {
        shared_state->partition_blocks = (PartitionBlockInfo *)shm_toc_lookup(toc,
            KMERSEARCH_KEY_PARTITION_BLOCKS, false);
        if (!shared_state->partition_blocks)
        {
            elog(ERROR, "kmersearch_analysis_worker: Failed to get partition blocks");
        }
    }

    ctx.total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;

    ctx.dna2_oid = TypenameGetTypid("dna2");
    ctx.dna4_oid = TypenameGetTypid("dna4");

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

        if (column_type_oid != ctx.dna2_oid && column_type_oid != ctx.dna4_oid)
        {
            elog(ERROR, "Column must be DNA2 or DNA4 type");
        }

        ctx.column_type_oid = column_type_oid;
    }

    /* Create temporary file for file hash table */
    snprintf(ctx.file_path, MAXPGPATH, "%s/pg_kmersearch_XXXXXX",
             shared_state->temp_dir_path);

    fd = mkstemp(ctx.file_path);
    if (fd < 0)
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not create temporary file \"%s\": %m", ctx.file_path)));
    }
    close(fd);
    unlink(ctx.file_path);

    /* Create file hash table based on total_bits */
    if (ctx.total_bits <= 16)
    {
        ctx.fht_ctx = kmersearch_fht16_create(ctx.file_path);
    }
    else if (ctx.total_bits <= 32)
    {
        ctx.fht_ctx = kmersearch_fht32_create(ctx.file_path, 0);
    }
    else
    {
        ctx.fht_ctx = kmersearch_fht64_create(ctx.file_path, 0);
    }

    ctx.batch_memory_context = AllocSetContextCreate(CurrentMemoryContext,
                                                     "KmerBatchMemoryContext",
                                                     ALLOCSET_DEFAULT_SIZES);

    {
        int tuples_per_page = 10;
        int pages_needed = (kmersearch_highfreq_analysis_batch_size / tuples_per_page) * 1.5;
        int ring_size_kb = Max(pages_needed * 8, 10 * 1024);
        int total_ring_size_kb = ring_size_kb * max_parallel_maintenance_workers;
        int shared_buffers_kb = NBuffers * 8;

        if (shared_buffers_kb < total_ring_size_kb * 2)
        {
            elog(DEBUG2, "Worker using normal buffer access: shared_buffers=%d KB < ring_requirement=%d KB",
                 shared_buffers_kb, total_ring_size_kb * 2);
            ctx.strategy = NULL;
        }
        else
        {
            elog(DEBUG2, "Worker using ring buffer strategy: size=%d KB for batch_size=%d, shared_buffers=%d KB",
                 ring_size_kb, kmersearch_highfreq_analysis_batch_size, shared_buffers_kb);
            ctx.strategy = GetAccessStrategyWithSize(BAS_BULKREAD, ring_size_kb);
        }
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
                BlockNumber global_block;

                global_block = pg_atomic_fetch_add_u32(&shared_state->next_global_block, 1);

                if (global_block < shared_state->total_blocks_all_partitions)
                {
                    PartitionBlockMapping mapping = kmersearch_map_global_to_partition_block(
                        global_block, shared_state);

                    kmersearch_process_block_with_batch(mapping.local_block_number,
                                                       mapping.partition_oid,
                                                       shared_state, &ctx);
                    blocks_processed++;
                    got_block = true;
                }
            }
            else
            {
                LWLockAcquire(&shared_state->mutex.lock, LW_EXCLUSIVE);

                if (shared_state->next_block < shared_state->total_blocks)
                {
                    current_block = shared_state->next_block;
                    shared_state->next_block++;
                    LWLockRelease(&shared_state->mutex.lock);

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

        kmersearch_flush_batch_to_fht(&ctx);

        total_rows = pg_atomic_add_fetch_u64(&shared_state->total_rows_processed, ctx.batch_count);
        batch_num = pg_atomic_add_fetch_u64(&shared_state->total_batches_committed, 1);

        ereport(INFO,
                (errmsg("Batch %lu completed: %lu / %lu rows processed of column %s in table %s (final)",
                        (unsigned long)batch_num, (unsigned long)total_rows,
                        (unsigned long)shared_state->total_rows,
                        shared_state->column_name, shared_state->table_name)));
    }

    if (ctx.batch_hash)
    {
        hash_destroy(ctx.batch_hash);
        ctx.batch_hash = NULL;
    }

    /* Close file hash table */
    if (ctx.fht_ctx)
    {
        if (ctx.total_bits <= 16)
        {
            kmersearch_fht16_close((FileHashTable16Context *)ctx.fht_ctx);
        }
        else if (ctx.total_bits <= 32)
        {
            kmersearch_fht32_close((FileHashTable32Context *)ctx.fht_ctx);
        }
        else
        {
            kmersearch_fht64_close((FileHashTable64Context *)ctx.fht_ctx);
        }
        ctx.fht_ctx = NULL;
    }

    if (ctx.batch_memory_context)
    {
        MemoryContextDelete(ctx.batch_memory_context);
        ctx.batch_memory_context = NULL;
    }

    if (ctx.strategy)
    {
        FreeAccessStrategy(ctx.strategy);
        ctx.strategy = NULL;
    }

    /* Register temporary file path for parent process */
    kmersearch_register_worker_temp_file(shared_state, ctx.file_path, worker_id);
}/*
 * Map global block number to partition and local block
 */
static PartitionBlockMapping
kmersearch_map_global_to_partition_block(BlockNumber global_block, 
                                        KmerAnalysisSharedState *state)
{
    PartitionBlockMapping result;
    int i;
    
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
 * Parallel merge worker entry point (file hash table based)
 */
PGDLLEXPORT void
kmersearch_parallel_merge_worker(dsm_segment *seg, shm_toc *toc)
{
    char *source_file;
    char *target_file;
    int worker_id;
    char *file_paths_base;
    int total_bits;

    worker_id = ParallelWorkerNumber;

    file_paths_base = (char *) shm_toc_lookup(toc, 0, false);
    if (!file_paths_base)
    {
        elog(ERROR, "Failed to get file paths from shared memory");
        return;
    }

    /* Get total_bits from shared memory (stored at end of file paths) */
    total_bits = *((int *)(file_paths_base + (max_parallel_maintenance_workers * 2 * MAXPGPATH)));

    source_file = file_paths_base + (worker_id * 2 * MAXPGPATH);
    target_file = file_paths_base + (worker_id * 2 * MAXPGPATH + MAXPGPATH);

    /* Merge based on key size */
    if (total_bits <= 16)
    {
        kmersearch_fht16_merge(target_file, source_file);
    }
    else if (total_bits <= 32)
    {
        kmersearch_fht32_merge(target_file, source_file);
    }
    else
    {
        kmersearch_fht64_merge(target_file, source_file);
    }

    elog(DEBUG1, "Merge worker %d: completed merge of %s into %s", worker_id, target_file, source_file);
}

/*
 * Parallel aggregation of temporary file hash tables
 */
static void
kmersearch_aggregate_temp_files_parallel(char file_paths[][MAXPGPATH], int num_files,
                                        const char *temp_dir_path, int total_bits)
{
    int current_files = num_files;
    int stage = 0;
    Size shm_size;

    PG_TRY();
    {
        while (current_files > 1)
        {
            int num_parallel = current_files / 2;
            ParallelContext *pcxt;
            shm_toc *toc;
            char *shm_pointer;

            stage++;
            ereport(INFO,
                    (errmsg("Merge stage %d: %d files, %d parallel workers",
                            stage, current_files, num_parallel)));

            EnterParallelMode();
            pcxt = CreateParallelContext("pg_kmersearch", "kmersearch_parallel_merge_worker", num_parallel);
            pcxt->nworkers = num_parallel;

            /* Include space for total_bits integer at the end */
            shm_size = MAXPGPATH * 2 * max_parallel_maintenance_workers + sizeof(int);
            shm_toc_estimate_chunk(&pcxt->estimator, shm_size);
            shm_toc_estimate_keys(&pcxt->estimator, 1);

            InitializeParallelDSM(pcxt);

            toc = pcxt->toc;

            shm_pointer = shm_toc_allocate(toc, shm_size);

            for (int i = 0; i < num_parallel; i++)
            {
                char *source_path = shm_pointer + (i * 2 * MAXPGPATH);
                char *target_path = shm_pointer + (i * 2 * MAXPGPATH + MAXPGPATH);

                strlcpy(source_path, file_paths[i * 2], MAXPGPATH);
                strlcpy(target_path, file_paths[i * 2 + 1], MAXPGPATH);
            }

            /* Store total_bits at end of shared memory */
            *((int *)(shm_pointer + (max_parallel_maintenance_workers * 2 * MAXPGPATH))) = total_bits;

            shm_toc_insert(toc, 0, shm_pointer);

            LaunchParallelWorkers(pcxt);

            WaitForParallelWorkersToFinish(pcxt);

            {
                int new_idx = 0;
                for (int i = 0; i < current_files; i += 2)
                {
                    if (i + 1 < current_files)
                    {
                        strlcpy(file_paths[new_idx], file_paths[i], MAXPGPATH);
                        new_idx++;
                    }
                    else
                    {
                        strlcpy(file_paths[new_idx], file_paths[i], MAXPGPATH);
                        new_idx++;
                    }
                }
                current_files = new_idx;
            }

            DestroyParallelContext(pcxt);
            ExitParallelMode();
        }
    }
    PG_CATCH();
    {
        if (temp_dir_path && temp_dir_path[0] != '\0')
        {
            PathNameDeleteTemporaryDir(temp_dir_path);
            elog(DEBUG1, "Deleted temporary directory on error: %s", temp_dir_path);
        }
        PG_RE_THROW();
    }
    PG_END_TRY();

    ereport(INFO,
            (errmsg("Merge completed after %d stages", stage)));
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
            char expected_prefix[32];
            
            /* Build expected prefix with current database OID */
            snprintf(expected_prefix, sizeof(expected_prefix), "pg_kmersearch_db%u_", MyDatabaseId);
            
            /* Check if it's a kmersearch temp directory for this database */
            if (strncmp(de->d_name, expected_prefix, strlen(expected_prefix)) == 0)
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
                
                /* Try to remove the directory after deleting files */
                if (rmdir(subdir_path) < 0)
                {
                    ereport(INFO,
                            (errmsg("could not remove directory \"%s\": %m (this is normal if still in use)", subdir_path)));
                }
                else
                {
                    ereport(INFO,
                            (errmsg("removed temporary directory \"%s\"", subdir_path)));
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
 * Flush batch hash table to file hash table
 */
static void
kmersearch_flush_batch_to_fht(FileHashWorkerContext *ctx)
{
    HASH_SEQ_STATUS status;
    void *entry;
    MemoryContext old_context;

    if (!ctx->fht_ctx)
    {
        ereport(ERROR,
                (errmsg("File hash table not properly initialized for worker")));
    }

    old_context = MemoryContextSwitchTo(ctx->batch_memory_context);

    hash_seq_init(&status, ctx->batch_hash);

    while ((entry = hash_seq_search(&status)) != NULL)
    {
        if (ctx->total_bits <= 16)
        {
            TempKmerFreqEntry16 *e = (TempKmerFreqEntry16 *)entry;
            kmersearch_fht16_add((FileHashTable16Context *)ctx->fht_ctx,
                                e->uintkey, e->appearance_nrow);
        }
        else if (ctx->total_bits <= 32)
        {
            TempKmerFreqEntry32 *e = (TempKmerFreqEntry32 *)entry;
            kmersearch_fht32_add((FileHashTable32Context *)ctx->fht_ctx,
                                e->uintkey, e->appearance_nrow);
        }
        else
        {
            TempKmerFreqEntry64 *e = (TempKmerFreqEntry64 *)entry;
            kmersearch_fht64_add((FileHashTable64Context *)ctx->fht_ctx,
                                e->uintkey, e->appearance_nrow);
        }
    }

    /* Flush to disk */
    if (ctx->total_bits <= 16)
    {
        kmersearch_fht16_flush((FileHashTable16Context *)ctx->fht_ctx);
    }
    else if (ctx->total_bits <= 32)
    {
        kmersearch_fht32_flush((FileHashTable32Context *)ctx->fht_ctx);
    }
    else
    {
        kmersearch_fht64_flush((FileHashTable64Context *)ctx->fht_ctx);
    }

    MemoryContextSwitchTo(old_context);
}

/*
 * Process a block with batch aggregation
 */
static void
kmersearch_process_block_with_batch(BlockNumber block,
                                   Oid table_oid,
                                   KmerAnalysisSharedState *shared_state,
                                   FileHashWorkerContext *ctx)
{
    Relation rel;
    TupleDesc tupdesc;
    Buffer buffer;
    Page page;
    OffsetNumber maxoff;
    bool batch_completed = false;

    /* Initialize batch hash if needed (start of new batch) */
    if (ctx->batch_hash == NULL)
    {
        HASHCTL hashctl;
        MemoryContext old_context;
        
        /* Switch to batch memory context */
        old_context = MemoryContextSwitchTo(ctx->batch_memory_context);
        
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
        
        /* Use HASH_CONTEXT to ensure hash table is created in batch memory context */
        hashctl.hcxt = ctx->batch_memory_context;
        
        ctx->batch_hash = hash_create("KmerBatchHash",
                                     kmersearch_highfreq_analysis_hashtable_size,
                                     &hashctl, HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
        
        /* Switch back to original context */
        MemoryContextSwitchTo(old_context);
    }
    
    /* Open table and read block */
    rel = table_open(table_oid, AccessShareLock);
    /* Cache the tuple descriptor to avoid repeated calls in the loop */
    tupdesc = RelationGetDescr(rel);
    
    /* Use the pre-allocated buffer strategy from context */
    if (ctx->strategy)
    {
        buffer = ReadBufferExtended(rel, MAIN_FORKNUM, block,
                                   RBM_NORMAL_NO_LOG, ctx->strategy);
    }
    else
    {
        buffer = ReadBuffer(rel, block);
    }
    LockBuffer(buffer, BUFFER_LOCK_SHARE);
    page = BufferGetPage(buffer);
    maxoff = PageGetMaxOffsetNumber(page);
    
    /* Process each tuple in the page */
    for (OffsetNumber offnum = FirstOffsetNumber;
         offnum <= maxoff;
         offnum = OffsetNumberNext(offnum))
    {
        MemoryContext old_context;
        ItemId itemid = PageGetItemId(page, offnum);
        HeapTupleData tuple;
        bool isnull;
        Datum datum;
        void *kmer_array = NULL;
        int kmer_count;
        
        if (!ItemIdIsNormal(itemid))
            continue;

        memset(&tuple, 0, sizeof(HeapTupleData));
        tuple.t_data = (HeapTupleHeader) PageGetItem(page, itemid);
        tuple.t_len = ItemIdGetLength(itemid);
        tuple.t_tableOid = table_oid;
        ItemPointerSet(&tuple.t_self, block, offnum);

        /* MVCC visibility check: skip tuples not visible to current snapshot */
        if (!HeapTupleSatisfiesVisibility(&tuple, GetActiveSnapshot(), buffer))
            continue;

        /* Switch to batch memory context BEFORE getting data to ensure all allocations happen there */
        old_context = MemoryContextSwitchTo(ctx->batch_memory_context);
        
        /* Get sequence data - now allocated in batch memory context */
        datum = heap_getattr(&tuple, shared_state->column_attnum,
                           tupdesc, &isnull);
        if (isnull)
        {
            /* Switch back to original context before continuing */
            MemoryContextSwitchTo(old_context);
            continue;
        }
        
        /* Extract k-mers based on data type (type OIDs are cached at worker start) */
        if (ctx->column_type_oid == ctx->dna2_oid)
        {
            VarBit *seq = DatumGetVarBitP(datum);
            kmersearch_extract_uintkey_from_dna2(seq, &kmer_array, &kmer_count);
            /* Free detoasted datum if it was copied */
            if ((Pointer)seq != DatumGetPointer(datum))
                pfree(seq);
        }
        else if (ctx->column_type_oid == ctx->dna4_oid)
        {
            VarBit *seq = DatumGetVarBitP(datum);
            kmersearch_extract_uintkey_from_dna4(seq, &kmer_array, &kmer_count);
            /* Free detoasted datum if it was copied */
            if ((Pointer)seq != DatumGetPointer(datum))
                pfree(seq);
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
                    freq_entry->appearance_nrow = 0;
                }
                freq_entry->appearance_nrow++;
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
                    freq_entry->appearance_nrow = 0;
                }
                freq_entry->appearance_nrow++;
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
                    freq_entry->appearance_nrow = 0;
                }
                freq_entry->appearance_nrow++;
            }
        }
        
        if (kmer_array)
            pfree(kmer_array);
        
        /* Switch back to original context */
        MemoryContextSwitchTo(old_context);
        
        ctx->batch_count++;
        
        /* Flush batch if size limit reached */
        if (ctx->batch_count >= kmersearch_highfreq_analysis_batch_size)
        {
            uint64 total_rows;
            uint64 batch_num;
            
            /* Release current buffer before clearing buffers */
            if (BufferIsValid(buffer))
            {
                UnlockReleaseBuffer(buffer);
                buffer = InvalidBuffer;
            }
            
            kmersearch_flush_batch_to_fht(ctx);
            
            /* Update shared progress counters */
            total_rows = pg_atomic_add_fetch_u64(&shared_state->total_rows_processed, ctx->batch_count);
            batch_num = pg_atomic_add_fetch_u64(&shared_state->total_batches_committed, 1);
            
            /* Report cumulative progress (deterministic) */
            ereport(INFO,
                    (errmsg("Batch %lu completed: %lu / %lu rows processed of column %s in table %s",
                            (unsigned long)batch_num, (unsigned long)total_rows,
                            (unsigned long)shared_state->total_rows,
                            shared_state->column_name, shared_state->table_name)));
            
            /* Destroy hash table for completed batch */
            if (ctx->batch_hash)
            {
                hash_destroy(ctx->batch_hash);
                ctx->batch_hash = NULL;
            }
            
            /* Reset memory context instead of delete/recreate */
            if (ctx->batch_memory_context)
            {
                MemoryContextReset(ctx->batch_memory_context);
            }
            
            /* Show memory usage after batch memory context reset */
            MemoryContextStats(TopMemoryContext);
            
            /* Re-acquire buffer after batch completion if needed */
            if (!BufferIsValid(buffer))
            {
                buffer = ReadBuffer(rel, block);
                LockBuffer(buffer, BUFFER_LOCK_SHARE);
                page = BufferGetPage(buffer);
            }

            /* Log memory statistics and trim unused memory if using GLIBC */
#ifdef __GLIBC__
            {
                struct mallinfo2 mi = mallinfo2();
                elog(DEBUG1, "glibc memory before trim: arena=%zu, ordblks=%zu, hblkhd=%zu, uordblks=%zu, fordblks=%zu",
                     (size_t)mi.arena, (size_t)mi.ordblks, (size_t)mi.hblkhd, 
                     (size_t)mi.uordblks, (size_t)mi.fordblks);
                
                /* Trim unused memory back to OS after MemoryContextReset */
                malloc_trim(0);
                
                mi = mallinfo2();
                elog(DEBUG1, "glibc memory after trim: arena=%zu, ordblks=%zu, hblkhd=%zu, uordblks=%zu, fordblks=%zu",
                     (size_t)mi.arena, (size_t)mi.ordblks, (size_t)mi.hblkhd,
                     (size_t)mi.uordblks, (size_t)mi.fordblks);
            }
#endif
            
            /* Create new hash table for next batch (reusing memory context) */
            {
                HASHCTL hashctl;
                MemoryContext batch_old_context;
                
                /* Switch to batch memory context */
                batch_old_context = MemoryContextSwitchTo(ctx->batch_memory_context);
                
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
                
                /* Use HASH_CONTEXT to ensure hash table is created in batch memory context */
                hashctl.hcxt = ctx->batch_memory_context;
                
                ctx->batch_hash = hash_create("KmerBatchHash",
                                             kmersearch_highfreq_analysis_hashtable_size,
                                             &hashctl, HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
                
                /* Switch back to original context */
                MemoryContextSwitchTo(batch_old_context);
            }
            
            ctx->batch_count = 0;
            batch_completed = true;
        }
    }

    /* Release buffer */
    if (BufferIsValid(buffer))
        UnlockReleaseBuffer(buffer);

    table_close(rel, AccessShareLock);

    /* Clear main and TOAST table buffers after batch completion */
    if (batch_completed)
    {
        Relation main_rel;
        Oid toast_oid;
        ForkNumber fork = MAIN_FORKNUM;
        BlockNumber clear_from_block = 0;  /* Clear all blocks */

        main_rel = table_open(table_oid, AccessShareLock);
        toast_oid = main_rel->rd_rel->reltoastrelid;

        /* Clear main table buffers */
        PG_TRY();
        {
            /* Flush dirty buffers first */
            FlushRelationBuffers(main_rel);

            /* DropRelationBuffers clears buffers from clear_from_block onwards */
            /* By specifying 0, we clear all blocks of the relation */
            DropRelationBuffers(RelationGetSmgr(main_rel),
                              &fork,
                              1,
                              &clear_from_block);
            elog(DEBUG1, "Cleared main table buffers after batch completion at block %u", block);
        }
        PG_CATCH();
        {
            /* Ignore errors - likely some buffers are pinned by other workers */
            FlushErrorState();
            elog(DEBUG2, "Could not clear main table buffers after batch (some pinned by other workers)");
        }
        PG_END_TRY();

        /* Clear main table index buffers - independently from main table buffers */
        {
            List *main_index_list;
            ListCell *lc;

            main_index_list = RelationGetIndexList(main_rel);
            foreach(lc, main_index_list)
            {
                Oid main_index_oid = lfirst_oid(lc);
                Relation main_index_rel;

                PG_TRY();
                {
                    main_index_rel = index_open(main_index_oid, AccessShareLock);

                    elog(DEBUG1, "Attempting to clear main table index buffers for OID %u", main_index_oid);

                    /* Flush dirty buffers for the index */
                    FlushRelationBuffers(main_index_rel);

                    /* Drop all buffers for the main table index */
                    DropRelationBuffers(RelationGetSmgr(main_index_rel),
                                      &fork,
                                      1,
                                      &clear_from_block);
                    elog(DEBUG1, "Successfully cleared main table index buffers");

                    index_close(main_index_rel, AccessShareLock);
                }
                PG_CATCH();
                {
                    /* Ignore errors from buffers pinned by other workers */
                    FlushErrorState();
                    elog(DEBUG2, "Could not clear main table index buffers for OID %u (pinned by other workers)", main_index_oid);
                }
                PG_END_TRY();
            }
            list_free(main_index_list);
        }

        /* Clear TOAST table buffers if present */
        if (OidIsValid(toast_oid))
        {
            Relation toast_rel;

            toast_rel = table_open(toast_oid, AccessShareLock);

            /* Clear TOAST table buffers */
            PG_TRY();
            {
                elog(DEBUG1, "Attempting to clear TOAST table buffers for OID %u after batch completion", toast_oid);

                /* Flush any dirty buffers first */
                FlushRelationBuffers(toast_rel);

                /* Drop all buffers for the TOAST table using SMgr */
                DropRelationBuffers(RelationGetSmgr(toast_rel),
                                  &fork,
                                  1,
                                  &clear_from_block);
                elog(DEBUG1, "Successfully cleared TOAST table buffers");
            }
            PG_CATCH();
            {
                /* Ignore errors from buffers pinned by other workers */
                FlushErrorState();
                elog(DEBUG2, "Could not clear TOAST table buffers (pinned by other workers)");
            }
            PG_END_TRY();

            /* Clear TOAST index buffers - independently from TOAST table buffers */
            {
                List *toast_index_list;
                ListCell *lc;

                toast_index_list = RelationGetIndexList(toast_rel);
                foreach(lc, toast_index_list)
                {
                    Oid toast_index_oid = lfirst_oid(lc);
                    Relation toast_index_rel;

                    PG_TRY();
                    {
                        toast_index_rel = index_open(toast_index_oid, AccessShareLock);

                        elog(DEBUG1, "Attempting to clear TOAST index buffers for OID %u", toast_index_oid);

                        /* Flush dirty buffers for the index */
                        FlushRelationBuffers(toast_index_rel);

                        /* Drop all buffers for the TOAST index */
                        DropRelationBuffers(RelationGetSmgr(toast_index_rel),
                                          &fork,
                                          1,
                                          &clear_from_block);
                        elog(DEBUG1, "Successfully cleared TOAST index buffers");

                        index_close(toast_index_rel, AccessShareLock);
                    }
                    PG_CATCH();
                    {
                        /* Ignore errors from buffers pinned by other workers */
                        FlushErrorState();
                        elog(DEBUG2, "Could not clear TOAST index buffers for OID %u (pinned by other workers)", toast_index_oid);
                    }
                    PG_END_TRY();
                }
                list_free(toast_index_list);
            }

            table_close(toast_rel, AccessShareLock);
        }

        table_close(main_rel, AccessShareLock);
    }
}

