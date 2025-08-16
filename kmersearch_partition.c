/*-------------------------------------------------------------------------
 *
 * kmersearch_partition.c
 *    Partitioning support functions for pg_kmersearch
 *
 * IDENTIFICATION
 *    pg_kmersearch/kmersearch_partition.c
 *
 *-------------------------------------------------------------------------
 */
#include "postgres.h"
#include "fmgr.h"
#include "access/genam.h"
#include "access/htup_details.h"
#include "access/table.h"
#include "catalog/indexing.h"
#include "catalog/namespace.h"
#include "catalog/pg_attribute.h"
#include "catalog/pg_class.h"
#include "catalog/pg_inherits.h"
#include "catalog/pg_type.h"
#include "commands/tablespace.h"
#include "executor/spi.h"
#include "funcapi.h"
#include "miscadmin.h"
#include "nodes/makefuncs.h"
#include "storage/ipc.h"
#include "storage/lmgr.h"
#include "utils/builtins.h"
#include "utils/lsyscache.h"
#include "utils/memutils.h"
#include "utils/rel.h"
#include "utils/snapmgr.h"
#include "utils/syscache.h"
#include "utils/timestamp.h"

#include "kmersearch.h"

/* Function declarations */
PG_FUNCTION_INFO_V1(kmersearch_partition_table);

/* Helper function prototypes */
static void validate_table_for_partitioning(Oid table_oid, char **dna_column_name, Oid *dna_column_type);
static int calculate_partition_batch_size(Oid table_oid);
static void create_partition_table(const char *temp_table_name, const char *table_name, 
                                   const char *dna_column_name, int partition_count, const char *tablespace_name);
static void migrate_data_in_batches(const char *table_name, const char *temp_table_name, Oid table_oid);
static void replace_table_with_partition(const char *table_name, const char *temp_table_name);


/*
 * kmersearch_partition_table
 *
 * Convert a non-partitioned table to a hash partitioned table based on DNA2/DNA4 column
 */
Datum
kmersearch_partition_table(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_PP(0);
    int32 partition_count = PG_GETARG_INT32(1);
    text *tablespace_name_text = PG_ARGISNULL(2) ? NULL : PG_GETARG_TEXT_PP(2);
    char *table_name;
    char *tablespace_name = NULL;
    Oid table_oid;
    Oid dna_column_type;
    char *dna_column_name = NULL;
    char temp_table_name[NAMEDATALEN];
    int ret;
    LOCKMODE lockmode = AccessExclusiveLock;
    
    /* Parameter validation */
    if (partition_count < 1)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("partition_count must be at least 1")));
    
    /* Get table name and OID */
    table_name = text_to_cstring(table_name_text);
    table_oid = RangeVarGetRelid(makeRangeVar(NULL, table_name, -1), lockmode, false);
    
    /* Get tablespace name if provided */
    if (tablespace_name_text != NULL)
        tablespace_name = text_to_cstring(tablespace_name_text);
    
    /* Validate table for partitioning */
    validate_table_for_partitioning(table_oid, &dna_column_name, &dna_column_type);
    
    /* Generate temporary table name */
    snprintf(temp_table_name, sizeof(temp_table_name), 
             "%s_part_%ld", table_name, (long)(GetCurrentTimestamp() / 1000));
    
    /* Connect to SPI */
    if ((ret = SPI_connect()) != SPI_OK_CONNECT)
        elog(ERROR, "SPI_connect failed: %s", SPI_result_code_string(ret));
    
    PG_TRY();
    {
        /* Create partition table */
        create_partition_table(temp_table_name, table_name, dna_column_name, partition_count, tablespace_name);
        
        /* Migrate data in batches */
        migrate_data_in_batches(table_name, temp_table_name, table_oid);
        
        /* Replace table with partition */
        replace_table_with_partition(table_name, temp_table_name);
        
        SPI_finish();
    }
    PG_CATCH();
    {
        /* Clean up temporary tables on error */
        StringInfoData query;
        int i;
        
        initStringInfo(&query);
        
        /* Try to drop partitions */
        for (i = 0; i < partition_count; i++)
        {
            resetStringInfo(&query);
            appendStringInfo(&query, "DROP TABLE IF EXISTS %s_%d", table_name, i);
            SPI_execute(query.data, false, 0);
        }
        
        /* Try to drop temporary parent table */
        resetStringInfo(&query);
        appendStringInfo(&query, "DROP TABLE IF EXISTS %s", temp_table_name);
        SPI_execute(query.data, false, 0);
        
        pfree(query.data);
        
        SPI_finish();
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    if (dna_column_name)
        pfree(dna_column_name);
    if (tablespace_name)
        pfree(tablespace_name);
    
    PG_RETURN_VOID();
}

/*
 * validate_table_for_partitioning
 *
 * Validate that the table meets requirements for partitioning
 */
static void
validate_table_for_partitioning(Oid table_oid, char **dna_column_name, Oid *dna_column_type)
{
    Relation rel;
    TupleDesc tupdesc;
    int natts;
    int dna_column_count = 0;
    int i;
    Oid dna2_type_oid;
    Oid dna4_type_oid;
    
    /* Open relation */
    rel = table_open(table_oid, AccessShareLock);
    
    /* Check if already partitioned */
    if (rel->rd_rel->relkind == RELKIND_PARTITIONED_TABLE)
    {
        table_close(rel, AccessShareLock);
        ereport(ERROR,
                (errcode(ERRCODE_WRONG_OBJECT_TYPE),
                 errmsg("table \"%s\" is already a partitioned table",
                        RelationGetRelationName(rel))));
    }
    
    /* Get DNA type OIDs */
    dna2_type_oid = get_dna2_type_oid();
    dna4_type_oid = get_dna4_type_oid();
    
    /* Check for DNA2/DNA4 columns */
    tupdesc = RelationGetDescr(rel);
    natts = tupdesc->natts;
    
    for (i = 0; i < natts; i++)
    {
        Form_pg_attribute attr = TupleDescAttr(tupdesc, i);
        
        if (attr->attisdropped)
            continue;
            
        if (attr->atttypid == dna2_type_oid || attr->atttypid == dna4_type_oid)
        {
            dna_column_count++;
            if (dna_column_count == 1)
            {
                *dna_column_name = pstrdup(NameStr(attr->attname));
                *dna_column_type = attr->atttypid;
            }
        }
    }
    
    table_close(rel, AccessShareLock);
    
    /* Validate exactly one DNA column */
    if (dna_column_count == 0)
        ereport(ERROR,
                (errcode(ERRCODE_WRONG_OBJECT_TYPE),
                 errmsg("table must have at least one DNA2 or DNA4 column")));
    else if (dna_column_count > 1)
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("table has %d DNA2/DNA4 columns, but exactly one is required",
                        dna_column_count)));
}

/*
 * calculate_partition_batch_size
 *
 * Calculate batch size based on maintenance_work_mem
 */
static int
calculate_partition_batch_size(Oid table_oid)
{
    int64 maintenance_work_mem_bytes = maintenance_work_mem * 1024L;
    int64 avg_row_size;
    int batch_size;
    HeapTuple tuple;
    Form_pg_class classForm;
    float4 avg_width;
    
    /* Get average row width from pg_class */
    tuple = SearchSysCache1(RELOID, ObjectIdGetDatum(table_oid));
    if (!HeapTupleIsValid(tuple))
        elog(ERROR, "cache lookup failed for relation %u", table_oid);
        
    classForm = (Form_pg_class) GETSTRUCT(tuple);
    avg_width = classForm->reltuples > 0 ? 
                (classForm->relpages * BLCKSZ) / classForm->reltuples : 1024;
    
    ReleaseSysCache(tuple);
    
    avg_row_size = (int64)avg_width;
    if (avg_row_size <= 0)
        avg_row_size = 1024;  /* Default 1KB */
    
    /* Use 1/4 of maintenance_work_mem for batch */
    batch_size = (maintenance_work_mem_bytes / 4) / avg_row_size;
    
    /* Limit batch size */
    if (batch_size < 1000)
        batch_size = 1000;
    else if (batch_size > 100000)
        batch_size = 100000;
        
    elog(NOTICE, "Using batch size %d based on maintenance_work_mem=%dMB",
         batch_size, maintenance_work_mem / 1024);
         
    return batch_size;
}

/*
 * create_partition_table
 *
 * Create the partitioned table structure
 */
static void
create_partition_table(const char *temp_table_name, const char *table_name, 
                       const char *dna_column_name, int partition_count, const char *tablespace_name)
{
    StringInfoData query;
    int ret;
    int i;
    Oid table_oid;
    Oid tablespace_oid;
    char *target_tablespace = NULL;
    
    initStringInfo(&query);
    
    /* Determine target tablespace */
    if (tablespace_name != NULL && strlen(tablespace_name) > 0)
    {
        /* Use explicitly provided tablespace */
        target_tablespace = (char *)tablespace_name;
    }
    else
    {
        /* Get tablespace of original table */
        table_oid = RangeVarGetRelid(makeRangeVar(NULL, (char *)table_name, -1), NoLock, false);
        tablespace_oid = get_rel_tablespace(table_oid);
        if (OidIsValid(tablespace_oid))
            target_tablespace = get_tablespace_name(tablespace_oid);
    }
    
    /* Create parent partitioned table */
    /* Note: We use INCLUDING DEFAULTS INCLUDING GENERATED INCLUDING IDENTITY INCLUDING STATISTICS
     * instead of INCLUDING ALL to avoid constraint issues with partitioning.
     * Constraints will be handled separately if needed. */
    appendStringInfo(&query,
        "CREATE TABLE %s (LIKE %s INCLUDING DEFAULTS INCLUDING GENERATED "
        "INCLUDING IDENTITY INCLUDING STATISTICS) PARTITION BY HASH (%s)",
        temp_table_name, table_name, dna_column_name);
        
    if (target_tablespace)
        appendStringInfo(&query, " TABLESPACE %s", target_tablespace);
        
    ret = SPI_execute(query.data, false, 0);
    if (ret != SPI_OK_UTILITY)
        elog(ERROR, "CREATE TABLE failed: %s", SPI_result_code_string(ret));
        
    /* Create partitions */
    for (i = 0; i < partition_count; i++)
    {
        resetStringInfo(&query);
        appendStringInfo(&query,
            "CREATE TABLE %s_%d PARTITION OF %s "
            "FOR VALUES WITH (modulus %d, remainder %d)",
            table_name, i, temp_table_name, partition_count, i);
            
        if (target_tablespace)
            appendStringInfo(&query, " TABLESPACE %s", target_tablespace);
            
        ret = SPI_execute(query.data, false, 0);
        if (ret != SPI_OK_UTILITY)
            elog(ERROR, "CREATE TABLE partition failed: %s", SPI_result_code_string(ret));
    }
    
    pfree(query.data);
    /* Only free if we allocated it ourselves */
    if (tablespace_name == NULL && target_tablespace != NULL)
        pfree(target_tablespace);
}

/*
 * perform_simple_migration
 *
 * Helper function to perform simple data migration without batching
 */
static uint64
perform_simple_migration(const char *table_name, const char *temp_table_name)
{
    StringInfoData query;
    int ret;
    uint64 rows_processed;
    
    initStringInfo(&query);
    
    appendStringInfo(&query,
        "INSERT INTO %s SELECT * FROM %s",
        temp_table_name, table_name);
        
    ret = SPI_execute(query.data, false, 0);
    if (ret != SPI_OK_INSERT)
        elog(ERROR, "Data migration failed: %s", SPI_result_code_string(ret));
        
    rows_processed = SPI_processed;
    pfree(query.data);
    
    return rows_processed;
}

/*
 * migrate_data_in_batches
 *
 * Migrate data from original table to partitioned table in batches
 */
static void
migrate_data_in_batches(const char *table_name, const char *temp_table_name, Oid table_oid)
{
    StringInfoData query;
    int ret;
    uint64 total_rows_processed = 0;
    int batch_size;
    SPITupleTable *tuptable;
    uint64 total_rows = 0;
    uint64 rows_migrated = 0;
    Oid seq_oid = InvalidOid;
    char *seq_column = NULL;
    bool has_suitable_index = false;
    
    initStringInfo(&query);
    
    /* Calculate batch size based on maintenance_work_mem */
    batch_size = calculate_partition_batch_size(table_oid);
    
    /* Check if table has a suitable index for batch processing */
    /* Look for primary key or unique index to enable ordered batching */
    appendStringInfo(&query,
        "SELECT a.attname "
        "FROM pg_index i "
        "JOIN pg_attribute a ON a.attrelid = i.indrelid AND a.attnum = ANY(i.indkey) "
        "WHERE i.indrelid = %u AND i.indisprimary "
        "ORDER BY array_position(i.indkey, a.attnum) "
        "LIMIT 1", table_oid);
    
    ret = SPI_execute(query.data, true, 1);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        bool isnull;
        Datum column_datum;
        
        tuptable = SPI_tuptable;
        column_datum = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isnull);
        if (!isnull)
        {
            seq_column = TextDatumGetCString(column_datum);
            has_suitable_index = true;
        }
    }
    
    if (!has_suitable_index)
    {
        /* Look for any unique index as fallback */
        resetStringInfo(&query);
        appendStringInfo(&query,
            "SELECT a.attname "
            "FROM pg_index i "
            "JOIN pg_attribute a ON a.attrelid = i.indrelid AND a.attnum = ANY(i.indkey) "
            "WHERE i.indrelid = %u AND i.indisunique AND NOT i.indisprimary "
            "ORDER BY i.indexrelid, array_position(i.indkey, a.attnum) "
            "LIMIT 1", table_oid);
        
        ret = SPI_execute(query.data, true, 1);
        if (ret == SPI_OK_SELECT && SPI_processed > 0)
        {
            bool isnull;
            Datum column_datum;
            
            tuptable = SPI_tuptable;
            column_datum = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                seq_column = TextDatumGetCString(column_datum);
                has_suitable_index = true;
            }
        }
    }
    
    /* Get total row count */
    resetStringInfo(&query);
    appendStringInfo(&query, "SELECT COUNT(*) FROM %s", table_name);
    ret = SPI_execute(query.data, true, 1);
    if (ret == SPI_OK_SELECT && SPI_processed == 1)
    {
        bool isnull;
        Datum count_datum;
        
        tuptable = SPI_tuptable;
        count_datum = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isnull);
        if (!isnull)
            total_rows = DatumGetInt64(count_datum);
    }
    
    if (has_suitable_index && total_rows > batch_size)
    {
        /* Batch processing with ordering by indexed column */
        Datum last_value = (Datum)0;
        bool first_batch = true;
        Oid column_type = InvalidOid;
        
        /* Get column type for proper comparison */
        resetStringInfo(&query);
        appendStringInfo(&query,
            "SELECT atttypid FROM pg_attribute "
            "WHERE attrelid = %u AND attname = '%s'",
            table_oid, seq_column);
        
        ret = SPI_execute(query.data, true, 1);
        if (ret == SPI_OK_SELECT && SPI_processed == 1)
        {
            bool isnull;
            Datum type_datum;
            
            tuptable = SPI_tuptable;
            type_datum = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isnull);
            if (!isnull)
                column_type = DatumGetObjectId(type_datum);
        }
        
        /* If we couldn't get column type, fall back to simple migration */
        if (!OidIsValid(column_type))
        {
            elog(WARNING, "Could not determine column type for %s, falling back to simple migration", seq_column);
            elog(NOTICE, "Attempting type resolution through syscache for column %s", seq_column);
            
            /* Try alternative method to get column type */
            AttrNumber attnum = get_attnum(table_oid, seq_column);
            if (attnum != InvalidAttrNumber)
            {
                column_type = get_atttype(table_oid, attnum);
            }
            
            /* If still invalid, use simple migration */
            if (!OidIsValid(column_type))
            {
                elog(NOTICE, "Using simple migration due to type resolution failure");
                total_rows_processed = perform_simple_migration(table_name, temp_table_name);
                
                /* Truncate original table */
                resetStringInfo(&query);
                appendStringInfo(&query, "TRUNCATE TABLE %s", table_name);
                ret = SPI_execute(query.data, false, 0);
                if (ret != SPI_OK_UTILITY)
                    elog(ERROR, "TRUNCATE TABLE failed: %s", SPI_result_code_string(ret));
                    
                /* Clean up and return early */
                pfree(query.data);
                if (seq_column)
                    pfree(seq_column);
                    
                elog(NOTICE, "Data migration completed: " UINT64_FORMAT " rows migrated",
                     total_rows_processed);
                return;
            }
        }
        
        elog(NOTICE, "Starting batch migration of " UINT64_FORMAT " rows (batch size: %d)",
             total_rows, batch_size);
        
        while (rows_migrated < total_rows)
        {
            resetStringInfo(&query);
            
            if (first_batch)
            {
                /* First batch - no WHERE clause needed */
                appendStringInfo(&query,
                    "INSERT INTO %s SELECT * FROM %s "
                    "ORDER BY %s LIMIT %d",
                    temp_table_name, table_name, seq_column, batch_size);
                first_batch = false;
            }
            else
            {
                /* Subsequent batches - use WHERE clause with last value */
                appendStringInfo(&query,
                    "INSERT INTO %s SELECT * FROM %s "
                    "WHERE %s > ",
                    temp_table_name, table_name, seq_column);
                
                /* Append last value based on type */
                if (column_type == INT4OID)
                    appendStringInfo(&query, "%d", DatumGetInt32(last_value));
                else if (column_type == INT8OID)
                    appendStringInfo(&query, INT64_FORMAT, DatumGetInt64(last_value));
                else
                    appendStringInfo(&query, "'%s'", TextDatumGetCString(last_value));
                
                appendStringInfo(&query, " ORDER BY %s LIMIT %d", seq_column, batch_size);
            }
            
            ret = SPI_execute(query.data, false, 0);
            if (ret != SPI_OK_INSERT)
                elog(ERROR, "Batch data migration failed: %s", SPI_result_code_string(ret));
            
            rows_migrated += SPI_processed;
            
            /* If we processed fewer rows than batch size, we're done */
            if (SPI_processed < batch_size)
                break;
            
            /* Get the last value for next batch */
            resetStringInfo(&query);
            appendStringInfo(&query,
                "SELECT MAX(%s) FROM ("
                "SELECT %s FROM %s ORDER BY %s LIMIT %d OFFSET " UINT64_FORMAT
                ") AS batch",
                seq_column, seq_column, table_name, seq_column, batch_size, rows_migrated - batch_size);
            
            ret = SPI_execute(query.data, true, 1);
            if (ret == SPI_OK_SELECT && SPI_processed == 1)
            {
                bool isnull;
                tuptable = SPI_tuptable;
                last_value = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isnull);
                if (isnull)
                    break;
            }
            
            /* Progress reporting */
            if (rows_migrated % (batch_size * 10) == 0)
            {
                elog(NOTICE, "Migration progress: " UINT64_FORMAT "/" UINT64_FORMAT " rows (%.1f%%)",
                     rows_migrated, total_rows, (double)rows_migrated * 100.0 / total_rows);
            }
        }
        
        total_rows_processed = rows_migrated;
    }
    else
    {
        /* Fall back to simple approach for small tables or tables without suitable index */
        elog(NOTICE, "Using simple migration for table with %lld rows", (long long)total_rows);
        total_rows_processed = perform_simple_migration(table_name, temp_table_name);
    }
    
    /* Now truncate the original table */
    resetStringInfo(&query);
    appendStringInfo(&query, "TRUNCATE TABLE %s", table_name);
    
    ret = SPI_execute(query.data, false, 0);
    if (ret != SPI_OK_UTILITY)
        elog(ERROR, "TRUNCATE TABLE failed: %s", SPI_result_code_string(ret));
    
    elog(NOTICE, "Data migration completed: " UINT64_FORMAT " rows migrated",
         total_rows_processed);
         
    pfree(query.data);
    if (seq_column)
        pfree(seq_column);
}

/*
 * replace_table_with_partition
 *
 * Replace original table with partitioned table
 */
static void
replace_table_with_partition(const char *table_name, const char *temp_table_name)
{
    StringInfoData query;
    int ret;
    SPITupleTable *tuptable;
    
    initStringInfo(&query);
    
    /* Verify original table is empty */
    appendStringInfo(&query, "SELECT COUNT(*) FROM %s", table_name);
    ret = SPI_execute(query.data, true, 0);
    if (ret != SPI_OK_SELECT || SPI_processed != 1)
        elog(ERROR, "Failed to verify table is empty");
        
    tuptable = SPI_tuptable;
    if (tuptable && tuptable->vals && tuptable->vals[0])
    {
        bool isnull;
        Datum count_datum = SPI_getbinval(tuptable->vals[0], tuptable->tupdesc, 1, &isnull);
        if (!isnull && DatumGetInt64(count_datum) > 0)
            elog(ERROR, "Original table is not empty after migration");
    }
    else
        elog(ERROR, "Failed to get count result");
        
    /* Drop original table with CASCADE to handle sequence dependencies */
    /* Log information about what will be cascaded */
    elog(NOTICE, "Dropping original table %s with CASCADE to handle dependent objects", table_name);
    elog(NOTICE, "Any dependent objects (sequences, views, etc.) will be dropped automatically");
    
    resetStringInfo(&query);
    
    /* First, identify dependent objects for user awareness */
    appendStringInfo(&query,
        "SELECT DISTINCT "
        "  dc.relname AS dependent_object, "
        "  CASE dc.relkind "
        "    WHEN 'S' THEN 'sequence' "
        "    WHEN 'v' THEN 'view' "
        "    WHEN 'm' THEN 'materialized view' "
        "    ELSE 'other' "
        "  END AS object_type "
        "FROM pg_depend d "
        "JOIN pg_class c ON d.refobjid = c.oid "
        "JOIN pg_class dc ON d.objid = dc.oid "
        "WHERE c.relname = %s AND d.deptype != 'i'", 
        quote_literal_cstr(table_name));
    
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        int i;
        SPITupleTable *tuptable = SPI_tuptable;
        
        for (i = 0; i < SPI_processed; i++)
        {
            bool name_isnull, type_isnull;
            char *obj_name_str;
            char *obj_type_str;
            
            obj_name_str = SPI_getvalue(tuptable->vals[i], tuptable->tupdesc, 1);
            obj_type_str = SPI_getvalue(tuptable->vals[i], tuptable->tupdesc, 2);
            
            if (obj_name_str && obj_type_str)
            {
                elog(NOTICE, "CASCADE will drop %s: %s", 
                     obj_type_str, obj_name_str);
            }
            
            if (obj_name_str)
                pfree(obj_name_str);
            if (obj_type_str)
                pfree(obj_type_str);
        }
    }
    
    /* Now perform the actual DROP with CASCADE */
    resetStringInfo(&query);
    appendStringInfo(&query, "DROP TABLE %s CASCADE", table_name);
    ret = SPI_execute(query.data, false, 0);
    if (ret != SPI_OK_UTILITY)
        elog(ERROR, "DROP TABLE failed: %s", SPI_result_code_string(ret));
        
    elog(NOTICE, "Original table %s dropped successfully", table_name);
        
    /* Rename partitioned table */
    resetStringInfo(&query);
    appendStringInfo(&query, "ALTER TABLE %s RENAME TO %s", temp_table_name, table_name);
    ret = SPI_execute(query.data, false, 0);
    if (ret != SPI_OK_UTILITY)
        elog(ERROR, "ALTER TABLE RENAME failed: %s", SPI_result_code_string(ret));
        
    pfree(query.data);
}








