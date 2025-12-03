# GIN Index Settings Validation and Selection Implementation Plan

## Overview

This document describes the implementation plan for validating GIN index settings at query time and automatically selecting the appropriate index when multiple GIN indexes exist on the same DNA2/DNA4 column with different configurations.

## Problem Statement

Currently, pg_kmersearch does not validate whether the GUC settings at query time match the settings used when the GIN index was built. This can lead to incorrect search results without any warning to the user.

## Goals

1. Save GIN index build settings to a metadata table when creating an index
2. Validate settings at query time using a planner hook
3. Set prohibitively high cost for indexes with mismatched settings
4. Allow the planner to automatically select the index with matching settings

## Implementation Components

### Component 1: Metadata Table Schema Update

**File**: `pg_kmersearch--1.0.sql`

Update `kmersearch_index_info` table to include `preclude_highfreq_kmer`:

```sql
CREATE TABLE kmersearch_index_info (
    index_oid oid PRIMARY KEY,
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    kmer_size integer NOT NULL,
    occur_bitlen integer NOT NULL,
    max_appearance_rate real NOT NULL,
    max_appearance_nrow integer NOT NULL,
    preclude_highfreq_kmer boolean NOT NULL,  -- NEW COLUMN
    total_nrow bigint NOT NULL,
    highfreq_kmer_count integer NOT NULL,
    created_at timestamp with time zone DEFAULT now()
);
```

### Component 2: Event Trigger for Index Creation

**File**: `pg_kmersearch--1.0.sql`

Create an event trigger function to capture GIN index creation and save settings:

```sql
CREATE FUNCTION kmersearch_on_index_create()
RETURNS event_trigger
LANGUAGE plpgsql
AS $$
DECLARE
    obj record;
    idx_oid oid;
    tbl_oid oid;
    col_name name;
    opclass_name name;
BEGIN
    FOR obj IN SELECT * FROM pg_event_trigger_ddl_commands()
        WHERE object_type = 'index' AND command_tag = 'CREATE INDEX'
    LOOP
        idx_oid := obj.objid;

        -- Check if this is a kmersearch GIN index
        SELECT
            i.indrelid,
            a.attname,
            oc.opcname
        INTO tbl_oid, col_name, opclass_name
        FROM pg_index i
        JOIN pg_class ic ON ic.oid = i.indexrelid
        JOIN pg_am am ON am.oid = ic.relam
        JOIN pg_attribute a ON a.attrelid = i.indrelid
            AND a.attnum = ANY(i.indkey)
        JOIN pg_opclass oc ON oc.oid = i.indclass[0]
        WHERE i.indexrelid = idx_oid
          AND am.amname = 'gin'
          AND oc.opcname LIKE 'kmersearch_%';

        IF FOUND THEN
            -- Insert current GUC settings into kmersearch_index_info
            INSERT INTO kmersearch_index_info (
                index_oid, table_oid, column_name,
                kmer_size, occur_bitlen,
                max_appearance_rate, max_appearance_nrow,
                preclude_highfreq_kmer,
                total_nrow, highfreq_kmer_count
            ) VALUES (
                idx_oid, tbl_oid, col_name,
                current_setting('kmersearch.kmer_size')::integer,
                current_setting('kmersearch.occur_bitlen')::integer,
                current_setting('kmersearch.max_appearance_rate')::real,
                current_setting('kmersearch.max_appearance_nrow')::integer,
                current_setting('kmersearch.preclude_highfreq_kmer')::boolean,
                0, 0  -- total_nrow and highfreq_kmer_count will be updated later
            )
            ON CONFLICT (index_oid) DO UPDATE SET
                kmer_size = EXCLUDED.kmer_size,
                occur_bitlen = EXCLUDED.occur_bitlen,
                max_appearance_rate = EXCLUDED.max_appearance_rate,
                max_appearance_nrow = EXCLUDED.max_appearance_nrow,
                preclude_highfreq_kmer = EXCLUDED.preclude_highfreq_kmer,
                created_at = now();
        END IF;
    END LOOP;
END;
$$;

CREATE EVENT TRIGGER kmersearch_index_create_trigger
ON ddl_command_end
WHEN TAG IN ('CREATE INDEX')
EXECUTE FUNCTION kmersearch_on_index_create();
```

### Component 3: Event Trigger for Index Drop

**File**: `pg_kmersearch--1.0.sql`

Create an event trigger to clean up metadata when index is dropped:

```sql
CREATE FUNCTION kmersearch_on_index_drop()
RETURNS event_trigger
LANGUAGE plpgsql
AS $$
DECLARE
    obj record;
BEGIN
    FOR obj IN SELECT * FROM pg_event_trigger_dropped_objects()
        WHERE object_type = 'index'
    LOOP
        DELETE FROM kmersearch_index_info WHERE index_oid = obj.objid;
        DELETE FROM kmersearch_gin_index_meta WHERE index_oid = obj.objid;
    END LOOP;
END;
$$;

CREATE EVENT TRIGGER kmersearch_index_drop_trigger
ON sql_drop
WHEN TAG IN ('DROP INDEX')
EXECUTE FUNCTION kmersearch_on_index_drop();
```

### Component 4: Planner Hook Implementation

**File**: `kmersearch_planner.c` (new file)

```c
#include "postgres.h"
#include "optimizer/pathnode.h"
#include "optimizer/paths.h"
#include "optimizer/planner.h"
#include "nodes/pathnodes.h"
#include "catalog/pg_class.h"
#include "catalog/pg_opclass.h"
#include "utils/syscache.h"
#include "utils/lsyscache.h"
#include "executor/spi.h"
#include "kmersearch.h"

/* Hook storage */
static set_rel_pathlist_hook_type prev_set_rel_pathlist_hook = NULL;

/* Structure to cache index settings lookup results */
typedef struct IndexSettingsCache
{
    Oid         index_oid;
    bool        is_kmersearch_index;
    bool        settings_found;
    int         kmer_size;
    int         occur_bitlen;
    float       max_appearance_rate;
    int         max_appearance_nrow;
    bool        preclude_highfreq_kmer;
} IndexSettingsCache;

/* Forward declarations */
static bool is_kmersearch_gin_index(Oid index_oid);
static bool get_index_settings(Oid index_oid, IndexSettingsCache *cache);
static bool check_settings_match(IndexSettingsCache *cache);
static void kmersearch_set_rel_pathlist(PlannerInfo *root, RelOptInfo *rel,
                                         Index rti, RangeTblEntry *rte);

/*
 * Check if an index is a kmersearch GIN index
 */
static bool
is_kmersearch_gin_index(Oid index_oid)
{
    HeapTuple   indexTuple;
    HeapTuple   classTuple;
    Form_pg_index indexForm;
    Form_pg_class classForm;
    Oid         amoid;
    Oid         opclass_oid;
    char       *opclass_name;
    bool        result = false;

    /* Get index tuple */
    indexTuple = SearchSysCache1(INDEXRELID, ObjectIdGetDatum(index_oid));
    if (!HeapTupleIsValid(indexTuple))
        return false;

    indexForm = (Form_pg_index) GETSTRUCT(indexTuple);

    /* Get class tuple to check access method */
    classTuple = SearchSysCache1(RELOID, ObjectIdGetDatum(index_oid));
    if (!HeapTupleIsValid(classTuple))
    {
        ReleaseSysCache(indexTuple);
        return false;
    }

    classForm = (Form_pg_class) GETSTRUCT(classTuple);
    amoid = classForm->relam;

    /* Check if it's a GIN index */
    if (amoid != GIN_AM_OID)  /* Need to get GIN_AM_OID properly */
    {
        ReleaseSysCache(classTuple);
        ReleaseSysCache(indexTuple);
        return false;
    }

    /* Check operator class name */
    opclass_oid = indexForm->indclass.values[0];
    opclass_name = get_opclass_name(opclass_oid, InvalidOid, false);

    if (opclass_name && strncmp(opclass_name, "kmersearch_", 11) == 0)
        result = true;

    ReleaseSysCache(classTuple);
    ReleaseSysCache(indexTuple);

    return result;
}

/*
 * Get index settings from kmersearch_index_info table
 */
static bool
get_index_settings(Oid index_oid, IndexSettingsCache *cache)
{
    int ret;
    bool found = false;
    StringInfoData query;

    cache->index_oid = index_oid;
    cache->is_kmersearch_index = false;
    cache->settings_found = false;

    /* First check if this is a kmersearch index */
    if (!is_kmersearch_gin_index(index_oid))
        return false;

    cache->is_kmersearch_index = true;

    /* Query settings from kmersearch_index_info */
    if (SPI_connect() != SPI_OK_CONNECT)
        return false;

    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT kmer_size, occur_bitlen, max_appearance_rate, "
        "max_appearance_nrow, preclude_highfreq_kmer "
        "FROM kmersearch_index_info WHERE index_oid = %u",
        index_oid);

    ret = SPI_execute(query.data, true, 1);

    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        bool isnull;

        cache->kmer_size = DatumGetInt32(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull));
        cache->occur_bitlen = DatumGetInt32(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull));
        cache->max_appearance_rate = DatumGetFloat4(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull));
        cache->max_appearance_nrow = DatumGetInt32(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 4, &isnull));
        cache->preclude_highfreq_kmer = DatumGetBool(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 5, &isnull));

        cache->settings_found = true;
        found = true;
    }

    pfree(query.data);
    SPI_finish();

    return found;
}

/*
 * Check if index settings match current GUC settings
 */
static bool
check_settings_match(IndexSettingsCache *cache)
{
    if (!cache->settings_found)
        return false;  /* No settings found, consider it a mismatch */

    /* Compare each setting */
    if (cache->kmer_size != kmersearch_kmer_size)
        return false;

    if (cache->occur_bitlen != kmersearch_occur_bitlen)
        return false;

    /* For floating point, use small epsilon for comparison */
    if (fabs(cache->max_appearance_rate - kmersearch_max_appearance_rate) > 0.0001)
        return false;

    if (cache->max_appearance_nrow != kmersearch_max_appearance_nrow)
        return false;

    if (cache->preclude_highfreq_kmer != kmersearch_preclude_highfreq_kmer)
        return false;

    return true;
}

/*
 * Planner hook to adjust costs for mismatched indexes
 */
static void
kmersearch_set_rel_pathlist(PlannerInfo *root, RelOptInfo *rel,
                             Index rti, RangeTblEntry *rte)
{
    ListCell *lc;

    /* Call previous hook first */
    if (prev_set_rel_pathlist_hook)
        prev_set_rel_pathlist_hook(root, rel, rti, rte);

    /* Only process base relations */
    if (rel->reloptkind != RELOPT_BASEREL)
        return;

    /* Iterate through all paths */
    foreach(lc, rel->pathlist)
    {
        Path *path = (Path *) lfirst(lc);

        /* Check for IndexPath (includes BitmapIndexScan used by GIN) */
        if (IsA(path, IndexPath))
        {
            IndexPath *ipath = (IndexPath *) path;
            IndexSettingsCache cache;

            /* Check if this is a kmersearch index with mismatched settings */
            if (get_index_settings(ipath->indexinfo->indexoid, &cache))
            {
                if (!check_settings_match(&cache))
                {
                    /* Set prohibitively high cost */
                    path->startup_cost = 1.0e10;
                    path->total_cost = 1.0e10;

                    ereport(DEBUG1,
                        (errmsg("kmersearch: disabling index %u due to settings mismatch",
                                ipath->indexinfo->indexoid),
                         errdetail("Index kmer_size=%d, current=%d; "
                                   "Index occur_bitlen=%d, current=%d",
                                   cache.kmer_size, kmersearch_kmer_size,
                                   cache.occur_bitlen, kmersearch_occur_bitlen)));
                }
            }
        }
        /* Also check BitmapHeapPath which may use GIN indexes */
        else if (IsA(path, BitmapHeapPath))
        {
            BitmapHeapPath *bhpath = (BitmapHeapPath *) path;

            /* Recursively check bitmap paths */
            if (IsA(bhpath->bitmapqual, IndexPath))
            {
                IndexPath *ipath = (IndexPath *) bhpath->bitmapqual;
                IndexSettingsCache cache;

                if (get_index_settings(ipath->indexinfo->indexoid, &cache))
                {
                    if (!check_settings_match(&cache))
                    {
                        path->startup_cost = 1.0e10;
                        path->total_cost = 1.0e10;
                    }
                }
            }
        }
    }
}

/*
 * Initialize planner hook
 */
void
kmersearch_planner_init(void)
{
    prev_set_rel_pathlist_hook = set_rel_pathlist_hook;
    set_rel_pathlist_hook = kmersearch_set_rel_pathlist;
}

/*
 * Cleanup planner hook
 */
void
kmersearch_planner_fini(void)
{
    set_rel_pathlist_hook = prev_set_rel_pathlist_hook;
}
```

### Component 5: Integration with _PG_init

**File**: `kmersearch.c`

Add hook initialization to `_PG_init()`:

```c
/* Add to includes */
extern void kmersearch_planner_init(void);

/* Add at end of _PG_init() */
void
_PG_init(void)
{
    /* ... existing code ... */

    /* Initialize planner hook for index settings validation */
    kmersearch_planner_init();
}
```

### Component 6: Header File Updates

**File**: `kmersearch.h`

Add declarations:

```c
/* Planner hook functions */
extern void kmersearch_planner_init(void);
extern void kmersearch_planner_fini(void);
```

### Component 7: Makefile Updates

**File**: `Makefile`

Add new source file:

```makefile
OBJS = kmersearch.o kmersearch_gin.o kmersearch_freq.o kmersearch_cache.o \
       kmersearch_simd.o kmersearch_partition.o kmersearch_planner.o
```

## Implementation Order

### Phase 1: Metadata Storage
1. Update `kmersearch_index_info` table schema (add `preclude_highfreq_kmer` column)
2. Create event trigger for `CREATE INDEX`
3. Create event trigger for `DROP INDEX`
4. Test that settings are correctly captured on index creation

### Phase 2: Planner Hook
1. Create `kmersearch_planner.c` with hook implementation
2. Add `is_kmersearch_gin_index()` function
3. Add `get_index_settings()` function using SPI
4. Add `check_settings_match()` function
5. Implement `kmersearch_set_rel_pathlist()` hook function
6. Integrate with `_PG_init()`

### Phase 3: Testing
1. Create test case: single index with matching settings (should use index)
2. Create test case: single index with mismatched settings (should use seq scan)
3. Create test case: multiple indexes, one matching (should select matching one)
4. Create test case: multiple indexes, none matching (should use seq scan)
5. Create test case: REINDEX after settings change
6. Verify EXPLAIN output shows expected behavior

### Phase 4: Edge Cases and Robustness
1. Handle case where `kmersearch_index_info` table doesn't exist
2. Handle case where index entry is missing from table
3. Handle concurrent index creation
4. Performance optimization: cache settings lookup results per query

## Testing Plan

### Test 1: Basic Settings Match
```sql
SET kmersearch.kmer_size = 8;
SET kmersearch.occur_bitlen = 4;
CREATE INDEX idx1 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

-- Should use idx1
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';
```

### Test 2: Settings Mismatch
```sql
SET kmersearch.kmer_size = 8;
CREATE INDEX idx1 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

SET kmersearch.kmer_size = 16;
-- Should NOT use idx1 (settings mismatch)
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';
```

### Test 3: Multiple Indexes
```sql
SET kmersearch.kmer_size = 8;
SET kmersearch.occur_bitlen = 4;
CREATE INDEX idx_k8 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

SET kmersearch.kmer_size = 16;
SET kmersearch.occur_bitlen = 8;
CREATE INDEX idx_k16 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int4);

-- Should use idx_k8
SET kmersearch.kmer_size = 8;
SET kmersearch.occur_bitlen = 4;
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';

-- Should use idx_k16
SET kmersearch.kmer_size = 16;
SET kmersearch.occur_bitlen = 8;
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';
```

## Notes and Considerations

### Performance Impact
- The planner hook adds overhead to every query planning
- SPI queries in the hook should be minimized
- Consider caching index settings in backend-local memory

### Compatibility
- Event triggers require PostgreSQL 9.3+
- `set_rel_pathlist_hook` is available in all supported PostgreSQL versions
- Need to handle different PostgreSQL version APIs

### Limitations
- Cannot validate settings for indexes created before this feature is implemented
- REINDEX does not fire CREATE INDEX event trigger (need separate handling)

### Future Enhancements
- Add a function to manually register existing index settings
- Add a view to show index settings vs current GUC settings
- Consider using index reloptions instead of separate table (more complex)
