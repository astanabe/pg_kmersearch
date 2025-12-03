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

        -- Check if this is a kmersearch GIN index by querying system catalogs
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

### Component 4: Event Trigger for REINDEX

**File**: `pg_kmersearch--1.0.sql`

REINDEX does not fire CREATE INDEX event trigger. Create a separate event trigger for REINDEX:

```sql
CREATE FUNCTION kmersearch_on_reindex()
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
        WHERE command_tag IN ('REINDEX')
    LOOP
        -- For REINDEX, objid may be table OID or index OID depending on command
        -- We need to find all kmersearch indexes and update their settings
        FOR idx_oid, tbl_oid, col_name, opclass_name IN
            SELECT
                i.indexrelid,
                i.indrelid,
                a.attname,
                oc.opcname
            FROM pg_index i
            JOIN pg_class ic ON ic.oid = i.indexrelid
            JOIN pg_am am ON am.oid = ic.relam
            JOIN pg_attribute a ON a.attrelid = i.indrelid
                AND a.attnum = ANY(i.indkey)
            JOIN pg_opclass oc ON oc.oid = i.indclass[0]
            WHERE am.amname = 'gin'
              AND oc.opcname LIKE 'kmersearch_%'
              AND (i.indexrelid = obj.objid OR i.indrelid = obj.objid)
        LOOP
            -- Update settings for this index
            UPDATE kmersearch_index_info SET
                kmer_size = current_setting('kmersearch.kmer_size')::integer,
                occur_bitlen = current_setting('kmersearch.occur_bitlen')::integer,
                max_appearance_rate = current_setting('kmersearch.max_appearance_rate')::real,
                max_appearance_nrow = current_setting('kmersearch.max_appearance_nrow')::integer,
                preclude_highfreq_kmer = current_setting('kmersearch.preclude_highfreq_kmer')::boolean,
                created_at = now()
            WHERE index_oid = idx_oid;

            -- If not exists, insert new record
            IF NOT FOUND THEN
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
                    0, 0
                );
            END IF;
        END LOOP;
    END LOOP;
END;
$$;

CREATE EVENT TRIGGER kmersearch_reindex_trigger
ON ddl_command_end
WHEN TAG IN ('REINDEX')
EXECUTE FUNCTION kmersearch_on_reindex();
```

### Component 5: Planner Hook Implementation

**File**: `kmersearch_planner.c` (new file)

**Required includes:**
```c
#include "postgres.h"
#include "fmgr.h"
#include "optimizer/pathnode.h"
#include "optimizer/paths.h"
#include "nodes/pathnodes.h"
#include "nodes/pg_list.h"
#include "catalog/pg_am_d.h"      /* For GIN_AM_OID */
#include "catalog/pg_class.h"
#include "catalog/pg_opclass.h"
#include "utils/syscache.h"
#include "utils/rel.h"
#include "executor/spi.h"
#include "utils/builtins.h"
#include "kmersearch.h"
#include <math.h>
```

**Note**: `get_opclass_name()` is a static function in `ruleutils.c` and cannot be called directly. Use SPI query to pg_opclass instead.

```c
/* Hook storage */
static set_rel_pathlist_hook_type prev_set_rel_pathlist_hook = NULL;

/* Structure to cache index settings lookup results */
typedef struct KmersearchIndexSettings
{
    Oid         index_oid;
    bool        is_kmersearch_index;
    bool        settings_found;
    int         kmer_size;
    int         occur_bitlen;
    float       max_appearance_rate;
    int         max_appearance_nrow;
    bool        preclude_highfreq_kmer;
} KmersearchIndexSettings;

/* Forward declarations */
static bool kmersearch_is_kmersearch_gin_index(Oid index_oid);
static bool kmersearch_get_index_settings(Oid index_oid, KmersearchIndexSettings *settings);
static bool kmersearch_check_settings_match(KmersearchIndexSettings *settings);
static void kmersearch_set_rel_pathlist(PlannerInfo *root, RelOptInfo *rel,
                                         Index rti, RangeTblEntry *rte);
static void kmersearch_adjust_bitmap_path_cost(Path *bitmapqual);

/*
 * Check if an index is a kmersearch GIN index using SPI
 *
 * Note: get_opclass_name() is static in ruleutils.c, so we use SPI query instead.
 */
static bool
kmersearch_is_kmersearch_gin_index(Oid index_oid)
{
    int ret;
    bool result = false;
    StringInfoData query;

    if (SPI_connect() != SPI_OK_CONNECT)
        return false;

    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT 1 FROM pg_index i "
        "JOIN pg_class ic ON ic.oid = i.indexrelid "
        "JOIN pg_am am ON am.oid = ic.relam "
        "JOIN pg_opclass oc ON oc.oid = i.indclass[0] "
        "WHERE i.indexrelid = %u "
        "  AND am.amname = 'gin' "
        "  AND oc.opcname LIKE 'kmersearch_%%'",
        index_oid);

    ret = SPI_execute(query.data, true, 1);

    if (ret == SPI_OK_SELECT && SPI_processed > 0)
        result = true;

    pfree(query.data);
    SPI_finish();

    return result;
}

/*
 * Get index settings from kmersearch_index_info table
 */
static bool
kmersearch_get_index_settings(Oid index_oid, KmersearchIndexSettings *settings)
{
    int ret;
    bool found = false;
    StringInfoData query;

    settings->index_oid = index_oid;
    settings->is_kmersearch_index = false;
    settings->settings_found = false;

    /* First check if this is a kmersearch index */
    if (!kmersearch_is_kmersearch_gin_index(index_oid))
        return false;

    settings->is_kmersearch_index = true;

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

        settings->kmer_size = DatumGetInt32(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull));
        settings->occur_bitlen = DatumGetInt32(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull));
        settings->max_appearance_rate = DatumGetFloat4(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull));
        settings->max_appearance_nrow = DatumGetInt32(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 4, &isnull));
        settings->preclude_highfreq_kmer = DatumGetBool(
            SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 5, &isnull));

        settings->settings_found = true;
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
kmersearch_check_settings_match(KmersearchIndexSettings *settings)
{
    if (!settings->settings_found)
        return false;  /* No settings found, consider it a mismatch */

    /* Compare each setting */
    if (settings->kmer_size != kmersearch_kmer_size)
        return false;

    if (settings->occur_bitlen != kmersearch_occur_bitlen)
        return false;

    /* For floating point, use small epsilon for comparison */
    if (fabs(settings->max_appearance_rate - kmersearch_max_appearance_rate) > 0.0001)
        return false;

    if (settings->max_appearance_nrow != kmersearch_max_appearance_nrow)
        return false;

    if (settings->preclude_highfreq_kmer != kmersearch_preclude_highfreq_kmer)
        return false;

    return true;
}

/*
 * Recursively adjust cost for bitmap paths containing mismatched indexes
 *
 * GIN indexes always use bitmap scans, and complex queries may have
 * BitmapOr or BitmapAnd paths combining multiple indexes.
 */
static void
kmersearch_adjust_bitmap_path_cost(Path *bitmapqual)
{
    if (bitmapqual == NULL)
        return;

    if (IsA(bitmapqual, IndexPath))
    {
        IndexPath *ipath = (IndexPath *) bitmapqual;
        KmersearchIndexSettings settings;

        if (kmersearch_get_index_settings(ipath->indexinfo->indexoid, &settings))
        {
            if (!kmersearch_check_settings_match(&settings))
            {
                /* Set prohibitively high cost */
                bitmapqual->startup_cost = 1.0e10;
                bitmapqual->total_cost = 1.0e10;

                ereport(DEBUG1,
                    (errmsg("kmersearch: disabling index %u due to settings mismatch",
                            ipath->indexinfo->indexoid),
                     errdetail("Index: kmer_size=%d, occur_bitlen=%d, max_appearance_rate=%.4f, "
                               "max_appearance_nrow=%d, preclude_highfreq_kmer=%s; "
                               "Current: kmer_size=%d, occur_bitlen=%d, max_appearance_rate=%.4f, "
                               "max_appearance_nrow=%d, preclude_highfreq_kmer=%s",
                               settings.kmer_size, settings.occur_bitlen,
                               settings.max_appearance_rate, settings.max_appearance_nrow,
                               settings.preclude_highfreq_kmer ? "true" : "false",
                               kmersearch_kmer_size, kmersearch_occur_bitlen,
                               kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow,
                               kmersearch_preclude_highfreq_kmer ? "true" : "false")));
            }
        }
    }
    else if (IsA(bitmapqual, BitmapOrPath))
    {
        BitmapOrPath *orpath = (BitmapOrPath *) bitmapqual;
        ListCell *lc;

        foreach(lc, orpath->bitmapquals)
        {
            kmersearch_adjust_bitmap_path_cost((Path *) lfirst(lc));
        }

        /* If any child has high cost, propagate it */
        foreach(lc, orpath->bitmapquals)
        {
            Path *child = (Path *) lfirst(lc);
            if (child->total_cost >= 1.0e10)
            {
                bitmapqual->startup_cost = 1.0e10;
                bitmapqual->total_cost = 1.0e10;
                break;
            }
        }
    }
    else if (IsA(bitmapqual, BitmapAndPath))
    {
        BitmapAndPath *andpath = (BitmapAndPath *) bitmapqual;
        ListCell *lc;

        foreach(lc, andpath->bitmapquals)
        {
            kmersearch_adjust_bitmap_path_cost((Path *) lfirst(lc));
        }

        /* If any child has high cost, propagate it */
        foreach(lc, andpath->bitmapquals)
        {
            Path *child = (Path *) lfirst(lc);
            if (child->total_cost >= 1.0e10)
            {
                bitmapqual->startup_cost = 1.0e10;
                bitmapqual->total_cost = 1.0e10;
                break;
            }
        }
    }
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

        /* GIN indexes use BitmapHeapPath */
        if (IsA(path, BitmapHeapPath))
        {
            BitmapHeapPath *bhpath = (BitmapHeapPath *) path;

            /* Recursively check and adjust bitmap paths */
            kmersearch_adjust_bitmap_path_cost(bhpath->bitmapqual);

            /* Propagate high cost to parent path */
            if (bhpath->bitmapqual->total_cost >= 1.0e10)
            {
                path->startup_cost = 1.0e10;
                path->total_cost = 1.0e10;
            }
        }
        /* Also handle direct IndexPath (though GIN typically uses bitmap) */
        else if (IsA(path, IndexPath))
        {
            IndexPath *ipath = (IndexPath *) path;
            KmersearchIndexSettings settings;

            if (kmersearch_get_index_settings(ipath->indexinfo->indexoid, &settings))
            {
                if (!kmersearch_check_settings_match(&settings))
                {
                    path->startup_cost = 1.0e10;
                    path->total_cost = 1.0e10;
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
 * Cleanup planner hook (called from _PG_fini if needed)
 */
void
kmersearch_planner_fini(void)
{
    set_rel_pathlist_hook = prev_set_rel_pathlist_hook;
}
```

### Component 6: Integration with _PG_init

**File**: `kmersearch.c`

Add hook initialization to `_PG_init()`:

```c
/* Add at end of _PG_init() */
void
_PG_init(void)
{
    /* ... existing code ... */

    /* Initialize planner hook for index settings validation */
    kmersearch_planner_init();
}
```

### Component 7: Header File Updates

**File**: `kmersearch.h`

Add declarations:

```c
/* Planner hook functions */
extern void kmersearch_planner_init(void);
extern void kmersearch_planner_fini(void);
```

### Component 8: Makefile Updates

**File**: `Makefile`

Add new source file to OBJS:

```makefile
OBJS = kmersearch.o kmersearch_gin.o kmersearch_datatype.o kmersearch_kmer.o \
       kmersearch_cache.o kmersearch_freq.o kmersearch_partition.o \
       kmersearch_util.o kmersearch_planner.o
```

## Implementation Order

### Phase 1: Metadata Storage
1. Update `kmersearch_index_info` table schema (add `preclude_highfreq_kmer` column)
2. Create event trigger for `CREATE INDEX`
3. Create event trigger for `DROP INDEX`
4. Create event trigger for `REINDEX`
5. Test that settings are correctly captured on index creation, drop, and reindex

### Phase 2: Planner Hook
1. Create `kmersearch_planner.c` with hook implementation
2. Add `kmersearch_is_kmersearch_gin_index()` function using SPI
3. Add `kmersearch_get_index_settings()` function using SPI
4. Add `kmersearch_check_settings_match()` function
5. Add `kmersearch_adjust_bitmap_path_cost()` function for recursive bitmap path handling
6. Implement `kmersearch_set_rel_pathlist()` hook function
7. Add `kmersearch_planner_init()` and `kmersearch_planner_fini()`
8. Integrate with `_PG_init()`
9. Update Makefile

### Phase 3: Testing
1. Create test case: single index with matching settings (should use index)
2. Create test case: single index with mismatched settings (should use seq scan)
3. Create test case: multiple indexes, one matching (should select matching one)
4. Create test case: multiple indexes, none matching (should use seq scan)
5. Create test case: REINDEX after settings change (should update metadata)
6. Verify EXPLAIN output shows expected behavior
7. Test with DEBUG1 log level to verify detailed messages

### Phase 4: Edge Cases and Robustness
1. Handle case where `kmersearch_index_info` table doesn't exist (SPI error handling)
2. Handle case where index entry is missing from table (treat as mismatch)
3. Handle concurrent index creation (ON CONFLICT clause handles this)
4. Performance optimization: consider caching settings lookup results per query planning session

## Testing Plan

### Test 1: Basic Settings Match
```sql
-- Settings: kmer_size=7, occur_bitlen=2 -> 16 bits -> int2
SET kmersearch.kmer_size = 7;
SET kmersearch.occur_bitlen = 2;
SET kmersearch.max_appearance_rate = 0.5;
SET kmersearch.max_appearance_nrow = 0;
SET kmersearch.preclude_highfreq_kmer = false;

CREATE INDEX idx1 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

-- Should use idx1
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';
```

### Test 2: Settings Mismatch
```sql
-- Create index with kmer_size=7
SET kmersearch.kmer_size = 7;
SET kmersearch.occur_bitlen = 2;
CREATE INDEX idx1 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

-- Change to different kmer_size
SET kmersearch.kmer_size = 6;
-- Should NOT use idx1 (settings mismatch), should use Seq Scan
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';
```

### Test 3: Multiple Indexes with Different Settings
```sql
-- Create index with kmer_size=7 (7*2+2=16 bits -> int2)
SET kmersearch.kmer_size = 7;
SET kmersearch.occur_bitlen = 2;
CREATE INDEX idx_k7 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

-- Create index with kmer_size=12 (12*2+4=28 bits -> int4)
SET kmersearch.kmer_size = 12;
SET kmersearch.occur_bitlen = 4;
CREATE INDEX idx_k12 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int4);

-- Should use idx_k7
SET kmersearch.kmer_size = 7;
SET kmersearch.occur_bitlen = 2;
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';

-- Should use idx_k12
SET kmersearch.kmer_size = 12;
SET kmersearch.occur_bitlen = 4;
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';
```

### Test 4: REINDEX Behavior
```sql
-- Create index with initial settings
SET kmersearch.kmer_size = 7;
SET kmersearch.occur_bitlen = 2;
CREATE INDEX idx1 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

-- Check stored settings
SELECT * FROM kmersearch_index_info WHERE index_oid = 'idx1'::regclass;

-- Change settings and REINDEX
SET kmersearch.kmer_size = 6;
SET kmersearch.occur_bitlen = 3;
REINDEX INDEX idx1;

-- Check that settings were updated
SELECT * FROM kmersearch_index_info WHERE index_oid = 'idx1'::regclass;
```

### Test 5: Verify DEBUG Messages
```sql
SET client_min_messages = DEBUG1;
SET kmersearch.kmer_size = 7;
SET kmersearch.occur_bitlen = 2;
CREATE INDEX idx1 ON seqs USING gin(seq kmersearch_dna2_gin_ops_int2);

SET kmersearch.kmer_size = 8;  -- Mismatch
-- Should see DEBUG1 message about disabling index
EXPLAIN SELECT * FROM seqs WHERE seq =% 'ACGTACGT';
```

## Notes and Considerations

### Performance Impact
- The planner hook adds overhead to every query planning
- SPI queries in the hook should be minimized
- Two SPI queries per kmersearch index: one to check if it's a kmersearch index, one to get settings
- Consider caching settings lookup results per query planning session in future optimization

### Compatibility
- Event triggers require PostgreSQL 9.3+
- `set_rel_pathlist_hook` is available in all supported PostgreSQL versions
- Tested with PostgreSQL 16

### Limitations
- Cannot validate settings for indexes created before this feature is implemented
- Workaround: Provide a function to manually register existing index settings
- SPI cannot be used during parallel workers (but planner runs in leader only)

### GIN Index Specifics
- GIN indexes always use Bitmap Heap Scan, not regular Index Scan
- Complex queries may combine multiple indexes with BitmapOr/BitmapAnd
- The implementation handles these cases recursively

### Future Enhancements
- Add `kmersearch_register_index_settings(index_name, ...)` function for manual registration
- Add `kmersearch_index_settings_view` to show index settings vs current GUC settings
- Cache settings lookup results in backend-local memory for better performance
- Consider using index reloptions instead of separate table (more integrated but complex)
