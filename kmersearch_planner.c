/*-------------------------------------------------------------------------
 *
 * kmersearch_planner.c
 *		Planner hook for validating GIN index settings and automatic
 *		index selection based on current GUC settings.
 *
 * This module implements a planner hook that checks whether the GUC
 * settings at query time match the settings used when GIN indexes
 * were built. If they don't match, the hook sets a prohibitively high
 * cost for those indexes, effectively disabling them and allowing
 * the planner to select an index with matching settings or fall back
 * to a sequential scan.
 *
 * Copyright (c) 2024, pg_kmersearch contributors
 *
 *-------------------------------------------------------------------------
 */
#include "kmersearch.h"

#include "optimizer/pathnode.h"
#include "optimizer/paths.h"
#include "optimizer/cost.h"
#include "nodes/pathnodes.h"
#include "catalog/pg_am.h"
#include "catalog/pg_am_d.h"
#include "catalog/pg_index.h"
#include "catalog/pg_opclass.h"
#include "utils/fmgroids.h"
#include "access/table.h"
#include "access/genam.h"

/* Hook storage */
static set_rel_pathlist_hook_type prev_set_rel_pathlist_hook = NULL;

/* Recursion guard */
static bool in_kmersearch_hook = false;

/* Structure to cache index settings lookup results */
typedef struct KmersearchIndexSettings
{
	Oid			index_oid;
	bool		is_kmersearch_index;
	bool		settings_found;
	int			kmer_size;
	int			occur_bitlen;
	float		max_appearance_rate;
	int			max_appearance_nrow;
	bool		preclude_highfreq_kmer;
} KmersearchIndexSettings;

/* Forward declarations */
static bool kmersearch_is_kmersearch_gin_index_direct(Oid index_oid);
static bool kmersearch_get_index_settings_direct(Oid index_oid, KmersearchIndexSettings *settings);
static bool kmersearch_check_settings_match(KmersearchIndexSettings *settings);
static void kmersearch_set_rel_pathlist(PlannerInfo *root, RelOptInfo *rel,
										Index rti, RangeTblEntry *rte);
static void kmersearch_adjust_bitmap_path_cost(Path *bitmapqual);
static IndexOptInfo *kmersearch_find_matching_index(RelOptInfo *rel, IndexPath *existing_ipath);
static void kmersearch_add_matching_index_path(PlannerInfo *root, RelOptInfo *rel,
											   BitmapHeapPath *existing_bhpath);

/*
 * Check if an index is a kmersearch GIN index using direct catalog access
 */
static bool
kmersearch_is_kmersearch_gin_index_direct(Oid index_oid)
{
	HeapTuple	indexTuple;
	HeapTuple	classTuple;
	HeapTuple	opclassTuple;
	Form_pg_index indexForm;
	Form_pg_class classForm;
	Form_pg_opclass opclassForm;
	Oid			amoid;
	Oid			opclassoid;
	char	   *opcname;
	bool		result = false;
	Datum		indclassDatum;
	bool		isnull;
	oidvector  *indclass;

	indexTuple = SearchSysCache1(INDEXRELID, ObjectIdGetDatum(index_oid));
	if (!HeapTupleIsValid(indexTuple))
		return false;

	indexForm = (Form_pg_index) GETSTRUCT(indexTuple);

	classTuple = SearchSysCache1(RELOID, ObjectIdGetDatum(index_oid));
	if (!HeapTupleIsValid(classTuple))
	{
		ReleaseSysCache(indexTuple);
		return false;
	}

	classForm = (Form_pg_class) GETSTRUCT(classTuple);
	amoid = classForm->relam;

	if (amoid != GIN_AM_OID)
	{
		ReleaseSysCache(classTuple);
		ReleaseSysCache(indexTuple);
		return false;
	}

	if (indexForm->indnatts < 1)
	{
		ReleaseSysCache(classTuple);
		ReleaseSysCache(indexTuple);
		return false;
	}

	indclassDatum = SysCacheGetAttr(INDEXRELID, indexTuple,
									Anum_pg_index_indclass, &isnull);
	if (isnull)
	{
		ReleaseSysCache(classTuple);
		ReleaseSysCache(indexTuple);
		return false;
	}

	indclass = (oidvector *) DatumGetPointer(indclassDatum);
	opclassoid = indclass->values[0];

	opclassTuple = SearchSysCache1(CLAOID, ObjectIdGetDatum(opclassoid));
	if (!HeapTupleIsValid(opclassTuple))
	{
		ReleaseSysCache(classTuple);
		ReleaseSysCache(indexTuple);
		return false;
	}

	opclassForm = (Form_pg_opclass) GETSTRUCT(opclassTuple);
	opcname = NameStr(opclassForm->opcname);

	if (strncmp(opcname, "kmersearch_", 11) == 0)
		result = true;

	ReleaseSysCache(opclassTuple);
	ReleaseSysCache(classTuple);
	ReleaseSysCache(indexTuple);

	return result;
}

/*
 * Get index settings from kmersearch_index_info table using direct heap scan
 * This avoids SPI which can cause issues in planner hooks
 */
static bool
kmersearch_get_index_settings_direct(Oid index_oid, KmersearchIndexSettings *settings)
{
	Oid			info_table_oid;
	Relation	info_rel;
	TableScanDesc scan;
	HeapTuple	tuple;
	TupleDesc	tupdesc;
	bool		found = false;
	Snapshot	snapshot;

	settings->index_oid = index_oid;
	settings->is_kmersearch_index = false;
	settings->settings_found = false;

	if (!kmersearch_is_kmersearch_gin_index_direct(index_oid))
		return false;

	settings->is_kmersearch_index = true;

	{
		Oid public_namespace = get_namespace_oid("public", true);
		if (!OidIsValid(public_namespace))
			return false;
		info_table_oid = get_relname_relid("kmersearch_index_info", public_namespace);
	}
	if (!OidIsValid(info_table_oid))
		return false;

	info_rel = table_open(info_table_oid, AccessShareLock);
	tupdesc = RelationGetDescr(info_rel);

	snapshot = GetActiveSnapshot();
	if (snapshot == NULL)
		snapshot = GetTransactionSnapshot();

	scan = table_beginscan(info_rel, snapshot, 0, NULL);

	while ((tuple = heap_getnext(scan, ForwardScanDirection)) != NULL)
	{
		bool	isnull;
		Datum	datum;
		Oid		stored_oid;

		datum = heap_getattr(tuple, 1, tupdesc, &isnull);
		if (isnull)
			continue;

		stored_oid = DatumGetObjectId(datum);
		if (stored_oid != index_oid)
			continue;

		datum = heap_getattr(tuple, 4, tupdesc, &isnull);
		settings->kmer_size = isnull ? 0 : DatumGetInt32(datum);

		datum = heap_getattr(tuple, 5, tupdesc, &isnull);
		settings->occur_bitlen = isnull ? 0 : DatumGetInt32(datum);

		datum = heap_getattr(tuple, 8, tupdesc, &isnull);
		settings->max_appearance_rate = isnull ? 0.0 : DatumGetFloat4(datum);

		datum = heap_getattr(tuple, 9, tupdesc, &isnull);
		settings->max_appearance_nrow = isnull ? 0 : DatumGetInt32(datum);

		datum = heap_getattr(tuple, 10, tupdesc, &isnull);
		settings->preclude_highfreq_kmer = isnull ? false : DatumGetBool(datum);

		settings->settings_found = true;
		found = true;
		break;
	}

	table_endscan(scan);
	table_close(info_rel, AccessShareLock);

	return found;
}

/*
 * Check if index settings match current GUC settings
 */
static bool
kmersearch_check_settings_match(KmersearchIndexSettings *settings)
{
	if (!settings->settings_found)
		return false;

	if (settings->kmer_size != kmersearch_kmer_size)
		return false;

	if (settings->occur_bitlen != kmersearch_occur_bitlen)
		return false;

	if (fabs(settings->max_appearance_rate - kmersearch_max_appearance_rate) > 0.0001)
		return false;

	if (settings->max_appearance_nrow != kmersearch_max_appearance_nrow)
		return false;

	if (settings->preclude_highfreq_kmer != kmersearch_preclude_highfreq_kmer)
		return false;

	return true;
}

/*
 * Find a kmersearch index with matching settings from rel->indexlist
 * Returns NULL if no matching index is found
 */
static IndexOptInfo *
kmersearch_find_matching_index(RelOptInfo *rel, IndexPath *existing_ipath)
{
	ListCell   *lc;
	Oid			existing_index_oid = existing_ipath->indexinfo->indexoid;

	foreach(lc, rel->indexlist)
	{
		IndexOptInfo *index = (IndexOptInfo *) lfirst(lc);
		KmersearchIndexSettings settings;

		if (index->indexoid == existing_index_oid)
			continue;

		if (!kmersearch_is_kmersearch_gin_index_direct(index->indexoid))
			continue;

		if (kmersearch_get_index_settings_direct(index->indexoid, &settings))
		{
			if (kmersearch_check_settings_match(&settings))
				return index;
		}
	}

	return NULL;
}

/*
 * Create and add a new BitmapHeapPath for a matching index
 */
static void
kmersearch_add_matching_index_path(PlannerInfo *root, RelOptInfo *rel,
								   BitmapHeapPath *existing_bhpath)
{
	IndexPath	   *existing_ipath;
	IndexOptInfo   *matching_index;
	IndexPath	   *new_ipath;
	BitmapHeapPath *new_bhpath;
	List		   *new_indexclauses;
	ListCell	   *lc;

	if (!IsA(existing_bhpath->bitmapqual, IndexPath))
		return;

	existing_ipath = (IndexPath *) existing_bhpath->bitmapqual;

	matching_index = kmersearch_find_matching_index(rel, existing_ipath);
	if (matching_index == NULL)
		return;

	new_indexclauses = NIL;
	foreach(lc, existing_ipath->indexclauses)
	{
		IndexClause *old_ic = (IndexClause *) lfirst(lc);
		IndexClause *new_ic = makeNode(IndexClause);

		new_ic->rinfo = old_ic->rinfo;
		new_ic->indexquals = list_copy(old_ic->indexquals);
		new_ic->lossy = old_ic->lossy;
		new_ic->indexcol = old_ic->indexcol;
		new_ic->indexcols = list_copy(old_ic->indexcols);

		new_indexclauses = lappend(new_indexclauses, new_ic);
	}

	new_ipath = create_index_path(root,
								  matching_index,
								  new_indexclauses,
								  NIL,
								  NIL,
								  NIL,
								  ForwardScanDirection,
								  false,
								  NULL,
								  1.0,
								  false);

	new_bhpath = create_bitmap_heap_path(root,
										 rel,
										 (Path *) new_ipath,
										 NULL,
										 1.0,
										 0);

	add_path(rel, (Path *) new_bhpath);
}

/*
 * Recursively adjust cost for bitmap paths containing mismatched indexes
 * Returns true if the path uses a mismatched kmersearch index
 */
static void
kmersearch_adjust_bitmap_path_cost(Path *bitmapqual)
{
	if (bitmapqual == NULL)
		return;

	if (IsA(bitmapqual, IndexPath))
	{
		IndexPath *ipath = (IndexPath *) bitmapqual;

		if (kmersearch_is_kmersearch_gin_index_direct(ipath->indexinfo->indexoid))
		{
			KmersearchIndexSettings settings;

			if (kmersearch_get_index_settings_direct(ipath->indexinfo->indexoid, &settings))
			{
				if (!kmersearch_check_settings_match(&settings))
				{
					bitmapqual->startup_cost = 1.0e10;
					bitmapqual->total_cost = 1.0e10;
				}
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
	ListCell   *lc;
	bool		found_mismatched_bitmap = false;
	BitmapHeapPath *mismatched_bhpath = NULL;

	if (prev_set_rel_pathlist_hook)
		prev_set_rel_pathlist_hook(root, rel, rti, rte);

	if (in_kmersearch_hook)
		return;

	if (rel->reloptkind != RELOPT_BASEREL)
		return;

	in_kmersearch_hook = true;

	PG_TRY();
	{
		foreach(lc, rel->pathlist)
		{
			Path *path = (Path *) lfirst(lc);

			if (IsA(path, BitmapHeapPath))
			{
				BitmapHeapPath *bhpath = (BitmapHeapPath *) path;

				if (IsA(bhpath->bitmapqual, IndexPath))
				{
					IndexPath *ipath = (IndexPath *) bhpath->bitmapqual;

					if (kmersearch_is_kmersearch_gin_index_direct(ipath->indexinfo->indexoid))
					{
						KmersearchIndexSettings settings;

						if (kmersearch_get_index_settings_direct(ipath->indexinfo->indexoid, &settings))
						{
							if (!kmersearch_check_settings_match(&settings))
							{
								found_mismatched_bitmap = true;
								mismatched_bhpath = bhpath;
							}
						}
					}
				}

				kmersearch_adjust_bitmap_path_cost(bhpath->bitmapqual);

				if (bhpath->bitmapqual->total_cost >= 1.0e10)
				{
					path->startup_cost = 1.0e10;
					path->total_cost = 1.0e10;
				}
			}
			else if (IsA(path, IndexPath))
			{
				IndexPath *ipath = (IndexPath *) path;

				if (kmersearch_is_kmersearch_gin_index_direct(ipath->indexinfo->indexoid))
				{
					KmersearchIndexSettings settings;

					if (kmersearch_get_index_settings_direct(ipath->indexinfo->indexoid, &settings))
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

		if (found_mismatched_bitmap && mismatched_bhpath != NULL)
		{
			kmersearch_add_matching_index_path(root, rel, mismatched_bhpath);
		}
	}
	PG_FINALLY();
	{
		in_kmersearch_hook = false;
	}
	PG_END_TRY();
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
