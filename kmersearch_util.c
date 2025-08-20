/*
 * kmersearch_util.c - Utility functions for pg_kmersearch
 *
 * This module contains general utility functions that don't fit
 * into other specific modules.
 */

#include "kmersearch.h"

/* Function to show build information */
PG_FUNCTION_INFO_V1(kmersearch_show_buildno);
Datum
kmersearch_show_buildno(PG_FUNCTION_ARGS)
{
    PG_RETURN_TEXT_P(cstring_to_text(KMERSEARCH_BUILD_VERSION));
}