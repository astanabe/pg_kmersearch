/*
 * kmersearch_util.c - Utility functions for pg_kmersearch
 *
 * This module contains general utility functions that don't fit
 * into other specific modules.
 */

#include "kmersearch.h"

/*
 * Get DNA2 type OID
 */
Oid
get_dna2_type_oid(void)
{
    return TypenameGetTypid("dna2");
}

/*
 * Get DNA4 type OID  
 */
Oid
get_dna4_type_oid(void)
{
    return TypenameGetTypid("dna4");
}