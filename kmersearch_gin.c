/*
 * kmersearch_gin.c - GIN index support functions for pg_kmersearch
 *
 * This module contains all GIN-related functionality including:
 * - extract_value functions for DNA2 and DNA4 types
 * - extract_query function for query processing  
 * - consistent function for index consistency checking
 * - compare_partial function for partial key comparison
 * - Supporting utility functions for k-mer extraction and processing
 */

#include "kmersearch.h"

PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna2_int2);
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna2_int4);
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna2_int8);
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna4_int2);
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna4_int4);
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna4_int8);
PG_FUNCTION_INFO_V1(kmersearch_extract_query_int2);
PG_FUNCTION_INFO_V1(kmersearch_extract_query_int4);
PG_FUNCTION_INFO_V1(kmersearch_extract_query_int8);
PG_FUNCTION_INFO_V1(kmersearch_consistent_int2);
PG_FUNCTION_INFO_V1(kmersearch_consistent_int4);
PG_FUNCTION_INFO_V1(kmersearch_consistent_int8);

static void check_operator_class_compatibility(const char *opclass_type);

/*
 * Check operator class compatibility with current GUC settings
 */
static void
check_operator_class_compatibility(const char *opclass_type)
{
    int total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;
    int storage_bits;
    const char *optimal_type;
    
    /* Determine storage size based on operator class type */
    if (strcmp(opclass_type, "int2") == 0)
        storage_bits = 16;
    else if (strcmp(opclass_type, "int4") == 0)
        storage_bits = 32;
    else if (strcmp(opclass_type, "int8") == 0)
        storage_bits = 64;
    else
        return; /* Unknown type, skip check */
    
    /* Determine optimal operator class */
    if (total_bits <= 16)
        optimal_type = "int2";
    else if (total_bits <= 32)
        optimal_type = "int4";
    else
        optimal_type = "int8";
    
    /* Check if storage is sufficient */
    if (total_bits > storage_bits)
    {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("operator class kmersearch_*_gin_ops_%s cannot store current configuration", opclass_type),
                 errdetail("Required bits: %d (kmer_size=%d * 2 + occur_bitlen=%d), Storage capacity: %d bits",
                          total_bits, kmersearch_kmer_size, kmersearch_occur_bitlen, storage_bits),
                 errhint("Use kmersearch_*_gin_ops_%s operator class for this configuration", optimal_type)));
    }
    
    /* Error if not using optimal storage */
    if (strcmp(opclass_type, optimal_type) != 0)
    {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("operator class kmersearch_*_gin_ops_%s is not optimal for current configuration", opclass_type),
                 errdetail("Current configuration requires %d bits (kmer_size=%d * 2 + occur_bitlen=%d)",
                          total_bits, kmersearch_kmer_size, kmersearch_occur_bitlen),
                 errhint("Use kmersearch_*_gin_ops_%s operator class for optimal performance and memory usage", optimal_type)));
    }
}

/*
 * Get index information from index OID
 */
bool
kmersearch_get_index_info(Oid index_oid, Oid *table_oid, char **column_name, int *k_size)
{
    int ret;
    bool found = false;
    StringInfoData query;
    
    /* Connect to SPI */
    if (SPI_connect() != SPI_OK_CONNECT)
        return false;
    
    /* Build query to get index information */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT table_oid, column_name, kmer_size FROM kmersearch_index_info "
        "WHERE index_oid = %u",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    if (ret != SPI_OK_SELECT)
        ereport(ERROR, (errmsg("SPI_execute failed with code %d", ret)));
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        Datum table_oid_datum, k_size_datum;
        bool isnull;
        
        table_oid_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull && table_oid)
            *table_oid = DatumGetObjectId(table_oid_datum);
        
        if (column_name)
        {
            char *col_name = SPI_getvalue(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2);
            if (col_name)
                *column_name = pstrdup(col_name);
        }
        
        k_size_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull);
        if (!isnull && k_size)
            *k_size = DatumGetInt32(k_size_datum);
        
        found = true;
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return found;
}


/*
 * Filter Datum array for indexing by removing high-frequency k-mers
 */
static Datum *
kmersearch_filter_datum_for_indexing(Datum *keys, int *nkeys, size_t key_size, int k_size)
{
    Datum *filtered_keys = NULL;
    int filtered_count = 0;
    int i;
    bool has_highfreq = false;

    if (!kmersearch_preclude_highfreq_kmer || keys == NULL || *nkeys == 0)
        return keys;

    for (i = 0; i < *nkeys; i++)
    {
        uint64 uintkey_val;

        if (key_size == sizeof(uint16))
            uintkey_val = (uint64)DatumGetInt16(keys[i]);
        else if (key_size == sizeof(uint32))
            uintkey_val = (uint64)DatumGetInt32(keys[i]);
        else
            uintkey_val = (uint64)DatumGetInt64(keys[i]);

        if (kmersearch_is_uintkey_highfreq(uintkey_val, k_size))
        {
            has_highfreq = true;
            break;
        }
    }

    if (!has_highfreq)
        return keys;

    filtered_keys = (Datum *)palloc(*nkeys * sizeof(Datum));

    for (i = 0; i < *nkeys; i++)
    {
        uint64 uintkey_val;

        if (key_size == sizeof(uint16))
            uintkey_val = (uint64)DatumGetInt16(keys[i]);
        else if (key_size == sizeof(uint32))
            uintkey_val = (uint64)DatumGetInt32(keys[i]);
        else
            uintkey_val = (uint64)DatumGetInt64(keys[i]);

        if (!kmersearch_is_uintkey_highfreq(uintkey_val, k_size))
            filtered_keys[filtered_count++] = keys[i];
    }

    pfree(keys);

    if (filtered_count == 0)
    {
        pfree(filtered_keys);
        *nkeys = 0;
        return NULL;
    }

    *nkeys = filtered_count;
    return filtered_keys;
}

/*
 * New uintkey-based GIN extract_value functions for DNA2
 */
Datum
kmersearch_extract_value_dna2_int2(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    Datum *keys = NULL;

    /* Check operator class compatibility */
    check_operator_class_compatibility("int2");

    /* Extract and convert to Datum array in one step */
    keys = kmersearch_extract_datum_from_dna2(dna, nkeys, sizeof(uint16));

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    keys = kmersearch_filter_datum_for_indexing(keys, nkeys, sizeof(uint16), kmersearch_kmer_size);

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    PG_RETURN_POINTER(keys);
}

Datum
kmersearch_extract_value_dna2_int4(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    Datum *keys = NULL;

    /* Check operator class compatibility */
    check_operator_class_compatibility("int4");

    /* Extract and convert to Datum array in one step */
    keys = kmersearch_extract_datum_from_dna2(dna, nkeys, sizeof(uint32));

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    keys = kmersearch_filter_datum_for_indexing(keys, nkeys, sizeof(uint32), kmersearch_kmer_size);

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    PG_RETURN_POINTER(keys);
}

Datum
kmersearch_extract_value_dna2_int8(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    Datum *keys = NULL;

    /* Check operator class compatibility */
    check_operator_class_compatibility("int8");

    /* Extract and convert to Datum array in one step */
    keys = kmersearch_extract_datum_from_dna2(dna, nkeys, sizeof(uint64));

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    keys = kmersearch_filter_datum_for_indexing(keys, nkeys, sizeof(uint64), kmersearch_kmer_size);

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    PG_RETURN_POINTER(keys);
}

/*
 * New uintkey-based GIN extract_value functions for DNA4
 */
Datum
kmersearch_extract_value_dna4_int2(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    Datum *keys = NULL;

    /* Check operator class compatibility */
    check_operator_class_compatibility("int2");

    /* Extract and convert to Datum array in one step */
    keys = kmersearch_extract_datum_from_dna4(dna, nkeys, sizeof(uint16));

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    keys = kmersearch_filter_datum_for_indexing(keys, nkeys, sizeof(uint16), kmersearch_kmer_size);

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    PG_RETURN_POINTER(keys);
}

Datum
kmersearch_extract_value_dna4_int4(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    Datum *keys = NULL;

    /* Check operator class compatibility */
    check_operator_class_compatibility("int4");

    /* Extract and convert to Datum array in one step */
    keys = kmersearch_extract_datum_from_dna4(dna, nkeys, sizeof(uint32));

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    keys = kmersearch_filter_datum_for_indexing(keys, nkeys, sizeof(uint32), kmersearch_kmer_size);

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    PG_RETURN_POINTER(keys);
}

Datum
kmersearch_extract_value_dna4_int8(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    Datum *keys = NULL;

    /* Check operator class compatibility */
    check_operator_class_compatibility("int8");

    /* Extract and convert to Datum array in one step */
    keys = kmersearch_extract_datum_from_dna4(dna, nkeys, sizeof(uint64));

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    keys = kmersearch_filter_datum_for_indexing(keys, nkeys, sizeof(uint64), kmersearch_kmer_size);

    if (keys == NULL || *nkeys == 0)
        PG_RETURN_POINTER(NULL);

    PG_RETURN_POINTER(keys);
}

/*
 * Filter uintkey array and set actual_min_score in cache
 * This function filters out high-frequency k-mers and caches the actual_min_score
 */
void *
kmersearch_filter_uintkey_and_set_actual_min_score(void *uintkey, int *nkeys, 
                                        const char *query_string, int k_size)
{
    void *filtered_keys = NULL;
    int original_nkeys = *nkeys;
    int filtered_count = 0;
    int i;
    bool has_highfreq = false;
    int total_bits;
    
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    
    if (!kmersearch_preclude_highfreq_kmer || uintkey == NULL || *nkeys == 0)
        return uintkey;
    
    /* First pass: check if there are any high-frequency k-mers */
    if (total_bits <= 16)
    {
        uint16 *keys = (uint16 *)uintkey;
        for (i = 0; i < *nkeys; i++)
        {
            if (kmersearch_is_uintkey_highfreq((uint64)keys[i], k_size))
            {
                has_highfreq = true;
                break;
            }
        }
    }
    else if (total_bits <= 32)
    {
        uint32 *keys = (uint32 *)uintkey;
        for (i = 0; i < *nkeys; i++)
        {
            if (kmersearch_is_uintkey_highfreq((uint64)keys[i], k_size))
            {
                has_highfreq = true;
                break;
            }
        }
    }
    else
    {
        uint64 *keys = (uint64 *)uintkey;
        for (i = 0; i < *nkeys; i++)
        {
            if (kmersearch_is_uintkey_highfreq(keys[i], k_size))
            {
                has_highfreq = true;
                break;
            }
        }
    }
    
    /* If no high-frequency k-mers, return original array */
    if (!has_highfreq)
    {
        kmersearch_get_cached_actual_min_score_uintkey(uintkey, *nkeys, k_size);
        return uintkey;
    }
    
    /* Allocate new array for filtered keys */
    if (total_bits <= 16)
    {
        uint16 *original = (uint16 *)uintkey;
        uint16 *filtered = (uint16 *)palloc(*nkeys * sizeof(uint16));
        
        for (i = 0; i < *nkeys; i++)
        {
            if (!kmersearch_is_uintkey_highfreq((uint64)original[i], k_size))
                filtered[filtered_count++] = original[i];
        }
        
        if (filtered_count > 0)
        {
            filtered_keys = repalloc(filtered, filtered_count * sizeof(uint16));
        }
        else
        {
            pfree(filtered);
        }
    }
    else if (total_bits <= 32)
    {
        uint32 *original = (uint32 *)uintkey;
        uint32 *filtered = (uint32 *)palloc(*nkeys * sizeof(uint32));
        
        for (i = 0; i < *nkeys; i++)
        {
            if (!kmersearch_is_uintkey_highfreq((uint64)original[i], k_size))
                filtered[filtered_count++] = original[i];
        }
        
        if (filtered_count > 0)
        {
            filtered_keys = repalloc(filtered, filtered_count * sizeof(uint32));
        }
        else
        {
            pfree(filtered);
        }
    }
    else
    {
        uint64 *original = (uint64 *)uintkey;
        uint64 *filtered = (uint64 *)palloc(*nkeys * sizeof(uint64));
        
        for (i = 0; i < *nkeys; i++)
        {
            if (!kmersearch_is_uintkey_highfreq(original[i], k_size))
                filtered[filtered_count++] = original[i];
        }
        
        if (filtered_count > 0)
        {
            filtered_keys = repalloc(filtered, filtered_count * sizeof(uint64));
        }
        else
        {
            pfree(filtered);
        }
    }
    
    /* Cache actual_min_score - this will be retrieved in consistent function */
    kmersearch_get_cached_actual_min_score_uintkey(filtered_keys ? filtered_keys : uintkey, 
                                        filtered_keys ? filtered_count : *nkeys, k_size);
    
    *nkeys = filtered_count;
    
    return filtered_keys;
}

/*
 * Check if a uintkey is high-frequency
 */
bool
kmersearch_is_uintkey_highfreq(uint64 uintkey, int k_size)
{
    /* Extract just the k-mer portion (remove occurrence bits) */
    uint64 kmer_only;
    bool is_highfreq = false;
    int ret;
    int total_bits;
    
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    
    /* For uintkey format, k-mer is in higher bits, occurrence in lower bits */
    kmer_only = uintkey >> kmersearch_occur_bitlen;
    
    /* Priority 1: Check in global cache (highest priority) */
    if (global_highfreq_cache.is_valid && global_highfreq_cache.highfreq_hash) {
        return kmersearch_lookup_uintkey_in_global_cache(kmer_only, NULL, NULL);
    }
    
    /* Priority 2: Check in parallel cache */
    if (kmersearch_is_parallel_highfreq_cache_loaded()) {
        return kmersearch_lookup_uintkey_in_parallel_cache(kmer_only, NULL, NULL);
    }
    
    /* Priority 3: Check kmersearch_highfreq_kmer table directly */
    ret = SPI_connect();
    if (ret == SPI_OK_CONNECT) {
        StringInfoData query;
        
        initStringInfo(&query);
        
        /* Build query based on total_bits */
        if (total_bits <= 16) {
            appendStringInfo(&query,
                "SELECT 1 FROM kmersearch_highfreq_kmer "
                "WHERE uintkey = %u "
                "LIMIT 1",
                (unsigned int)kmer_only);
        } else if (total_bits <= 32) {
            appendStringInfo(&query,
                "SELECT 1 FROM kmersearch_highfreq_kmer "
                "WHERE uintkey = %u "
                "LIMIT 1", 
                (unsigned int)kmer_only);
        } else {
            appendStringInfo(&query,
                "SELECT 1 FROM kmersearch_highfreq_kmer "
                "WHERE uintkey = %lu "
                "LIMIT 1",
                kmer_only);
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
 * New uintkey-based GIN extract_query functions
 */
Datum
kmersearch_extract_query_int2(PG_FUNCTION_ARGS)
{
    Datum query = PG_GETARG_DATUM(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    StrategyNumber strategy = PG_GETARG_UINT16(2);
    bool **pmatch = (bool **) PG_GETARG_POINTER(3);
    Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4);
    bool **nullFlags = (bool **) PG_GETARG_POINTER(5);
    int32 *searchMode = (int32 *) PG_GETARG_POINTER(6);
    
    text *query_text = DatumGetTextP(query);
    char *query_string = text_to_cstring(query_text);
    Datum *keys = NULL;
    void *uintkey = NULL;
    uint16 *uint16_keys;
    int i;
    
    /* Use cached query-kmer extraction */
    uintkey = kmersearch_get_cached_query_uintkey(query_string, kmersearch_kmer_size, nkeys);
    
    /* Filter high-frequency k-mers and cache actual_min_score */
    if (uintkey != NULL && *nkeys > 0)
    {
        uintkey = kmersearch_filter_uintkey_and_set_actual_min_score(uintkey, nkeys, 
                                                           query_string, kmersearch_kmer_size);
    }
    
    if (uintkey == NULL || *nkeys == 0)
    {
        PG_RETURN_POINTER(NULL);
    }
    
    /* Convert to Datum array */
    uint16_keys = (uint16 *)uintkey;
    keys = (Datum *) palloc(*nkeys * sizeof(Datum));
    for (i = 0; i < *nkeys; i++)
    {
        keys[i] = Int16GetDatum(uint16_keys[i]);
    }
    
    *searchMode = GIN_SEARCH_MODE_DEFAULT;
    PG_RETURN_POINTER(keys);
}

Datum
kmersearch_extract_query_int4(PG_FUNCTION_ARGS)
{
    Datum query = PG_GETARG_DATUM(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    StrategyNumber strategy = PG_GETARG_UINT16(2);
    bool **pmatch = (bool **) PG_GETARG_POINTER(3);
    Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4);
    bool **nullFlags = (bool **) PG_GETARG_POINTER(5);
    int32 *searchMode = (int32 *) PG_GETARG_POINTER(6);
    
    text *query_text = DatumGetTextP(query);
    char *query_string = text_to_cstring(query_text);
    Datum *keys = NULL;
    void *uintkey = NULL;
    uint32 *uint32_keys;
    int i;
    
    /* Use cached query-kmer extraction */
    uintkey = kmersearch_get_cached_query_uintkey(query_string, kmersearch_kmer_size, nkeys);
    
    /* Filter high-frequency k-mers and cache actual_min_score */
    if (uintkey != NULL && *nkeys > 0)
    {
        uintkey = kmersearch_filter_uintkey_and_set_actual_min_score(uintkey, nkeys, 
                                                           query_string, kmersearch_kmer_size);
    }
    
    if (uintkey == NULL || *nkeys == 0)
    {
        PG_RETURN_POINTER(NULL);
    }
    
    /* Convert to Datum array */
    uint32_keys = (uint32 *)uintkey;
    keys = (Datum *) palloc(*nkeys * sizeof(Datum));
    for (i = 0; i < *nkeys; i++)
    {
        keys[i] = Int32GetDatum(uint32_keys[i]);
    }
    
    *searchMode = GIN_SEARCH_MODE_DEFAULT;
    PG_RETURN_POINTER(keys);
}

Datum
kmersearch_extract_query_int8(PG_FUNCTION_ARGS)
{
    Datum query = PG_GETARG_DATUM(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    StrategyNumber strategy = PG_GETARG_UINT16(2);
    bool **pmatch = (bool **) PG_GETARG_POINTER(3);
    Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4);
    bool **nullFlags = (bool **) PG_GETARG_POINTER(5);
    int32 *searchMode = (int32 *) PG_GETARG_POINTER(6);
    
    text *query_text = DatumGetTextP(query);
    char *query_string = text_to_cstring(query_text);
    Datum *keys = NULL;
    void *uintkey = NULL;
    uint64 *uint64_keys;
    int i;
    
    /* Use cached query-kmer extraction */
    uintkey = kmersearch_get_cached_query_uintkey(query_string, kmersearch_kmer_size, nkeys);
    
    /* Filter high-frequency k-mers and cache actual_min_score */
    if (uintkey != NULL && *nkeys > 0)
    {
        uintkey = kmersearch_filter_uintkey_and_set_actual_min_score(uintkey, nkeys, 
                                                           query_string, kmersearch_kmer_size);
    }
    
    if (uintkey == NULL || *nkeys == 0)
    {
        PG_RETURN_POINTER(NULL);
    }
    
    /* Convert to Datum array */
    uint64_keys = (uint64 *)uintkey;
    keys = (Datum *) palloc(*nkeys * sizeof(Datum));
    for (i = 0; i < *nkeys; i++)
    {
        keys[i] = Int64GetDatum(uint64_keys[i]);
    }
    
    *searchMode = GIN_SEARCH_MODE_DEFAULT;
    PG_RETURN_POINTER(keys);
}

/*
 * New uintkey-based GIN consistent functions
 */
Datum
kmersearch_consistent_int2(PG_FUNCTION_ARGS)
{
    bool *check = (bool *) PG_GETARG_POINTER(0);
    StrategyNumber strategy = PG_GETARG_UINT16(1);
    text *query_text = PG_GETARG_TEXT_P(2);
    int32 nkeys = PG_GETARG_INT32(3);
    Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4);
    bool *recheck = (bool *) PG_GETARG_POINTER(5);
    Datum *queryKeys = (Datum *) PG_GETARG_POINTER(6);
    bool *nullFlags = (bool *) PG_GETARG_POINTER(7);
    
    int shared_count = 0;
    int actual_min_score;
    int i;
    
    *recheck = false;
    
    /* Count matches */
    for (i = 0; i < nkeys; i++)
    {
        if (check[i])
            shared_count++;
    }
    
    /* Get cached actual_min_score */
    actual_min_score = kmersearch_get_cached_actual_min_score_datum_int2(queryKeys, nkeys);
    
    PG_RETURN_BOOL(shared_count >= actual_min_score);
}

Datum
kmersearch_consistent_int4(PG_FUNCTION_ARGS)
{
    bool *check = (bool *) PG_GETARG_POINTER(0);
    StrategyNumber strategy = PG_GETARG_UINT16(1);
    text *query_text = PG_GETARG_TEXT_P(2);
    int32 nkeys = PG_GETARG_INT32(3);
    Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4);
    bool *recheck = (bool *) PG_GETARG_POINTER(5);
    Datum *queryKeys = (Datum *) PG_GETARG_POINTER(6);
    bool *nullFlags = (bool *) PG_GETARG_POINTER(7);
    
    int shared_count = 0;
    int actual_min_score;
    int i;
    
    *recheck = false;
    
    /* Count matches */
    for (i = 0; i < nkeys; i++)
    {
        if (check[i])
            shared_count++;
    }
    
    /* Get cached actual_min_score */
    actual_min_score = kmersearch_get_cached_actual_min_score_datum_int4(queryKeys, nkeys);
    
    PG_RETURN_BOOL(shared_count >= actual_min_score);
}

Datum
kmersearch_consistent_int8(PG_FUNCTION_ARGS)
{
    bool *check = (bool *) PG_GETARG_POINTER(0);
    StrategyNumber strategy = PG_GETARG_UINT16(1);
    text *query_text = PG_GETARG_TEXT_P(2);
    int32 nkeys = PG_GETARG_INT32(3);
    Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4);
    bool *recheck = (bool *) PG_GETARG_POINTER(5);
    Datum *queryKeys = (Datum *) PG_GETARG_POINTER(6);
    bool *nullFlags = (bool *) PG_GETARG_POINTER(7);
    
    int shared_count = 0;
    int actual_min_score;
    int i;
    
    *recheck = false;
    
    /* Count matches */
    for (i = 0; i < nkeys; i++)
    {
        if (check[i])
            shared_count++;
    }
    
    /* Get cached actual_min_score */
    actual_min_score = kmersearch_get_cached_actual_min_score_datum_int8(queryKeys, nkeys);
    
    PG_RETURN_BOOL(shared_count >= actual_min_score);
}

