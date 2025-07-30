/*
 * kmersearch_cache.c - Cache management functions for pg_kmersearch
 *
 * This module contains all cache-related functionality including:
 * - Query pattern cache for storing parsed k-mer patterns
 * - Actual min score cache for threshold calculations
 * - Cache manager creation, lookup, eviction, and cleanup functions
 */

#include "kmersearch.h"

/* PostgreSQL function info declarations for cache functions */
PG_FUNCTION_INFO_V1(kmersearch_query_pattern_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_query_pattern_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_load);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_free_all);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_load);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_free_all);

/* Global high-frequency k-mer cache */
HighfreqKmerCache global_highfreq_cache = {0};

/* Global testing variable for dshash usage (not exposed to users) */
bool kmersearch_force_use_parallel_highfreq_kmer_cache = false;

/* Global parallel cache state */
ParallelHighfreqKmerCache *parallel_highfreq_cache = NULL;
dsm_segment *parallel_cache_segment = NULL;
dsa_area *parallel_cache_dsa = NULL;
dshash_table *parallel_cache_hash = NULL;

/*
 * Forward declarations for internal functions
 */

/* Identity hash functions for k-mer values */
static uint32 kmersearch_uint16_identity_hash(const void *key, size_t keysize, void *arg);
static uint32 kmersearch_uint32_identity_hash(const void *key, size_t keysize, void *arg);

/* Query pattern cache functions */
static void init_query_pattern_cache_manager(QueryPatternCacheManager **manager);
static uint64 generate_query_pattern_cache_key(const char *query_string, int k_size);
static QueryPatternCacheEntry *lookup_query_pattern_cache_entry(QueryPatternCacheManager *manager, const char *query_string, int k_size);
static void store_query_pattern_cache_entry(QueryPatternCacheManager *manager, uint64 hash_key, const char *query_string, int k_size, VarBit **kmers, int kmer_count);
static void lru_touch_query_pattern_cache(QueryPatternCacheManager *manager, QueryPatternCacheEntry *entry);
static void lru_evict_oldest_query_pattern_cache(QueryPatternCacheManager *manager);
void free_query_pattern_cache_manager(QueryPatternCacheManager **manager);

/* Actual min score cache functions */
static void create_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);
void free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);

/* High-frequency k-mer cache functions */
void kmersearch_highfreq_kmer_cache_init(void);
bool kmersearch_highfreq_kmer_cache_load_internal(Oid table_oid, const char *column_name, int k_value);
void kmersearch_highfreq_kmer_cache_free_internal(void);
bool kmersearch_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value);


/* High-level cache functions */

/* External global variables (defined in kmersearch.c) */

/*
 * Query Pattern Cache Functions
 */

/*
 * Initialize query pattern cache manager
 */
static void
init_query_pattern_cache_manager(QueryPatternCacheManager **manager)
{
    if (*manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        HASHCTL hash_ctl;
        
        /* Allocate manager in TopMemoryContext */
        *manager = (QueryPatternCacheManager *) palloc0(sizeof(QueryPatternCacheManager));
        
        /* Create query pattern cache context under TopMemoryContext */
        (*manager)->query_pattern_cache_context = AllocSetContextCreate(TopMemoryContext,
                                                                        "QueryPatternCache",
                                                                        ALLOCSET_DEFAULT_SIZES);
        
        /* Initialize query pattern cache parameters */
        (*manager)->max_entries = kmersearch_query_pattern_cache_max_entries;
        (*manager)->current_entries = 0;
        (*manager)->hits = 0;
        (*manager)->misses = 0;
        (*manager)->lru_head = NULL;
        (*manager)->lru_tail = NULL;
        
        /* Create hash table */
        MemSet(&hash_ctl, 0, sizeof(hash_ctl));
        hash_ctl.keysize = sizeof(uint64);
        hash_ctl.entrysize = sizeof(QueryPatternCacheEntry);
        hash_ctl.hcxt = (*manager)->query_pattern_cache_context;
        
        (*manager)->hash_table = hash_create("QueryPatternCache", 256, &hash_ctl,
                                           HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
        
        MemoryContextSwitchTo(old_context);
    }
}

/*
 * Generate cache key for query pattern
 */
static uint64
generate_query_pattern_cache_key(const char *query_string, int k_size)
{
    uint64 query_hash, k_hash;
    
    /* Hash query string */
    query_hash = hash_any_extended((unsigned char*)query_string, strlen(query_string), 0);
    
    /* Hash k_size with different seed */
    k_hash = hash_any_extended((unsigned char*)&k_size, sizeof(k_size), 1);
    
    /* Combine hashes */
    return query_hash ^ (k_hash << 1);
}

/*
 * Move entry to head of LRU chain (most recently used)
 */
static void
lru_touch_query_pattern_cache(QueryPatternCacheManager *manager, QueryPatternCacheEntry *entry)
{
    if (entry == manager->lru_head)
        return;  /* Already at head */
    
    /* Remove from current position */
    if (entry->prev)
        entry->prev->next = entry->next;
    else
        manager->lru_tail = entry->next;
    
    if (entry->next)
        entry->next->prev = entry->prev;
    else
        manager->lru_head = entry->prev;
    
    /* Add to head */
    entry->prev = NULL;
    entry->next = manager->lru_head;
    if (manager->lru_head)
        manager->lru_head->prev = entry;
    else
        manager->lru_tail = entry;
    manager->lru_head = entry;
}

/*
 * Evict oldest entry from query pattern cache
 */
static void
lru_evict_oldest_query_pattern_cache(QueryPatternCacheManager *manager)
{
    QueryPatternCacheEntry *tail = manager->lru_tail;
    bool found;
    int i;
    
    if (!tail)
        return;
    
    /* Remove from hash table */
    hash_search(manager->hash_table, &tail->hash_key, HASH_REMOVE, &found);
    
    /* Remove from LRU chain */
    if (tail->prev)
        tail->prev->next = NULL;
    else
        manager->lru_head = NULL;
    manager->lru_tail = tail->prev;
    
    /* Free allocated memory */
    if (tail->query_string_copy)
        pfree(tail->query_string_copy);
    if (tail->extracted_kmers)
    {
        for (i = 0; i < tail->kmer_count; i++)
        {
            if (tail->extracted_kmers[i])
                pfree(tail->extracted_kmers[i]);
        }
        pfree(tail->extracted_kmers);
    }
    
    manager->current_entries--;
}

/*
 * Look up query pattern cache entry
 */
static QueryPatternCacheEntry *
lookup_query_pattern_cache_entry(QueryPatternCacheManager *manager, const char *query_string, int k_size)
{
    uint64 hash_key = generate_query_pattern_cache_key(query_string, k_size);
    QueryPatternCacheEntry *entry;
    bool found;
    
    entry = (QueryPatternCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_FIND, &found);
    
    if (found && entry && strcmp(entry->query_string_copy, query_string) == 0 && entry->kmer_size == k_size)
    {
        /* Cache hit - move to head of LRU */
        lru_touch_query_pattern_cache(manager, entry);
        manager->hits++;
        return entry;
    }
    
    return NULL;  /* Cache miss */
}

/*
 * Store entry in query pattern cache
 */
static void
store_query_pattern_cache_entry(QueryPatternCacheManager *manager, uint64 hash_key, 
                               const char *query_string, int k_size, VarBit **kmers, int kmer_count)
{
    QueryPatternCacheEntry *entry;
    MemoryContext old_context;
    bool found;
    int i;
    
    /* Evict oldest entries if cache is full */
    while (manager->current_entries >= manager->max_entries)
        lru_evict_oldest_query_pattern_cache(manager);
    
    old_context = MemoryContextSwitchTo(manager->query_pattern_cache_context);
    
    /* Create new entry */
    entry = (QueryPatternCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_ENTER, &found);
    if (!found)
    {
        entry->hash_key = hash_key;
        entry->query_string_copy = pstrdup(query_string);
        entry->kmer_size = k_size;
        entry->kmer_count = kmer_count;
        
        /* Copy k-mers */
        entry->extracted_kmers = (VarBit **) palloc(kmer_count * sizeof(VarBit *));
        for (i = 0; i < kmer_count; i++)
        {
            entry->extracted_kmers[i] = (VarBit *) palloc(VARSIZE(kmers[i]));
            memcpy(entry->extracted_kmers[i], kmers[i], VARSIZE(kmers[i]));
        }
        
        /* Add to LRU chain */
        entry->prev = NULL;
        entry->next = manager->lru_head;
        if (manager->lru_head)
            manager->lru_head->prev = entry;
        else
            manager->lru_tail = entry;
        manager->lru_head = entry;
        
        manager->current_entries++;
    }
    
    MemoryContextSwitchTo(old_context);
}

/*
 * Get cached query k-mers or extract and cache them
 */
VarBit **
get_cached_query_kmer(const char *query_string, int k_size, int *nkeys)
{
    QueryPatternCacheEntry *cache_entry;
    VarBit **result_kmers = NULL;
    VarBit **extracted_kmers = NULL;
    uint64 hash_key;
    int i;
    MemoryContext old_context;
    
    *nkeys = 0;
    
    /* Initialize query pattern cache manager if not already done */
    if (query_pattern_cache_manager == NULL)
    {
        old_context = MemoryContextSwitchTo(TopMemoryContext);
        init_query_pattern_cache_manager(&query_pattern_cache_manager);
        MemoryContextSwitchTo(old_context);
    }
    
    /* Try to find in cache first */
    cache_entry = lookup_query_pattern_cache_entry(query_pattern_cache_manager, query_string, k_size);
    if (cache_entry != NULL)
    {
        /* Cache hit - return pointers to cached k-mers directly */
        *nkeys = cache_entry->kmer_count;
        return cache_entry->extracted_kmers;
    }
    
    /* Cache miss - extract k-mers and store in cache */
    query_pattern_cache_manager->misses++;
    extracted_kmers = kmersearch_extract_query_ngram_key2(query_string, k_size, nkeys);
    
    if (extracted_kmers != NULL && *nkeys > 0)
    {
        /* Store in cache */
        hash_key = generate_query_pattern_cache_key(query_string, k_size);
        store_query_pattern_cache_entry(query_pattern_cache_manager, hash_key, 
                                       query_string, k_size, extracted_kmers, *nkeys);
        
        /* Find the cache entry we just stored and return its k-mers */
        cache_entry = lookup_query_pattern_cache_entry(query_pattern_cache_manager, query_string, k_size);
        if (cache_entry != NULL) {
            /* Free the original extracted k-mers (cache has its own copy) */
            for (i = 0; i < *nkeys; i++)
            {
                pfree(extracted_kmers[i]);
            }
            pfree(extracted_kmers);
            
            return cache_entry->extracted_kmers;
        }
    }
    
    return extracted_kmers;
}

/*
 * Free query pattern cache manager
 */
void
free_query_pattern_cache_manager(QueryPatternCacheManager **manager)
{
    if (*manager)
    {
        /* Delete the query pattern cache context, which will free all allocated memory */
        if ((*manager)->query_pattern_cache_context)
            MemoryContextDelete((*manager)->query_pattern_cache_context);
        
        /* Free the manager itself (allocated in TopMemoryContext) */
        pfree(*manager);
        *manager = NULL;
    }
}

/*
 * Actual Min Score Cache Functions
 */

/*
 * Create actual min score cache manager
 */
static void
create_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager)
{
    HASHCTL hash_ctl;
    
    /* Allocate manager in current context (caller should have set appropriate context) */
    *manager = (ActualMinScoreCacheManager *) palloc0(sizeof(ActualMinScoreCacheManager));
    
    /* Create actual min score cache context under current context */
    (*manager)->cache_context = AllocSetContextCreate(CurrentMemoryContext,
                                                   "ActualMinScoreCache",
                                                   ALLOCSET_DEFAULT_SIZES);
    
    /* Initialize parameters */
    (*manager)->hits = 0;
    (*manager)->misses = 0;
    (*manager)->max_entries = kmersearch_actual_min_score_cache_max_entries;
    (*manager)->current_entries = 0;
    
    /* Create hash table */
    MemSet(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64);  /* Hash value as key */
    hash_ctl.entrysize = sizeof(ActualMinScoreCacheEntry);
    hash_ctl.hcxt = (*manager)->cache_context;
    
    (*manager)->cache_hash = hash_create("ActualMinScoreCache", 256, &hash_ctl,
                                     HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
}

/*
 * Free actual min score cache manager
 */
void
free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager)
{
    if (*manager)
    {
        /* Delete the actual min score cache context, which will free all allocated memory */
        if ((*manager)->cache_context)
            MemoryContextDelete((*manager)->cache_context);
        *manager = NULL;
    }
}

/*
 * Calculate actual minimum score considering thresholds
 */
int
calculate_actual_min_score(VarBit **query_keys, int nkeys, int query_total_kmers)
{
    int base_min_score;
    int highfreq_count;
    int actual_min_score;
    int absolute_min = kmersearch_min_score;
    int relative_min = 0;
    
    /* Validate input parameters */
    if (query_keys == NULL) {
        return kmersearch_min_score;
    }
    
    /* Calculate base minimum score (maximum of absolute and relative) */
    
    if (query_total_kmers > 0) {
        double relative_threshold = kmersearch_min_shared_ngram_key_rate * query_total_kmers;
        relative_min = (int)ceil(relative_threshold);
    }
    
    base_min_score = (absolute_min > relative_min) ? absolute_min : relative_min;
    
    /* If high-frequency k-mer filtering is enabled, subtract high-frequency k-mer count */
    if (kmersearch_is_highfreq_filtering_enabled()) {
        /* Debug: Log before calling kmersearch_count_highfreq_kmer_in_query */
        elog(DEBUG1, "calculate_actual_min_score: calling kmersearch_count_highfreq_kmer_in_query with nkeys = %d", nkeys);
        highfreq_count = kmersearch_count_highfreq_kmer_in_query(query_keys, nkeys);
        actual_min_score = base_min_score - highfreq_count;
        
        /* Ensure non-negative value */
        if (actual_min_score < 0) {
            actual_min_score = 0;
        }
    } else {
        actual_min_score = base_min_score;
    }
    
    return actual_min_score;
}

/*
 * Get cached actual min score using TopMemoryContext cache (global)
 * Modified to create cache key from filtered keys while calculating with original keys
 */
int
get_cached_actual_min_score(VarBit **query_keys, int nkeys)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash;
    bool found;
    MemoryContext old_context;
    int actual_min_score;
    VarBit **filtered_keys = NULL;
    int filtered_nkeys = nkeys;
    int i;
    
    /* Validate input parameters */
    if (query_keys == NULL) {
        return calculate_actual_min_score(query_keys, nkeys, nkeys);
    }
    
    /* Validate query_keys array elements */
    for (i = 0; i < nkeys; i++) {
        if (query_keys[i] == NULL) {
            return calculate_actual_min_score(query_keys, nkeys, nkeys);
        }
    }
    
    /* Create cache manager in TopMemoryContext if not exists */
    if (actual_min_score_cache_manager == NULL)
    {
        old_context = MemoryContextSwitchTo(TopMemoryContext);
        
        PG_TRY();
        {
            create_actual_min_score_cache_manager(&actual_min_score_cache_manager);
        }
        PG_CATCH();
        {
            MemoryContextSwitchTo(old_context);
            /* Fallback to direct calculation if cache creation fails */
            return calculate_actual_min_score(query_keys, nkeys, nkeys);
        }
        PG_END_TRY();
        
        MemoryContextSwitchTo(old_context);
    }
    
    /* Step 1: Filter high-frequency k-mers for cache key generation */
    if (kmersearch_preclude_highfreq_kmer && kmersearch_is_highfreq_filtering_enabled()) {
        filtered_keys = (VarBit **) palloc(nkeys * sizeof(VarBit *));
        filtered_nkeys = 0;
        
        for (i = 0; i < nkeys; i++) {
            if (!kmersearch_is_kmer_highfreq(query_keys[i])) {
                filtered_keys[filtered_nkeys++] = query_keys[i];
            }
        }
        
        /* Calculate hash using filtered keys */
        query_hash = 0;
        for (i = 0; i < filtered_nkeys; i++) {
            uint64 kmer_hash;
            kmer_hash = hash_any_extended((unsigned char *)VARBITS(filtered_keys[i]), 
                                          VARBITBYTES(filtered_keys[i]), query_hash);
            query_hash = kmer_hash;
        }
    } else {
        /* No filtering - use original keys for hash */
        query_hash = 0;
        for (i = 0; i < nkeys; i++) {
            uint64 kmer_hash;
            
            /* Validate VarBit structure before accessing */
            if (query_keys[i] == NULL) {
                elog(ERROR, "get_cached_actual_min_score: query_keys[%d] is NULL", i);
            }
            
            kmer_hash = hash_any_extended((unsigned char *)VARBITS(query_keys[i]), 
                                          VARBITBYTES(query_keys[i]), query_hash);
            query_hash = kmer_hash;
        }
    }
    
    
    /* Look up in hash table */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, &found);
    
    
    if (found) {
        actual_min_score_cache_manager->hits++;
        if (filtered_keys) pfree(filtered_keys);
        return cache_entry->actual_min_score;
    }
    
    /* Not found - calculate and cache */
    actual_min_score_cache_manager->misses++;
    /* Note: Calculate with original keys (includes high-frequency k-mers) */
    actual_min_score = calculate_actual_min_score(query_keys, nkeys, nkeys);
    
    /* Add to cache if not at capacity */
    if (actual_min_score_cache_manager->current_entries < actual_min_score_cache_manager->max_entries)
    {
        old_context = MemoryContextSwitchTo(actual_min_score_cache_manager->cache_context);
        
        PG_TRY();
        {
            cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                                  &query_hash, HASH_ENTER, &found);
            
            if (cache_entry != NULL && !found) {
                cache_entry->query_hash = query_hash;
                cache_entry->actual_min_score = actual_min_score;
                actual_min_score_cache_manager->current_entries++;
            }
        }
        PG_CATCH();
        {
            /* Ignore cache storage errors and continue */
        }
        PG_END_TRY();
        
        MemoryContextSwitchTo(old_context);
    }
    
    if (filtered_keys) pfree(filtered_keys);
    
    return actual_min_score;
}

/*
 * Get cached actual_min_score or error if not found
 * For use in kmersearch_consistent where cache must already be populated
 * Note: query_keys should already be filtered (high-frequency k-mers removed)
 */
int
get_cached_actual_min_score_or_error(VarBit **query_keys, int nkeys)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash;
    int i;
    
    /* Cache must be initialized */
    if (actual_min_score_cache_manager == NULL) {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("actual_min_score cache not initialized"),
                 errhint("This should not happen. Extract_query should have initialized the cache.")));
    }
    
    /* Calculate hash from query_keys (already filtered) */
    query_hash = 0;
    for (i = 0; i < nkeys; i++) {
        uint64 kmer_hash;
        
        if (query_keys[i] == NULL) {
            elog(ERROR, "get_cached_actual_min_score_or_error: query_keys[%d] is NULL", i);
        }
        
        kmer_hash = hash_any_extended((unsigned char *)VARBITS(query_keys[i]), 
                                      VARBITBYTES(query_keys[i]), query_hash);
        query_hash = kmer_hash;
    }
    
    /* Look up in cache */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, NULL);
    
    if (cache_entry == NULL) {
        /* Cache miss should not happen - extract_query should have cached it */
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("actual_min_score not found in cache"),
                 errdetail("Query hash: %lu, nkeys: %d", query_hash, nkeys),
                 errhint("This should not happen. Extract_query should have cached this value.")));
    }
    
    actual_min_score_cache_manager->hits++;
    
    return cache_entry->actual_min_score;
}








/*
 * Query pattern cache statistics function
 */
Datum
kmersearch_query_pattern_cache_stats(PG_FUNCTION_ARGS)
{
    TupleDesc tupdesc;
    Datum values[4];
    bool nulls[4] = {false};
    HeapTuple tuple;
    
    /* Build tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("function returning record called in context that cannot accept a set")));
    
    /* Query pattern cache statistics */
    if (query_pattern_cache_manager)
    {
        values[0] = Int64GetDatum(query_pattern_cache_manager->hits);
        values[1] = Int64GetDatum(query_pattern_cache_manager->misses);
        values[2] = Int32GetDatum(query_pattern_cache_manager->current_entries);
        values[3] = Int32GetDatum(query_pattern_cache_manager->max_entries);
    }
    else
    {
        /* Cache manager not initialized */
        values[0] = Int64GetDatum(0);  /* hits */
        values[1] = Int64GetDatum(0);  /* misses */
        values[2] = Int32GetDatum(0);  /* current_entries */
        values[3] = Int32GetDatum(0);  /* max_entries */
    }
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}

/*
 * Query pattern cache free function
 */
Datum
kmersearch_query_pattern_cache_free(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    /* Count entries before freeing */
    if (query_pattern_cache_manager)
        freed_entries += query_pattern_cache_manager->current_entries;
    
    /* Free query pattern cache manager (uses TopMemoryContext - needs manual cleanup) */
    /* DNA2/DNA4 cache managers are now local and automatically freed with QueryContext */
    free_query_pattern_cache_manager(&query_pattern_cache_manager);
    
    PG_RETURN_INT32(freed_entries);
}

/*
 * Actual min score cache statistics function
 */
Datum
kmersearch_actual_min_score_cache_stats(PG_FUNCTION_ARGS)
{
    TupleDesc tupdesc;
    Datum values[4];
    bool nulls[4] = {false};
    HeapTuple tuple;
    
    /* Build tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("function returning record called in context that cannot accept a set")));
    
    /* Actual min score cache statistics */
    if (actual_min_score_cache_manager)
    {
        values[0] = Int64GetDatum(actual_min_score_cache_manager->hits);
        values[1] = Int64GetDatum(actual_min_score_cache_manager->misses);
        values[2] = Int32GetDatum(actual_min_score_cache_manager->current_entries);
        values[3] = Int32GetDatum(actual_min_score_cache_manager->max_entries);
    }
    else
    {
        /* Cache manager not initialized */
        values[0] = Int64GetDatum(0);  /* hits */
        values[1] = Int64GetDatum(0);  /* misses */
        values[2] = Int32GetDatum(0);  /* current_entries */
        values[3] = Int32GetDatum(0);  /* max_entries */
    }
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}

/*
 * Free actual min score cache
 */
Datum
kmersearch_actual_min_score_cache_free(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    /* Count entries before freeing */
    if (actual_min_score_cache_manager)
        freed_entries = actual_min_score_cache_manager->current_entries;
    
    /* Free actual min score cache manager (uses TopMemoryContext - needs manual cleanup) */
    free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
    PG_RETURN_INT32(freed_entries);
}



/*
 * High-frequency k-mer cache management functions implementation
 */
void
kmersearch_highfreq_kmer_cache_init(void)
{
    MemoryContext old_context;
    
    /* Switch to TopMemoryContext */
    old_context = MemoryContextSwitchTo(TopMemoryContext);
    
    /* Initialize cache structure */
    memset(&global_highfreq_cache, 0, sizeof(HighfreqKmerCache));
    global_highfreq_cache.is_valid = false;
    memset(&global_highfreq_cache.current_cache_key, 0, sizeof(HighfreqCacheKey));
    global_highfreq_cache.current_cache_key.table_oid = InvalidOid;
    
    /* Create dedicated memory context for high-frequency k-mer cache */
    global_highfreq_cache.cache_context = AllocSetContextCreate(TopMemoryContext,
                                                                "HighfreqKmerCache",
                                                                ALLOCSET_DEFAULT_SIZES);
    
    global_highfreq_cache.highfreq_hash = NULL;
    global_highfreq_cache.highfreq_kmers = NULL;
    global_highfreq_cache.highfreq_count = 0;
    
    MemoryContextSwitchTo(old_context);
}

bool
kmersearch_highfreq_kmer_cache_load_internal(Oid table_oid, const char *column_name, int k_value)
{
    MemoryContext old_context;
    VarBit **highfreq_kmers;
    int highfreq_count;
    VarBit **cache_kmers;
    int i;
    HASHCTL hash_ctl;
    int total_inserted;
    int batch_num;
    int offset;
    
    if (!column_name || k_value <= 0) {
        return false;
    }
    
    /* Initialize cache if not already done */
    if (global_highfreq_cache.cache_context == NULL) {
        kmersearch_highfreq_kmer_cache_init();
    }
    
    /* Validate current GUC settings against metadata table */
    if (!kmersearch_validate_guc_against_metadata(table_oid, column_name, k_value)) {
        return false;
    }
    
    /* Clear existing cache if valid */
    if (global_highfreq_cache.is_valid) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Existing cache is valid, clearing")));
        kmersearch_highfreq_kmer_cache_free_internal();
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache clearing completed")));
    } else {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: No existing valid cache to clear")));
    }
    
    /* Reinitialize cache context if it was freed */
    if (global_highfreq_cache.cache_context == NULL) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache context was freed, reinitializing")));
        kmersearch_highfreq_kmer_cache_init();
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache reinitialization completed")));
    } else {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache context still exists")));
    }
    
    
    /* Count total k-mers first for hash table size initialization */
    {
        StringInfoData count_query;
        char *escaped_column_name;
        char *cp;
        const char *sp;
        int ret;
        
        /* Connect to SPI for counting */
        if (SPI_connect() != SPI_OK_CONNECT) {
            ereport(ERROR, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: SPI_connect failed for counting")));
        }
        
        /* Escape column name to prevent SQL injection */
        escaped_column_name = SPI_palloc(strlen(column_name) * 2 + 1);
        cp = escaped_column_name;
        sp = column_name;
        
        while (*sp) {
            if (*sp == '\'') {
                *cp++ = '\'';
                *cp++ = '\'';
            } else {
                *cp++ = *sp;
            }
            sp++;
        }
        *cp = '\0';
        
        /* Build count query */
        initStringInfo(&count_query);
        appendStringInfo(&count_query,
            "SELECT COUNT(DISTINCT hkm.kmer2_as_uint) FROM kmersearch_highfreq_kmer hkm "
            "WHERE hkm.table_oid = %u "
            "AND hkm.column_name = '%s' "
            "AND EXISTS ("
            "    SELECT 1 FROM kmersearch_highfreq_kmer_meta hkm_meta "
            "    WHERE hkm_meta.table_oid = %u "
            "    AND hkm_meta.column_name = '%s' "
            "    AND hkm_meta.kmer_size = %d"
            ")",
            table_oid, escaped_column_name, table_oid, escaped_column_name, k_value);
        
        /* Execute count query */
        ret = SPI_execute(count_query.data, true, 0);
        
        if (ret == SPI_OK_SELECT && SPI_processed > 0) {
            bool isnull;
            Datum count_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull) {
                highfreq_count = DatumGetInt64(count_datum);
            } else {
                highfreq_count = 0;
            }
        } else {
            highfreq_count = 0;
        }
        
        /* Cleanup */
        pfree(count_query.data);
        pfree(escaped_column_name);
        SPI_finish();
    }
    
    if (highfreq_count <= 0) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: No high-frequency k-mers found, cache remains invalid")));
        return false;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Found %d total high-frequency k-mers, will load in batches of %d", 
                           highfreq_count, kmersearch_highfreq_kmer_cache_load_batch_size)));
    
    /* Switch to cache context for cache storage */
    old_context = MemoryContextSwitchTo(global_highfreq_cache.cache_context);
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Switched to cache context successfully")));
    
    /* Store in cache - build cache key */
    global_highfreq_cache.current_cache_key.table_oid = table_oid;
    global_highfreq_cache.current_cache_key.column_name_hash = hash_any((unsigned char*)column_name, strlen(column_name));
    global_highfreq_cache.current_cache_key.kmer_size = k_value;
    global_highfreq_cache.current_cache_key.occur_bitlen = kmersearch_occur_bitlen;
    global_highfreq_cache.current_cache_key.max_appearance_rate = kmersearch_max_appearance_rate;
    global_highfreq_cache.current_cache_key.max_appearance_nrow = kmersearch_max_appearance_nrow;
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Building cache key and initializing hash table")));
    
    /* Initialize hash table in cache context */
    MemSet(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64);  /* Use hash value as key */
    hash_ctl.entrysize = sizeof(HighfreqKmerHashEntry);
    hash_ctl.hash = tag_hash;
    hash_ctl.hcxt = CurrentMemoryContext;
    
    global_highfreq_cache.highfreq_hash = hash_create("HighfreqKmerHash",
                                                      highfreq_count,
                                                      &hash_ctl,
                                                      HASH_ELEM | HASH_FUNCTION | HASH_CONTEXT);
    
    if (!global_highfreq_cache.highfreq_hash) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Hash table creation failed")));
        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(global_highfreq_cache.cache_context);
        global_highfreq_cache.cache_context = NULL;
        global_highfreq_cache.is_valid = false;
        return false;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Hash table created successfully, starting batch population")));
    
    /* Populate hash table with k-mers using batch processing */
    total_inserted = 0;
    batch_num = 0;
    offset = 0;
    
    /* Process k-mers in batches to reduce memory usage */
    while (offset < highfreq_count) {
        uint64 *batch_kmers;
        int batch_count;
        int current_batch_limit;
        
        /* Ensure we don't exceed the total count */
        if (offset + kmersearch_highfreq_kmer_cache_load_batch_size > highfreq_count) {
            current_batch_limit = highfreq_count - offset;
        } else {
            current_batch_limit = kmersearch_highfreq_kmer_cache_load_batch_size;
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Processing batch %d (offset=%d, batch_size=%d)", 
                               batch_num, offset, current_batch_limit)));
        
        /* Get batch of k-mers from table using LIMIT/OFFSET */
        {
            StringInfoData query;
            char *escaped_column_name;
            char *cp;
            const char *sp;
            int ret;
            
            /* Switch to parent context for SPI operations */
            MemoryContextSwitchTo(old_context);
            
            /* Connect to SPI */
            if (SPI_connect() != SPI_OK_CONNECT) {
                ereport(ERROR, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: SPI_connect failed for batch %d", batch_num)));
            }
            
            /* Escape column name to prevent SQL injection */
            escaped_column_name = SPI_palloc(strlen(column_name) * 2 + 1);
            cp = escaped_column_name;
            sp = column_name;
            
            while (*sp) {
                if (*sp == '\'') {
                    *cp++ = '\'';
                    *cp++ = '\'';
                } else {
                    *cp++ = *sp;
                }
                sp++;
            }
            *cp = '\0';
            
            /* Build query with LIMIT/OFFSET for batch processing */
            initStringInfo(&query);
            appendStringInfo(&query,
                "SELECT DISTINCT hkm.kmer2_as_uint FROM kmersearch_highfreq_kmer hkm "
                "WHERE hkm.table_oid = %u "
                "AND hkm.column_name = '%s' "
                "AND EXISTS ("
                "    SELECT 1 FROM kmersearch_highfreq_kmer_meta hkm_meta "
                "    WHERE hkm_meta.table_oid = %u "
                "    AND hkm_meta.column_name = '%s' "
                "    AND hkm_meta.kmer_size = %d"
                ") "
                "ORDER BY hkm.kmer2_as_uint "
                "LIMIT %d OFFSET %d",
                table_oid, escaped_column_name, table_oid, escaped_column_name, k_value, current_batch_limit, offset);
            
            /* Execute batch query */
            ret = SPI_execute(query.data, true, 0);
            
            if (ret == SPI_OK_SELECT && SPI_processed > 0) {
                batch_count = SPI_processed;
                batch_kmers = (uint64 *) palloc(batch_count * sizeof(uint64));
                
                for (i = 0; i < batch_count; i++) {
                    bool isnull;
                    Datum kmer_datum;
                    
                    kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
                    if (!isnull) {
                        /* Convert based on k_value size */
                        if (k_value <= 8) {
                            batch_kmers[i] = (uint64)DatumGetInt16(kmer_datum);
                        } else if (k_value <= 16) {
                            batch_kmers[i] = (uint64)DatumGetInt32(kmer_datum);
                        } else {
                            batch_kmers[i] = DatumGetInt64(kmer_datum);
                        }
                    } else {
                        /* This should never happen - kmersearch_highfreq_kmer should not contain NULL values */
                        ereport(ERROR, (errmsg("Unexpected NULL k-mer value in kmersearch_highfreq_kmer table")));
                    }
                }
            } else {
                batch_count = 0;
                batch_kmers = NULL;
            }
            
            /* Cleanup */
            pfree(query.data);
            pfree(escaped_column_name);
            SPI_finish();
            
            /* Switch back to cache context */
            MemoryContextSwitchTo(global_highfreq_cache.cache_context);
        }
        
        if (!batch_kmers || batch_count <= 0) {
            ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: No more k-mers in batch %d", batch_num)));
            break;
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Inserting %d k-mers from batch %d into hash table", 
                               batch_count, batch_num)));
        
        /* Insert batch k-mers into hash table */
        for (i = 0; i < batch_count; i++) {
            uint64 kmer2_as_uint;
            HighfreqKmerHashEntry *entry;
            bool found;
            
            /* All k-mer values are valid, including 0 (which represents "AAAA") */
            
            /* Use kmer2_as_uint value directly */
            kmer2_as_uint = batch_kmers[i];
            
            ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Inserting k-mer %d/%d from batch %d with kmer2_as_uint %lu", 
                                   i + 1, batch_count, batch_num, kmer2_as_uint)));
            
            /* Insert into hash table using kmer2_as_uint as both key and value */
            entry = (HighfreqKmerHashEntry *) hash_search(global_highfreq_cache.highfreq_hash,
                                                         (void *) &kmer2_as_uint,
                                                         HASH_ENTER,
                                                         &found);
            
            if (entry && !found) {
                entry->kmer_key = NULL;  /* No longer storing VarBit */
                entry->hash_value = kmer2_as_uint;
                total_inserted++;
                ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Successfully inserted k-mer %d (total: %d)", 
                                       i + 1, total_inserted)));
            } else if (found) {
                ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: K-mer %d already exists in hash table", i + 1)));
            } else {
                ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Failed to insert k-mer %d", i + 1)));
            }
        }
        
        /* Clean up batch memory in parent context */
        MemoryContextSwitchTo(old_context);
        if (batch_kmers) {
            /* uint64 arrays don't need individual pfree calls */
            pfree(batch_kmers);
        }
        MemoryContextSwitchTo(global_highfreq_cache.cache_context);
        
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Completed batch %d, total inserted so far: %d", 
                               batch_num, total_inserted)));
        
        /* Move to next batch */
        offset += batch_count;
        batch_num++;
        
        /* Break if we got fewer results than requested (end of data) */
        if (batch_count < current_batch_limit) {
            ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Reached end of data (batch_count=%d < batch_size=%d)", 
                                   batch_count, current_batch_limit)));
            break;
        }
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Completed all batches, total inserted: %d/%d", 
                           total_inserted, highfreq_count)));
    
    /* Set cache metadata */
    global_highfreq_cache.highfreq_kmers = NULL;  /* We don't store the array anymore */
    global_highfreq_cache.highfreq_count = total_inserted;
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Hash table population %s", 
                           global_highfreq_cache.highfreq_hash ? "successful" : "failed")));
    
    
    if (global_highfreq_cache.highfreq_hash) {
        global_highfreq_cache.is_valid = true;
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache loading successful, cache is now valid")));
    } else {
        /* Hash table creation failed, clean up by deleting the context */
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Hash table creation failed, cleaning up cache context")));
        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(global_highfreq_cache.cache_context);
        global_highfreq_cache.cache_context = NULL;
        global_highfreq_cache.is_valid = false;
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache cleanup completed, returning false")));
        return false;
    }
    
    MemoryContextSwitchTo(old_context);
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Function completed successfully, cache is valid=%s", 
                           global_highfreq_cache.is_valid ? "true" : "false")));
    
    return global_highfreq_cache.is_valid;
}

void
kmersearch_highfreq_kmer_cache_free_internal(void)
{
    if (!global_highfreq_cache.is_valid) {
        return;
    }
    
    /* Delete the entire cache context, which frees all allocated memory */
    if (global_highfreq_cache.cache_context) {
        MemoryContextDelete(global_highfreq_cache.cache_context);
        global_highfreq_cache.cache_context = NULL;
    }
    
    /* Reset cache state */
    global_highfreq_cache.is_valid = false;
    memset(&global_highfreq_cache.current_cache_key, 0, sizeof(HighfreqCacheKey));
    global_highfreq_cache.current_cache_key.table_oid = InvalidOid;
    global_highfreq_cache.highfreq_count = 0;
    global_highfreq_cache.highfreq_hash = NULL;
    global_highfreq_cache.highfreq_kmers = NULL;
}

bool
kmersearch_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value)
{
    /* Build expected cache key */
    HighfreqCacheKey expected_key;
    expected_key.table_oid = table_oid;
    expected_key.column_name_hash = hash_any((unsigned char*)column_name, strlen(column_name));
    expected_key.kmer_size = k_value;
    expected_key.occur_bitlen = kmersearch_occur_bitlen;
    expected_key.max_appearance_rate = kmersearch_max_appearance_rate;
    expected_key.max_appearance_nrow = kmersearch_max_appearance_nrow;
    
    /* Compare cache keys */
    return (global_highfreq_cache.is_valid &&
            memcmp(&global_highfreq_cache.current_cache_key, &expected_key, sizeof(HighfreqCacheKey)) == 0);
}

/*
 * Check if global_highfreq_cache is loaded
 */
bool
kmersearch_is_global_highfreq_cache_loaded(void)
{
    return (global_highfreq_cache.is_valid && 
            global_highfreq_cache.highfreq_count > 0);
}

/*
 * Validate that the cache key matches the specified table and column
 */
bool
kmersearch_validate_cache_key_match(Oid table_oid, const char *column_name)
{
    HighfreqCacheKey expected_key;
    bool matches;
    
    if (!global_highfreq_cache.is_valid) {
        ereport(DEBUG1, (errmsg("kmersearch_validate_cache_key_match: global cache is not valid")));
        return false;
    }
    
    /* Build expected cache key using current GUC values */
    expected_key.table_oid = table_oid;
    expected_key.column_name_hash = hash_any((unsigned char*)column_name, strlen(column_name));
    expected_key.kmer_size = kmersearch_kmer_size;
    expected_key.occur_bitlen = kmersearch_occur_bitlen;
    expected_key.max_appearance_rate = kmersearch_max_appearance_rate;
    expected_key.max_appearance_nrow = kmersearch_max_appearance_nrow;
    
    /* Compare with current cache key */
    matches = (memcmp(&global_highfreq_cache.current_cache_key, &expected_key, sizeof(HighfreqCacheKey)) == 0);
    
    if (!matches) {
        ereport(DEBUG1, (errmsg("kmersearch_validate_cache_key_match: cache key mismatch - "
                               "expected table_oid=%u, kmer_size=%d, occur_bitlen=%d, max_appearance_rate=%f, max_appearance_nrow=%d "
                               "but cache has table_oid=%u, kmer_size=%d, occur_bitlen=%d, max_appearance_rate=%f, max_appearance_nrow=%d",
                               expected_key.table_oid, expected_key.kmer_size, expected_key.occur_bitlen,
                               expected_key.max_appearance_rate, expected_key.max_appearance_nrow,
                               global_highfreq_cache.current_cache_key.table_oid,
                               global_highfreq_cache.current_cache_key.kmer_size,
                               global_highfreq_cache.current_cache_key.occur_bitlen,
                               global_highfreq_cache.current_cache_key.max_appearance_rate,
                               global_highfreq_cache.current_cache_key.max_appearance_nrow)));
    }
    
    return matches;
}

/*
 * Validate that the parallel cache key matches the specified table and column
 */
bool
kmersearch_validate_parallel_cache_key_match(Oid table_oid, const char *column_name)
{
    HighfreqCacheKey expected_key;
    bool matches;
    
    if (parallel_highfreq_cache == NULL || !parallel_highfreq_cache->is_initialized) {
        ereport(DEBUG1, (errmsg("kmersearch_validate_parallel_cache_key_match: parallel cache is not initialized")));
        return false;
    }
    
    /* Build expected cache key using current GUC values */
    expected_key.table_oid = table_oid;
    expected_key.column_name_hash = hash_any((unsigned char*)column_name, strlen(column_name));
    expected_key.kmer_size = kmersearch_kmer_size;
    expected_key.occur_bitlen = kmersearch_occur_bitlen;
    expected_key.max_appearance_rate = kmersearch_max_appearance_rate;
    expected_key.max_appearance_nrow = kmersearch_max_appearance_nrow;
    
    /* Compare with current parallel cache key */
    matches = (memcmp(&parallel_highfreq_cache->cache_key, &expected_key, sizeof(HighfreqCacheKey)) == 0);
    
    if (!matches) {
        ereport(DEBUG1, (errmsg("kmersearch_validate_parallel_cache_key_match: cache key mismatch - "
                               "expected table_oid=%u, kmer_size=%d, occur_bitlen=%d, max_appearance_rate=%f, max_appearance_nrow=%d "
                               "but parallel cache has table_oid=%u, kmer_size=%d, occur_bitlen=%d, max_appearance_rate=%f, max_appearance_nrow=%d",
                               expected_key.table_oid, expected_key.kmer_size, expected_key.occur_bitlen,
                               expected_key.max_appearance_rate, expected_key.max_appearance_nrow,
                               parallel_highfreq_cache->cache_key.table_oid,
                               parallel_highfreq_cache->cache_key.kmer_size,
                               parallel_highfreq_cache->cache_key.occur_bitlen,
                               parallel_highfreq_cache->cache_key.max_appearance_rate,
                               parallel_highfreq_cache->cache_key.max_appearance_nrow)));
    }
    
    return matches;
}

/*
 * Lookup k-mer in global_highfreq_cache
 */
bool
kmersearch_lookup_in_global_cache(VarBit *kmer_key)
{
    uint64 hash = kmersearch_ngram_key_to_hash(kmer_key);
    bool found;
    
    if (!global_highfreq_cache.is_valid || global_highfreq_cache.highfreq_count == 0)
        return false;
    
    /* Search for the hash in the global high-frequency k-mer cache */
    hash_search(global_highfreq_cache.highfreq_hash, &hash, HASH_FIND, &found);
    
    return found;
}

/*
 * SQL-accessible high-frequency cache load function
 */
Datum
kmersearch_highfreq_kmer_cache_load(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_P(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    
    char *table_name = text_to_cstring(table_name_text);
    char *column_name = text_to_cstring(column_name_text);
    bool success;
    Oid table_oid;
    
    /* Get table OID from table name */
    table_oid = RelnameGetRelid(table_name);
    if (!OidIsValid(table_oid))
    {
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("relation \"%s\" does not exist", table_name)));
    }
    
    success = kmersearch_highfreq_kmer_cache_load_internal(table_oid, column_name, kmersearch_kmer_size);
    
    PG_RETURN_BOOL(success);
}

/*
 * SQL-accessible high-frequency cache free function
 */
Datum
kmersearch_highfreq_kmer_cache_free(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_P(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    
    char *table_name = text_to_cstring(table_name_text);
    char *column_name = text_to_cstring(column_name_text);
    int freed_entries = 0;
    
    /* Get table OID from table name */
    Oid table_oid = RelnameGetRelid(table_name);
    if (!OidIsValid(table_oid))
    {
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("relation \"%s\" does not exist", table_name)));
    }
    
    /* Validate cache key matches table/column before freeing */
    if (!kmersearch_validate_cache_key_match(table_oid, column_name)) {
        ereport(WARNING,
                (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                 errmsg("cache key mismatch for table \"%s\" column \"%s\"", table_name, column_name),
                 errhint("The cache was not loaded for this table/column combination, or was loaded with different parameters.")));
        PG_RETURN_INT32(0);
    }
    
    /* Count entries before freeing */
    if (global_highfreq_cache.is_valid)
        freed_entries = global_highfreq_cache.highfreq_count;
    
    ereport(LOG, (errmsg("kmersearch_highfreq_kmer_cache_free: freeing %d entries for table \"%s\" column \"%s\"", 
                         freed_entries, table_name, column_name)));
    
    /* Free the cache */
    kmersearch_highfreq_kmer_cache_free_internal();
    
    PG_RETURN_INT32(freed_entries);
}

/*
 * SQL-accessible high-frequency cache free function without parameters
 * For backwards compatibility with test cases
 */
Datum
kmersearch_highfreq_kmer_cache_free_all(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    if (global_highfreq_cache.is_valid)
        freed_entries = 1;
    
    /* Free the cache */
    kmersearch_highfreq_kmer_cache_free_internal();
    
    PG_RETURN_INT32(freed_entries);
}

/*
 * Validate current GUC settings against metadata table values
 * Returns true if validation passes, false otherwise
 */
bool
kmersearch_validate_guc_against_metadata(Oid table_oid, const char *column_name, int k_value)
{
    int ret;
    StringInfoData query;
    bool validation_passed = true;
    /* Connect to SPI */
    if (SPI_connect() != SPI_OK_CONNECT) {
        ereport(ERROR, (errmsg("kmersearch_validate_guc_against_metadata: SPI_connect failed")));
    }
    
    /* Build query to get metadata from high-frequency k-mer metadata table */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT kmer_size, occur_bitlen, max_appearance_rate, max_appearance_nrow "
        "FROM kmersearch_highfreq_kmer_meta "
        "WHERE table_oid = %u AND column_name = '%s' AND kmer_size = %d",
        table_oid, column_name, k_value);
    
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        bool isnull;
        Datum kmer_size_datum, occur_bitlen_datum, max_appearance_rate_datum, max_appearance_nrow_datum;
        int stored_kmer_size, stored_occur_bitlen, stored_max_appearance_nrow;
        float stored_max_appearance_rate;
        
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Found high-frequency k-mer metadata record, extracting values")));
        
        /* Get and validate kmer_size */
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Getting kmer_size value")));
        kmer_size_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull)
        {
            stored_kmer_size = DatumGetInt32(kmer_size_datum);
            ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: stored_kmer_size=%d, current=%d",
                                   stored_kmer_size, k_value)));
            if (stored_kmer_size != k_value)
            {
                ereport(ERROR,
                        (errcode(ERRCODE_CONFIG_FILE_ERROR),
                         errmsg("GUC validation failed: kmersearch.kmer_size mismatch"),
                         errdetail("Current setting: %d, Required by metadata: %d",
                                 k_value, stored_kmer_size),
                         errhint("Set kmersearch.kmer_size = %d to match the metadata configuration.",
                                stored_kmer_size)));
                validation_passed = false;
            }
        }
        
        /* Get stored metadata values */
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Getting occur_bitlen value")));
        occur_bitlen_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull);
        if (!isnull)
        {
            stored_occur_bitlen = DatumGetInt32(occur_bitlen_datum);
            ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: stored_occur_bitlen=%d, current=%d",
                                   stored_occur_bitlen, kmersearch_occur_bitlen)));
            if (stored_occur_bitlen != kmersearch_occur_bitlen)
            {
                ereport(ERROR,
                        (errcode(ERRCODE_CONFIG_FILE_ERROR),
                         errmsg("GUC validation failed: kmersearch.occur_bitlen mismatch"),
                         errdetail("Current setting: %d, Required by metadata: %d",
                                 kmersearch_occur_bitlen, stored_occur_bitlen),
                         errhint("Set kmersearch.occur_bitlen = %d before loading cache.",
                                stored_occur_bitlen)));
                validation_passed = false;
            }
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Getting max_appearance_rate value")));
        max_appearance_rate_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull);
        if (!isnull)
        {
            stored_max_appearance_rate = DatumGetFloat4(max_appearance_rate_datum);
            ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: stored_max_appearance_rate=%.4f, current=%.4f",
                                   stored_max_appearance_rate, kmersearch_max_appearance_rate)));
            if (fabs(stored_max_appearance_rate - kmersearch_max_appearance_rate) > 0.0001)
            {
                ereport(ERROR,
                        (errcode(ERRCODE_CONFIG_FILE_ERROR),
                         errmsg("GUC validation failed: kmersearch.max_appearance_rate mismatch"),
                         errdetail("Current setting: %.4f, Required by metadata: %.4f",
                                 kmersearch_max_appearance_rate, stored_max_appearance_rate),
                         errhint("Set kmersearch.max_appearance_rate = %.4f before loading cache.",
                                stored_max_appearance_rate)));
                validation_passed = false;
            }
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Getting max_appearance_nrow value")));
        max_appearance_nrow_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 4, &isnull);
        if (!isnull)
        {
            stored_max_appearance_nrow = DatumGetInt32(max_appearance_nrow_datum);
            ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: stored_max_appearance_nrow=%d, current=%d",
                                   stored_max_appearance_nrow, kmersearch_max_appearance_nrow)));
            if (stored_max_appearance_nrow != kmersearch_max_appearance_nrow)
            {
                ereport(ERROR,
                        (errcode(ERRCODE_CONFIG_FILE_ERROR),
                         errmsg("GUC validation failed: kmersearch.max_appearance_nrow mismatch"),
                         errdetail("Current setting: %d, Required by metadata: %d",
                                 kmersearch_max_appearance_nrow, stored_max_appearance_nrow),
                         errhint("Set kmersearch.max_appearance_nrow = %d before loading cache.",
                                stored_max_appearance_nrow)));
                validation_passed = false;
            }
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: All metadata values validated successfully")));
    }
    else
    {
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: No metadata found or query failed")));
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("No metadata found for table_oid=%u, column_name='%s', kmer_size=%d",
                       table_oid, column_name, k_value),
                 errhint("Run kmersearch_perform_highfreq_analysis() first to create metadata.")));
        validation_passed = false;
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    
    return validation_passed;
}

HTAB *
kmersearch_create_highfreq_hash_from_array(VarBit **kmers, int nkeys)
{
    HTAB *hash_table;
    HASHCTL hash_ctl;
    int i;
    
    
    if (!kmers || nkeys <= 0) {
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Invalid parameters, kmers=%p, nkeys=%d",
                               kmers, nkeys)));
        return NULL;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Parameters validated, setting up hash table")));
    
    /* Set up hash table using hash value as key */
    MemSet(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64);  /* Use hash value as key */
    hash_ctl.entrysize = sizeof(HighfreqKmerHashEntry);
    hash_ctl.hash = tag_hash;
    hash_ctl.hcxt = CurrentMemoryContext;
    
    
    hash_table = hash_create("HighfreqKmerHash",
                            nkeys,
                            &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_CONTEXT);
    
    if (!hash_table) {
        return NULL;
    }
    
    
    /* Add each k-mer to the hash table */
    for (i = 0; i < nkeys; i++)
    {
        HighfreqKmerHashEntry *entry;
        bool found;
        uint64 hash_value;
        unsigned char *bits_ptr;
        int bytes_len;
        
        
        /* Check if the k-mer pointer itself is valid before any access */
        
        if (!kmers[i]) {
            continue;
        }
        
        
        
        /* Validate VarBit data before hash calculation */
        if (VARSIZE(kmers[i]) < VARHDRSZ) {
            continue;
        }
        
        
        /* Check VARBITS pointer */
        bits_ptr = (unsigned char *) VARBITS(kmers[i]);
        if (!bits_ptr) {
            continue;
        }
        
        
        /* Calculate bytes length manually for ngram_key2 (more reliable than VARBITBYTES) */
        {
            int bit_length = VARBITLEN(kmers[i]);
            bytes_len = (bit_length + 7) / 8;  /* Round up to next byte */
        }
        
        if (bytes_len <= 0 || bytes_len > 1000) {  /* Reasonable upper limit */
            continue;
        }
        
        
        /* Calculate hash value for this ngram_key2 (kmer2 + occurrence bits) */
        hash_value = DatumGetUInt64(hash_any(bits_ptr, bytes_len));
        
        
        entry = (HighfreqKmerHashEntry *) hash_search(hash_table,
                                                     (void *) &hash_value,
                                                     HASH_ENTER,
                                                     &found);
        
        
        if (entry && !found)
        {
            entry->kmer_key = kmers[i];
            entry->hash_value = hash_value;
        }
    }
    
    
    return hash_table;
}

/*
 * Convert a VarBit k-mer key to hash value
 * This is used for high-frequency k-mer cache lookups
 */
uint64
kmersearch_ngram_key_to_hash(VarBit *ngram_key)
{
    if (!ngram_key)
        return 0;
    
    /* Hash the VarBit content */
    return DatumGetUInt64(hash_any((unsigned char *) VARBITS(ngram_key), VARBITBYTES(ngram_key)));
}

/*
 * Parallel high-frequency k-mer cache load function
 */
Datum
kmersearch_parallel_highfreq_kmer_cache_load(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_P(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    
    char *table_name = text_to_cstring(table_name_text);
    char *column_name = text_to_cstring(column_name_text);
    bool result;
    Oid table_oid;
    
    /* Get table OID from table name */
    table_oid = RelnameGetRelid(table_name);
    if (!OidIsValid(table_oid))
    {
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("relation \"%s\" does not exist", table_name)));
    }
    
    /* Initialize parallel cache if not already done */
    if (parallel_highfreq_cache == NULL) {
        kmersearch_parallel_highfreq_kmer_cache_init();
    }
    
    /* Load cache data into DSM */
    result = kmersearch_parallel_highfreq_kmer_cache_load_internal(table_oid, column_name, kmersearch_kmer_size);
    
    /* Register cleanup function for process exit */
    if (result) {
        on_proc_exit(kmersearch_parallel_cache_cleanup_on_exit, 0);
    }
    
    pfree(column_name);
    PG_RETURN_BOOL(result);
}

/*
 * Parallel high-frequency k-mer cache free function
 */
Datum
kmersearch_parallel_highfreq_kmer_cache_free(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_P(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    
    char *table_name = text_to_cstring(table_name_text);
    char *column_name = text_to_cstring(column_name_text);
    int32 freed_entries = 0;
    
    /* Get table OID from table name */
    Oid table_oid = RelnameGetRelid(table_name);
    if (!OidIsValid(table_oid))
    {
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("relation \"%s\" does not exist", table_name)));
    }
    
    ereport(LOG, (errmsg("kmersearch_parallel_highfreq_kmer_cache_free: Starting function call for table %s, column %s", table_name, column_name)));
    
    /* Validate cache key matches table/column before freeing */
    if (!kmersearch_validate_parallel_cache_key_match(table_oid, column_name)) {
        ereport(WARNING,
                (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                 errmsg("parallel cache key mismatch for table \"%s\" column \"%s\"", table_name, column_name),
                 errhint("The parallel cache was not loaded for this table/column combination, or was loaded with different parameters.")));
        PG_RETURN_INT32(0);
    }
    
    /* Get the actual number of entries from the cache */
    if (parallel_highfreq_cache != NULL && parallel_highfreq_cache->is_initialized) {
        freed_entries = parallel_highfreq_cache->num_entries;
        ereport(LOG, (errmsg("kmersearch_parallel_highfreq_kmer_cache_free: Found %d entries to free", freed_entries)));
    } else {
        ereport(LOG, (errmsg("kmersearch_parallel_highfreq_kmer_cache_free: parallel_highfreq_cache is NULL or not initialized")));
        freed_entries = 0;
    }
    
    /* Free parallel cache */
    kmersearch_parallel_highfreq_kmer_cache_free_internal();
    
    
    PG_RETURN_INT32(freed_entries);
}

/*
 * SQL-accessible parallel high-frequency cache free function without parameters
 * For backwards compatibility with test cases
 */
Datum
kmersearch_parallel_highfreq_kmer_cache_free_all(PG_FUNCTION_ARGS)
{
    int32 freed_entries = 0;
    
    if (parallel_highfreq_cache != NULL && parallel_highfreq_cache->is_initialized) {
        freed_entries = parallel_highfreq_cache->num_entries;
    }
    
    /* Free parallel cache */
    kmersearch_parallel_highfreq_kmer_cache_free_internal();
    
    PG_RETURN_INT32(freed_entries);
}

/*
 * GUC hook function for query pattern cache max entries changes
 */
void
kmersearch_query_pattern_cache_max_entries_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Clear query pattern cache to recreate with new size limit */
    if (query_pattern_cache_manager)
        free_query_pattern_cache_manager(&query_pattern_cache_manager);
}

/* Additional global variable for parallel cache exit callback */
bool parallel_cache_exit_callback_registered = false;

/*
 * Internal cleanup function for parallel cache resources
 */
static void
kmersearch_parallel_cache_cleanup_internal(void)
{
    MemoryContext oldcontext;
    
    ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Starting cleanup")));
    
    /* Prevent double cleanup - check if already cleaned up */
    if (parallel_cache_hash == NULL && parallel_cache_dsa == NULL && parallel_cache_segment == NULL) {
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Already cleaned up, skipping")));
        return;
    }
    
    /* Use TopMemoryContext for dshash operations */
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    /* Step 1: Handle dshash table cleanup based on process type */
    if (parallel_cache_hash != NULL) {
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: parallel_cache_hash=%p", parallel_cache_hash)));
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: parallel_cache_dsa=%p", parallel_cache_dsa)));
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: parallel_cache_segment=%p", parallel_cache_segment)));
        
        /* Check if DSA and DSM are still valid before operating on dshash */
        if (parallel_cache_dsa != NULL && parallel_cache_segment != NULL) {
            if (!IsParallelWorker()) {
                /* Main process: destroy the hash table */
                ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Main process destroying dshash table")));
                dshash_destroy(parallel_cache_hash);
            } else {
                /* Parallel worker: detach from hash table */
                ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Parallel worker detaching from dshash table")));
                dshash_detach(parallel_cache_hash);
            }
        } else {
            ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: DSA or DSM already invalid, just detaching")));
            /* DSA/DSM already destroyed, just detach without destroy */
            dshash_detach(parallel_cache_hash);
        }
        parallel_cache_hash = NULL;
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: dshash table cleanup completed")));
    } else {
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: No dshash table to cleanup")));
    }
    
    /* Switch back to original context before DSA/DSM operations */
    MemoryContextSwitchTo(oldcontext);
    
    /* Step 2: Handle DSA area cleanup */
    if (parallel_cache_dsa != NULL) {
        if (!IsParallelWorker()) {
            /* Main process: unpin and detach */
            ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Main process unpinning and detaching DSA area")));
            PG_TRY();
            {
                /* Unpin DSA area */
                dsa_unpin(parallel_cache_dsa);
                /* Detach from DSA area */
                dsa_detach(parallel_cache_dsa);
                ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: DSA area unpinned and detached successfully")));
            }
            PG_CATCH();
            {
                ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: DSA area cleanup failed, but continuing")));
                FlushErrorState();
            }
            PG_END_TRY();
        } else {
            /* Parallel worker: detach only */
            ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Parallel worker detaching from DSA area")));
            dsa_detach(parallel_cache_dsa);
            ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: DSA area detached successfully")));
        }
        parallel_cache_dsa = NULL;
    } else {
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: No DSA area to cleanup")));
    }
    
    /* Step 3: Handle DSM segment properly */
    if (parallel_cache_segment != NULL) {
        if (!IsParallelWorker()) {
            /* Main process: unpin and detach */
            ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Main process unpinning and detaching DSM segment")));
            PG_TRY();
            {
                /* Get DSM handle before detaching */
                dsm_handle handle = dsm_segment_handle(parallel_cache_segment);
                /* Unpin mapping first */
                dsm_unpin_mapping(parallel_cache_segment);
                /* Detach from the DSM segment */
                dsm_detach(parallel_cache_segment);
                /* Unpin the DSM segment to allow cleanup */
                dsm_unpin_segment(handle);
                ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: DSM segment detached and unpinned successfully")));
            }
            PG_CATCH();
            {
                ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: DSM segment cleanup failed, but continuing")));
                FlushErrorState();
            }
            PG_END_TRY();
        } else {
            /* Parallel worker: detach only */
            ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: Parallel worker detaching from DSM segment")));
            dsm_detach(parallel_cache_segment);
            ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: DSM segment detached successfully")));
        }
        parallel_cache_segment = NULL;
    } else {
        ereport(LOG, (errmsg("kmersearch_parallel_cache_cleanup_internal: No DSM segment to handle")));
    }
    
    /* Reset cache pointer and callback flag */
    parallel_highfreq_cache = NULL;
    parallel_cache_exit_callback_registered = false;
}

/*
 * Exit callback for DSM cleanup
 */
static void
dshash_cache_cleanup_callback(int code, Datum arg)
{
    ereport(LOG, (errmsg("dshash_cache_cleanup_callback: Starting cleanup on process exit")));
    
    kmersearch_parallel_cache_cleanup_internal();
    
    ereport(LOG, (errmsg("dshash_cache_cleanup_callback: Cleanup completed")));
}

/*
 * Initialize parallel high-frequency k-mer cache
 */
void
kmersearch_parallel_highfreq_kmer_cache_init(void)
{
    /* Initialize parallel cache state */
    parallel_highfreq_cache = NULL;
    parallel_cache_segment = NULL;
    parallel_cache_hash = NULL;
}

/*
 * Load data into parallel high-frequency k-mer cache
 */
bool
kmersearch_parallel_highfreq_kmer_cache_load_internal(Oid table_oid, const char *column_name, int k_value)
{
    dshash_parameters params;
    Size segment_size;
    ParallelHighfreqKmerCacheEntry *entry;
    bool found;
    int i;
    Size cache_struct_size;
    Size entries_size;
    Size dsa_min_size;
    Size dshash_overhead;
    char *dsa_start;
    Size dsa_size;
    MemoryContext oldcontext;
    int total_kmer_count = 0;
    int total_inserted = 0;
    int batch_num = 0;
    int offset = 0;
    VarBit **count_kmers;
    
    if (!column_name || k_value <= 0)
        return false;
    
    /* Validate current GUC settings against metadata table */
    if (!kmersearch_validate_guc_against_metadata(table_oid, column_name, k_value))
        return false;
    
    /* Check if cache is already loaded for this table */
    if (parallel_cache_segment != NULL && 
        parallel_cache_dsa != NULL && 
        parallel_cache_hash != NULL &&
        parallel_highfreq_cache != NULL) {
        
        /* Verify cache data is still valid */
        if (parallel_highfreq_cache->is_initialized &&
            parallel_highfreq_cache->cache_key.table_oid == table_oid &&
            parallel_highfreq_cache->cache_key.kmer_size == k_value) {
            ereport(LOG, (errmsg("dshash_cache_load: Cache already loaded for table %u, k=%d", 
                                table_oid, k_value)));
            return true;
        } else {
            /* Cache exists but is for different table/k_value, clean it up */
            ereport(LOG, (errmsg("dshash_cache_load: Cache exists but for different table/k_value, cleaning up")));
            kmersearch_parallel_highfreq_kmer_cache_free_internal();
        }
    }
    
    /* Count total k-mers first for DSM segment size calculation */
    {
        StringInfoData count_query;
        char *escaped_column_name;
        char *cp;
        const char *sp;
        int ret;
        
        /* Connect to SPI for counting */
        if (SPI_connect() != SPI_OK_CONNECT) {
            ereport(ERROR, (errmsg("kmersearch_parallel_highfreq_kmer_cache_load_internal: SPI_connect failed for counting")));
        }
        
        /* Escape column name to prevent SQL injection */
        escaped_column_name = SPI_palloc(strlen(column_name) * 2 + 1);
        cp = escaped_column_name;
        sp = column_name;
        
        while (*sp) {
            if (*sp == '\'') {
                *cp++ = '\'';
                *cp++ = '\'';
            } else {
                *cp++ = *sp;
            }
            sp++;
        }
        *cp = '\0';
        
        /* Build count query */
        initStringInfo(&count_query);
        appendStringInfo(&count_query,
            "SELECT COUNT(DISTINCT hkm.kmer2_as_uint) FROM kmersearch_highfreq_kmer hkm "
            "WHERE hkm.table_oid = %u "
            "AND hkm.column_name = '%s' "
            "AND EXISTS ("
            "    SELECT 1 FROM kmersearch_highfreq_kmer_meta hkm_meta "
            "    WHERE hkm_meta.table_oid = %u "
            "    AND hkm_meta.column_name = '%s' "
            "    AND hkm_meta.kmer_size = %d"
            ")",
            table_oid, escaped_column_name, table_oid, escaped_column_name, k_value);
        
        /* Execute count query */
        ret = SPI_execute(count_query.data, true, 0);
        
        if (ret == SPI_OK_SELECT && SPI_processed > 0) {
            bool isnull;
            Datum count_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull) {
                total_kmer_count = DatumGetInt64(count_datum);
            } else {
                total_kmer_count = 0;
            }
        } else {
            total_kmer_count = 0;
        }
        
        /* Cleanup */
        pfree(count_query.data);
        pfree(escaped_column_name);
        SPI_finish();
    }
    
    if (total_kmer_count <= 0) {
        /* No high-frequency k-mers found */
        ereport(LOG, (errmsg("dshash_cache_load: No high-frequency k-mers found")));
        return false;
    }
    
    ereport(LOG, (errmsg("dshash_cache_load: Found %d total high-frequency k-mers, will load in batches of %d", 
                        total_kmer_count, kmersearch_highfreq_kmer_cache_load_batch_size)));
    
    /* Calculate required segment size using total count */
    
    cache_struct_size = MAXALIGN(sizeof(ParallelHighfreqKmerCache));
    entries_size = total_kmer_count * sizeof(ParallelHighfreqKmerCacheEntry);
    dsa_min_size = 8192; /* Minimum DSA area size */
    dshash_overhead = MAXALIGN(512); /* Extra space for dshash overhead */
    
    /* Total size = cache structure + DSA area (at least 8192) + entries + overhead */
    segment_size = cache_struct_size + dsa_min_size + entries_size + dshash_overhead;
    
    /* Ensure total size is reasonable */
    if (segment_size < 16384) {
        segment_size = 16384; /* At least 16KB total */
    }
    
    /* Create DSM segment */
    ereport(LOG, (errmsg("dshash_cache_load: Creating DSM segment of size %zu", segment_size)));
    
    parallel_cache_segment = dsm_create(segment_size, 0);
    if (!parallel_cache_segment) {
        ereport(ERROR,
                (errcode(ERRCODE_OUT_OF_MEMORY),
                 errmsg("failed to create DSM segment for parallel cache")));
    }
    
    /* Pin the DSM segment to prevent automatic cleanup when query ends */
    dsm_pin_segment(parallel_cache_segment);
    dsm_pin_mapping(parallel_cache_segment);
    
    ereport(LOG, (errmsg("dshash_cache_load: DSM segment created successfully")));
    
    /* Initialize parallel cache structure in DSM */
    parallel_highfreq_cache = (ParallelHighfreqKmerCache *) dsm_segment_address(parallel_cache_segment);
    
    /* Initialize complete cache key */
    parallel_highfreq_cache->cache_key.table_oid = table_oid;
    parallel_highfreq_cache->cache_key.column_name_hash = hash_any((unsigned char*)column_name, strlen(column_name));
    parallel_highfreq_cache->cache_key.kmer_size = k_value;
    parallel_highfreq_cache->cache_key.occur_bitlen = kmersearch_occur_bitlen;
    parallel_highfreq_cache->cache_key.max_appearance_rate = kmersearch_max_appearance_rate;
    parallel_highfreq_cache->cache_key.max_appearance_nrow = kmersearch_max_appearance_nrow;
    
    /* num_entries will be set after batch processing */
    parallel_highfreq_cache->segment_size = segment_size;
    parallel_highfreq_cache->dsm_handle = dsm_segment_handle(parallel_cache_segment);
    
    /* Create DSA area from DSM segment */
    /* Skip the ParallelHighfreqKmerCache structure to avoid overlap */
    
    dsa_start = (char *) dsm_segment_address(parallel_cache_segment) + cache_struct_size;
    dsa_size = segment_size - cache_struct_size;
    
    /* Ensure DSA area is at least 8192 bytes */
    if (dsa_size < 8192) {
        ereport(ERROR,
                (errcode(ERRCODE_OUT_OF_MEMORY),
                 errmsg("DSA area size %zu is too small, need at least 8192 bytes", dsa_size)));
    }
    
    /* Use TopMemoryContext for persistent dshash objects */
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    parallel_cache_dsa = dsa_create_in_place(dsa_start,
                                             dsa_size,
                                             LWTRANCHE_PARALLEL_QUERY_DSA,
                                             parallel_cache_segment);
    
    if (!parallel_cache_dsa) {
        MemoryContextSwitchTo(oldcontext);
        dsm_detach(parallel_cache_segment);
        parallel_cache_segment = NULL;
        ereport(ERROR,
                (errcode(ERRCODE_OUT_OF_MEMORY),
                 errmsg("failed to create DSA area for parallel cache")));
    }
    
    /* Pin the DSA area to prevent unexpected cleanup */
    dsa_pin(parallel_cache_dsa);
    dsa_pin_mapping(parallel_cache_dsa);
    
    /* Set up dshash parameters based on k-mer size */
    memset(&params, 0, sizeof(params));
    params.compare_function = dshash_memcmp;
    params.tranche_id = LWTRANCHE_KMERSEARCH_CACHE;
    
    if (k_value <= 8) {
        params.key_size = sizeof(uint16);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry16);
        params.hash_function = kmersearch_uint16_identity_hash;
    } else if (k_value <= 16) {
        params.key_size = sizeof(uint32);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry32);
        params.hash_function = kmersearch_uint32_identity_hash;
    } else {
        params.key_size = sizeof(uint64);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry64);
        params.hash_function = dshash_memhash;
    }
    
    /* Create dshash table */
    ereport(LOG, (errmsg("dshash_cache_load: Creating dshash table with %d entries", total_kmer_count)));
    
    parallel_cache_hash = dshash_create(parallel_cache_dsa, &params, NULL);
    if (!parallel_cache_hash) {
        dsa_detach(parallel_cache_dsa);
        parallel_cache_dsa = NULL;
        dsm_detach(parallel_cache_segment);
        parallel_cache_segment = NULL;
        MemoryContextSwitchTo(oldcontext);
        ereport(ERROR,
                (errcode(ERRCODE_OUT_OF_MEMORY),
                 errmsg("failed to create dshash table for parallel cache")));
    }
    
    ereport(LOG, (errmsg("dshash_cache_load: dshash table created successfully")));
    
    /* Store the dshash table handle */
    parallel_highfreq_cache->hash_handle = dshash_get_hash_table_handle(parallel_cache_hash);
    parallel_highfreq_cache->is_initialized = true;
    
    /* Populate the hash table with high-frequency k-mers using batch processing */
    ereport(LOG, (errmsg("dshash_cache_load: Starting batch population of %d k-mers", total_kmer_count)));
    
    /* Process k-mers in batches to reduce memory usage */
    while (offset < total_kmer_count) {
        uint64 *batch_kmers;
        int batch_count;
        int current_batch_limit;
        
        /* Ensure we don't exceed the total count */
        if (offset + kmersearch_highfreq_kmer_cache_load_batch_size > total_kmer_count) {
            current_batch_limit = total_kmer_count - offset;
        } else {
            current_batch_limit = kmersearch_highfreq_kmer_cache_load_batch_size;
        }
        
        ereport(LOG, (errmsg("dshash_cache_load: Processing batch %d (offset=%d, batch_size=%d)", 
                             batch_num, offset, current_batch_limit)));
        
        /* Get batch of k-mers from table using LIMIT/OFFSET */
        {
            StringInfoData query;
            char *escaped_column_name;
            char *cp;
            const char *sp;
            int ret;
            
            /* Connect to SPI */
            if (SPI_connect() != SPI_OK_CONNECT) {
                ereport(ERROR, (errmsg("dshash_cache_load: SPI_connect failed for batch %d", batch_num)));
            }
            
            /* Escape column name to prevent SQL injection */
            escaped_column_name = SPI_palloc(strlen(column_name) * 2 + 1);
            cp = escaped_column_name;
            sp = column_name;
            
            while (*sp) {
                if (*sp == '\'') {
                    *cp++ = '\'';
                    *cp++ = '\'';
                } else {
                    *cp++ = *sp;
                }
                sp++;
            }
            *cp = '\0';
            
            /* Build query with LIMIT/OFFSET for batch processing */
            initStringInfo(&query);
            appendStringInfo(&query,
                "SELECT DISTINCT hkm.kmer2_as_uint FROM kmersearch_highfreq_kmer hkm "
                "WHERE hkm.table_oid = %u "
                "AND hkm.column_name = '%s' "
                "AND EXISTS ("
                "    SELECT 1 FROM kmersearch_highfreq_kmer_meta hkm_meta "
                "    WHERE hkm_meta.table_oid = %u "
                "    AND hkm_meta.column_name = '%s' "
                "    AND hkm_meta.kmer_size = %d"
                ") "
                "ORDER BY hkm.kmer2_as_uint "
                "LIMIT %d OFFSET %d",
                table_oid, escaped_column_name, table_oid, escaped_column_name, k_value, current_batch_limit, offset);
            
            /* Execute batch query */
            ret = SPI_execute(query.data, true, 0);
            
            if (ret == SPI_OK_SELECT && SPI_processed > 0) {
                batch_count = SPI_processed;
                batch_kmers = (uint64 *) palloc(batch_count * sizeof(uint64));
                
                for (i = 0; i < batch_count; i++) {
                    bool isnull;
                    Datum kmer_datum;
                    
                    kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
                    if (!isnull) {
                        /* Convert based on k_value size */
                        if (k_value <= 8) {
                            batch_kmers[i] = (uint64)DatumGetInt16(kmer_datum);
                        } else if (k_value <= 16) {
                            batch_kmers[i] = (uint64)DatumGetInt32(kmer_datum);
                        } else {
                            batch_kmers[i] = DatumGetInt64(kmer_datum);
                        }
                    } else {
                        /* This should never happen - kmersearch_highfreq_kmer should not contain NULL values */
                        ereport(ERROR, (errmsg("Unexpected NULL k-mer value in kmersearch_highfreq_kmer table")));
                    }
                }
            } else {
                batch_count = 0;
                batch_kmers = NULL;
            }
            
            /* Cleanup */
            pfree(query.data);
            pfree(escaped_column_name);
            SPI_finish();
        }
        
        if (!batch_kmers || batch_count <= 0) {
            ereport(LOG, (errmsg("dshash_cache_load: No more k-mers in batch %d", batch_num)));
            break;
        }
        
        ereport(LOG, (errmsg("dshash_cache_load: Inserting %d k-mers from batch %d into dshash", 
                             batch_count, batch_num)));
        
        /* Insert batch k-mers into dshash */
        for (i = 0; i < batch_count; i++) {
            uint64 kmer_value;
            void *key_ptr;
            uint16 key16;
            uint32 key32;
            uint64 key64;
            
            /* All k-mer values from the database are valid */
            
            /* All k-mer values are valid, including 0 (which represents "AAAA") */
            
            /* Use kmer2_as_uint value directly */
            kmer_value = batch_kmers[i];
            
            /* Prepare key based on k-mer size */
            if (k_value <= 8) {
                key16 = (uint16)kmer_value;
                key_ptr = &key16;
            } else if (k_value <= 16) {
                key32 = (uint32)kmer_value;
                key_ptr = &key32;
            } else {
                key64 = kmer_value;
                key_ptr = &key64;
            }
            
            ereport(DEBUG1, (errmsg("dshash_cache_load: Inserting k-mer %d/%d from batch %d with value %lu", 
                                   i + 1, batch_count, batch_num, kmer_value)));
            
            /* Insert into dshash table with error handling */
            PG_TRY();
            {
                if (k_value <= 8) {
                    ParallelHighfreqKmerCacheEntry16 *entry16;
                    entry16 = (ParallelHighfreqKmerCacheEntry16 *) dshash_find_or_insert(parallel_cache_hash, 
                                                                                        key_ptr, &found);
                    if (entry16) {
                        entry16->kmer2_as_uint = key16;
                        entry16->frequency_count = 1; /* Mark as high-frequency */
                        dshash_release_lock(parallel_cache_hash, entry16);
                        total_inserted++;
                    }
                } else if (k_value <= 16) {
                    ParallelHighfreqKmerCacheEntry32 *entry32;
                    entry32 = (ParallelHighfreqKmerCacheEntry32 *) dshash_find_or_insert(parallel_cache_hash, 
                                                                                        key_ptr, &found);
                    if (entry32) {
                        entry32->kmer2_as_uint = key32;
                        entry32->frequency_count = 1; /* Mark as high-frequency */
                        dshash_release_lock(parallel_cache_hash, entry32);
                        total_inserted++;
                    }
                } else {
                    ParallelHighfreqKmerCacheEntry64 *entry64;
                    entry64 = (ParallelHighfreqKmerCacheEntry64 *) dshash_find_or_insert(parallel_cache_hash, 
                                                                                        key_ptr, &found);
                    if (entry64) {
                        entry64->kmer2_as_uint = key64;
                        entry64->frequency_count = 1; /* Mark as high-frequency */
                        dshash_release_lock(parallel_cache_hash, entry64);
                        total_inserted++;
                    }
                }
                ereport(DEBUG1, (errmsg("dshash_cache_load: Successfully inserted k-mer %d (total: %d)", 
                                       i + 1, total_inserted)));
            }
            PG_CATCH();
            {
                /* Error handling is done by PG_RE_THROW */
                ereport(ERROR,
                        (errcode(ERRCODE_INTERNAL_ERROR),
                         errmsg("Failed to insert k-mer into dshash table at batch %d index %d", batch_num, i)));
            }
            PG_END_TRY();
        }
        
        /* Clean up batch memory */
        if (batch_kmers) {
            /* uint64 arrays don't need individual pfree calls */
            pfree(batch_kmers);
        }
        
        ereport(LOG, (errmsg("dshash_cache_load: Completed batch %d, total inserted so far: %d", 
                             batch_num, total_inserted)));
        
        /* Move to next batch */
        offset += batch_count;
        batch_num++;
        
        /* Break if we got fewer results than requested (end of data) */
        if (batch_count < current_batch_limit) {
            ereport(LOG, (errmsg("dshash_cache_load: Reached end of data (batch_count=%d < batch_size=%d)", 
                                 batch_count, current_batch_limit)));
            break;
        }
    }
    
    ereport(LOG, (errmsg("dshash_cache_load: Completed all batches, total inserted: %d/%d", 
                         total_inserted, total_kmer_count)));
    
    /* Set the actual number of entries that were successfully inserted */
    parallel_highfreq_cache->num_entries = total_inserted;
    parallel_highfreq_cache->is_initialized = true;
    
    ereport(LOG, (errmsg("dshash_cache_load: Batch processing completed, no additional cleanup needed")));
    
    /* Switch back to original context */
    MemoryContextSwitchTo(oldcontext);
    
    /* Register exit callback for proper DSM cleanup on process exit */
    if (!parallel_cache_exit_callback_registered) {
        ereport(LOG, (errmsg("dshash_cache_load: Registering exit callback for DSM cleanup")));
        on_shmem_exit(dshash_cache_cleanup_callback, 0);
        parallel_cache_exit_callback_registered = true;
    } else {
        ereport(LOG, (errmsg("dshash_cache_load: Exit callback already registered")));
    }
    
    return true;
}

/*
 * Free parallel high-frequency k-mer cache
 */
void
kmersearch_parallel_highfreq_kmer_cache_free_internal(void)
{
    ereport(LOG, (errmsg("kmersearch_parallel_highfreq_kmer_cache_free_internal: Starting parallel cache free")));
    
    /* Use the unified cleanup function for proper resource destruction */
    kmersearch_parallel_cache_cleanup_internal();
    
    ereport(LOG, (errmsg("kmersearch_parallel_highfreq_kmer_cache_free_internal: Parallel cache free completed")));
}

/*
 * Check if parallel high-frequency k-mer cache is valid
 */
bool
kmersearch_parallel_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value)
{
    /* Check if parallel cache exists and is valid */
    if (parallel_highfreq_cache == NULL || !parallel_highfreq_cache->is_initialized)
        return false;
    
    /* Check if cache matches the requested parameters */
    if (parallel_highfreq_cache->cache_key.table_oid != table_oid || 
        parallel_highfreq_cache->cache_key.kmer_size != k_value)
        return false;
    
    return true;
}

/*
 * Lookup entry in parallel high-frequency k-mer cache
 */
bool
kmersearch_parallel_cache_lookup(uint64 kmer_hash)
{
    ParallelHighfreqKmerCacheEntry *entry;
    bool found = false;
    MemoryContext oldcontext;
    
    /* Return false if parallel cache is not initialized */
    if (parallel_cache_hash == NULL)
        return false;
    
    /* Switch to TopMemoryContext for dshash operations */
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    /* Lookup in dshash table */
    entry = (ParallelHighfreqKmerCacheEntry *) dshash_find(parallel_cache_hash, &kmer_hash, false);
    
    if (entry != NULL) {
        found = true;
        /* Must release lock after dshash_find() */
        dshash_release_lock(parallel_cache_hash, entry);
    }
    
    MemoryContextSwitchTo(oldcontext);
    
    return found;
}

/*
 * Attach to existing parallel cache from worker process
 */
static bool
kmersearch_parallel_cache_attach(dsm_handle handle)
{
    dshash_parameters params;
    MemoryContext oldcontext;
    Size cache_struct_size;
    char *dsa_start;
    bool success;
    
    /* Use TopMemoryContext for persistent dshash objects */
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    /* Attach to DSM segment */
    parallel_cache_segment = dsm_attach(handle);
    if (!parallel_cache_segment) {
        MemoryContextSwitchTo(oldcontext);
        return false;
    }
    
    /* Get parallel cache structure from DSM */
    parallel_highfreq_cache = (ParallelHighfreqKmerCache *) dsm_segment_address(parallel_cache_segment);
    
    if (!parallel_highfreq_cache->is_initialized) {
        MemoryContextSwitchTo(oldcontext);
        return false;
    }
    
    /* Set up dshash parameters based on k-mer size from cache */
    memset(&params, 0, sizeof(params));
    params.compare_function = dshash_memcmp;
    params.tranche_id = LWTRANCHE_KMERSEARCH_CACHE;
    
    if (parallel_highfreq_cache->cache_key.kmer_size <= 8) {
        params.key_size = sizeof(uint16);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry16);
        params.hash_function = kmersearch_uint16_identity_hash;
    } else if (parallel_highfreq_cache->cache_key.kmer_size <= 16) {
        params.key_size = sizeof(uint32);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry32);
        params.hash_function = kmersearch_uint32_identity_hash;
    } else {
        params.key_size = sizeof(uint64);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry64);
        params.hash_function = dshash_memhash;
    }
    
    /* Attach to DSA area */
    cache_struct_size = MAXALIGN(sizeof(ParallelHighfreqKmerCache));
    dsa_start = (char *) dsm_segment_address(parallel_cache_segment) + cache_struct_size;
    parallel_cache_dsa = dsa_attach_in_place(dsa_start,
                                             parallel_cache_segment);
    if (!parallel_cache_dsa) {
        MemoryContextSwitchTo(oldcontext);
        return false;
    }
    
    /* Pin the DSA mapping to prevent unexpected cleanup */
    dsa_pin_mapping(parallel_cache_dsa);
    
    /* Attach to dshash table using stored handle */
    parallel_cache_hash = dshash_attach(parallel_cache_dsa, &params, 
                                       parallel_highfreq_cache->hash_handle, NULL);
    
    success = (parallel_cache_hash != NULL);
    
    MemoryContextSwitchTo(oldcontext);
    
    return success;
}


/*
 * Cleanup function for parallel cache on process exit
 */
void
kmersearch_parallel_cache_cleanup_on_exit(int code, Datum arg)
{
    ereport(LOG, (errmsg("parallel_cache_cleanup_on_exit: Starting cleanup, code=%d", code)));
    
    /* Clean up parallel cache resources */
    if (parallel_cache_hash != NULL || parallel_cache_dsa != NULL || parallel_cache_segment != NULL) {
        ereport(LOG, (errmsg("parallel_cache_cleanup_on_exit: Cleaning up resources")));
        kmersearch_parallel_highfreq_kmer_cache_free_internal();
    } else {
        ereport(LOG, (errmsg("parallel_cache_cleanup_on_exit: No resources to clean up")));
    }
}

/*
 * Functions moved from kmersearch_freq.c for better modular organization
 */


/*
 * Check if parallel_highfreq_cache is loaded
 */
bool
kmersearch_is_parallel_highfreq_cache_loaded(void)
{
    return (parallel_highfreq_cache != NULL && 
            parallel_highfreq_cache->is_initialized &&
            parallel_highfreq_cache->num_entries > 0);
}

/*
 * Lookup k-mer in parallel_highfreq_cache  
 */
bool
kmersearch_lookup_in_parallel_cache(VarBit *kmer_key)
{
    MemoryContext oldcontext;
    uint64 kmer_hash;
    void *entry = NULL;
    bool found = false;
    void *key_ptr;
    uint16 key16;
    uint32 key32;
    uint64 key64;
    int k_value;
    
    /* Basic validation checks */
    if (!parallel_highfreq_cache || !parallel_highfreq_cache->is_initialized || 
        parallel_highfreq_cache->num_entries == 0)
        return false;
    
    if (!parallel_cache_hash)
        return false;
    
    /* Get k-mer size from cache */
    k_value = parallel_highfreq_cache->cache_key.kmer_size;
    
    /* Switch to TopMemoryContext for dshash operations */
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    PG_TRY();
    {
        /* Calculate hash using same logic as global cache */
        kmer_hash = kmersearch_ngram_key_to_hash(kmer_key);
        
        /* Prepare key based on k-mer size */
        if (k_value <= 8) {
            key16 = (uint16)kmer_hash;
            key_ptr = &key16;
        } else if (k_value <= 16) {
            key32 = (uint32)kmer_hash;
            key_ptr = &key32;
        } else {
            key64 = kmer_hash;
            key_ptr = &key64;
        }
        
        /* Lookup in dshash table */
        entry = dshash_find(parallel_cache_hash, key_ptr, false);
        
        if (entry != NULL) {
            found = true;
            /* Must release lock after dshash_find() */
            dshash_release_lock(parallel_cache_hash, entry);
        }
    }
    PG_CATCH();
    {
        /* Ensure lock is released even in error cases */
        if (entry)
            dshash_release_lock(parallel_cache_hash, entry);
        MemoryContextSwitchTo(oldcontext);
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    MemoryContextSwitchTo(oldcontext);
    return found;
}

/*
 * Check if uint k-mer exists in global high-frequency cache
 */
bool
kmersearch_lookup_kmer2_as_uint_in_global_cache(uint64 kmer2_as_uint, const char *table_name, const char *column_name)
{
    bool found;
    
    if (!global_highfreq_cache.is_valid || global_highfreq_cache.highfreq_count == 0)
        return false;
    
    hash_search(global_highfreq_cache.highfreq_hash, &kmer2_as_uint, HASH_FIND, &found);
    
    return found;
}

/*
 * Check if uint k-mer exists in parallel high-frequency cache
 */
bool
kmersearch_lookup_kmer2_as_uint_in_parallel_cache(uint64 kmer2_as_uint, const char *table_name, const char *column_name)
{
    MemoryContext oldcontext;
    void *entry = NULL;
    bool found = false;
    void *key_ptr;
    uint16 key16;
    uint32 key32;
    uint64 key64;
    int k_value;
    
    if (!parallel_highfreq_cache || !parallel_highfreq_cache->is_initialized || 
        parallel_highfreq_cache->num_entries == 0)
        return false;
    
    if (!parallel_cache_hash)
        return false;
    
    /* Get k-mer size from cache */
    k_value = parallel_highfreq_cache->cache_key.kmer_size;
    
    /* Prepare key based on k-mer size */
    if (k_value <= 8) {
        key16 = (uint16)kmer2_as_uint;
        key_ptr = &key16;
    } else if (k_value <= 16) {
        key32 = (uint32)kmer2_as_uint;
        key_ptr = &key32;
    } else {
        key64 = kmer2_as_uint;
        key_ptr = &key64;
    }
    
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    PG_TRY();
    {
        entry = dshash_find(parallel_cache_hash, key_ptr, false);
        
        if (entry != NULL) {
            found = true;
            dshash_release_lock(parallel_cache_hash, entry);
        }
    }
    PG_CATCH();
    {
        if (entry)
            dshash_release_lock(parallel_cache_hash, entry);
        MemoryContextSwitchTo(oldcontext);
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    MemoryContextSwitchTo(oldcontext);
    return found;
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
