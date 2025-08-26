/*
 * kmersearch_cache.c - Cache management functions for pg_kmersearch
 *
 * This module contains all cache-related functionality including:
 * - Query-kmer cache for storing parsed k-mer patterns
 * - Actual min score cache for threshold calculations
 * - Cache manager creation, lookup, eviction, and cleanup functions
 */

#include "kmersearch.h"

PG_FUNCTION_INFO_V1(kmersearch_query_kmer_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_query_kmer_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_load);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_free_all);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_load);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_free_all);

HighfreqKmerCache global_highfreq_cache = {0};

bool kmersearch_force_use_parallel_highfreq_kmer_cache = false;

ParallelHighfreqKmerCache *parallel_highfreq_cache = NULL;
dsm_segment *parallel_cache_segment = NULL;
dsa_area *parallel_cache_dsa = NULL;
dshash_table *parallel_cache_hash = NULL;

static uint32 kmersearch_uint16_identity_hash(const void *key, size_t keysize, void *arg);
static uint32 kmersearch_uint32_identity_hash(const void *key, size_t keysize, void *arg);

static void init_query_kmer_cache_manager(QueryKmerCacheManager **manager);
static uint64 generate_query_kmer_cache_key(const char *query_string, int k_size);
static QueryKmerCacheEntry *lookup_query_kmer_cache_entry(QueryKmerCacheManager *manager, const char *query_string, int k_size);
static void store_query_kmer_cache_entry(QueryKmerCacheManager *manager, uint64 hash_key, const char *query_string, int k_size, void *uintkeys, int kmer_count);
static void lru_touch_query_kmer_cache(QueryKmerCacheManager *manager, QueryKmerCacheEntry *entry);
static void lru_evict_oldest_query_kmer_cache(QueryKmerCacheManager *manager);
void kmersearch_free_query_kmer_cache_manager(QueryKmerCacheManager **manager);

static void create_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);
void kmersearch_free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);
static int calculate_actual_min_score_from_uintkey(void *uintkey, int nkeys, int k_size);

void kmersearch_highfreq_kmer_cache_init(void);
bool kmersearch_highfreq_kmer_cache_load_internal(Oid table_oid, const char *column_name, int k_value);
void kmersearch_highfreq_kmer_cache_free_internal(void);
bool kmersearch_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value);

/*
 * Initialize query-kmer cache manager
 */
static void
init_query_kmer_cache_manager(QueryKmerCacheManager **manager)
{
    if (*manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        HASHCTL hash_ctl;
        
        /* Allocate manager in TopMemoryContext */
        *manager = (QueryKmerCacheManager *) palloc0(sizeof(QueryKmerCacheManager));
        
        /* Create query-kmer cache context under TopMemoryContext */
        (*manager)->query_kmer_cache_context = AllocSetContextCreate(TopMemoryContext,
                                                                        "QueryKmerCache",
                                                                        ALLOCSET_DEFAULT_SIZES);
        
        /* Initialize query-kmer cache parameters */
        (*manager)->max_entries = kmersearch_query_kmer_cache_max_entries;
        (*manager)->current_entries = 0;
        (*manager)->hits = 0;
        (*manager)->misses = 0;
        (*manager)->lru_head = NULL;
        (*manager)->lru_tail = NULL;
        
        /* Create hash table */
        MemSet(&hash_ctl, 0, sizeof(hash_ctl));
        hash_ctl.keysize = sizeof(uint64);
        hash_ctl.entrysize = sizeof(QueryKmerCacheEntry);
        hash_ctl.hcxt = (*manager)->query_kmer_cache_context;
        
        (*manager)->hash_table = hash_create("QueryKmerCache", 256, &hash_ctl,
                                           HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
        
        MemoryContextSwitchTo(old_context);
    }
}

/*
 * Generate cache key for query-kmer
 */
static uint64
generate_query_kmer_cache_key(const char *query_string, int k_size)
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
lru_touch_query_kmer_cache(QueryKmerCacheManager *manager, QueryKmerCacheEntry *entry)
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
 * Evict oldest entry from query-kmer cache
 */
static void
lru_evict_oldest_query_kmer_cache(QueryKmerCacheManager *manager)
{
    QueryKmerCacheEntry *tail = manager->lru_tail;
    bool found;
    
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
    if (tail->extracted_uintkey)
        pfree(tail->extracted_uintkey);
    
    manager->current_entries--;
}

/*
 * Look up query-kmer cache entry
 */
static QueryKmerCacheEntry *
lookup_query_kmer_cache_entry(QueryKmerCacheManager *manager, const char *query_string, int k_size)
{
    uint64 hash_key = generate_query_kmer_cache_key(query_string, k_size);
    QueryKmerCacheEntry *entry;
    bool found;
    
    entry = (QueryKmerCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_FIND, &found);
    
    if (found && entry && strcmp(entry->query_string_copy, query_string) == 0 && entry->kmer_size == k_size)
    {
        /* Cache hit - move to head of LRU */
        lru_touch_query_kmer_cache(manager, entry);
        manager->hits++;
        return entry;
    }
    
    return NULL;  /* Cache miss */
}

/*
 * Store entry in query-kmer cache
 */
static void
store_query_kmer_cache_entry(QueryKmerCacheManager *manager, uint64 hash_key, 
                               const char *query_string, int k_size, void *uintkeys, int kmer_count)
{
    QueryKmerCacheEntry *entry;
    MemoryContext old_context;
    bool found;
    size_t uintkey_size;
    int total_bits;
    
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    
    /* Evict oldest entries if cache is full */
    while (manager->current_entries >= manager->max_entries)
        lru_evict_oldest_query_kmer_cache(manager);
    
    old_context = MemoryContextSwitchTo(manager->query_kmer_cache_context);
    
    /* Create new entry */
    entry = (QueryKmerCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_ENTER, &found);
    if (!found)
    {
        entry->hash_key = hash_key;
        entry->query_string_copy = pstrdup(query_string);
        entry->kmer_size = k_size;
        entry->kmer_count = kmer_count;
        
        /* Determine size of uintkey array based on total_bits */
        if (total_bits <= 16)
            uintkey_size = kmer_count * sizeof(uint16);
        else if (total_bits <= 32)
            uintkey_size = kmer_count * sizeof(uint32);
        else
            uintkey_size = kmer_count * sizeof(uint64);
        
        /* Copy uintkey array */
        entry->extracted_uintkey = palloc(uintkey_size);
        memcpy(entry->extracted_uintkey, uintkeys, uintkey_size);
        
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
 * Get cached query uintkeys or extract and cache them
 */
void *
kmersearch_get_cached_query_uintkey(const char *query_string, int k_size, int *nkeys)
{
    QueryKmerCacheEntry *cache_entry;
    void *extracted_uintkeys = NULL;
    uint64 hash_key;
    MemoryContext old_context;
    
    *nkeys = 0;
    
    /* Initialize query-kmer cache manager if not already done */
    if (query_kmer_cache_manager == NULL)
    {
        old_context = MemoryContextSwitchTo(TopMemoryContext);
        init_query_kmer_cache_manager(&query_kmer_cache_manager);
        MemoryContextSwitchTo(old_context);
    }
    
    /* Try to find in cache first */
    cache_entry = lookup_query_kmer_cache_entry(query_kmer_cache_manager, query_string, k_size);
    if (cache_entry != NULL)
    {
        /* Cache hit - return pointer to cached uintkeys directly */
        *nkeys = cache_entry->kmer_count;
        return cache_entry->extracted_uintkey;
    }
    
    /* Cache miss - extract uintkeys and store in cache */
    query_kmer_cache_manager->misses++;
    kmersearch_extract_uintkey_from_text(query_string, &extracted_uintkeys, nkeys);
    
    if (extracted_uintkeys != NULL && *nkeys > 0)
    {
        /* Store in cache */
        hash_key = generate_query_kmer_cache_key(query_string, k_size);
        store_query_kmer_cache_entry(query_kmer_cache_manager, hash_key, 
                                       query_string, k_size, extracted_uintkeys, *nkeys);
        
        /* Find the cache entry we just stored and return its uintkeys */
        cache_entry = lookup_query_kmer_cache_entry(query_kmer_cache_manager, query_string, k_size);
        if (cache_entry != NULL) {
            /* Free the original extracted uintkeys (cache has its own copy) */
            pfree(extracted_uintkeys);
            
            return cache_entry->extracted_uintkey;
        }
    }
    
    return extracted_uintkeys;
}

/*
 * Free query-kmer cache manager
 */
void
kmersearch_free_query_kmer_cache_manager(QueryKmerCacheManager **manager)
{
    if (*manager)
    {
        /* Delete the query-kmer cache context, which will free all allocated memory */
        if ((*manager)->query_kmer_cache_context)
            MemoryContextDelete((*manager)->query_kmer_cache_context);
        
        /* Free the manager itself (allocated in TopMemoryContext) */
        pfree(*manager);
        *manager = NULL;
    }
}

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
kmersearch_free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager)
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
 * Calculate actual min score from uintkey array
 * Helper function for cache miss case
 */
static int
calculate_actual_min_score_from_uintkey(void *uintkey, int nkeys, int k_size)
{
    int base_min_score;
    int highfreq_count = 0;
    int actual_min_score;
    int relative_min = 0;
    int query_total_kmers = nkeys;
    int total_bits;
    
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    
    /* Calculate base minimum score (maximum of absolute and relative) */
    if (query_total_kmers > 0)
    {
        relative_min = (int)ceil(kmersearch_min_shared_kmer_rate * query_total_kmers);
    }
    
    base_min_score = (kmersearch_min_score > relative_min) ? kmersearch_min_score : relative_min;
    
    /* If high-frequency k-mer filtering is enabled, subtract high-frequency k-mer count */
    if (kmersearch_is_highfreq_filtering_enabled())
    {
        /* Count high-frequency k-mers */
        if (total_bits <= 16)
        {
            uint16 *keys = (uint16 *)uintkey;
            for (int i = 0; i < nkeys; i++)
            {
                if (kmersearch_is_uintkey_highfreq((uint64)keys[i], k_size))
                    highfreq_count++;
            }
        }
        else if (total_bits <= 32)
        {
            uint32 *keys = (uint32 *)uintkey;
            for (int i = 0; i < nkeys; i++)
            {
                if (kmersearch_is_uintkey_highfreq((uint64)keys[i], k_size))
                    highfreq_count++;
            }
        }
        else
        {
            uint64 *keys = (uint64 *)uintkey;
            for (int i = 0; i < nkeys; i++)
            {
                if (kmersearch_is_uintkey_highfreq(keys[i], k_size))
                    highfreq_count++;
            }
        }
        
        actual_min_score = base_min_score - highfreq_count;
        
        /* Ensure minimum value of 1 */
        if (actual_min_score < 1)
        {
            actual_min_score = 1;
        }
    }
    else
    {
        actual_min_score = base_min_score;
    }
    
    return actual_min_score;
}

/*
 * Get cached actual_min_score from uintkey array
 * For use with new uintkey-based extraction
 */
int
kmersearch_get_cached_actual_min_score_uintkey(void *uintkey, int nkeys, int k_size)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash = 0;
    MemoryContext old_context;
    int total_bits;
    
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    
    /* Initialize cache manager if needed */
    if (actual_min_score_cache_manager == NULL)
    {
        old_context = MemoryContextSwitchTo(TopMemoryContext);
        create_actual_min_score_cache_manager(&actual_min_score_cache_manager);
        MemoryContextSwitchTo(old_context);
    }
    
    /* Calculate hash based on total_bits to determine uintkey type */
    if (total_bits <= 16)
    {
        uint16 *keys = (uint16 *)uintkey;
        for (int i = 0; i < nkeys; i++)
            query_hash = query_hash * 31 + keys[i];
    }
    else if (total_bits <= 32)
    {
        uint32 *keys = (uint32 *)uintkey;
        for (int i = 0; i < nkeys; i++)
            query_hash = query_hash * 31 + keys[i];
    }
    else
    {
        uint64 *keys = (uint64 *)uintkey;
        for (int i = 0; i < nkeys; i++)
            query_hash = query_hash * 31 + keys[i];
    }
    
    /* Look up in cache */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, NULL);
    
    if (cache_entry == NULL)
    {
        /* Cache miss - calculate and store */
        int actual_min_score = calculate_actual_min_score_from_uintkey(uintkey, nkeys, k_size);
        bool found;
        
        old_context = MemoryContextSwitchTo(actual_min_score_cache_manager->cache_context);
        cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                              &query_hash, HASH_ENTER, &found);
        if (cache_entry != NULL && !found) {
            cache_entry->query_hash = query_hash;
            cache_entry->actual_min_score = actual_min_score;
            actual_min_score_cache_manager->current_entries++;
        }
        MemoryContextSwitchTo(old_context);
        
        actual_min_score_cache_manager->misses++;
    }
    else
    {
        actual_min_score_cache_manager->hits++;
    }
    
    return cache_entry->actual_min_score;
}

/*
 * Get cached actual_min_score from int2 Datum array
 */
int
kmersearch_get_cached_actual_min_score_datum_int2(Datum *queryKeys, int nkeys)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash = 0;
    
    /* Initialize cache manager if needed */
    if (actual_min_score_cache_manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        create_actual_min_score_cache_manager(&actual_min_score_cache_manager);
        MemoryContextSwitchTo(old_context);
    }
    
    /* Calculate hash from uint16 values - must match extract_query */
    for (int i = 0; i < nkeys; i++)
    {
        uint16 key = (uint16)DatumGetInt16(queryKeys[i]);
        query_hash = query_hash * 31 + key;
    }
    
    /* Look up in cache */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, NULL);
    
    if (cache_entry == NULL)
    {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("actual_min_score not found in cache for int2"),
                 errdetail("Query hash: %lu, nkeys: %d", query_hash, nkeys)));
    }
    
    actual_min_score_cache_manager->hits++;
    return cache_entry->actual_min_score;
}

/*
 * Get cached actual_min_score from int4 Datum array
 */
int
kmersearch_get_cached_actual_min_score_datum_int4(Datum *queryKeys, int nkeys)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash = 0;
    
    /* Initialize cache manager if needed */
    if (actual_min_score_cache_manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        create_actual_min_score_cache_manager(&actual_min_score_cache_manager);
        MemoryContextSwitchTo(old_context);
    }
    
    /* Calculate hash from uint32 values - must match extract_query */
    for (int i = 0; i < nkeys; i++)
    {
        uint32 key = (uint32)DatumGetInt32(queryKeys[i]);
        query_hash = query_hash * 31 + key;
    }
    
    /* Look up in cache */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, NULL);
    
    if (cache_entry == NULL)
    {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("actual_min_score not found in cache for int4"),
                 errdetail("Query hash: %lu, nkeys: %d", query_hash, nkeys)));
    }
    
    actual_min_score_cache_manager->hits++;
    return cache_entry->actual_min_score;
}

/*
 * Get cached actual_min_score from int8 Datum array
 */
int
kmersearch_get_cached_actual_min_score_datum_int8(Datum *queryKeys, int nkeys)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash = 0;
    
    /* Initialize cache manager if needed */
    if (actual_min_score_cache_manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        create_actual_min_score_cache_manager(&actual_min_score_cache_manager);
        MemoryContextSwitchTo(old_context);
    }
    
    /* Calculate hash from uint64 values - must match extract_query */
    for (int i = 0; i < nkeys; i++)
    {
        uint64 key = (uint64)DatumGetInt64(queryKeys[i]);
        query_hash = query_hash * 31 + key;
    }
    
    /* Look up in cache */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, NULL);
    
    if (cache_entry == NULL)
    {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("actual_min_score not found in cache for int8"),
                 errdetail("Query hash: %lu, nkeys: %d", query_hash, nkeys)));
    }
    
    actual_min_score_cache_manager->hits++;
    return cache_entry->actual_min_score;
}

/*
 * Query-kmer cache statistics function
 */
Datum
kmersearch_query_kmer_cache_stats(PG_FUNCTION_ARGS)
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
    
    /* Query-kmer cache statistics */
    if (query_kmer_cache_manager)
    {
        values[0] = Int64GetDatum(query_kmer_cache_manager->hits);
        values[1] = Int64GetDatum(query_kmer_cache_manager->misses);
        values[2] = Int32GetDatum(query_kmer_cache_manager->current_entries);
        values[3] = Int32GetDatum(query_kmer_cache_manager->max_entries);
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
 * Query-kmer cache free function
 */
Datum
kmersearch_query_kmer_cache_free(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    /* Count entries before freeing */
    if (query_kmer_cache_manager)
        freed_entries += query_kmer_cache_manager->current_entries;
    
    /* Free query-kmer cache manager (uses TopMemoryContext - needs manual cleanup) */
    /* DNA2/DNA4 cache managers are now local and automatically freed with QueryContext */
    kmersearch_free_query_kmer_cache_manager(&query_kmer_cache_manager);
    
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
    kmersearch_free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
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
    int total_bits;
    
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
        kmersearch_highfreq_kmer_cache_free_internal();
    } else {
    }
    
    /* Check cache context state and investigate why it might be NULL */
    if (global_highfreq_cache.cache_context == NULL) {
        /* Log detailed information about why the context might be freed */
        ereport(WARNING, 
                (errmsg("Cache context unexpectedly freed"),
                 errdetail("Cache was %svalid, had %d entries",
                          global_highfreq_cache.is_valid ? "" : "not ",
                          global_highfreq_cache.highfreq_count),
                 errhint("This may indicate a memory management issue or improper cache cleanup sequence")));
        
        /* Reinitialize the cache with additional logging */
        kmersearch_highfreq_kmer_cache_init();
        
        /* Verify reinitialization was successful */
        if (global_highfreq_cache.cache_context == NULL) {
            ereport(ERROR,
                    (errcode(ERRCODE_INTERNAL_ERROR),
                     errmsg("Failed to reinitialize cache context"),
                     errhint("Check PostgreSQL memory settings and available memory")));
        }
    } else {
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
            "SELECT COUNT(DISTINCT hkm.uintkey) FROM kmersearch_highfreq_kmer hkm "
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
        return false;
    }
    
    /* Switch to cache context for cache storage */
    old_context = MemoryContextSwitchTo(global_highfreq_cache.cache_context);
    
    /* Store in cache - build cache key */
    global_highfreq_cache.current_cache_key.table_oid = table_oid;
    global_highfreq_cache.current_cache_key.column_name_hash = hash_any((unsigned char*)column_name, strlen(column_name));
    global_highfreq_cache.current_cache_key.kmer_size = k_value;
    global_highfreq_cache.current_cache_key.occur_bitlen = kmersearch_occur_bitlen;
    global_highfreq_cache.current_cache_key.max_appearance_rate = kmersearch_max_appearance_rate;
    global_highfreq_cache.current_cache_key.max_appearance_nrow = kmersearch_max_appearance_nrow;
    
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
        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(global_highfreq_cache.cache_context);
        global_highfreq_cache.cache_context = NULL;
        global_highfreq_cache.is_valid = false;
        return false;
    }
    
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
                "SELECT DISTINCT hkm.uintkey FROM kmersearch_highfreq_kmer hkm "
                "WHERE hkm.table_oid = %u "
                "AND hkm.column_name = '%s' "
                "AND EXISTS ("
                "    SELECT 1 FROM kmersearch_highfreq_kmer_meta hkm_meta "
                "    WHERE hkm_meta.table_oid = %u "
                "    AND hkm_meta.column_name = '%s' "
                "    AND hkm_meta.kmer_size = %d"
                ") "
                "ORDER BY hkm.uintkey "
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
                        /* Convert based on total bit length (k-mer2 bits + occurrence bits) */
                        total_bits = k_value * 2 + kmersearch_occur_bitlen;
                        if (total_bits <= 16) {
                            batch_kmers[i] = (uint64)DatumGetInt16(kmer_datum);
                        } else if (total_bits <= 32) {
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
            break;
        }
        
        /* Insert batch k-mers into hash table */
        for (i = 0; i < batch_count; i++) {
            uint64 uintkey;
            HighfreqKmerHashEntry *entry;
            bool found;
            
            /* All k-mer values are valid, including 0 (which represents "AAAA") */
            
            /* Use uintkey value directly */
            uintkey = batch_kmers[i];
            
            /* Insert into hash table using uintkey as both key and value */
            entry = (HighfreqKmerHashEntry *) hash_search(global_highfreq_cache.highfreq_hash,
                                                         (void *) &uintkey,
                                                         HASH_ENTER,
                                                         &found);
            
            if (entry && !found) {
                entry->kmer_key = NULL;  /* No longer storing VarBit */
                entry->hash_value = uintkey;
                total_inserted++;
            } else if (found) {
            } else {
            }
        }
        
        /* Clean up batch memory in parent context */
        MemoryContextSwitchTo(old_context);
        if (batch_kmers) {
            /* uint64 arrays don't need individual pfree calls */
            pfree(batch_kmers);
        }
        MemoryContextSwitchTo(global_highfreq_cache.cache_context);
        
        /* Move to next batch */
        offset += batch_count;
        batch_num++;
        
        /* Break if we got fewer results than requested (end of data) */
        if (batch_count < current_batch_limit) {
            break;
        }
    }
    
    /* Set cache metadata */
    global_highfreq_cache.highfreq_kmers = NULL;  /* We don't store the array anymore */
    global_highfreq_cache.highfreq_count = total_inserted;
    
    
    if (global_highfreq_cache.highfreq_hash) {
        global_highfreq_cache.is_valid = true;
    } else {
        /* Hash table creation failed, clean up by deleting the context */
        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(global_highfreq_cache.cache_context);
        global_highfreq_cache.cache_context = NULL;
        global_highfreq_cache.is_valid = false;
        return false;
    }
    
    MemoryContextSwitchTo(old_context);
    
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
    }
    
    return matches;
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
        
        /* Get and validate kmer_size */
        kmer_size_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull)
        {
            stored_kmer_size = DatumGetInt32(kmer_size_datum);
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
        occur_bitlen_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull);
        if (!isnull)
        {
            stored_occur_bitlen = DatumGetInt32(occur_bitlen_datum);
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
        max_appearance_rate_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull);
        if (!isnull)
        {
            stored_max_appearance_rate = DatumGetFloat4(max_appearance_rate_datum);
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
        max_appearance_nrow_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 4, &isnull);
        if (!isnull)
        {
            stored_max_appearance_nrow = DatumGetInt32(max_appearance_nrow_datum);
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
    }
    else
    {
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
    } else {
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
 * GUC hook function for query-kmer cache max entries changes
 */
void
kmersearch_query_kmer_cache_max_entries_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Clear query-kmer cache to recreate with new size limit */
    if (query_kmer_cache_manager)
        kmersearch_free_query_kmer_cache_manager(&query_kmer_cache_manager);
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
    
    /* Prevent double cleanup - check if already cleaned up */
    if (parallel_cache_hash == NULL && parallel_cache_dsa == NULL && parallel_cache_segment == NULL) {
        return;
    }
    
    /* Use TopMemoryContext for dshash operations */
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    
    /* Step 1: Handle dshash table cleanup based on process type */
    if (parallel_cache_hash != NULL) {
        
        /* Check if DSA and DSM are still valid before operating on dshash */
        if (parallel_cache_dsa != NULL && parallel_cache_segment != NULL) {
            if (!IsParallelWorker()) {
                /* Main process: destroy the hash table */
                dshash_destroy(parallel_cache_hash);
            } else {
                /* Parallel worker: detach from hash table */
                dshash_detach(parallel_cache_hash);
            }
        } else {
            /* DSA/DSM already destroyed, just detach without destroy */
            dshash_detach(parallel_cache_hash);
        }
        parallel_cache_hash = NULL;
    } else {
    }
    
    /* Switch back to original context before DSA/DSM operations */
    MemoryContextSwitchTo(oldcontext);
    
    /* Step 2: Handle DSA area cleanup */
    if (parallel_cache_dsa != NULL) {
        if (!IsParallelWorker()) {
            /* Main process: unpin and detach */
            PG_TRY();
            {
                /* Unpin DSA area */
                dsa_unpin(parallel_cache_dsa);
                /* Detach from DSA area */
                dsa_detach(parallel_cache_dsa);
            }
            PG_CATCH();
            {
                FlushErrorState();
            }
            PG_END_TRY();
        } else {
            /* Parallel worker: detach only */
            dsa_detach(parallel_cache_dsa);
        }
        parallel_cache_dsa = NULL;
    } else {
    }
    
    /* Step 3: Handle DSM segment properly */
    if (parallel_cache_segment != NULL) {
        if (!IsParallelWorker()) {
            /* Main process: unpin and detach */
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
            }
            PG_CATCH();
            {
                FlushErrorState();
            }
            PG_END_TRY();
        } else {
            /* Parallel worker: detach only */
            dsm_detach(parallel_cache_segment);
        }
        parallel_cache_segment = NULL;
    } else {
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
    
    kmersearch_parallel_cache_cleanup_internal();
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
    int total_bits;
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
            parallel_highfreq_cache->cache_key.column_name_hash == hash_any((unsigned char*)column_name, strlen(column_name)) &&
            parallel_highfreq_cache->cache_key.kmer_size == k_value &&
            parallel_highfreq_cache->cache_key.occur_bitlen == kmersearch_occur_bitlen &&
            fabs(parallel_highfreq_cache->cache_key.max_appearance_rate - kmersearch_max_appearance_rate) < 0.0001 &&
            parallel_highfreq_cache->cache_key.max_appearance_nrow == kmersearch_max_appearance_nrow) {
            /* All parameters match (table, column, and GUC variables) - return true */
            return true;
        } else {
            /* Cache exists but GUC variables don't match - keep cache and return false */
            return false;
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
            "SELECT COUNT(DISTINCT hkm.uintkey) FROM kmersearch_highfreq_kmer hkm "
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
        return false;
    }
    
    /* Calculate required segment size using total count */
    
    cache_struct_size = MAXALIGN(sizeof(ParallelHighfreqKmerCache));
    
    /* Calculate entry size based on total_bits (kmer_size * 2 + occur_bitlen) */
    total_bits = k_value * 2 + kmersearch_occur_bitlen;
    if (total_bits <= 16) {
        entries_size = total_kmer_count * sizeof(ParallelHighfreqKmerCacheEntry16);
    } else if (total_bits <= 32) {
        entries_size = total_kmer_count * sizeof(ParallelHighfreqKmerCacheEntry32);
    } else {
        entries_size = total_kmer_count * sizeof(ParallelHighfreqKmerCacheEntry64);
    }
    
    dsa_min_size = 8192; /* Minimum DSA area size */
    dshash_overhead = MAXALIGN(512); /* Extra space for dshash overhead */
    
    /* Total size = cache structure + DSA area (at least 8192) + entries + overhead */
    segment_size = cache_struct_size + dsa_min_size + entries_size + dshash_overhead;
    
    /* Ensure minimum segment size for DSM/DSA requirements */
    if (segment_size < 16384) {
        segment_size = 16384; /* At least 16KB total */
    }
    
    elog(DEBUG1, "DSM segment size calculation: cache_struct=%zu, dsa_min=%zu, entries=%zu, overhead=%zu, total=%zu",
         cache_struct_size, dsa_min_size, entries_size, dshash_overhead, segment_size);
    
    /* Create DSM segment */
    
    parallel_cache_segment = dsm_create(segment_size, 0);
    if (!parallel_cache_segment) {
        ereport(ERROR,
                (errcode(ERRCODE_OUT_OF_MEMORY),
                 errmsg("failed to create DSM segment for parallel cache")));
    }
    
    /* Pin the DSM segment to prevent automatic cleanup when query ends */
    dsm_pin_segment(parallel_cache_segment);
    dsm_pin_mapping(parallel_cache_segment);
    
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
    
    /* Set up dshash parameters based on total bit length (k-mer2 bits + occurrence bits) */
    memset(&params, 0, sizeof(params));
    params.compare_function = dshash_memcmp;
    params.tranche_id = LWTRANCHE_KMERSEARCH_CACHE;
    
    total_bits = k_value * 2 + kmersearch_occur_bitlen;
    if (total_bits <= 16) {
        params.key_size = sizeof(uint16);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry16);
        params.hash_function = kmersearch_uint16_identity_hash;
    } else if (total_bits <= 32) {
        params.key_size = sizeof(uint32);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry32);
        params.hash_function = kmersearch_uint32_identity_hash;
    } else {
        params.key_size = sizeof(uint64);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry64);
        params.hash_function = dshash_memhash;
    }
    
    /* Create dshash table */
    
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
    
    /* Store the dshash table handle */
    parallel_highfreq_cache->hash_handle = dshash_get_hash_table_handle(parallel_cache_hash);
    parallel_highfreq_cache->is_initialized = true;
    
    /* Populate the hash table with high-frequency k-mers using batch processing */
    
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
                "SELECT DISTINCT hkm.uintkey FROM kmersearch_highfreq_kmer hkm "
                "WHERE hkm.table_oid = %u "
                "AND hkm.column_name = '%s' "
                "AND EXISTS ("
                "    SELECT 1 FROM kmersearch_highfreq_kmer_meta hkm_meta "
                "    WHERE hkm_meta.table_oid = %u "
                "    AND hkm_meta.column_name = '%s' "
                "    AND hkm_meta.kmer_size = %d"
                ") "
                "ORDER BY hkm.uintkey "
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
                        /* Convert based on total bit length (k-mer2 bits + occurrence bits) */
                        total_bits = k_value * 2 + kmersearch_occur_bitlen;
                        if (total_bits <= 16) {
                            batch_kmers[i] = (uint64)DatumGetInt16(kmer_datum);
                        } else if (total_bits <= 32) {
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
            break;
        }
        
        /* Insert batch k-mers into dshash */
        for (i = 0; i < batch_count; i++) {
            uint64 kmer_value;
            void *key_ptr;
            uint16 key16;
            uint32 key32;
            uint64 key64;
            
            /* All k-mer values from the database are valid */
            
            /* All k-mer values are valid, including 0 (which represents "AAAA") */
            
            /* Use uintkey value directly */
            kmer_value = batch_kmers[i];
            
            /* Prepare key based on total bit length (k-mer2 bits + occurrence bits) */
            if (total_bits <= 16) {
                key16 = (uint16)kmer_value;
                key_ptr = &key16;
            } else if (total_bits <= 32) {
                key32 = (uint32)kmer_value;
                key_ptr = &key32;
            } else {
                key64 = kmer_value;
                key_ptr = &key64;
            }
            
            /* Insert into dshash table with error handling */
            PG_TRY();
            {
                if (total_bits <= 16) {
                    ParallelHighfreqKmerCacheEntry16 *entry16;
                    entry16 = (ParallelHighfreqKmerCacheEntry16 *) dshash_find_or_insert(parallel_cache_hash, 
                                                                                        key_ptr, &found);
                    if (entry16) {
                        entry16->uintkey = key16;
                        entry16->frequency_count = 1; /* Mark as high-frequency */
                        dshash_release_lock(parallel_cache_hash, entry16);
                        total_inserted++;
                    }
                } else if (total_bits <= 32) {
                    ParallelHighfreqKmerCacheEntry32 *entry32;
                    entry32 = (ParallelHighfreqKmerCacheEntry32 *) dshash_find_or_insert(parallel_cache_hash, 
                                                                                        key_ptr, &found);
                    if (entry32) {
                        entry32->uintkey = key32;
                        entry32->frequency_count = 1; /* Mark as high-frequency */
                        dshash_release_lock(parallel_cache_hash, entry32);
                        total_inserted++;
                    }
                } else {
                    ParallelHighfreqKmerCacheEntry64 *entry64;
                    entry64 = (ParallelHighfreqKmerCacheEntry64 *) dshash_find_or_insert(parallel_cache_hash, 
                                                                                        key_ptr, &found);
                    if (entry64) {
                        entry64->uintkey = key64;
                        entry64->frequency_count = 1; /* Mark as high-frequency */
                        dshash_release_lock(parallel_cache_hash, entry64);
                        total_inserted++;
                    }
                }
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
        
        /* Move to next batch */
        offset += batch_count;
        batch_num++;
        
        /* Break if we got fewer results than requested (end of data) */
        if (batch_count < current_batch_limit) {
            break;
        }
    }
    
    /* Set the actual number of entries that were successfully inserted */
    parallel_highfreq_cache->num_entries = total_inserted;
    parallel_highfreq_cache->is_initialized = true;
    
    /* Switch back to original context */
    MemoryContextSwitchTo(oldcontext);
    
    /* Register exit callback for proper DSM cleanup on process exit */
    if (!parallel_cache_exit_callback_registered) {
        on_shmem_exit(dshash_cache_cleanup_callback, 0);
        parallel_cache_exit_callback_registered = true;
    } else {
    }
    
    return true;
}

/*
 * Free parallel high-frequency k-mer cache
 */
void
kmersearch_parallel_highfreq_kmer_cache_free_internal(void)
{
    
    /* Use the unified cleanup function for proper resource destruction */
    kmersearch_parallel_cache_cleanup_internal();
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
    int total_bits;
    
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
    
    /* Set up dshash parameters based on total_bits from cache */
    /* Calculate total bits to determine key type */
    total_bits = parallel_highfreq_cache->cache_key.kmer_size * 2 + 
                 parallel_highfreq_cache->cache_key.occur_bitlen;
    
    memset(&params, 0, sizeof(params));
    params.compare_function = dshash_memcmp;
    params.tranche_id = LWTRANCHE_KMERSEARCH_CACHE;
    
    if (total_bits <= 16) {
        /* Use uint16 for keys when total_bits <= 16 */
        params.key_size = sizeof(uint16);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry16);
        params.hash_function = kmersearch_uint16_identity_hash;
    } else if (total_bits <= 32) {
        /* Use uint32 for keys when total_bits <= 32 */
        params.key_size = sizeof(uint32);
        params.entry_size = sizeof(ParallelHighfreqKmerCacheEntry32);
        params.hash_function = kmersearch_uint32_identity_hash;
    } else {
        /* Use uint64 for keys when total_bits <= 64 */
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
    
    /* Clean up parallel cache resources */
    if (parallel_cache_hash != NULL || parallel_cache_dsa != NULL || parallel_cache_segment != NULL) {
        kmersearch_parallel_highfreq_kmer_cache_free_internal();
    } else {
    }
}

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
 * Check if uintkey exists in global high-frequency cache
 */
bool
kmersearch_lookup_uintkey_in_global_cache(uint64 uintkey, const char *table_name, const char *column_name)
{
    bool found;
    
    if (!global_highfreq_cache.is_valid || global_highfreq_cache.highfreq_count == 0)
        return false;
    
    hash_search(global_highfreq_cache.highfreq_hash, &uintkey, HASH_FIND, &found);
    
    return found;
}

/*
 * Check if uintkey exists in parallel high-frequency cache
 */
bool
kmersearch_lookup_uintkey_in_parallel_cache(uint64 uintkey, const char *table_name, const char *column_name)
{
    MemoryContext oldcontext;
    void *entry = NULL;
    bool found = false;
    void *key_ptr;
    uint16 key16;
    uint32 key32;
    uint64 key64;
    int k_value;
    int total_bits;
    
    if (!parallel_highfreq_cache || !parallel_highfreq_cache->is_initialized || 
        parallel_highfreq_cache->num_entries == 0)
        return false;
    
    if (!parallel_cache_hash)
        return false;
    
    /* Get k-mer size from cache */
    k_value = parallel_highfreq_cache->cache_key.kmer_size;
    
    /* Calculate total bits needed for type selection */
    total_bits = k_value * 2 + kmersearch_occur_bitlen;
    
    /* Prepare key based on total bits */
    if (total_bits <= 16) {
        key16 = (uint16)uintkey;
        key_ptr = &key16;
    } else if (total_bits <= 32) {
        key32 = (uint32)uintkey;
        key_ptr = &key32;
    } else {
        key64 = uintkey;
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
