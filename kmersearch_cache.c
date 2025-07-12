/*
 * kmersearch_cache.c - Cache management functions for pg_kmersearch
 *
 * This module contains all cache-related functionality including:
 * - Query pattern cache for storing parsed k-mer patterns
 * - Actual min score cache for threshold calculations
 * - Rawscore cache for computed similarity scores
 * - Cache manager creation, lookup, eviction, and cleanup functions
 */

#include "kmersearch.h"

/* PostgreSQL function info declarations for cache functions */
PG_FUNCTION_INFO_V1(kmersearch_rawscore_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_rawscore_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_query_pattern_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_query_pattern_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_load);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmer_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_load);
PG_FUNCTION_INFO_V1(kmersearch_parallel_highfreq_kmer_cache_free);

/* Global high-frequency k-mer cache */
HighfreqKmerCache global_highfreq_cache = {0};

/* Global testing variable for dshash usage (not exposed to users) */
bool kmersearch_force_use_dshash = false;

/* Global parallel cache state */
ParallelHighfreqKmerCache *parallel_highfreq_cache = NULL;
dsm_segment *parallel_cache_segment = NULL;
dsa_area *parallel_cache_dsa = NULL;
dshash_table *parallel_cache_hash = NULL;

/*
 * Forward declarations for internal functions
 */

/* Query pattern cache functions */
static void init_query_pattern_cache_manager(QueryPatternCacheManager **manager);
static uint64 generate_query_pattern_cache_key(const char *query_string, int k_size);
static QueryPatternCacheEntry *lookup_query_pattern_cache_entry(QueryPatternCacheManager *manager, const char *query_string, int k_size);
static void store_query_pattern_cache_entry(QueryPatternCacheManager *manager, uint64 hash_key, const char *query_string, int k_size, VarBit **kmers, int kmer_count);
static void lru_touch_query_pattern_cache(QueryPatternCacheManager *manager, QueryPatternCacheEntry *entry);
static void lru_evict_oldest_query_pattern_cache(QueryPatternCacheManager *manager);
static void free_query_pattern_cache_manager(QueryPatternCacheManager **manager);

/* Actual min score cache functions */
static void create_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);
static void free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);

/* Rawscore cache functions */
static RawscoreCacheManager *create_rawscore_cache_manager(const char *name);
static void free_rawscore_cache_manager(RawscoreCacheManager **manager);
static uint64 generate_cache_key(VarBit *sequence, const char *query_string);
static bool sequences_equal(VarBit *a, VarBit *b);
static RawscoreCacheEntry *lookup_rawscore_cache_entry(RawscoreCacheManager *manager, VarBit *sequence, const char *query_string);
static void store_rawscore_cache_entry(RawscoreCacheManager *manager, uint64 hash_key, VarBit *sequence, VarBit **query_keys, const char *query_string, KmerMatchResult result);

/* High-frequency k-mer cache functions */
void kmersearch_highfreq_kmer_cache_init(void);
bool kmersearch_highfreq_kmer_cache_load_internal(Oid table_oid, const char *column_name, int k_value);
void kmersearch_highfreq_kmer_cache_free_internal(void);
bool kmersearch_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value);

/* Rawscore cache heap management */
static void rawscore_heap_swap(RawscoreCacheManager *manager, int i, int j);
static void rawscore_heap_bubble_up(RawscoreCacheManager *manager, int index);
static void rawscore_heap_bubble_down(RawscoreCacheManager *manager, int index);
static void rawscore_heap_insert(RawscoreCacheManager *manager, RawscoreCacheEntry *entry);
static void rawscore_heap_remove(RawscoreCacheManager *manager, RawscoreCacheEntry *entry);
static void rawscore_heap_evict_lowest_score(RawscoreCacheManager *manager);

/* High-level cache functions */

/* External global variables (defined in kmersearch.c) */
extern ActualMinScoreCacheManager *actual_min_score_cache_manager;
extern QueryPatternCacheManager *query_pattern_cache_manager;
extern RawscoreCacheManager *rawscore_cache_manager;

extern int kmersearch_actual_min_score_cache_max_entries;
extern int kmersearch_query_pattern_cache_max_entries;
extern int kmersearch_rawscore_cache_max_entries;

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
        /* Cache hit - copy cached k-mers to current memory context */
        *nkeys = cache_entry->kmer_count;
        if (*nkeys > 0)
        {
            result_kmers = (VarBit **) palloc(*nkeys * sizeof(VarBit *));
            for (i = 0; i < *nkeys; i++)
            {
                /* Copy each k-mer from cache context to current context */
                result_kmers[i] = (VarBit *) palloc(VARSIZE(cache_entry->extracted_kmers[i]));
                memcpy(result_kmers[i], cache_entry->extracted_kmers[i], VARSIZE(cache_entry->extracted_kmers[i]));
            }
        }
        return result_kmers;
    }
    
    /* Cache miss - extract k-mers and store in cache */
    query_pattern_cache_manager->misses++;
    extracted_kmers = kmersearch_extract_query_kmer_with_degenerate(query_string, k_size, nkeys);
    
    if (extracted_kmers != NULL && *nkeys > 0)
    {
        /* Store in cache */
        hash_key = generate_query_pattern_cache_key(query_string, k_size);
        store_query_pattern_cache_entry(query_pattern_cache_manager, hash_key, 
                                       query_string, k_size, extracted_kmers, *nkeys);
        
        /* Return the extracted k-mers */
        result_kmers = extracted_kmers;
    }
    
    return result_kmers;
}

/*
 * Free query pattern cache manager
 */
static void
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
static void
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
    int absolute_min;
    int relative_min;
    
    elog(LOG, "calculate_actual_min_score: Started with nkeys=%d, query_total_kmers=%d", nkeys, query_total_kmers);
    
    /* Validate input parameters */
    if (query_keys == NULL) {
        elog(LOG, "calculate_actual_min_score: query_keys is NULL, returning default min score");
        return kmersearch_min_score;
    }
    
    elog(LOG, "calculate_actual_min_score: About to call kmersearch_get_adjusted_min_score");
    /* Use adjusted minimum score that considers high-frequency k-mer filtering */
    absolute_min = kmersearch_get_adjusted_min_score(query_keys, nkeys);
    elog(LOG, "calculate_actual_min_score: kmersearch_get_adjusted_min_score returned %d", absolute_min);
    
    if (query_total_kmers > 0) {
        /* Calculate minimum score from relative threshold */
        /* kmersearch_min_shared_ngram_key_rate <= shared_count / query_total_kmers */
        /* shared_count >= kmersearch_min_shared_ngram_key_rate * query_total_kmers */
        double relative_threshold = kmersearch_min_shared_ngram_key_rate * query_total_kmers;
        relative_min = (int)ceil(relative_threshold);
    } else {
        relative_min = 0;
    }
    
    /* Return the maximum of absolute and relative minimums */
    return (absolute_min > relative_min) ? absolute_min : relative_min;
}

/*
 * Get cached actual min score using TopMemoryContext cache (global)
 */
int
get_cached_actual_min_score(VarBit **query_keys, int nkeys)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash;
    bool found;
    MemoryContext old_context;
    int actual_min_score;
    int i;
    
    elog(LOG, "get_cached_actual_min_score: Started with nkeys=%d", nkeys);
    
    /* Validate input parameters */
    if (query_keys == NULL) {
        elog(LOG, "get_cached_actual_min_score: query_keys is NULL, returning fallback score");
        return calculate_actual_min_score(query_keys, nkeys, nkeys);
    }
    
    /* Validate query_keys array elements */
    for (i = 0; i < nkeys; i++) {
        if (query_keys[i] == NULL) {
            elog(LOG, "get_cached_actual_min_score: query_keys[%d] is NULL, returning fallback score", i);
            return calculate_actual_min_score(query_keys, nkeys, nkeys);
        }
        elog(LOG, "get_cached_actual_min_score: query_keys[%d] pointer valid: %p", i, query_keys[i]);
    }
    
    /* Create cache manager in TopMemoryContext if not exists */
    if (actual_min_score_cache_manager == NULL)
    {
        elog(LOG, "get_cached_actual_min_score: Creating cache manager");
        old_context = MemoryContextSwitchTo(TopMemoryContext);
        
        PG_TRY();
        {
            create_actual_min_score_cache_manager(&actual_min_score_cache_manager);
        }
        PG_CATCH();
        {
            MemoryContextSwitchTo(old_context);
            elog(LOG, "get_cached_actual_min_score: Cache creation failed, falling back to direct calculation");
            /* Fallback to direct calculation if cache creation fails */
            return calculate_actual_min_score(query_keys, nkeys, nkeys);
        }
        PG_END_TRY();
        
        MemoryContextSwitchTo(old_context);
        elog(LOG, "get_cached_actual_min_score: Cache manager created successfully");
    }
    
    elog(LOG, "get_cached_actual_min_score: About to calculate hash value");
    
    /* Calculate hash value for query keys content (not pointers) */
    query_hash = 0;
    for (i = 0; i < nkeys; i++) {
        uint64 kmer_hash;
        elog(LOG, "get_cached_actual_min_score: Hashing k-mer %d", i+1);
        kmer_hash = hash_any_extended((unsigned char *)VARBITS(query_keys[i]), 
                                      VARBITBYTES(query_keys[i]), query_hash);
        query_hash = kmer_hash;
        elog(LOG, "get_cached_actual_min_score: K-mer %d hash calculated: %lu", i+1, kmer_hash);
    }
    
    elog(LOG, "get_cached_actual_min_score: Hash calculation completed, query_hash=%lu", query_hash);
    
    /* Look up in hash table */
    elog(LOG, "get_cached_actual_min_score: About to search hash table");
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, &found);
    
    elog(LOG, "get_cached_actual_min_score: Hash search completed, found=%d", found ? 1 : 0);
    
    if (found) {
        elog(LOG, "get_cached_actual_min_score: Cache hit, returning cached score=%d", cache_entry->actual_min_score);
        actual_min_score_cache_manager->hits++;
        return cache_entry->actual_min_score;
    }
    
    /* Not found - calculate and cache */
    elog(LOG, "get_cached_actual_min_score: Cache miss, calculating actual min score");
    actual_min_score_cache_manager->misses++;
    actual_min_score = calculate_actual_min_score(query_keys, nkeys, nkeys);
    elog(LOG, "get_cached_actual_min_score: Calculated actual_min_score=%d", actual_min_score);
    
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
    
    elog(LOG, "get_cached_actual_min_score: Function completed successfully, returning actual_min_score=%d", actual_min_score);
    return actual_min_score;
}

/* ================================
 * RAWSCORE CACHE IMPLEMENTATION
 * ================================ */

/* Global rawscore cache statistics */
static struct {
    uint64 dna2_hits;
    uint64 dna2_misses;
    int dna2_current_entries;
    int dna2_max_entries;
    uint64 dna4_hits;
    uint64 dna4_misses;
    int dna4_current_entries;
    int dna4_max_entries;
} rawscore_cache_stats = {0};

/*
 * GUC hook function for rawscore cache max entries changes
 */
void
kmersearch_rawscore_cache_max_entries_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Clear rawscore cache to recreate with new size limit */
    if (rawscore_cache_manager)
        free_rawscore_cache_manager(&rawscore_cache_manager);
}

/*
 * Generate cache key from sequence and query string
 */
static uint64
generate_cache_key(VarBit *sequence, const char *query_string)
{
    uint64 seq_hash, query_hash, kmer_hash;
    
    /* Safety checks */
    if (sequence == NULL || query_string == NULL) {
        elog(WARNING, "generate_cache_key: NULL parameter detected");
        return 0;
    }
    
    /* Validate VarBit structure */
    if (VARBITLEN(sequence) == 0 || VARBITBYTES(sequence) == 0) {
        elog(WARNING, "generate_cache_key: Invalid VarBit structure");
        return 0;
    }
    
    /* Hash sequence data */
    seq_hash = hash_any_extended(VARBITS(sequence), VARBITBYTES(sequence), 0);
    
    /* Hash query string with different seed */
    query_hash = hash_any_extended((unsigned char*)query_string, strlen(query_string), 1);
    
    /* Hash k-mer size to ensure different cache entries for different k values */
    kmer_hash = hash_any_extended((unsigned char*)&kmersearch_kmer_size, sizeof(int), 2);
    
    /* Combine hashes */
    return seq_hash ^ (query_hash << 1) ^ (kmer_hash << 2);
}

/*
 * Check if two sequences are exactly equal
 */
static bool
sequences_equal(VarBit *a, VarBit *b)
{
    if (VARBITLEN(a) != VARBITLEN(b))
        return false;
    if (VARSIZE(a) != VARSIZE(b))
        return false;
    return memcmp(VARBITS(a), VARBITS(b), VARBITBYTES(a)) == 0;
}

/*
 * Helper functions for min-heap operations
 */
static void
rawscore_heap_swap(RawscoreCacheManager *manager, int i, int j)
{
    RawscoreCacheEntry *temp = manager->min_heap[i];
    manager->min_heap[i] = manager->min_heap[j];
    manager->min_heap[j] = temp;
    
    /* Update heap indices */
    manager->min_heap[i]->heap_index = i;
    manager->min_heap[j]->heap_index = j;
}

static void
rawscore_heap_bubble_up(RawscoreCacheManager *manager, int index)
{
    int parent;
    
    while (index > 0)
    {
        parent = (index - 1) / 2;
        if (manager->min_heap[index]->result.shared_count >= manager->min_heap[parent]->result.shared_count)
            break;
        
        rawscore_heap_swap(manager, index, parent);
        index = parent;
    }
}

static void
rawscore_heap_bubble_down(RawscoreCacheManager *manager, int index)
{
    int smallest = index;
    int left = 2 * index + 1;
    int right = 2 * index + 2;
    
    if (left < manager->heap_size && 
        manager->min_heap[left]->result.shared_count < manager->min_heap[smallest]->result.shared_count)
        smallest = left;
    
    if (right < manager->heap_size && 
        manager->min_heap[right]->result.shared_count < manager->min_heap[smallest]->result.shared_count)
        smallest = right;
    
    if (smallest != index)
    {
        rawscore_heap_swap(manager, index, smallest);
        rawscore_heap_bubble_down(manager, smallest);
    }
}

/*
 * Insert entry into min-heap
 */
static void
rawscore_heap_insert(RawscoreCacheManager *manager, RawscoreCacheEntry *entry)
{
    if (manager->heap_size >= manager->max_entries)
        return;  /* Heap is full */
    
    /* Add to end of heap */
    manager->min_heap[manager->heap_size] = entry;
    entry->heap_index = manager->heap_size;
    manager->heap_size++;
    
    /* Bubble up to maintain heap property */
    rawscore_heap_bubble_up(manager, entry->heap_index);
}

/*
 * Remove entry from min-heap
 */
static void
rawscore_heap_remove(RawscoreCacheManager *manager, RawscoreCacheEntry *entry)
{
    int index = entry->heap_index;
    
    if (index < 0 || index >= manager->heap_size)
        return;  /* Entry not in heap */
    
    /* Move last element to this position */
    manager->heap_size--;
    if (index != manager->heap_size)
    {
        manager->min_heap[index] = manager->min_heap[manager->heap_size];
        manager->min_heap[index]->heap_index = index;
        
        /* Restore heap property */
        rawscore_heap_bubble_up(manager, index);
        rawscore_heap_bubble_down(manager, index);
    }
    
    entry->heap_index = -1;
}

/*
 * Evict entry with lowest rawscore
 */
static void
rawscore_heap_evict_lowest_score(RawscoreCacheManager *manager)
{
    RawscoreCacheEntry *lowest;
    bool found;
    
    if (manager->heap_size == 0)
        return;
    
    /* Root of min-heap has lowest score */
    lowest = manager->min_heap[0];
    
    /* Remove from heap */
    rawscore_heap_remove(manager, lowest);
    
    /* Remove from hash table */
    hash_search(manager->hash_table, &lowest->hash_key, HASH_REMOVE, &found);
    
    manager->current_entries--;
}

/*
 * Create rawscore cache manager
 */
static RawscoreCacheManager *
create_rawscore_cache_manager(const char *name)
{
    RawscoreCacheManager *manager;
    HASHCTL hash_ctl;
    Size heap_size;
    
    /* Allocate manager in current context (should be TopMemoryContext) */
    manager = (RawscoreCacheManager *) palloc0(sizeof(RawscoreCacheManager));
    
    /* Create dedicated memory context for rawscore cache */
    manager->cache_context = AllocSetContextCreate(CurrentMemoryContext,
                                                   "RawscoreCache",
                                                   ALLOCSET_DEFAULT_SIZES);
    
    /* Initialize rawscore cache parameters */
    manager->max_entries = kmersearch_rawscore_cache_max_entries;
    manager->current_entries = 0;
    manager->hits = 0;
    manager->misses = 0;
    manager->heap_size = 0;
    
    /* Allocate min-heap array in dedicated context */
    heap_size = (Size) manager->max_entries * sizeof(RawscoreCacheEntry *);
    manager->min_heap = (RawscoreCacheEntry **) MemoryContextAllocZero(manager->cache_context, heap_size);
    
    /* Create hash table using dedicated context */
    MemSet(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64);
    hash_ctl.entrysize = sizeof(RawscoreCacheEntry);
    hash_ctl.hcxt = manager->cache_context;
    
    manager->hash_table = hash_create(name, 1024, &hash_ctl,
                                     HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
    
    return manager;
}

/*
 * Free rawscore cache manager
 */
static void
free_rawscore_cache_manager(RawscoreCacheManager **manager)
{
    if (*manager)
    {
        /* Delete the dedicated cache context, which will free all allocated memory */
        if ((*manager)->cache_context)
            MemoryContextDelete((*manager)->cache_context);
        
        /* Free the manager itself (allocated in parent context) */
        pfree(*manager);
        *manager = NULL;
    }
}

/*
 * Lookup rawscore cache entry
 */
static RawscoreCacheEntry *
lookup_rawscore_cache_entry(RawscoreCacheManager *manager, VarBit *sequence, const char *query_string)
{
    uint64 hash_key;
    RawscoreCacheEntry *entry;
    bool found;
    
    /* Safety checks */
    if (manager == NULL || manager->hash_table == NULL || sequence == NULL || query_string == NULL) {
        elog(LOG, "lookup_cache_entry: NULL parameter detected");
        return NULL;
    }
    
    elog(LOG, "lookup_cache_entry: Looking up cache for query '%s'", query_string);
    hash_key = generate_cache_key(sequence, query_string);
    elog(LOG, "lookup_cache_entry: Generated hash key %lu", hash_key);
    
    /* Check if hash key generation failed */
    if (hash_key == 0) {
        elog(LOG, "lookup_cache_entry: Invalid hash key, skipping lookup");
        return NULL;
    }
    
    entry = (RawscoreCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_FIND, &found);
    elog(LOG, "lookup_cache_entry: Hash search completed, found=%d", found);
    
    if (found && entry != NULL && entry->sequence_copy != NULL && entry->query_string_copy != NULL &&
        sequences_equal(entry->sequence_copy, sequence) &&
        strcmp(entry->query_string_copy, query_string) == 0)
    {
        /* Cache hit - no need to update position in score-based cache */
        elog(LOG, "lookup_cache_entry: Cache hit found");
        manager->hits++;
        return entry;
    }
    
    elog(LOG, "lookup_cache_entry: Cache miss");
    return NULL;  /* Cache miss */
}

/*
 * Store rawscore entry in cache
 */
static void
store_rawscore_cache_entry(RawscoreCacheManager *manager, uint64 hash_key, VarBit *sequence, 
                          VarBit **query_keys, const char *query_string, KmerMatchResult result)
{
    RawscoreCacheEntry *entry;
    MemoryContext old_context;
    bool found;
    int actual_min_score;
    
    /* Safety checks */
    if (manager == NULL || manager->hash_table == NULL || sequence == NULL || query_string == NULL || hash_key == 0) {
        elog(WARNING, "store_rawscore_cache_entry: Invalid parameters, skipping cache storage");
        return;
    }
    
    /* Check if this result is worth caching based on actual minimum score */
    actual_min_score = get_cached_actual_min_score(query_keys, result.query_nkeys);
    if (result.shared_count < actual_min_score)
    {
        /* Don't cache results that can never match the condition */
        return;
    }
    
    /* Check if we need to evict */
    while (manager->current_entries >= manager->max_entries)
        rawscore_heap_evict_lowest_score(manager);
    
    old_context = MemoryContextSwitchTo(manager->cache_context);
    
    /* Create new entry */
    entry = (RawscoreCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_ENTER, &found);
    
    if (!found)
    {
        /* New entry */
        entry->hash_key = hash_key;
        entry->sequence_copy = (VarBit *) palloc(VARSIZE(sequence));
        memcpy(entry->sequence_copy, sequence, VARSIZE(sequence));
        entry->query_string_copy = pstrdup(query_string);
        entry->result = result;
        entry->heap_index = -1;
        
        /* Add to min-heap */
        rawscore_heap_insert(manager, entry);
        manager->current_entries++;
    }
    
    MemoryContextSwitchTo(old_context);
}

/*
 * Get cached rawscore result for DNA2
 */
KmerMatchResult
get_cached_rawscore_dna2(VarBit *sequence, const char *query_string)
{
    KmerMatchResult result;
    RawscoreCacheEntry *cache_entry;
    uint64 cache_key;
    VarBit **query_keys;
    int i;
    
    /* Create global cache manager if not exists */
    if (rawscore_cache_manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        rawscore_cache_manager = create_rawscore_cache_manager("GlobalRawscoreCache");
        MemoryContextSwitchTo(old_context);
        /* Update global statistics with max entries for both DNA2 and DNA4 */
        rawscore_cache_stats.dna2_max_entries = rawscore_cache_manager->max_entries;
        rawscore_cache_stats.dna4_max_entries = rawscore_cache_manager->max_entries;
    }
    
    /* Try cache lookup first */
    cache_entry = lookup_rawscore_cache_entry(rawscore_cache_manager, sequence, query_string);
    if (cache_entry != NULL)
    {
        /* Cache hit - return cached result */
        rawscore_cache_manager->hits++;
        rawscore_cache_stats.dna2_hits++;
        return cache_entry->result;
    }
    
    /* Cache miss - calculate result */
    rawscore_cache_manager->misses++;
    rawscore_cache_stats.dna2_misses++;
    result = kmersearch_calculate_kmer_match_and_score_dna2(sequence, query_string);
    
    /* Store result in cache if valid */
    if (result.valid)
    {
        cache_key = generate_cache_key(sequence, query_string);
        if (cache_key != 0)  /* Only proceed if cache key is valid */
        {
            query_keys = kmersearch_extract_kmer_from_query(query_string, kmersearch_kmer_size, &result.query_nkeys);
            if (query_keys != NULL)
        {
            store_rawscore_cache_entry(rawscore_cache_manager, cache_key, sequence, query_keys, query_string, result);
            /* Update global current entries count for both DNA2 and DNA4 since they share the same cache */
            rawscore_cache_stats.dna2_current_entries = rawscore_cache_manager->current_entries;
            rawscore_cache_stats.dna4_current_entries = rawscore_cache_manager->current_entries;
            /* Free query_keys */
            for (i = 0; i < result.query_nkeys; i++)
            {
                if (query_keys[i])
                    pfree(query_keys[i]);
            }
            pfree(query_keys);
            }
        }
    }
    
    return result;
}

/*
 * Get cached rawscore result for DNA4
 */
KmerMatchResult
get_cached_rawscore_dna4(VarBit *sequence, const char *query_string)
{
    KmerMatchResult result;
    RawscoreCacheEntry *cache_entry;
    uint64 cache_key;
    VarBit **query_keys;
    int i;
    
    /* Create global cache manager if not exists (shared with DNA2) */
    if (rawscore_cache_manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        rawscore_cache_manager = create_rawscore_cache_manager("GlobalRawscoreCache");
        MemoryContextSwitchTo(old_context);
        /* Update global statistics with max entries for both DNA2 and DNA4 */
        rawscore_cache_stats.dna2_max_entries = rawscore_cache_manager->max_entries;
        rawscore_cache_stats.dna4_max_entries = rawscore_cache_manager->max_entries;
    }
    
    /* Try cache lookup first */
    cache_entry = lookup_rawscore_cache_entry(rawscore_cache_manager, sequence, query_string);
    if (cache_entry != NULL)
    {
        /* Cache hit - return cached result */
        rawscore_cache_manager->hits++;
        rawscore_cache_stats.dna4_hits++;
        return cache_entry->result;
    }
    
    /* Cache miss - calculate result */
    rawscore_cache_manager->misses++;
    rawscore_cache_stats.dna4_misses++;
    result = kmersearch_calculate_kmer_match_and_score_dna4(sequence, query_string);
    
    /* Store result in cache if valid */
    if (result.valid)
    {
        cache_key = generate_cache_key(sequence, query_string);
        if (cache_key != 0)  /* Only proceed if cache key is valid */
        {
            query_keys = kmersearch_extract_kmer_from_query(query_string, kmersearch_kmer_size, &result.query_nkeys);
            if (query_keys != NULL)
        {
            store_rawscore_cache_entry(rawscore_cache_manager, cache_key, sequence, query_keys, query_string, result);
            /* Update global current entries count for both DNA2 and DNA4 since they share the same cache */
            rawscore_cache_stats.dna2_current_entries = rawscore_cache_manager->current_entries;
            rawscore_cache_stats.dna4_current_entries = rawscore_cache_manager->current_entries;
            /* Free query_keys */
            for (i = 0; i < result.query_nkeys; i++)
            {
                if (query_keys[i])
                    pfree(query_keys[i]);
            }
            pfree(query_keys);
            }
        }
    }
    
    return result;
}

/*
 * SQL interface function for rawscore cache statistics
 */
Datum
kmersearch_rawscore_cache_stats(PG_FUNCTION_ARGS)
{
    TupleDesc tupdesc;
    Datum values[8];
    bool nulls[8] = {false};
    HeapTuple tuple;
    
    /* Build tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("function returning record called in context that cannot accept a set")));
    
    /* DNA2 cache statistics */
    values[0] = Int64GetDatum(rawscore_cache_stats.dna2_hits);  /* hits */
    values[1] = Int64GetDatum(rawscore_cache_stats.dna2_misses);  /* misses */
    values[2] = Int32GetDatum(rawscore_cache_stats.dna2_current_entries);  /* current_entries */
    values[3] = Int32GetDatum(rawscore_cache_stats.dna2_max_entries);  /* max_entries */
    
    /* DNA4 cache statistics */
    values[4] = Int64GetDatum(rawscore_cache_stats.dna4_hits);  /* hits */
    values[5] = Int64GetDatum(rawscore_cache_stats.dna4_misses);  /* misses */
    values[6] = Int32GetDatum(rawscore_cache_stats.dna4_current_entries);  /* current_entries */
    values[7] = Int32GetDatum(rawscore_cache_stats.dna4_max_entries);  /* max_entries */
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}

/*
 * SQL interface function for freeing rawscore cache
 */
Datum
kmersearch_rawscore_cache_free(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    /* Count current entries before clearing */
    if (rawscore_cache_manager)
        freed_entries = rawscore_cache_manager->current_entries;
    
    /* Free the global rawscore cache manager completely */
    if (rawscore_cache_manager)
        free_rawscore_cache_manager(&rawscore_cache_manager);
    
    /* Reset all rawscore cache statistics */
    rawscore_cache_stats.dna2_hits = 0;
    rawscore_cache_stats.dna2_misses = 0;
    rawscore_cache_stats.dna2_current_entries = 0;
    rawscore_cache_stats.dna2_max_entries = 0;
    rawscore_cache_stats.dna4_hits = 0;
    rawscore_cache_stats.dna4_misses = 0;
    rawscore_cache_stats.dna4_current_entries = 0;
    rawscore_cache_stats.dna4_max_entries = 0;
    
    PG_RETURN_INT32(freed_entries);
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
 * Internal function to free query pattern cache (without PostgreSQL function wrapper)
 */
void
kmersearch_free_query_pattern_cache_internal(void)
{
    free_query_pattern_cache_manager(&query_pattern_cache_manager);
}

/*
 * Internal function to free actual min score cache (without PostgreSQL function wrapper)
 */
void
kmersearch_free_actual_min_score_cache_internal(void)
{
    free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
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
    global_highfreq_cache.current_table_oid = InvalidOid;
    global_highfreq_cache.current_column_name = NULL;
    global_highfreq_cache.current_kmer_size = 0;
    
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
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Started with table_oid=%u, column_name=%s, k_value=%d",
                           table_oid, column_name ? column_name : "NULL", k_value)));
    
    if (!column_name || k_value <= 0) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Invalid parameters, returning false")));
        return false;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Parameters validated successfully")));
    
    /* Initialize cache if not already done */
    if (global_highfreq_cache.cache_context == NULL) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache context is NULL, initializing")));
        kmersearch_highfreq_kmer_cache_init();
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache initialization completed")));
    } else {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cache context already exists")));
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to validate GUC settings")));
    
    /* Validate current GUC settings against metadata table */
    if (!kmersearch_validate_guc_against_metadata(table_oid, column_name, k_value)) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: GUC validation failed, returning false")));
        return false;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: GUC validation passed")));
    
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
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to get high-frequency k-mers from table")));
    
    /* Get high-frequency k-mers list */
    highfreq_kmers = kmersearch_get_highfreq_kmer_from_table(table_oid, column_name, k_value, &highfreq_count);
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Retrieved %d high-frequency k-mers",
                           highfreq_count)));
    
    if (!highfreq_kmers || highfreq_count <= 0) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: No high-frequency k-mers found, cache remains invalid")));
        return false;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to switch to cache context")));
    
    /* Switch to cache context for cache storage */
    old_context = MemoryContextSwitchTo(global_highfreq_cache.cache_context);
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Switched to cache context successfully")));
    
    /* Store in cache */
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to store table_oid in cache")));
    global_highfreq_cache.current_table_oid = table_oid;
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to store column_name in cache")));
    global_highfreq_cache.current_column_name = pstrdup(column_name);  /* Cache context copy */
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to store k_value in cache")));
    global_highfreq_cache.current_kmer_size = k_value;
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to copy highfreq_kmers array to cache context")));
    
    /* Copy the k-mer array to cache context to prevent memory context issues */
    cache_kmers = (VarBit **) palloc(highfreq_count * sizeof(VarBit *));
    for (i = 0; i < highfreq_count; i++) {
        if (highfreq_kmers[i]) {
            /* Copy each k-mer to the cache context */
            cache_kmers[i] = DatumGetVarBitPCopy(PointerGetDatum(highfreq_kmers[i]));
            ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Copied k-mer %d to cache context", i+1)));
        } else {
            cache_kmers[i] = NULL;
            ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: K-mer %d is NULL, copied as NULL", i+1)));
        }
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to store copied highfreq_kmers array in cache")));
    global_highfreq_cache.highfreq_kmers = cache_kmers;
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to store highfreq_count in cache")));
    global_highfreq_cache.highfreq_count = highfreq_count;
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to create hash table from array")));
    
    /* Create hash table in cache context */
    global_highfreq_cache.highfreq_hash = kmersearch_create_highfreq_hash_from_array(cache_kmers, highfreq_count);
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Hash table creation completed, checking result")));
    
    if (global_highfreq_cache.highfreq_hash) {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Hash table created successfully, setting cache as valid")));
        global_highfreq_cache.is_valid = true;
    } else {
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Hash table creation failed, cleaning up")));
        /* Hash table creation failed, clean up by deleting the context */
        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(global_highfreq_cache.cache_context);
        global_highfreq_cache.cache_context = NULL;
        global_highfreq_cache.is_valid = false;
        ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Cleanup completed, returning false")));
        return false;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: About to switch back to old context")));
    MemoryContextSwitchTo(old_context);
    
    ereport(DEBUG1, (errmsg("kmersearch_highfreq_kmer_cache_load_internal: Function completed successfully, cache is_valid=%d",
                           global_highfreq_cache.is_valid)));
    
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
    global_highfreq_cache.current_table_oid = InvalidOid;
    global_highfreq_cache.current_kmer_size = 0;
    global_highfreq_cache.highfreq_count = 0;
    global_highfreq_cache.highfreq_hash = NULL;
    global_highfreq_cache.highfreq_kmers = NULL;
    global_highfreq_cache.current_column_name = NULL;
}

bool
kmersearch_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value)
{
    return (global_highfreq_cache.is_valid &&
            global_highfreq_cache.current_table_oid == table_oid &&
            global_highfreq_cache.current_kmer_size == k_value &&
            global_highfreq_cache.current_column_name &&
            column_name &&
            strcmp(global_highfreq_cache.current_column_name, column_name) == 0);
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
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k_value = PG_GETARG_INT32(2);
    char *column_name = text_to_cstring(column_name_text);
    bool success;
    
    success = kmersearch_highfreq_kmer_cache_load_internal(table_oid, column_name, k_value);
    
    PG_RETURN_BOOL(success);
}

/*
 * SQL-accessible high-frequency cache free function
 */
Datum
kmersearch_highfreq_kmer_cache_free(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    /* Count entries before freeing */
    if (global_highfreq_cache.is_valid)
        freed_entries = global_highfreq_cache.highfreq_count;
    
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
    
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Started with table_oid=%u, column_name=%s, k_value=%d",
                           table_oid, column_name ? column_name : "NULL", k_value)));
    
    /* Connect to SPI */
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: About to connect to SPI")));
    if (SPI_connect() != SPI_OK_CONNECT) {
        ereport(ERROR, (errmsg("kmersearch_validate_guc_against_metadata: SPI_connect failed")));
    }
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: SPI connected successfully")));
    
    /* Build query to get metadata */
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Building metadata query")));
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT occur_bitlen, max_appearance_rate, max_appearance_nrow "
        "FROM kmersearch_highfreq_kmer_meta "
        "WHERE table_oid = %u AND column_name = '%s' AND k_value = %d",
        table_oid, column_name, k_value);
    
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Query built: %s", query.data)));
    
    /* Execute query */
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: About to execute query")));
    ret = SPI_execute(query.data, true, 1);
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Query executed, ret=%d, processed=%zu", ret, SPI_processed)));
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        bool isnull;
        Datum occur_bitlen_datum, max_appearance_rate_datum, max_appearance_nrow_datum;
        int stored_occur_bitlen, stored_max_appearance_nrow;
        float stored_max_appearance_rate;
        
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Found metadata record, extracting values")));
        
        /* Get stored metadata values */
        ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Getting occur_bitlen value")));
        occur_bitlen_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
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
        max_appearance_rate_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull);
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
        max_appearance_nrow_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull);
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
                 errmsg("No metadata found for table_oid=%u, column_name='%s', k_value=%d",
                       table_oid, column_name, k_value),
                 errhint("Run kmersearch_analyze_table() first to create metadata.")));
        validation_passed = false;
    }
    
    /* Cleanup */
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: About to cleanup and finish SPI")));
    pfree(query.data);
    SPI_finish();
    
    ereport(DEBUG1, (errmsg("kmersearch_validate_guc_against_metadata: Function completed, validation_passed=%d", validation_passed)));
    
    return validation_passed;
}

VarBit **
kmersearch_get_highfreq_kmer_from_table(Oid table_oid, const char *column_name, int k, int *nkeys)
{
    VarBit **result = NULL;
    int ret;
    StringInfoData query;
    int i;
    
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Started with table_oid=%u, column_name=%s, k=%d",
                           table_oid, column_name ? column_name : "NULL", k)));
    
    if (!nkeys) {
        ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: nkeys is NULL, returning NULL")));
        return NULL;
    }
    
    *nkeys = 0;
    
    /* Connect to SPI */
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: About to connect to SPI")));
    if (SPI_connect() != SPI_OK_CONNECT) {
        ereport(ERROR, (errmsg("kmersearch_get_highfreq_kmer_from_table: SPI_connect failed")));
    }
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: SPI connected successfully")));
    
    /* Build query to get highly frequent k-mers */
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Building query to get high-frequency ngram_key2")));
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT DISTINCT hkm.ngram_key FROM kmersearch_highfreq_kmer hkm "
        "WHERE hkm.index_oid IN ("
        "    SELECT indexrelid FROM pg_stat_user_indexes pui "
        "    JOIN pg_class pc ON pui.relid = pc.oid "
        "    WHERE pc.oid = %u "
        "    AND EXISTS ("
        "        SELECT 1 FROM kmersearch_highfreq_kmer_meta hkm_meta "
        "        WHERE hkm_meta.table_oid = %u "
        "        AND hkm_meta.column_name = '%s' "
        "        AND hkm_meta.k_value = %d"
        "    )"
        ") "
        "ORDER BY hkm.ngram_key",
        table_oid, table_oid, column_name, k);
    
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Query built: %s", query.data)));
    
    /* Execute query */
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: About to execute query")));
    ret = SPI_execute(query.data, true, 0);
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Query executed, ret=%d, processed=%zu", ret, SPI_processed)));
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        *nkeys = SPI_processed;
        ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Found %d high-frequency k-mers, allocating result array", *nkeys)));
        result = (VarBit **) palloc(*nkeys * sizeof(VarBit *));
        
        ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Starting to process k-mer results")));
        for (i = 0; i < *nkeys; i++)
        {
            bool isnull;
            Datum kmer_datum;
            
            ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Processing k-mer %d/%d", i+1, *nkeys)));
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                /* Copy the varbit value */
                result[i] = DatumGetVarBitPCopy(kmer_datum);
                ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Successfully copied k-mer %d", i+1)));
            }
            else
            {
                result[i] = NULL;
                ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: K-mer %d is NULL", i+1)));
            }
        }
        ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Finished processing all k-mers")));
    }
    else
    {
        ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: No high-frequency k-mers found or query failed")));
    }
    
    /* Cleanup */
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: About to cleanup and finish SPI")));
    pfree(query.data);
    SPI_finish();
    
    ereport(DEBUG1, (errmsg("kmersearch_get_highfreq_kmer_from_table: Function completed, returning %d k-mers", *nkeys)));
    
    return result;
}

HTAB *
kmersearch_create_highfreq_hash_from_array(VarBit **kmers, int nkeys)
{
    HTAB *hash_table;
    HASHCTL hash_ctl;
    int i;
    
    ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Started with nkeys=%d", nkeys)));
    
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
    
    ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Hash control structure set up, about to create hash table")));
    
    hash_table = hash_create("HighfreqKmerHash",
                            nkeys,
                            &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_CONTEXT);
    
    if (!hash_table) {
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Hash table creation failed")));
        return NULL;
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Hash table created successfully, adding %d k-mers", nkeys)));
    
    /* Add each k-mer to the hash table */
    for (i = 0; i < nkeys; i++)
    {
        HighfreqKmerHashEntry *entry;
        bool found;
        uint64 hash_value;
        unsigned char *bits_ptr;
        int bytes_len;
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Processing k-mer %d/%d", i+1, nkeys)));
        
        /* Check if the k-mer pointer itself is valid before any access */
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Checking k-mer %d pointer: %p", i+1, kmers[i])));
        
        if (!kmers[i]) {
            ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d is NULL, skipping", i+1)));
            continue;
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d pointer is valid", i+1)));
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Calculating hash value for ngram_key2 %d", i+1)));
        
        /* Validate VarBit data before hash calculation */
        if (VARSIZE(kmers[i]) < VARHDRSZ) {
            ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d has invalid size %u, skipping", i+1, VARSIZE(kmers[i]))));
            continue;
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d size validation passed, size=%u", i+1, VARSIZE(kmers[i]))));
        
        /* Check VARBITS pointer */
        bits_ptr = (unsigned char *) VARBITS(kmers[i]);
        if (!bits_ptr) {
            ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d VARBITS returned NULL, skipping", i+1)));
            continue;
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d VARBITS pointer valid: %p", i+1, bits_ptr)));
        
        /* Calculate bytes length manually for ngram_key2 (more reliable than VARBITBYTES) */
        {
            int bit_length = VARBITLEN(kmers[i]);
            bytes_len = (bit_length + 7) / 8;  /* Round up to next byte */
        }
        
        if (bytes_len <= 0 || bytes_len > 1000) {  /* Reasonable upper limit */
            ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d has invalid bytes length %d, skipping", i+1, bytes_len)));
            continue;
        }
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d bytes length calculated: %d", i+1, bytes_len)));
        
        /* Calculate hash value for this ngram_key2 (kmer2 + occurrence bits) */
        hash_value = DatumGetUInt64(hash_any(bits_ptr, bytes_len));
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Hash value calculated for k-mer %d: %lu", i+1, hash_value)));
        
        entry = (HighfreqKmerHashEntry *) hash_search(hash_table,
                                                     (void *) &hash_value,
                                                     HASH_ENTER,
                                                     &found);
        
        ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Hash search completed for k-mer %d, entry=%p, found=%d", i+1, entry, found)));
        
        if (entry && !found)
        {
            entry->kmer_key = kmers[i];
            entry->hash_value = hash_value;
            ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Successfully added k-mer %d to hash table", i+1)));
        } else if (found) {
            ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: K-mer %d already exists in hash table", i+1)));
        } else {
            ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Failed to add k-mer %d to hash table", i+1)));
        }
    }
    
    ereport(DEBUG1, (errmsg("kmersearch_create_highfreq_hash_from_array: Function completed successfully, hash table created")));
    
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
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_PP(1);
    int k_value = PG_GETARG_INT32(2);
    char *column_name = text_to_cstring(column_name_text);
    bool result;
    
    /* Initialize parallel cache if not already done */
    if (parallel_highfreq_cache == NULL) {
        kmersearch_parallel_highfreq_kmer_cache_init();
    }
    
    /* Load cache data into DSM */
    result = kmersearch_parallel_highfreq_kmer_cache_load_internal(table_oid, column_name, k_value);
    
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
    int32 freed_entries = 0;
    
    ereport(LOG, (errmsg("kmersearch_parallel_highfreq_kmer_cache_free: Starting function call")));
    
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
    
    ereport(LOG, (errmsg("kmersearch_parallel_highfreq_kmer_cache_free: Function completed, returning %d", freed_entries)));
    
    PG_RETURN_INT32(freed_entries);
}