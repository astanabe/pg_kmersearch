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

/* PostgreSQL function info declarations for cache functions will be moved here later */

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
static VarBit **get_cached_query_kmer(const char *query_string, int k_size, int *nkeys);
static void free_query_pattern_cache_manager(QueryPatternCacheManager **manager);

/* Actual min score cache functions */
static void create_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);
static void free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);
static int get_cached_actual_min_score(VarBit **query_keys, int nkeys);

/* Rawscore cache functions */
static RawscoreCacheManager *create_rawscore_cache_manager(const char *name);
static void free_rawscore_cache_manager(RawscoreCacheManager **manager);
static uint64 generate_cache_key(VarBit *sequence, const char *query_string);
static bool sequences_equal(VarBit *a, VarBit *b);
static RawscoreCacheEntry *lookup_rawscore_cache_entry(RawscoreCacheManager *manager, VarBit *sequence, const char *query_string);
static void store_rawscore_cache_entry(RawscoreCacheManager *manager, uint64 hash_key, VarBit *sequence, VarBit **query_keys, const char *query_string, KmerMatchResult result);

/* Rawscore cache heap management */
static void rawscore_heap_swap(RawscoreCacheManager *manager, int i, int j);
static void rawscore_heap_bubble_up(RawscoreCacheManager *manager, int index);
static void rawscore_heap_bubble_down(RawscoreCacheManager *manager, int index);
static void rawscore_heap_insert(RawscoreCacheManager *manager, RawscoreCacheEntry *entry);
static void rawscore_heap_remove(RawscoreCacheManager *manager, RawscoreCacheEntry *entry);
static void rawscore_heap_evict_lowest_score(RawscoreCacheManager *manager);

/* High-level cache functions */
static KmerMatchResult get_cached_rawscore_dna2(VarBit *sequence, const char *query_string);
static KmerMatchResult get_cached_rawscore_dna4(VarBit *sequence, const char *query_string);

/* External global variables (defined in pg_kmersearch.c) */
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
static VarBit **
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
    
    /* Use adjusted minimum score that considers high-frequency k-mer filtering */
    absolute_min = kmersearch_get_adjusted_min_score(query_keys, nkeys);
    
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
static int
get_cached_actual_min_score(VarBit **query_keys, int nkeys)
{
    ActualMinScoreCacheEntry *cache_entry;
    uint64 query_hash;
    bool found;
    MemoryContext old_context;
    int actual_min_score;
    
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
    
    /* Calculate hash value for query keys */
    query_hash = hash_any_extended((unsigned char *)query_keys, nkeys * sizeof(VarBit *), 0);
    
    /* Look up in hash table */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, &found);
    
    if (found) {
        actual_min_score_cache_manager->hits++;
        return cache_entry->actual_min_score;
    }
    
    /* Not found - calculate and cache */
    actual_min_score_cache_manager->misses++;
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
    
    return actual_min_score;
}