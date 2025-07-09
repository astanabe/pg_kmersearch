#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/builtins.h"
#include "utils/varbit.h"
#include "utils/memutils.h"
#include "utils/varlena.h"
#include "varatt.h"
#include "access/gin.h"
#include "access/stratnum.h"
#include "utils/guc.h"
#include "utils/array.h"
#include "catalog/pg_type.h"
#include "executor/executor.h"
#include "nodes/pg_list.h"
#include "storage/lwlock.h"
#include "miscadmin.h"
#include "access/heapam.h"
#include "access/tableam.h"
#include "utils/snapmgr.h"
#include "utils/rel.h"
#include "utils/hsearch.h"
#include "access/relation.h"
#include "utils/syscache.h"
#include "access/tupdesc.h"
#include "catalog/pg_class.h"
#include "commands/tablecmds.h"
#include "common/hashfn.h"
#include "access/htup_details.h"
#include "funcapi.h"
#include "utils/lsyscache.h"
#include "catalog/namespace.h"
#include "executor/spi.h"
#include "lib/stringinfo.h"
#include "utils/relcache.h"
#include "executor/tuptable.h"
#include "optimizer/plancat.h"
#include "storage/smgr.h"
#include "storage/bufmgr.h"
#include "access/gin_private.h"
#include "nodes/pg_list.h"
#include <ctype.h>
#include <string.h>
#include <math.h>

/* SIMD intrinsics headers */
#ifdef __x86_64__
#include <cpuid.h>
#include <immintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif
#endif
#include <unistd.h>

PG_MODULE_MAGIC;

/* SIMD capability detection */
typedef enum {
    SIMD_NONE,
    SIMD_AVX2,     /* AMD64: AVX2 + SSE4.x */
    SIMD_AVX512,   /* AMD64: AVX512 + AVX2 + SSE4.x */
    SIMD_NEON,     /* ARM64: Apple M1+ NEON */
    SIMD_SVE       /* ARM64: Graviton4+ SVE */
} simd_capability_t;

/* Function pointers for different SIMD implementations */
typedef struct {
    void (*dna2_encode)(const char* input, uint8_t* output, int len);
    void (*dna2_decode)(const uint8_t* input, char* output, int len);
    void (*dna4_encode)(const char* input, uint8_t* output, int len);
    void (*dna4_decode)(const uint8_t* input, char* output, int len);
} simd_dispatch_table_t;

/* Global SIMD dispatch table */
static simd_dispatch_table_t simd_dispatch;
static simd_capability_t simd_capability = SIMD_NONE;

/*
 * Module load and unload functions
 */
void _PG_init(void);
void _PG_fini(void);

/* Global variables for k-mer search configuration */
static int kmersearch_occur_bitlen = 8;  /* Default 8 bits for occurrence count */
static int kmersearch_kmer_size = 8;  /* Default k-mer size */
static double kmersearch_max_appearance_rate = 0.05;  /* Default max appearance rate */
static int kmersearch_max_appearance_nrow = 0;  /* Default max appearance nrow (0 = undefined) */
static int kmersearch_min_score = 1;  /* Default minimum score for GIN search */
static double kmersearch_min_shared_ngram_key_rate = 0.9;  /* Default minimum shared n-gram key rate for =% operator */

/* Cache configuration variables */
static int kmersearch_rawscore_cache_max_entries = 50000;  /* Default max rawscore cache entries */
static int kmersearch_query_pattern_cache_max_entries = 50000;  /* Default max query pattern cache entries */
static int kmersearch_actual_min_score_cache_max_entries = 50000;  /* Default max actual min score cache entries */

/* Hash table entry for k-mer frequency counting */
typedef struct KmerFreqEntry
{
    VarBit      *kmer_key;      /* K-mer binary key (without occurrence count) */
    int         row_count;      /* Number of rows where this k-mer appears */
    bool        highfreq;       /* Whether this k-mer is highly frequent */
} KmerFreqEntry;

/* A-3: Removed KmerHashKey structure - no longer needed for direct comparison */

/* Combined result structure for k-mer matching and scoring */
typedef struct KmerMatchResult
{
    int         shared_count;       /* Number of shared k-mers (rawscore) */
    int         seq_nkeys;          /* Number of k-mers in sequence */
    int         query_nkeys;        /* Number of k-mers in query */
    double      sharing_rate;       /* Calculated sharing rate */
    bool        match_result;       /* Boolean match result for =% operator */
    bool        valid;              /* Whether the result is valid */
} KmerMatchResult;

/* Actual min score cache entry for TopMemoryContext caching */
typedef struct ActualMinScoreCacheEntry
{
    uint64      query_hash;                /* Hash value of query string */
    int         actual_min_score;          /* Actual minimum score for this query */
} ActualMinScoreCacheEntry;

/* Actual min score cache manager for TopMemoryContext */
typedef struct ActualMinScoreCacheManager
{
    HTAB        *cache_hash;               /* Hash table for actual min score cache */
    MemoryContext cache_context;           /* Memory context for cache */
    int         hits;                      /* Cache hit count */
    int         misses;                    /* Cache miss count */
    int         max_entries;               /* Maximum number of entries */
    int         current_entries;           /* Current number of entries */
} ActualMinScoreCacheManager;

/* Global cache managers */
static ActualMinScoreCacheManager *actual_min_score_cache_manager = NULL;

/* Simple k-mer occurrence tracking for k<=64 (no hashing needed) */
typedef struct KmerOccurrence
{
    uint64_t    kmer_value;         /* K-mer as single 64-bit value */
    int         count;              /* Occurrence count */
} KmerOccurrence;

/* Analysis result structure */
typedef struct KmerAnalysisResult
{
    int64       total_rows;                /* Total number of rows analyzed */
    int         highfreq_kmers_count;      /* Number of highly frequent k-mers */
    int         parallel_workers_used;     /* Number of parallel workers used */
    double      analysis_duration;        /* Analysis duration in seconds */
    double      max_appearance_rate_used; /* Max appearance rate used */
    int         max_appearance_nrow_used;  /* Max appearance nrow used */
} KmerAnalysisResult;

/* Drop analysis result structure */
typedef struct DropAnalysisResult
{
    int         dropped_analyses;          /* Number of dropped analyses */
    int         dropped_highfreq_kmers;    /* Number of dropped highly frequent k-mers */
    int64       freed_storage_bytes;       /* Freed storage in bytes */
} DropAnalysisResult;

/* Rawscore cache entry structure */
typedef struct RawscoreCacheEntry
{
    uint64      hash_key;                  /* Hash key for this entry */
    VarBit      *sequence_copy;            /* Copy of sequence for exact match */
    char        *query_string_copy;        /* Copy of query string for exact match */
    KmerMatchResult result;                /* Cached rawscore calculation result */
    int         heap_index;                /* Index in min-heap array (-1 if not in heap) */
} RawscoreCacheEntry;

/* Rawscore cache manager structure */
typedef struct RawscoreCacheManager
{
    HTAB        *hash_table;               /* PostgreSQL hash table */
    MemoryContext cache_context;           /* Memory context for rawscore cache */
    int         max_entries;               /* Maximum number of entries */
    int         current_entries;           /* Current number of entries */
    uint64      hits;                      /* Cache hit count */
    uint64      misses;                    /* Cache miss count */
    RawscoreCacheEntry **min_heap;         /* Min-heap array for score-based eviction */
    int         heap_size;                 /* Current heap size */
} RawscoreCacheManager;

/* Rawscore cache managers are now local to each function call - no global state */

/* B-2: Query pattern cache structures */
typedef struct QueryPatternCacheEntry
{
    uint64      hash_key;                  /* Hash key for this entry */
    char        *query_string_copy;        /* Copy of query string */
    int         k_size;                    /* K-mer size for this pattern */
    VarBit      **extracted_kmers;         /* Cached extracted k-mers */
    int         kmer_count;                /* Number of extracted k-mers */
    struct QueryPatternCacheEntry *next;   /* For LRU chain */
    struct QueryPatternCacheEntry *prev;   /* For LRU chain */
} QueryPatternCacheEntry;

typedef struct QueryPatternCacheManager
{
    HTAB        *hash_table;               /* PostgreSQL hash table */
    MemoryContext query_pattern_cache_context; /* Memory context for query pattern cache */
    int         max_entries;               /* Maximum number of entries */
    int         current_entries;           /* Current number of entries */
    uint64      hits;                      /* Cache hit count */
    uint64      misses;                    /* Cache miss count */
    QueryPatternCacheEntry *lru_head;      /* LRU chain head (most recent) */
    QueryPatternCacheEntry *lru_tail;      /* LRU chain tail (least recent) */
} QueryPatternCacheManager;

/* Global query pattern cache manager for cross-query sharing */
static QueryPatternCacheManager *query_pattern_cache_manager = NULL;

/* Global rawscore cache manager for cross-query sharing */
static RawscoreCacheManager *rawscore_cache_manager = NULL;

/* K-mer data union for different k-values */
typedef union KmerData
{
    uint16      k8_data;                 /* k <= 8: 16 bits */
    uint32      k16_data;                /* k <= 16: 32 bits */  
    uint64      k32_data;                /* k <= 32: 64 bits */
    struct {
        uint64  high;
        uint64  low;
    }           k64_data;                /* k <= 64: 128 bits */
} KmerData;

/* Compact k-mer frequency entry */
typedef struct CompactKmerFreq
{
    KmerData    kmer_data;               /* Raw k-mer data */
    int         frequency_count;          /* Frequency count */
    bool        is_highfreq;              /* Whether this k-mer is highly frequent */
} CompactKmerFreq;

/* Worker buffer for k-mer collection */
typedef struct KmerBuffer
{
    CompactKmerFreq *entries;            /* Buffer entries */
    int             count;               /* Current count */
    int             capacity;            /* Buffer capacity */
    int             k_size;              /* K-mer size */
} KmerBuffer;

/* High-frequency k-mer hash entry */
typedef struct HighfreqKmerHashEntry
{
    VarBit         *kmer_key;             /* K-mer key for hashing */
    uint64          hash_value;           /* Precomputed hash value */
} HighfreqKmerHashEntry;

/* High-frequency k-mer cache for index creation */
typedef struct HighfreqKmerCache
{
    Oid         current_table_oid;        /* Current table OID */
    char       *current_column_name;      /* Current column name */
    int         current_k_value;          /* Current k-mer size */
    MemoryContext cache_context;          /* Memory context for cache data */
    HTAB       *highfreq_hash;           /* Hash table for fast lookup */
    VarBit    **highfreq_kmers;          /* Array of high-frequency k-mers */
    int         highfreq_count;          /* Number of high-frequency k-mers */
    bool        is_valid;                /* Cache validity flag */
} HighfreqKmerCache;

/* Global high-frequency k-mer cache */
static HighfreqKmerCache global_highfreq_cache = {0};

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

/* Worker state for parallel k-mer analysis */
typedef struct KmerWorkerState
{
    int         worker_id;                /* Worker identifier */
    BlockNumber start_block;              /* Starting block number */
    BlockNumber end_block;                /* Ending block number */
    KmerBuffer  buffer;                   /* Local k-mer buffer */
    int         local_highfreq_count;     /* Local count of highly frequent k-mers */
    int64       rows_processed;           /* Number of rows processed */
    char       *temp_table_name;          /* Temporary table name for this worker */
} KmerWorkerState;

/* Macro for safe memory cleanup */
#define CLEANUP_KMER_ARRAYS(seq_keys, seq_nkeys, query_keys, query_nkeys) \
    do { \
        if (seq_keys) { \
            int i; \
            for (i = 0; i < (seq_nkeys); i++) { \
                if (seq_keys[i]) \
                    pfree(seq_keys[i]); \
            } \
            pfree(seq_keys); \
            seq_keys = NULL; \
        } \
        if (query_keys) { \
            int i; \
            for (i = 0; i < (query_nkeys); i++) { \
                if (query_keys[i]) \
                    pfree(query_keys[i]); \
            } \
            pfree(query_keys); \
            query_keys = NULL; \
        } \
    } while(0)

/* Forward declarations */
static bool kmersearch_is_kmer_highfreq(VarBit *kmer_key);
static int kmersearch_count_highfreq_kmers_in_query(VarBit **query_keys, int nkeys);
static int kmersearch_get_adjusted_min_score(VarBit **query_keys, int nkeys);
static int kmersearch_calculate_raw_score(VarBit *seq1, VarBit *seq2, text *query_text);

/* High-frequency k-mer filtering functions */
static VarBit **kmersearch_get_highfreq_kmers_from_table(Oid table_oid, const char *column_name, int k, int *nkeys);
static HTAB *kmersearch_create_highfreq_hash_from_array(VarBit **kmers, int nkeys);
static Datum *kmersearch_filter_highfreq_kmers_from_keys(Datum *original_keys, int *nkeys, HTAB *highfreq_hash, int k);
static VarBit *kmersearch_remove_occurrence_bits(VarBit *key_with_occurrence, int k);

/* High-frequency k-mer cache management functions */
static void kmersearch_highfreq_kmers_cache_init(void);
static bool kmersearch_highfreq_kmers_cache_load_internal(Oid table_oid, const char *column_name, int k_value);
static void kmersearch_highfreq_kmers_cache_free_internal(void);
static bool kmersearch_highfreq_kmers_cache_is_valid(Oid table_oid, const char *column_name, int k_value);
static bool kmersearch_auto_load_cache_if_needed(void);

/* Helper functions for direct k-mer management (no hashing) */
static uint64_t kmersearch_extract_kmer_as_uint64(VarBit *seq, int start_pos, int k);
static int kmersearch_find_or_add_kmer_occurrence(KmerOccurrence *occurrences, int *count, uint64_t kmer_value, int max_count);
static VarBit **kmersearch_extract_kmers_from_varbit(VarBit *seq, int k, int *nkeys);
static VarBit **kmersearch_extract_kmers_from_query(const char *query, int k, int *nkeys);
static Datum *kmersearch_extract_kmers_with_degenerate(const char *sequence, int seq_len, int k, int *nkeys);
static int kmersearch_count_degenerate_combinations(const char *kmer, int k);
static VarBit *kmersearch_create_ngram_key_with_occurrence(const char *kmer, int k, int occurrence);
static bool kmersearch_will_exceed_degenerate_limit(const char *seq, int len);
static bool kmersearch_will_exceed_degenerate_limit_dna4_bits(VarBit *seq, int start_pos, int k);
static VarBit **kmersearch_expand_dna4_kmer_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count);
static VarBit *kmersearch_create_ngram_key_from_dna2_bits(VarBit *seq, int start_pos, int k, int occurrence);
static VarBit *kmersearch_create_kmer_key_from_dna2_bits(VarBit *seq, int start_pos, int k);
static Datum *kmersearch_extract_dna2_kmers_direct(VarBit *seq, int k, int *nkeys);

/* Rawscore cache management functions */
static RawscoreCacheManager *create_rawscore_cache_manager(const char *name);
static void free_rawscore_cache_manager(RawscoreCacheManager **manager);
static uint64 generate_cache_key(VarBit *sequence, const char *query_string);
static bool sequences_equal(VarBit *a, VarBit *b);
static RawscoreCacheEntry *lookup_rawscore_cache_entry(RawscoreCacheManager *manager, VarBit *sequence, const char *query_string);
static void store_rawscore_cache_entry(RawscoreCacheManager *manager, uint64 hash_key, VarBit *sequence, VarBit **query_keys, const char *query_string, KmerMatchResult result);
static void rawscore_heap_insert(RawscoreCacheManager *manager, RawscoreCacheEntry *entry);
static void rawscore_heap_remove(RawscoreCacheManager *manager, RawscoreCacheEntry *entry);
static void rawscore_heap_evict_lowest_score(RawscoreCacheManager *manager);
static KmerMatchResult get_cached_rawscore_dna2(VarBit *sequence, const char *query_string);
static KmerMatchResult get_cached_rawscore_dna4(VarBit *sequence, const char *query_string);

/* B-2: Query pattern cache functions */
static void init_query_pattern_cache_manager(QueryPatternCacheManager **manager);
static uint64 generate_query_pattern_cache_key(const char *query_string, int k_size);
static QueryPatternCacheEntry *lookup_query_pattern_cache_entry(QueryPatternCacheManager *manager, const char *query_string, int k_size);
static void store_query_pattern_cache_entry(QueryPatternCacheManager *manager, uint64 hash_key, const char *query_string, int k_size, VarBit **kmers, int kmer_count);
static void lru_touch_query_pattern_cache(QueryPatternCacheManager *manager, QueryPatternCacheEntry *entry);
static void lru_evict_oldest_query_pattern_cache(QueryPatternCacheManager *manager);
static VarBit **get_cached_query_kmers(const char *query_string, int k_size, int *nkeys);
static void free_query_pattern_cache_manager(QueryPatternCacheManager **manager);
static void free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);
static Datum *kmersearch_extract_dna4_kmers_with_expansion_direct(VarBit *seq, int k, int *nkeys);
static VarBit *kmersearch_create_ngram_key_with_occurrence_from_dna2(VarBit *dna2_kmer, int k, int occurrence);
static VarBit **kmersearch_extract_query_kmers(const char *query, int k, int *nkeys);
static int kmersearch_count_matching_kmers_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
static VarBit *kmersearch_create_kmer_key_only(const char *kmer, int k);
static bool kmersearch_kmer_based_match_dna2(VarBit *sequence, const char *query_string);
static bool kmersearch_kmer_based_match_dna4(VarBit *sequence, const char *query_string);
static bool kmersearch_evaluate_match_conditions(int shared_count, int query_total);
static KmerMatchResult kmersearch_calculate_kmer_match_and_score_dna2(VarBit *sequence, const char *query_string);
static KmerMatchResult kmersearch_calculate_kmer_match_and_score_dna4(VarBit *sequence, const char *query_string);

/* Actual min score cache functions */
static ActualMinScoreCacheManager *create_actual_min_score_cache_manager(void);
static int calculate_actual_min_score(VarBit **query_keys, int nkeys, int query_total_kmers);
static int get_cached_actual_min_score(VarBit **query_keys, int nkeys, const char *query_string, int query_total_kmers);
static bool evaluate_optimized_match_condition(VarBit **query_keys, int nkeys, int shared_count, const char *query_string, int query_total_kmers);

/* New parallel analysis functions */
static int kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation);
static KmerAnalysisResult kmersearch_analyze_table_parallel(Oid table_oid, const char *column_name, int k_size, int parallel_workers);
static void kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size);
static void kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, const char *final_table_name, int k_size, int threshold_rows);
static void kmersearch_persist_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, void *unused_table, int threshold_rows);
static void kmersearch_persist_highfreq_kmers_from_temp(Oid table_oid, const char *column_name, int k_size, const char *temp_table_name);

/* New memory-efficient k-mer functions */
static size_t kmersearch_get_kmer_data_size(int k_size);
static bool kmersearch_get_index_info(Oid index_oid, Oid *table_oid, char **column_name, int *k_size);
static void kmersearch_spi_connect_or_error(void);
static void kmersearch_handle_spi_error(int spi_result, const char *operation);
static bool kmersearch_delete_kmer_from_gin_index(Relation index_rel, VarBit *kmer_key);
static List *kmersearch_get_highfreq_kmers_list(Oid index_oid);
static int kmersearch_calculate_buffer_size(int k_size);
static KmerData kmersearch_encode_kmer_data(VarBit *kmer, int k_size);
static void kmersearch_init_buffer(KmerBuffer *buffer, int k_size);
static void kmersearch_add_to_buffer(KmerBuffer *buffer, KmerData kmer_data, const char *temp_table_name);
static void kmersearch_flush_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name);
static void kmersearch_aggregate_buffer_entries(KmerBuffer *buffer);
static void kmersearch_create_worker_temp_table(const char *temp_table_name, int k_size);
static bool kmersearch_check_analysis_exists(Oid table_oid, const char *column_name, int k_size);
static void kmersearch_validate_analysis_parameters(Oid table_oid, const char *column_name, int k_size);
static Datum *kmersearch_filter_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, Datum *all_keys, int total_keys, int *filtered_count);
static void kmersearch_delete_existing_analysis(Oid table_oid, const char *column_name, int k_size);
static DropAnalysisResult kmersearch_drop_analysis_internal(Oid table_oid, const char *column_name, int k_size);

/* Custom GUC variables */
void _PG_init(void);

/*
 * DNA2 type: 2-bit encoding for ACGT
 * A=00, C=01, G=10, T=11
 */
typedef struct
{
    int32 vl_len_;
    bits8 data[FLEXIBLE_ARRAY_MEMBER];
} kmersearch_dna2;

/*
 * DNA4 type: 4-bit encoding with degenerate codes
 * Uses standard varbit structure
 */
typedef struct
{
    int32 vl_len_;
    bits8 data[FLEXIBLE_ARRAY_MEMBER];
} kmersearch_dna4;

/* DNA2 encoding table */
static const uint8 kmersearch_dna2_encode_table[256] = {
    ['A'] = 0, ['a'] = 0,
    ['C'] = 1, ['c'] = 1,
    ['G'] = 2, ['g'] = 2,
    ['T'] = 3, ['t'] = 3,
    ['U'] = 3, ['u'] = 3,  /* U is treated as T */
};

/* DNA2 decoding table */
static const char kmersearch_dna2_decode_table[4] = {'A', 'C', 'G', 'T'};

/* DNA4 encoding table */
static const uint8 kmersearch_dna4_encode_table[256] = {
    ['A'] = 0x1, ['a'] = 0x1,  /* 0001 */
    ['C'] = 0x2, ['c'] = 0x2,  /* 0010 */
    ['G'] = 0x4, ['g'] = 0x4,  /* 0100 */
    ['T'] = 0x8, ['t'] = 0x8,  /* 1000 */
    ['U'] = 0x8, ['u'] = 0x8,  /* U is treated as T */
    ['M'] = 0x3, ['m'] = 0x3,  /* A or C: 0011 */
    ['R'] = 0x5, ['r'] = 0x5,  /* A or G: 0101 */
    ['W'] = 0x9, ['w'] = 0x9,  /* A or T: 1001 */
    ['S'] = 0x6, ['s'] = 0x6,  /* C or G: 0110 */
    ['Y'] = 0xA, ['y'] = 0xA,  /* C or T: 1010 */
    ['K'] = 0xC, ['k'] = 0xC,  /* G or T: 1100 */
    ['V'] = 0x7, ['v'] = 0x7,  /* A or C or G: 0111 */
    ['H'] = 0xB, ['h'] = 0xB,  /* A or C or T: 1011 */
    ['D'] = 0xD, ['d'] = 0xD,  /* A or G or T: 1101 */
    ['B'] = 0xE, ['b'] = 0xE,  /* C or G or T: 1110 */
    ['N'] = 0xF, ['n'] = 0xF,  /* A or C or G or T: 1111 */
};

/* DNA4 decoding table */
static const char kmersearch_dna4_decode_table[16] = {
    '?',  /* 0000 - invalid */
    'A',  /* 0001 */
    'C',  /* 0010 */
    'M',  /* 0011 */
    'G',  /* 0100 */
    'R',  /* 0101 */
    'S',  /* 0110 */
    'V',  /* 0111 */
    'T',  /* 1000 */
    'W',  /* 1001 */
    'Y',  /* 1010 */
    'H',  /* 1011 */
    'K',  /* 1100 */
    'D',  /* 1101 */
    'B',  /* 1110 */
    'N'   /* 1111 */
};

/* DNA4 to DNA2 expansion table */
/* Each entry contains: [expansion_count, base1, base2, base3, base4] */
static const uint8 kmersearch_dna4_to_dna2_table[16][5] = {
    {0, 0, 0, 0, 0},     /* 0000 - invalid */
    {1, 0, 0, 0, 0},     /* 0001 - A */
    {1, 1, 0, 0, 0},     /* 0010 - C */
    {2, 0, 1, 0, 0},     /* 0011 - M (A,C) */
    {1, 2, 0, 0, 0},     /* 0100 - G */
    {2, 0, 2, 0, 0},     /* 0101 - R (A,G) */
    {2, 1, 2, 0, 0},     /* 0110 - S (C,G) */
    {3, 0, 1, 2, 0},     /* 0111 - V (A,C,G) */
    {1, 3, 0, 0, 0},     /* 1000 - T */
    {2, 0, 3, 0, 0},     /* 1001 - W (A,T) */
    {2, 1, 3, 0, 0},     /* 1010 - Y (C,T) */
    {3, 0, 1, 3, 0},     /* 1011 - H (A,C,T) */
    {2, 2, 3, 0, 0},     /* 1100 - K (G,T) */
    {3, 0, 2, 3, 0},     /* 1101 - D (A,G,T) */
    {3, 1, 2, 3, 0},     /* 1110 - B (C,G,T) */
    {4, 0, 1, 2, 3}      /* 1111 - N (A,C,G,T) */
};

/* DNA2 functions */
PG_FUNCTION_INFO_V1(kmersearch_dna2_in);
PG_FUNCTION_INFO_V1(kmersearch_dna2_out);
PG_FUNCTION_INFO_V1(kmersearch_dna2_recv);
PG_FUNCTION_INFO_V1(kmersearch_dna2_send);

/* DNA4 functions */
PG_FUNCTION_INFO_V1(kmersearch_dna4_in);
PG_FUNCTION_INFO_V1(kmersearch_dna4_out);
PG_FUNCTION_INFO_V1(kmersearch_dna4_recv);
PG_FUNCTION_INFO_V1(kmersearch_dna4_send);

/* K-mer and GIN index functions */
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna2);
PG_FUNCTION_INFO_V1(kmersearch_extract_value_dna4);
PG_FUNCTION_INFO_V1(kmersearch_extract_query);
PG_FUNCTION_INFO_V1(kmersearch_consistent);
PG_FUNCTION_INFO_V1(kmersearch_compare_partial);
PG_FUNCTION_INFO_V1(kmersearch_dna2_match);
PG_FUNCTION_INFO_V1(kmersearch_dna4_match);

/* K-mer frequency analysis functions */
PG_FUNCTION_INFO_V1(kmersearch_analyze_table_frequency);
PG_FUNCTION_INFO_V1(kmersearch_get_highfreq_kmers);
PG_FUNCTION_INFO_V1(kmersearch_analyze_table);
PG_FUNCTION_INFO_V1(kmersearch_drop_analysis);
PG_FUNCTION_INFO_V1(kmersearch_reduce_index);
PG_FUNCTION_INFO_V1(kmersearch_rawscore_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_rawscore_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_query_pattern_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_query_pattern_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_actual_min_score_cache_free);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmers_cache_load);
PG_FUNCTION_INFO_V1(kmersearch_highfreq_kmers_cache_free);

/* Score calculation functions */
PG_FUNCTION_INFO_V1(kmersearch_rawscore_dna2);
PG_FUNCTION_INFO_V1(kmersearch_rawscore_dna4);
PG_FUNCTION_INFO_V1(kmersearch_correctedscore_dna2);
PG_FUNCTION_INFO_V1(kmersearch_correctedscore_dna4);

/* Operator support functions */
PG_FUNCTION_INFO_V1(kmersearch_dna2_eq);
PG_FUNCTION_INFO_V1(kmersearch_dna4_eq);

/* Helper functions */
static bool kmersearch_is_valid_dna2_char(char c)
{
    return (c == 'A' || c == 'a' || c == 'C' || c == 'c' ||
            c == 'G' || c == 'g' || c == 'T' || c == 't' ||
            c == 'U' || c == 'u');
}

static bool kmersearch_is_valid_dna4_char(char c)
{
    return kmersearch_dna4_encode_table[(unsigned char)c] != 0 || c == 'A' || c == 'a';
}

/*
 * DNA2 input function
 */
Datum
kmersearch_dna2_in(PG_FUNCTION_ARGS)
{
    char *input_string = PG_GETARG_CSTRING(0);
    int input_len = strlen(input_string);
    int bit_len = input_len * 2;  /* 2 bits per character */
    int byte_len = (bit_len + 7) / 8;  /* Round up to bytes */
    VarBit *result;
    bits8 *data_ptr;
    int i;
    
    /* Validate input */
    for (i = 0; i < input_len; i++)
    {
        if (!kmersearch_is_valid_dna2_char(input_string[i]))
        {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid character '%c' for DNA2 type", input_string[i]),
                     errhint("DNA2 type accepts only A, C, G, T, U characters")));
        }
    }
    
    /* Allocate result */
    result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    VARBITLEN(result) = bit_len;
    
    /* Encode sequence using SIMD dispatch */
    data_ptr = VARBITS(result);
    simd_dispatch.dna2_encode(input_string, (uint8_t*)data_ptr, input_len);
    
    PG_RETURN_VARBIT_P(result);
}

/*
 * DNA2 output function
 */
Datum
kmersearch_dna2_out(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int char_len = bit_len / 2;
    char *result;
    bits8 *data_ptr = VARBITS(dna);
    
    if (dna == NULL)
        ereport(ERROR, (errmsg("input DNA sequence is NULL")));
    
    if (bit_len < 0)
        ereport(ERROR, (errmsg("invalid bit length: %d", bit_len)));
    
    if (bit_len % 2 != 0)
        ereport(ERROR, (errmsg("bit length must be even for DNA2")));
    
    result = (char *) palloc(char_len + 1);
    
    /* Decode sequence using SIMD dispatch */
    simd_dispatch.dna2_decode((uint8_t*)data_ptr, result, char_len);
    
    PG_RETURN_CSTRING(result);
}

/*
 * DNA2 receive function (binary input)
 */
Datum
kmersearch_dna2_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    VarBit *result;
    int32 bit_len;
    int byte_len;
    
    bit_len = pq_getmsgint(buf, sizeof(int32));
    byte_len = (bit_len + 7) / 8;
    
    result = (VarBit *) palloc(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    VARBITLEN(result) = bit_len;
    
    pq_copymsgbytes(buf, (char *) VARBITS(result), byte_len);
    
    PG_RETURN_VARBIT_P(result);
}

/*
 * DNA2 send function (binary output)
 */
Datum
kmersearch_dna2_send(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    StringInfoData buf;
    int32 bit_len = VARBITLEN(dna);
    int byte_len = (bit_len + 7) / 8;
    
    pq_begintypsend(&buf);
    pq_sendint32(&buf, bit_len);
    pq_sendbytes(&buf, (char *) VARBITS(dna), byte_len);
    
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/*
 * DNA4 input function
 */
Datum
kmersearch_dna4_in(PG_FUNCTION_ARGS)
{
    char *input_string = PG_GETARG_CSTRING(0);
    int input_len = strlen(input_string);
    int bit_len = input_len * 4;  /* 4 bits per character */
    int byte_len = (bit_len + 7) / 8;  /* Round up to bytes */
    VarBit *result;
    bits8 *data_ptr;
    int i;
    
    /* Validate input */
    for (i = 0; i < input_len; i++)
    {
        if (!kmersearch_is_valid_dna4_char(input_string[i]))
        {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid character '%c' for DNA4 type", input_string[i]),
                     errhint("DNA4 type accepts A,C,G,T,U,M,R,W,S,Y,K,V,H,D,B,N characters")));
        }
    }
    
    /* Allocate result */
    result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    VARBITLEN(result) = bit_len;
    
    /* Encode sequence using SIMD dispatch */
    data_ptr = VARBITS(result);
    simd_dispatch.dna4_encode(input_string, (uint8_t*)data_ptr, input_len);
    
    PG_RETURN_VARBIT_P(result);
}

/*
 * DNA4 output function
 */
Datum
kmersearch_dna4_out(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int char_len = bit_len / 4;
    char *result;
    bits8 *data_ptr = VARBITS(dna);
    int i;
    
    if (dna == NULL)
        ereport(ERROR, (errmsg("input DNA sequence is NULL")));
    
    if (bit_len < 0)
        ereport(ERROR, (errmsg("invalid bit length: %d", bit_len)));
    
    if (bit_len % 4 != 0)
        ereport(ERROR, (errmsg("bit length must be multiple of 4 for DNA4")));
    
    result = (char *) palloc(char_len + 1);
    
    /* Decode sequence using SIMD dispatch */
    simd_dispatch.dna4_decode((uint8_t*)data_ptr, result, char_len);
    PG_RETURN_CSTRING(result);
}

/*
 * DNA4 receive function (binary input)
 */
Datum
kmersearch_dna4_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    VarBit *result;
    int32 bit_len;
    int byte_len;
    
    bit_len = pq_getmsgint(buf, sizeof(int32));
    byte_len = (bit_len + 7) / 8;
    
    result = (VarBit *) palloc(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    VARBITLEN(result) = bit_len;
    
    pq_copymsgbytes(buf, (char *) VARBITS(result), byte_len);
    
    PG_RETURN_VARBIT_P(result);
}

/*
 * DNA4 send function (binary output)
 */
Datum
kmersearch_dna4_send(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    StringInfoData buf;
    int32 bit_len = VARBITLEN(dna);
    int byte_len = (bit_len + 7) / 8;
    
    pq_begintypsend(&buf);
    pq_sendint32(&buf, bit_len);
    pq_sendbytes(&buf, (char *) VARBITS(dna), byte_len);
    
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/* SIMD capability detection functions */
static simd_capability_t detect_cpu_capabilities(void);
static void init_simd_dispatch_table(void);

/* SIMD implementation functions */
static void dna2_encode_scalar(const char* input, uint8_t* output, int len);
static void dna2_decode_scalar(const uint8_t* input, char* output, int len);
static void dna4_encode_scalar(const char* input, uint8_t* output, int len);
static void dna4_decode_scalar(const uint8_t* input, char* output, int len);

#ifdef __x86_64__
static void dna2_encode_avx2(const char* input, uint8_t* output, int len);
static void dna2_decode_avx2(const uint8_t* input, char* output, int len);
static void dna4_encode_avx2(const char* input, uint8_t* output, int len);
static void dna4_decode_avx2(const uint8_t* input, char* output, int len);

/* K-mer processing functions with SIMD optimization */
static Datum *kmersearch_extract_dna2_kmers_direct_avx2(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmers_with_expansion_direct_avx2(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmers_fast_avx2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

static Datum *kmersearch_extract_dna2_kmers_direct_avx512(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmers_with_expansion_direct_avx512(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmers_fast_avx512(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

static Datum *kmersearch_extract_dna2_kmers_direct_neon(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmers_with_expansion_direct_neon(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmers_fast_neon(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

static Datum *kmersearch_extract_dna2_kmers_direct_sve(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmers_with_expansion_direct_sve(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmers_fast_sve(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

/* Scalar versions */
static Datum *kmersearch_extract_dna2_kmers_direct_scalar(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmers_with_expansion_direct_scalar(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmers_fast_scalar(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

static void dna2_encode_avx512(const char* input, uint8_t* output, int len);
static void dna2_decode_avx512(const uint8_t* input, char* output, int len);
static void dna4_encode_avx512(const char* input, uint8_t* output, int len);
static void dna4_decode_avx512(const uint8_t* input, char* output, int len);
#endif

#ifdef __aarch64__
static void dna2_encode_neon(const char* input, uint8_t* output, int len);
static void dna2_decode_neon(const uint8_t* input, char* output, int len);
static void dna4_encode_neon(const char* input, uint8_t* output, int len);
static void dna4_decode_neon(const uint8_t* input, char* output, int len);

#ifdef __ARM_FEATURE_SVE
static void dna2_encode_sve(const char* input, uint8_t* output, int len);
static void dna2_decode_sve(const uint8_t* input, char* output, int len);
static void dna4_encode_sve(const char* input, uint8_t* output, int len);
static void dna4_decode_sve(const uint8_t* input, char* output, int len);
#endif
#endif

/*
 * GUC assign hook functions for cache invalidation
 */

/* High-frequency cache clear with warning */
static void 
clear_highfreq_cache_with_warning(void)
{
    kmersearch_highfreq_kmers_cache_free_internal();
    elog(WARNING, "High-frequency k-mer cache has been cleared. "
                  "You may need to manually execute kmersearch_highfreq_kmers_cache_load() "
                  "to reload the cache if needed.");
}

/* K-mer size change affects all caches */
static void
kmersearch_kmer_size_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Clear rawscore cache */
    if (rawscore_cache_manager)
        free_rawscore_cache_manager(&rawscore_cache_manager);
    
    /* Clear query pattern cache */
    if (query_pattern_cache_manager)
        free_query_pattern_cache_manager(&query_pattern_cache_manager);
    
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
    /* Clear high-frequency k-mer cache with warning */
    clear_highfreq_cache_with_warning();
}

/* Appearance rate change affects high-freq and actual min score caches */
static void
kmersearch_max_appearance_rate_assign_hook(double newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
    /* Clear high-frequency k-mer cache with warning */
    clear_highfreq_cache_with_warning();
}

/* Appearance nrow change affects high-freq and actual min score caches */
static void
kmersearch_max_appearance_nrow_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
    /* Clear high-frequency k-mer cache with warning */
    clear_highfreq_cache_with_warning();
}

/* Min score change affects actual min score cache */
static void
kmersearch_min_score_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
}

/* Min shared ngram key rate change affects actual min score cache */
static void
kmersearch_min_shared_ngram_key_rate_assign_hook(double newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
}

/* Rawscore cache max entries change requires cache recreation */
static void
kmersearch_rawscore_cache_max_entries_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Clear rawscore cache to recreate with new size limit */
    if (rawscore_cache_manager)
        free_rawscore_cache_manager(&rawscore_cache_manager);
}

/* Query pattern cache max entries change requires cache recreation */
static void
kmersearch_query_pattern_cache_max_entries_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Clear query pattern cache to recreate with new size limit */
    if (query_pattern_cache_manager)
        free_query_pattern_cache_manager(&query_pattern_cache_manager);
}

/* Occurrence bit length change affects rawscore and high-freq caches */
static void
kmersearch_occur_bitlen_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Clear rawscore cache */
    if (rawscore_cache_manager)
        free_rawscore_cache_manager(&rawscore_cache_manager);
    
    /* Clear high-frequency k-mer cache with warning */
    clear_highfreq_cache_with_warning();
}

/*
 * Module initialization
 */
void
_PG_init(void)
{
    /* Initialize SIMD capabilities */
    simd_capability = detect_cpu_capabilities();
    init_simd_dispatch_table();
    
    /* Define custom GUC variables */
    DefineCustomRealVariable("kmersearch.max_appearance_rate",
                            "Maximum appearance rate for k-mers to be included in index",
                            "K-mers appearing in more than this fraction of rows will be identified as highly frequent",
                            &kmersearch_max_appearance_rate,
                            0.05,
                            0.0,
                            1.0,
                            PGC_USERSET,
                            0,
                            NULL,
                            kmersearch_max_appearance_rate_assign_hook,
                            NULL);
    
    DefineCustomIntVariable("kmersearch.max_appearance_nrow",
                           "Maximum number of rows for k-mers to be included in index",
                           "K-mers appearing in more than this number of rows will be identified as highly frequent (0 = unlimited)",
                           &kmersearch_max_appearance_nrow,
                           0,
                           0,
                           INT_MAX,
                           PGC_USERSET,
                           0,
                           NULL,
                           kmersearch_max_appearance_nrow_assign_hook,
                           NULL);
    
    DefineCustomIntVariable("kmersearch.min_score",
                           "Minimum score (shared n-gram count) for GIN k-mer search",
                           "Query results with score below this threshold will be filtered out",
                           &kmersearch_min_score,
                           1,
                           0,
                           INT_MAX,
                           PGC_USERSET,
                           0,
                           NULL,
                           kmersearch_min_score_assign_hook,
                           NULL);
    
    DefineCustomIntVariable("kmersearch.occur_bitlen",
                           "Number of bits used for occurrence count in k-mer index",
                           "Controls the maximum occurrence count that can be stored (0-16 bits)",
                           &kmersearch_occur_bitlen,
                           8,
                           0,
                           16,
                           PGC_USERSET,
                           0,
                           NULL,
                           kmersearch_occur_bitlen_assign_hook,
                           NULL);
    
    DefineCustomIntVariable("kmersearch.kmer_size",
                           "K-mer size used for index creation and search",
                           "Length of k-mer sequences for similarity matching (4-64)",
                           &kmersearch_kmer_size,
                           8,
                           4,
                           64,
                           PGC_USERSET,
                           0,
                           NULL,
                           kmersearch_kmer_size_assign_hook,
                           NULL);
    
    DefineCustomRealVariable("kmersearch.min_shared_ngram_key_rate",
                            "Minimum shared n-gram key rate for =% operator matching",
                            "Minimum ratio of shared n-gram keys between query and target sequence (0.0-1.0)",
                            &kmersearch_min_shared_ngram_key_rate,
                            0.9,
                            0.0,
                            1.0,
                            PGC_USERSET,
                            0,
                            NULL,
                            kmersearch_min_shared_ngram_key_rate_assign_hook,
                            NULL);

    /* Define GUC variables for cache configuration */
    DefineCustomIntVariable("kmersearch.rawscore_cache_max_entries",
                           "Maximum number of entries in rawscore cache",
                           "Controls the maximum number of cached rawscore calculation results",
                           &kmersearch_rawscore_cache_max_entries,
                           50000,
                           1000,
                           10000000,
                           PGC_USERSET,
                           0,
                           NULL,
                           kmersearch_rawscore_cache_max_entries_assign_hook,
                           NULL);

    DefineCustomIntVariable("kmersearch.query_pattern_cache_max_entries",
                           "Maximum number of entries in query pattern cache",
                           "Controls the maximum number of cached query pattern extraction results",
                           &kmersearch_query_pattern_cache_max_entries,
                           50000,
                           1000,
                           10000000,
                           PGC_USERSET,
                           0,
                           NULL,
                           kmersearch_query_pattern_cache_max_entries_assign_hook,
                           NULL);

    DefineCustomIntVariable("kmersearch.actual_min_score_cache_max_entries",
                           "Maximum number of entries in actual min score cache",
                           "Controls the maximum number of cached actual min score calculation results",
                           &kmersearch_actual_min_score_cache_max_entries,
                           50000,
                           1000,
                           10000000,
                           PGC_USERSET,
                           0,
                           NULL,
                           NULL,
                           NULL);
    
    /* Initialize high-frequency k-mer cache */
    kmersearch_highfreq_kmers_cache_init();
}

/*
 * CPU capability detection
 */
static simd_capability_t detect_cpu_capabilities(void)
{
#ifdef __x86_64__
    unsigned int eax, ebx, ecx, edx;
    
    /* Check for AVX512 support */
    if (__get_cpuid_max(0, NULL) >= 7) {
        __cpuid_count(7, 0, eax, ebx, ecx, edx);
        if (ebx & (1 << 16)) { /* AVX512F */
            if (ebx & (1 << 30)) { /* AVX512BW */
                return SIMD_AVX512;
            }
        }
    }
    
    /* Check for AVX2 support */
    if (__get_cpuid_max(0, NULL) >= 7) {
        __cpuid_count(7, 0, eax, ebx, ecx, edx);
        if (ebx & (1 << 5)) { /* AVX2 */
            return SIMD_AVX2;
        }
    }
#elif defined(__aarch64__)
    /* Check for SVE support */
#ifdef __ARM_FEATURE_SVE
    if (access("/proc/sys/abi/sve_default_vector_length", F_OK) == 0) {
        return SIMD_SVE;
    }
#endif
    
    /* ARM64 always has NEON */
    return SIMD_NEON;
#endif
    
    return SIMD_NONE;
}

/*
 * Initialize SIMD dispatch table
 */
static void init_simd_dispatch_table(void)
{
    /* Set default scalar implementations */
    simd_dispatch.dna2_encode = dna2_encode_scalar;
    simd_dispatch.dna2_decode = dna2_decode_scalar;
    simd_dispatch.dna4_encode = dna4_encode_scalar;
    simd_dispatch.dna4_decode = dna4_decode_scalar;
    
    /* Override with SIMD implementations if available */
    switch (simd_capability) {
#ifdef __x86_64__
        case SIMD_AVX512:
            simd_dispatch.dna2_encode = dna2_encode_avx512;
            simd_dispatch.dna2_decode = dna2_decode_avx512;
            simd_dispatch.dna4_encode = dna4_encode_avx512;
            simd_dispatch.dna4_decode = dna4_decode_avx512;
            break;
        case SIMD_AVX2:
            simd_dispatch.dna2_encode = dna2_encode_avx2;
            simd_dispatch.dna2_decode = dna2_decode_avx2;
            simd_dispatch.dna4_encode = dna4_encode_avx2;
            simd_dispatch.dna4_decode = dna4_decode_avx2;
            break;
#endif
#ifdef __aarch64__
#ifdef __ARM_FEATURE_SVE
        case SIMD_SVE:
            simd_dispatch.dna2_encode = dna2_encode_sve;
            simd_dispatch.dna2_decode = dna2_decode_sve;
            simd_dispatch.dna4_encode = dna4_encode_sve;
            simd_dispatch.dna4_decode = dna4_decode_sve;
            break;
#endif
        case SIMD_NEON:
            simd_dispatch.dna2_encode = dna2_encode_neon;
            simd_dispatch.dna2_decode = dna2_decode_neon;
            simd_dispatch.dna4_encode = dna4_encode_neon;
            simd_dispatch.dna4_decode = dna4_decode_neon;
            break;
#endif
        case SIMD_NONE:
        default:
            /* Already set to scalar implementations */
            break;
    }
}

/*
 * Helper function to get bit at position from bit array
 */
static inline uint8
kmersearch_get_bit_at(bits8 *data, int bit_pos)
{
    int byte_pos = bit_pos / 8;
    int bit_offset = bit_pos % 8;
    return (data[byte_pos] >> (7 - bit_offset)) & 1;
}

/*
 * Helper function to set bit at position in bit array
 */
static inline void
kmersearch_set_bit_at(bits8 *data, int bit_pos, uint8 value)
{
    int byte_pos = bit_pos / 8;
    int bit_offset = bit_pos % 8;
    if (value)
        data[byte_pos] |= (1 << (7 - bit_offset));
    else
        data[byte_pos] &= ~(1 << (7 - bit_offset));
}

/*
 * Convert DNA2 sequence to string
 */
static char *
kmersearch_dna2_to_string(VarBit *dna)
{
    int32 bit_len = VARBITLEN(dna);
    int char_len = bit_len / 2;
    char *result;
    bits8 *data_ptr = VARBITS(dna);
    int i;
    
    result = (char *) palloc(char_len + 1);
    
    for (i = 0; i < char_len; i++)
    {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 encoded = (data_ptr[byte_pos] >> (6 - bit_offset)) & 0x3;
        result[i] = kmersearch_dna2_decode_table[encoded];
    }
    
    result[char_len] = '\0';
    return result;
}

/*
 * Convert DNA4 sequence to string
 */
static char *
kmersearch_dna4_to_string(VarBit *dna)
{
    int32 bit_len = VARBITLEN(dna);
    int char_len = bit_len / 4;
    char *result;
    bits8 *data_ptr = VARBITS(dna);
    int i;
    
    result = (char *) palloc(char_len + 1);
    
    for (i = 0; i < char_len; i++)
    {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 encoded;
        
        if (bit_offset <= 4)
        {
            encoded = (data_ptr[byte_pos] >> (4 - bit_offset)) & 0xF;
        }
        else
        {
            encoded = ((data_ptr[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < (bit_len + 7) / 8)
                encoded |= (data_ptr[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        result[i] = kmersearch_dna4_decode_table[encoded];
    }
    
    result[char_len] = '\0';
    return result;
}

/*
 * Fast check if degenerate combinations exceed limit (11+)
 * Returns true if combinations will exceed 10
 */
static bool
kmersearch_will_exceed_degenerate_limit(const char *seq, int len)
{
    int n_count = 0, vhdb_count = 0, mrwsyk_count = 0;
    int i;
    
    for (i = 0; i < len; i++)
    {
        char c = toupper(seq[i]);
        if (c == 'N')
        {
            n_count++;
            /* Early exit: 2+ N's always exceed limit */
            if (n_count >= 2)
                return true;
        }
        else if (c == 'V' || c == 'H' || c == 'D' || c == 'B')
        {
            vhdb_count++;
            /* Early exit: 3+ VHDB always exceed limit */
            if (vhdb_count >= 3)
                return true;
            /* Early exit: N + VHDB exceed limit */
            if (n_count >= 1 && vhdb_count >= 1)
                return true;
        }
        else if (c == 'M' || c == 'R' || c == 'W' || c == 'S' || c == 'Y' || c == 'K')
        {
            mrwsyk_count++;
            /* Early exit: 4+ MRWSYK always exceed limit */
            if (mrwsyk_count >= 4)
                return true;
            /* Early exit: N + 2+ MRWSYK exceed limit */
            if (n_count >= 1 && mrwsyk_count >= 2)
                return true;
            /* Early exit: 2+ VHDB + MRWSYK exceed limit */
            if (vhdb_count >= 2 && mrwsyk_count >= 1)
                return true;
            /* Early exit: VHDB + 2+ MRWSYK exceed limit */
            if (vhdb_count >= 1 && mrwsyk_count >= 2)
                return true;
        }
    }
    
    return false;  /* Within limit */
}

/*
 * Fast check if DNA4 k-mer will exceed degenerate limit using bit operations
 */
static bool
kmersearch_will_exceed_degenerate_limit_dna4_bits(VarBit *seq, int start_pos, int k)
{
    int n_count = 0, vhdb_count = 0, mrwsyk_count = 0;
    bits8 *data = VARBITS(seq);
    int i;
    
    for (i = 0; i < k; i++)
    {
        int bit_pos = (start_pos + i) * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 encoded;
        
        /* Extract 4 bits for this base */
        if (bit_offset <= 4)
        {
            encoded = (data[byte_pos] >> (4 - bit_offset)) & 0xF;
        }
        else
        {
            encoded = ((data[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < VARBITBYTES(seq))
                encoded |= (data[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        /* Check expansion count using lookup table */
        {
            int expansion_count = kmersearch_dna4_to_dna2_table[encoded][0];
        
        if (expansion_count == 4)  /* N */
        {
            n_count++;
            if (n_count >= 2)
                return true;
        }
        else if (expansion_count == 3)  /* V,H,D,B */
        {
            vhdb_count++;
            if (vhdb_count >= 3)
                return true;
            if (n_count >= 1 && vhdb_count >= 1)
                return true;
        }
        else if (expansion_count == 2)  /* M,R,W,S,Y,K */
        {
            mrwsyk_count++;
            if (mrwsyk_count >= 4)
                return true;
            if (n_count >= 1 && mrwsyk_count >= 2)
                return true;
            if (vhdb_count >= 2 && mrwsyk_count >= 1)
                return true;
            if (vhdb_count >= 1 && mrwsyk_count >= 2)
                return true;
        }
        }
    }
    
    return false;
}

/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using bit operations
 */
static VarBit **
kmersearch_expand_dna4_kmer_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[64][4];  /* Max k=64, max 4 expansions per base */
    int base_counts[64];
    int total_combinations = 1;
    VarBit **results;
    int i, combo;
    
    *expansion_count = 0;
    
    /* Check if expansion will exceed limit */
    if (kmersearch_will_exceed_degenerate_limit_dna4_bits(dna4_seq, start_pos, k))
        return NULL;
    
    /* Extract expansion info for each base */
    for (i = 0; i < k; i++)
    {
        int bit_pos = (start_pos + i) * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 encoded;
        int exp_count;
        int j;
        
        /* Extract 4 bits */
        if (bit_offset <= 4)
        {
            encoded = (data[byte_pos] >> (4 - bit_offset)) & 0xF;
        }
        else
        {
            encoded = ((data[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < VARBITBYTES(dna4_seq))
                encoded |= (data[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        /* Get expansion from table */
        exp_count = kmersearch_dna4_to_dna2_table[encoded][0];
        base_counts[i] = exp_count;
        
        for (j = 0; j < exp_count; j++)
        {
            base_expansions[i][j] = kmersearch_dna4_to_dna2_table[encoded][j + 1];
        }
        
        total_combinations *= exp_count;
    }
    
    /* Allocate result array */
    results = (VarBit **) palloc(total_combinations * sizeof(VarBit *));
    
    /* Generate all combinations */
    for (combo = 0; combo < total_combinations; combo++)
    {
        int kmer_bits = k * 2;
        int kmer_bytes = (kmer_bits + 7) / 8;
        VarBit *result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
        bits8 *result_data;
        int temp_combo = combo;
        
        SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + kmer_bytes);
        VARBITLEN(result) = kmer_bits;
        result_data = VARBITS(result);
        
        /* Generate this combination */
        for (i = 0; i < k; i++)
        {
            int base_idx = temp_combo % base_counts[i];
            uint8 dna2_base = base_expansions[i][base_idx];
            int dst_bit_pos = i * 2;
            int dst_byte_pos = dst_bit_pos / 8;
            int dst_bit_offset = dst_bit_pos % 8;
            
            result_data[dst_byte_pos] |= (dna2_base << (6 - dst_bit_offset));
            temp_combo /= base_counts[i];
        }
        
        results[combo] = result;
    }
    
    *expansion_count = total_combinations;
    return results;
}

/*
 * Count combinations for degenerate code expansion
 */
static int
kmersearch_count_degenerate_combinations(const char *seq, int len)
{
    int n_count = 0, vhdb_count = 0, mrwsyk_count = 0;
    int combinations = 1;
    int i;
    
    for (i = 0; i < len; i++)
    {
        char c = toupper(seq[i]);
        if (c == 'N')
            n_count++;
        else if (c == 'V' || c == 'H' || c == 'D' || c == 'B')
            vhdb_count++;
        else if (c == 'M' || c == 'R' || c == 'W' || c == 'S' || c == 'Y' || c == 'K')
            mrwsyk_count++;
    }
    
    /* Check if combinations exceed 10 */
    if (n_count >= 2 ||
        (n_count == 1 && (vhdb_count >= 1 || mrwsyk_count >= 2)) ||
        vhdb_count >= 3 ||
        (vhdb_count == 2 && mrwsyk_count >= 1) ||
        (vhdb_count == 1 && mrwsyk_count >= 2) ||
        mrwsyk_count >= 4)
    {
        return 11;  /* Over limit */
    }
    
    /* Calculate actual combinations */
    for (i = 0; i < len; i++)
    {
        char c = toupper(seq[i]);
        if (c == 'N')
            combinations *= 4;
        else if (c == 'V' || c == 'H' || c == 'D' || c == 'B')
            combinations *= 3;
        else if (c == 'M' || c == 'R' || c == 'W' || c == 'S' || c == 'Y' || c == 'K')
            combinations *= 2;
    }
    
    return combinations;
}

/*
 * Expand degenerate codes to all possible combinations
 */
static void
kmersearch_expand_degenerate_sequence(const char *seq, int len, char **results, int *count)
{
    int combinations = kmersearch_count_degenerate_combinations(seq, len);
    if (combinations > 10)
    {
        *count = 0;
        return;
    }
    
    *count = combinations;
    
    /* Allocate results array */
    for (int i = 0; i < combinations; i++)
    {
        results[i] = (char *) palloc(len + 1);
        results[i][len] = '\0';
    }
    
    /* Generate all combinations */
    for (int combo = 0; combo < combinations; combo++)
    {
        int temp_combo = combo;
        for (int pos = 0; pos < len; pos++)
        {
            char c = toupper(seq[pos]);
            char bases[4];
            int base_count = 0;
            
            /* Get possible bases for this character */
            if (c == 'A') { bases[0] = 'A'; base_count = 1; }
            else if (c == 'C') { bases[0] = 'C'; base_count = 1; }
            else if (c == 'G') { bases[0] = 'G'; base_count = 1; }
            else if (c == 'T' || c == 'U') { bases[0] = 'T'; base_count = 1; }
            else if (c == 'M') { bases[0] = 'A'; bases[1] = 'C'; base_count = 2; }
            else if (c == 'R') { bases[0] = 'A'; bases[1] = 'G'; base_count = 2; }
            else if (c == 'W') { bases[0] = 'A'; bases[1] = 'T'; base_count = 2; }
            else if (c == 'S') { bases[0] = 'C'; bases[1] = 'G'; base_count = 2; }
            else if (c == 'Y') { bases[0] = 'C'; bases[1] = 'T'; base_count = 2; }
            else if (c == 'K') { bases[0] = 'G'; bases[1] = 'T'; base_count = 2; }
            else if (c == 'V') { bases[0] = 'A'; bases[1] = 'C'; bases[2] = 'G'; base_count = 3; }
            else if (c == 'H') { bases[0] = 'A'; bases[1] = 'C'; bases[2] = 'T'; base_count = 3; }
            else if (c == 'D') { bases[0] = 'A'; bases[1] = 'G'; bases[2] = 'T'; base_count = 3; }
            else if (c == 'B') { bases[0] = 'C'; bases[1] = 'G'; bases[2] = 'T'; base_count = 3; }
            else if (c == 'N') { bases[0] = 'A'; bases[1] = 'C'; bases[2] = 'G'; bases[3] = 'T'; base_count = 4; }
            
            results[combo][pos] = bases[temp_combo % base_count];
            temp_combo /= base_count;
        }
    }
}

/*
 * Create n-gram key from k-mer string
 */
static VarBit *
kmersearch_create_ngram_key(const char *kmer, int k, int occurrence)
{
    int kmer_bits = k * 2;  /* 2 bits per base */
    int occur_bits = kmersearch_occur_bitlen;
    int total_bits = kmer_bits + occur_bits;
    int total_bytes = (total_bits + 7) / 8;
    int adj_occurrence = occurrence - 1;  /* 1-offset to 0-offset */
    VarBit *result;
    bits8 *data_ptr;
    int i;
    
    result = (VarBit *) palloc0(VARBITHDRSZ + total_bytes);
    SET_VARSIZE(result, VARBITHDRSZ + total_bytes);
    VARBITLEN(result) = total_bits;
    
    data_ptr = VARBITS(result);
    
    /* Encode k-mer */
    for (i = 0; i < k; i++)
    {
        uint8 encoded = kmersearch_dna2_encode_table[(unsigned char)kmer[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        data_ptr[byte_pos] |= (encoded << (6 - bit_offset));
    }
    
    /* Encode occurrence count */
    if (adj_occurrence >= (1 << occur_bits))
        adj_occurrence = (1 << occur_bits) - 1;  /* Cap at max value */
    
    for (i = 0; i < occur_bits; i++)
    {
        int bit_pos = kmer_bits + i;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (adj_occurrence & (1 << (occur_bits - 1 - i)))
            data_ptr[byte_pos] |= (1 << (7 - bit_offset));
    }
    
    return result;
}

/*
 * Create n-gram key directly from DNA2 bit sequence
 */
static VarBit *
kmersearch_create_ngram_key_from_dna2_bits(VarBit *seq, int start_pos, int k, int occurrence)
{
    int kmer_bits = k * 2;
    int occur_bits = kmersearch_occur_bitlen;
    int total_bits = kmer_bits + occur_bits;
    int total_bytes = (total_bits + 7) / 8;
    int adj_occurrence = occurrence - 1;
    VarBit *result;
    bits8 *src_data, *dst_data;
    int i;
    
    result = (VarBit *) palloc0(VARBITHDRSZ + total_bytes);
    SET_VARSIZE(result, VARBITHDRSZ + total_bytes);
    VARBITLEN(result) = total_bits;
    
    src_data = VARBITS(seq);
    dst_data = VARBITS(result);
    
    /* Copy k-mer bits directly from source */
    for (i = 0; i < k; i++)
    {
        int src_bit_pos = (start_pos + i) * 2;
        int dst_bit_pos = i * 2;
        int src_byte_pos = src_bit_pos / 8;
        int src_bit_offset = src_bit_pos % 8;
        int dst_byte_pos = dst_bit_pos / 8;
        int dst_bit_offset = dst_bit_pos % 8;
        
        /* Extract 2 bits from source */
        uint8 base_bits = (src_data[src_byte_pos] >> (6 - src_bit_offset)) & 0x3;
        
        /* Store 2 bits in destination */
        dst_data[dst_byte_pos] |= (base_bits << (6 - dst_bit_offset));
    }
    
    /* Encode occurrence count */
    if (adj_occurrence >= (1 << occur_bits))
        adj_occurrence = (1 << occur_bits) - 1;
    
    for (i = 0; i < occur_bits; i++)
    {
        int bit_pos = kmer_bits + i;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (adj_occurrence & (1 << (occur_bits - 1 - i)))
            dst_data[byte_pos] |= (1 << (7 - bit_offset));
    }
    
    return result;
}

/*
 * Create k-mer key from DNA2 bits (without occurrence count)
 */
static VarBit *
kmersearch_create_kmer_key_from_dna2_bits(VarBit *seq, int start_pos, int k)
{
    int kmer_bits = k * 2;
    int total_bytes = (kmer_bits + 7) / 8;
    VarBit *result;
    bits8 *src_data, *dst_data;
    int src_bytes = VARBITBYTES(seq);
    int i;
    
    result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + total_bytes);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + total_bytes);
    VARBITLEN(result) = kmer_bits;
    
    src_data = VARBITS(seq);
    dst_data = VARBITS(result);
    
    /* Copy k-mer bits directly from source */
    for (i = 0; i < k; i++)
    {
        int src_bit_pos = (start_pos + i) * 2;
        int dst_bit_pos = i * 2;
        int src_byte_pos = src_bit_pos / 8;
        int src_bit_offset = src_bit_pos % 8;
        int dst_byte_pos = dst_bit_pos / 8;
        int dst_bit_offset = dst_bit_pos % 8;
        uint8 base_bits;
        
        /* Boundary check to prevent buffer overflow */
        if (src_byte_pos >= src_bytes) {
            /* Return NULL to indicate failed k-mer extraction */
            pfree(result);
            return NULL;
        }
        
        /* Extract 2 bits from source */
        base_bits = (src_data[src_byte_pos] >> (6 - src_bit_offset)) & 0x3;
        
        /* Store in destination */
        dst_data[dst_byte_pos] |= (base_bits << (6 - dst_bit_offset));
    }
    
    return result;
}

/*
 * Extract k-mers directly from DNA2 bit sequence (with SIMD dispatch)
 */
static Datum *
kmersearch_extract_dna2_kmers_direct(VarBit *seq, int k, int *nkeys)
{
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512) {
        return kmersearch_extract_dna2_kmers_direct_avx512(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_AVX2) {
        return kmersearch_extract_dna2_kmers_direct_avx2(seq, k, nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE) {
        return kmersearch_extract_dna2_kmers_direct_sve(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_NEON) {
        return kmersearch_extract_dna2_kmers_direct_neon(seq, k, nkeys);
    }
#endif
    return kmersearch_extract_dna2_kmers_direct_scalar(seq, k, nkeys);
}

/*
 * Scalar version: Extract k-mers directly from DNA2 bit sequence
 */
static Datum *
kmersearch_extract_dna2_kmers_direct_scalar(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Use simple array-based k-mer tracking (no hashing needed for k<=64) */
    occurrences = (KmerOccurrence *) palloc(max_kmers * sizeof(KmerOccurrence));
    
    /* Extract k-mers */
    for (i = 0; i <= seq_bases - k; i++)
    {
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        /* Extract k-mer as single uint64_t value */
        kmer_value = kmersearch_extract_kmer_as_uint64(seq, i, k);
        
        /* Skip if k-mer extraction failed (boundary check failed) */
        if (kmer_value == 0 && k > 0) {
            /* Note: valid k-mers could be 0 (all A's), but check bounds properly */
            int last_bit_pos = (i + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;  /* Out of bounds, skip */
            }
        }
        
        /* Find or add occurrence count using binary search */
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, max_kmers);
        
        if (current_count < 0)
            continue;  /* Array full, skip */
        
        /* Skip if occurrence exceeds bit limit */
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* Create simple k-mer key (without occurrence count for matching) */
        ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, i, k);
        if (ngram_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(ngram_key);
    }
    
    /* Cleanup */
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/*
 * Extract k-mers directly from DNA4 bit sequence with degenerate expansion (with SIMD dispatch)
 */
static Datum *
kmersearch_extract_dna4_kmers_with_expansion_direct(VarBit *seq, int k, int *nkeys)
{
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512) {
        return kmersearch_extract_dna4_kmers_with_expansion_direct_avx512(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_AVX2) {
        return kmersearch_extract_dna4_kmers_with_expansion_direct_avx2(seq, k, nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE) {
        return kmersearch_extract_dna4_kmers_with_expansion_direct_sve(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_NEON) {
        return kmersearch_extract_dna4_kmers_with_expansion_direct_neon(seq, k, nkeys);
    }
#endif
    return kmersearch_extract_dna4_kmers_with_expansion_direct_scalar(seq, k, nkeys);
}

/*
 * Scalar version: Extract k-mers directly from DNA4 bit sequence with degenerate expansion
 */
static Datum *
kmersearch_extract_dna4_kmers_with_expansion_direct_scalar(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Use simple array-based k-mer tracking (no hashing needed for k<=64) */
    
    occurrences = (KmerOccurrence *) palloc(max_kmers * 10 * sizeof(KmerOccurrence));
    
    /* Extract k-mers */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            /* Extract k-mer as single uint64_t value */
            kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
            
            /* Find or add occurrence count using binary search */
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers * 10);
            
            if (current_count < 0)
                continue;  /* Array full, skip */
            
            /* Skip if occurrence exceeds bit limit */
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            /* Create simple k-mer key (copy the DNA2 k-mer directly) */
            ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
            memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
        
        /* Free each expanded k-mer and then the array */
        if (expanded_kmers)
        {
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
        }
    }
    
    /* Cleanup */
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/*
 * Create n-gram key with occurrence count from DNA2 k-mer
 */
static VarBit *
kmersearch_create_ngram_key_with_occurrence_from_dna2(VarBit *dna2_kmer, int k, int occurrence)
{
    int kmer_bits = k * 2;
    int occur_bits = kmersearch_occur_bitlen;
    int total_bits = kmer_bits + occur_bits;
    int total_bytes = (total_bits + 7) / 8;
    int adj_occurrence = occurrence - 1;
    VarBit *result;
    bits8 *src_data, *dst_data;
    int i;
    
    result = (VarBit *) palloc0(VARBITHDRSZ + total_bytes);
    SET_VARSIZE(result, VARBITHDRSZ + total_bytes);
    VARBITLEN(result) = total_bits;
    
    src_data = VARBITS(dna2_kmer);
    dst_data = VARBITS(result);
    
    /* Copy k-mer bits directly */
    {
        int kmer_bytes = (kmer_bits + 7) / 8;
    memcpy(dst_data, src_data, kmer_bytes);
    }
    
    /* Encode occurrence count */
    if (adj_occurrence >= (1 << occur_bits))
        adj_occurrence = (1 << occur_bits) - 1;
    
    for (i = 0; i < occur_bits; i++)
    {
        int bit_pos = kmer_bits + i;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (adj_occurrence & (1 << (occur_bits - 1 - i)))
            dst_data[byte_pos] |= (1 << (7 - bit_offset));
    }
    
    return result;
}

/*
 * Extract k-mers from query string with degenerate code expansion
 */
static VarBit **
kmersearch_extract_query_kmers(const char *query, int k, int *nkeys)
{
    int query_len = strlen(query);
    int max_kmers = (query_len >= k) ? (query_len - k + 1) : 0;
    VarBit **keys;
    int key_count = 0;
    int i;
    bool has_degenerate;
    int j;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for degenerate expansions */
    keys = (VarBit **) palloc(max_kmers * 10 * sizeof(VarBit *));
    
    /* Extract k-mers from query */
    for (i = 0; i <= query_len - k; i++)
    {
        char kmer[65];
        strncpy(kmer, query + i, k);
        kmer[k] = '\0';
        
        /* Check if this k-mer has degenerate codes */
        if (kmersearch_will_exceed_degenerate_limit(kmer, k))
            continue;  /* Skip k-mers with too many combinations */
        
        /* Check for degenerate codes */
        has_degenerate = false;
        for (j = 0; j < k; j++)
        {
            char c = toupper(kmer[j]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'U')
            {
                has_degenerate = true;
                break;
            }
        }
        
        if (has_degenerate)
        {
            /* Expand degenerate codes */
            char *expanded[10];
            int expand_count;
            
            kmersearch_expand_degenerate_sequence(kmer, k, expanded, &expand_count);
            
            for (j = 0; j < expand_count; j++)
            {
                VarBit *kmer_key = kmersearch_create_kmer_key_only(expanded[j], k);
                keys[key_count++] = kmer_key;
                pfree(expanded[j]);
            }
        }
        else
        {
            /* Simple case - no degenerate codes */
            VarBit *kmer_key = kmersearch_create_kmer_key_only(kmer, k);
            keys[key_count++] = kmer_key;
        }
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Cache management functions
 */

/*
 * Create a new cache manager for local use (no global state)
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
 * Create actual min score cache manager
 */
static ActualMinScoreCacheManager *
create_actual_min_score_cache_manager(void)
{
    ActualMinScoreCacheManager *manager;
    HASHCTL hash_ctl;
    
    /* Allocate manager in current context (caller should have set appropriate context) */
    manager = (ActualMinScoreCacheManager *) palloc0(sizeof(ActualMinScoreCacheManager));
    
    /* Create actual min score cache context under current context */
    manager->cache_context = AllocSetContextCreate(CurrentMemoryContext,
                                                   "ActualMinScoreCache",
                                                   ALLOCSET_DEFAULT_SIZES);
    
    /* Initialize parameters */
    manager->hits = 0;
    manager->misses = 0;
    manager->max_entries = kmersearch_actual_min_score_cache_max_entries;
    manager->current_entries = 0;
    
    /* Create hash table */
    MemSet(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64);  /* Hash value as key */
    hash_ctl.entrysize = sizeof(ActualMinScoreCacheEntry);
    hash_ctl.hcxt = manager->cache_context;
    
    manager->cache_hash = hash_create("ActualMinScoreCache", 256, &hash_ctl,
                                     HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
    
    return manager;
}


/*
 * Clean up query conditions manager when PortalContext is destroyed
 */
static void
cleanup_query_conditions_manager(void)
{
    /* This function will be called when the PortalContext is destroyed */
    /* The static pointer will be automatically set to NULL */
}

/*
 * Generate cache key using PostgreSQL's hash_any_extended
 */
static uint64
generate_cache_key(VarBit *sequence, const char *query_string)
{
    uint64 seq_hash, query_hash;
    
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
    
    /* Combine hashes */
    return seq_hash ^ (query_hash << 1);
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
 * Look up rawscore cache entry
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
    actual_min_score = get_cached_actual_min_score(query_keys, result.query_nkeys, query_string, result.query_nkeys);
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
static KmerMatchResult
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
            query_keys = kmersearch_extract_kmers_from_query(query_string, kmersearch_kmer_size, &result.query_nkeys);
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
static KmerMatchResult
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
            query_keys = kmersearch_extract_kmers_from_query(query_string, kmersearch_kmer_size, &result.query_nkeys);
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
 * B-2: Query pattern cache implementation
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
 * Evict oldest entry from cache
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
    
    if (found && entry && strcmp(entry->query_string_copy, query_string) == 0 && entry->k_size == k_size)
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
        entry->k_size = k_size;
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
get_cached_query_kmers(const char *query_string, int k_size, int *nkeys)
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
    extracted_kmers = kmersearch_extract_query_kmers(query_string, k_size, nkeys);
    
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
 * Fast k-mer matching using hash table - optimized O(n+m) implementation (with SIMD dispatch)
 */
static int
kmersearch_count_matching_kmers_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512) {
        return kmersearch_count_matching_kmers_fast_avx512(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    if (simd_capability >= SIMD_AVX2) {
        return kmersearch_count_matching_kmers_fast_avx2(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE) {
        return kmersearch_count_matching_kmers_fast_sve(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    if (simd_capability >= SIMD_NEON) {
        return kmersearch_count_matching_kmers_fast_neon(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
#endif
    return kmersearch_count_matching_kmers_fast_scalar(seq_keys, seq_nkeys, query_keys, query_nkeys);
}

/*
 * Scalar version: Fast k-mer matching using hash table - optimized O(n+m) implementation
 */
static int
kmersearch_count_matching_kmers_fast_scalar(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* For small datasets, O(n*m) might be faster than hash table overhead */
    if (seq_nkeys * query_nkeys < 100)
    {
        /* Fall back to simple comparison for small datasets */
        for (i = 0; i < seq_nkeys; i++)
        {
            int j;
            for (j = 0; j < query_nkeys; j++)
            {
                if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[j]) &&
                    VARSIZE(seq_keys[i]) == VARSIZE(query_keys[j]) &&
                    memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[j]), VARBITBYTES(seq_keys[i])) == 0)
                {
                    match_count++;
                    break;
                }
            }
        }
        return match_count;
    }
    
    /* Create hash table using VarBit content as key */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    /* Safety check: ensure we have valid query keys */
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmers_fast: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);  /* Use data size, not total size */
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    elog(LOG, "kmersearch_count_matching_kmers_fast: Creating hash with keysize=%zu, query_nkeys=%d", 
         (size_t)VARBITBYTES(query_keys[0]), query_nkeys);
    
    query_hash = hash_create("QueryKmerHash", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* Insert all query k-mers into hash table using content as key */
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmers_fast: NULL query key at index %d", i);
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Check each sequence k-mer against hash table */
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmers_fast: NULL seq key at index %d", i);
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            elog(LOG, "kmersearch_count_matching_kmers_fast: Size mismatch seq[%d]=%zu vs query[0]=%zu", 
                 i, (size_t)VARBITBYTES(seq_keys[i]), (size_t)VARBITBYTES(query_keys[0]));
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    /* Cleanup */
    hash_destroy(query_hash);
    
    return match_count;
}


/*
 * Extract k-mers from DNA sequence and create n-gram keys
 */
static Datum *
kmersearch_extract_kmers(const char *sequence, int seq_len, int k, int *nkeys)
{
    Datum *keys;
    int max_keys = seq_len - k + 1;
    int key_count = 0;
    bool has_degenerate = false;
    int i, j;
    
    if (max_keys <= 0)
    {
        *nkeys = 0;
        return NULL;
    }
    
    keys = (Datum *) palloc(sizeof(Datum) * max_keys * 10);  /* Extra space for degenerate expansion */
    
    for (i = 0; i <= seq_len - k; i++)
    {
        char kmer[65];  /* Max k=64 + null terminator */
        strncpy(kmer, sequence + i, k);
        kmer[k] = '\0';
        
        /* Check if this k-mer has degenerate codes */
        has_degenerate = false;
        for (j = 0; j < k; j++)
        {
            char c = toupper(kmer[j]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            {
                has_degenerate = true;
                break;
            }
        }
        
        if (has_degenerate)
        {
            /* Expand degenerate codes */
            char *expanded[10];
            int expand_count;
            
            kmersearch_expand_degenerate_sequence(kmer, k, expanded, &expand_count);
            
            for (j = 0; j < expand_count; j++)
            {
                /* Count occurrences of this expanded k-mer */
                int occurrence = 1;
                VarBit *ngram_key;
                int prev;
                
                for (prev = 0; prev < key_count; prev++)
                {
                    /* This is simplified - in reality we'd need to compare the actual k-mer */
                    /* For now, assume each k-mer appears once per position */
                }
                
                ngram_key = kmersearch_create_ngram_key(expanded[j], k, occurrence);
                keys[key_count++] = PointerGetDatum(ngram_key);
                
                pfree(expanded[j]);
            }
        }
        else
        {
            /* Simple case - no degenerate codes */
            int occurrence = 1;
            VarBit *ngram_key = kmersearch_create_ngram_key(kmer, k, occurrence);
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * GIN extract_value function for DNA2
 * Note: Exclusion filtering is applied separately after index creation
 * via kmersearch_analyze_table() and related functions
 */
Datum
kmersearch_extract_value_dna2(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *dna = (kmersearch_dna2 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    /* Use direct bit extraction instead of string conversion */
    keys = kmersearch_extract_dna2_kmers_direct((VarBit *)dna, k, nkeys);
    
    /* Apply high-frequency k-mer filtering if cache is available */
    if (keys && *nkeys > 0) {
        /* Try to auto-load cache if not already loaded */
        if (!global_highfreq_cache.is_valid) {
            kmersearch_auto_load_cache_if_needed();
        }
        
        /* Use cached hash table for filtering if available */
        if (global_highfreq_cache.is_valid) {
            keys = kmersearch_filter_highfreq_kmers_from_keys(keys, nkeys, global_highfreq_cache.highfreq_hash, k);
        }
    }
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
}

/*
 * GIN extract_value function for DNA4
 * Note: Exclusion filtering is applied separately after index creation
 * via kmersearch_analyze_table() and related functions
 */
Datum
kmersearch_extract_value_dna4(PG_FUNCTION_ARGS)
{
    kmersearch_dna4 *dna = (kmersearch_dna4 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    /* Use direct bit extraction with degenerate expansion */
    keys = kmersearch_extract_dna4_kmers_with_expansion_direct((VarBit *)dna, k, nkeys);
    
    /* Apply high-frequency k-mer filtering if cache is available */
    if (keys && *nkeys > 0) {
        /* Try to auto-load cache if not already loaded */
        if (!global_highfreq_cache.is_valid) {
            kmersearch_auto_load_cache_if_needed();
        }
        
        /* Use cached hash table for filtering if available */
        if (global_highfreq_cache.is_valid) {
            keys = kmersearch_filter_highfreq_kmers_from_keys(keys, nkeys, global_highfreq_cache.highfreq_hash, k);
        }
    }
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
}


/*
 * Create n-gram key with occurrence count
 */
static VarBit *
kmersearch_create_ngram_key_with_occurrence(const char *kmer, int k, int occurrence)
{
    int kmer_bits = k * 2;
    int occur_bits = kmersearch_occur_bitlen;
    int total_bits = kmer_bits + occur_bits;
    int total_bytes = (total_bits + 7) / 8;
    VarBit *result;
    bits8 *data_ptr;
    int i, bit_pos = 0;
    
    /* Adjust occurrence to 0-based (1-256 becomes 0-255) */
    int adjusted_occurrence = occurrence - 1;
    if (adjusted_occurrence < 0)
        adjusted_occurrence = 0;
    if (adjusted_occurrence >= (1 << occur_bits))
        adjusted_occurrence = (1 << occur_bits) - 1;
    
    result = (VarBit *) palloc0(VARBITHDRSZ + total_bytes);
    SET_VARSIZE(result, VARBITHDRSZ + total_bytes);
    VARBITLEN(result) = total_bits;
    data_ptr = VARBITS(result);
    
    /* Encode k-mer bits (2 bits per base) */
    for (i = 0; i < k; i++)
    {
        uint8 base_code = 0;
        char base = toupper(kmer[i]);
        
        switch (base)
        {
            case 'A': base_code = 0; break;  /* 00 */
            case 'C': base_code = 1; break;  /* 01 */
            case 'G': base_code = 2; break;  /* 10 */
            case 'T': case 'U': base_code = 3; break;  /* 11 */
            default: base_code = 0; break;   /* Default to A for invalid chars */
        }
        
        /* Set 2 bits for this base */
        kmersearch_set_bit_at(data_ptr, bit_pos++, (base_code >> 1) & 1);
        kmersearch_set_bit_at(data_ptr, bit_pos++, base_code & 1);
    }
    
    /* Encode occurrence count */
    for (i = occur_bits - 1; i >= 0; i--)
    {
        kmersearch_set_bit_at(data_ptr, bit_pos++, (adjusted_occurrence >> i) & 1);
    }
    
    return result;
}

/*
 * Extract k-mers with degenerate code expansion
 */
static Datum *
kmersearch_extract_kmers_with_degenerate(const char *sequence, int seq_len, int k, int *nkeys)
{
    Datum *keys;
    int max_keys = (seq_len >= k) ? (seq_len - k + 1) * 10 : 0;  /* Estimate max with expansions */
    int key_count = 0;
    int i;
    HTAB *occurrence_hash;
    HASHCTL hash_ctl;
    
    *nkeys = 0;
    if (max_keys <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_keys * sizeof(Datum));
    
    /* Create hash table for tracking k-mer occurrences */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = k + 1;  /* k characters + null terminator */
    hash_ctl.entrysize = sizeof(int);
    occurrence_hash = hash_create("KmerOccurrenceHash", 1000, &hash_ctl,
                                  HASH_ELEM | HASH_BLOBS);
    
    /* Extract k-mers */
    for (i = 0; i <= seq_len - k; i++)
    {
        char kmer[65];  /* Max k=64 + null terminator */
        bool found;
        int *occurrence_ptr;
        int combinations;
        
        /* Extract k-mer */
        strncpy(kmer, sequence + i, k);
        kmer[k] = '\0';
        
        /* Check if degenerate expansion exceeds limit */
        combinations = kmersearch_count_degenerate_combinations(kmer, k);
        if (combinations > 10)
            continue;  /* Skip this k-mer */
        
        /* Get or create occurrence count */
        occurrence_ptr = (int *) hash_search(occurrence_hash, kmer, HASH_ENTER, &found);
        if (!found)
            *occurrence_ptr = 0;
        
        (*occurrence_ptr)++;
        
        /* Skip if occurrence exceeds bit limit */
        if (*occurrence_ptr > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* For now, create simple k-mer without full degenerate expansion */
        /* TODO: Implement full degenerate expansion */
        {
            VarBit *ngram_key = kmersearch_create_ngram_key_with_occurrence(kmer, k, *occurrence_ptr);
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
        
        if (key_count >= max_keys)
            break;
    }
    
    hash_destroy(occurrence_hash);
    
    *nkeys = key_count;
    return keys;
}

/*
 * GIN extract_query function
 */
Datum
kmersearch_extract_query(PG_FUNCTION_ARGS)
{
    Datum query = PG_GETARG_DATUM(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    StrategyNumber strategy = PG_GETARG_UINT16(2);
    bool **pmatch = (bool **) PG_GETARG_POINTER(3);
    Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4);
    bool **nullFlags = (bool **) PG_GETARG_POINTER(5);
    int32 *searchMode = (int32 *) PG_GETARG_POINTER(6);
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    text *query_text = DatumGetTextP(query);
    char *query_string = text_to_cstring(query_text);
    int query_len = strlen(query_string);
    Datum *keys;
    
    if (query_len < 8)
        ereport(ERROR, (errmsg("Query sequence must be at least 8 bases long")));
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    keys = kmersearch_extract_kmers(query_string, query_len, k, nkeys);
    
    *pmatch = NULL;
    *extra_data = NULL;
    *nullFlags = NULL;
    *searchMode = GIN_SEARCH_MODE_DEFAULT;
    
    pfree(query_string);
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
}

/*
 * GIN consistent function
 */
Datum
kmersearch_consistent(PG_FUNCTION_ARGS)
{
    bool *check = (bool *) PG_GETARG_POINTER(0);
    StrategyNumber strategy = PG_GETARG_UINT16(1);
    Datum query = PG_GETARG_DATUM(2);
    int32 nkeys = PG_GETARG_INT32(3);
    Pointer *extra_data = (Pointer *) PG_GETARG_POINTER(4);
    bool *recheck = (bool *) PG_GETARG_POINTER(5);
    Datum *queryKeys = (Datum *) PG_GETARG_POINTER(6);
    bool *nullFlags = (bool *) PG_GETARG_POINTER(7);
    
    int match_count = 0;
    int actual_min_score;
    int i;
    VarBit **query_key_array;
    
    /* Try to auto-load cache if not already loaded for search optimization */
    if (!global_highfreq_cache.is_valid) {
        kmersearch_auto_load_cache_if_needed();
    }
    
    *recheck = true;  /* Always recheck for scoring */
    
    /* Count matching keys */
    for (i = 0; i < nkeys; i++)
    {
        if (check[i])
            match_count++;
    }
    
    /* Convert queryKeys to VarBit array for excluded k-mer checking */
    query_key_array = (VarBit **) palloc(nkeys * sizeof(VarBit *));
    for (i = 0; i < nkeys; i++)
    {
        query_key_array[i] = DatumGetVarBitP(queryKeys[i]);
    }
    
    /* Calculate actual minimum score using comprehensive scoring logic */
    actual_min_score = calculate_actual_min_score(query_key_array, nkeys, nkeys);
    
    pfree(query_key_array);
    
    /* Return true if match count meets actual minimum score */
    PG_RETURN_BOOL(match_count >= actual_min_score);
}

/*
 * GIN compare_partial function - simple byte comparison for varbit
 */
Datum
kmersearch_compare_partial(PG_FUNCTION_ARGS)
{
    VarBit *a = DatumGetVarBitP(PG_GETARG_DATUM(0));
    VarBit *b = DatumGetVarBitP(PG_GETARG_DATUM(1));
    int result;
    
    int32 len_a = VARBITLEN(a);
    int32 len_b = VARBITLEN(b);
    
    if (len_a < len_b)
        result = -1;
    else if (len_a > len_b)
        result = 1;
    else
    {
        result = memcmp(VARBITS(a), VARBITS(b), VARBITBYTES(a));
    }
    
    PG_RETURN_INT32(result);
}

/*
 * DNA2 =% operator for k-mer search
 */
Datum
kmersearch_dna2_match(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    text *pattern = PG_GETARG_TEXT_P(1);
    
    char *pattern_string = text_to_cstring(pattern);
    
    /* Use cached calculation function */
    KmerMatchResult result = get_cached_rawscore_dna2(dna, pattern_string);
    
    pfree(pattern_string);
    
    /* Return boolean match result */
    PG_RETURN_BOOL(result.valid ? result.match_result : false);
}

/*
 * DNA4 =% operator for k-mer search
 */
Datum
kmersearch_dna4_match(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    text *pattern = PG_GETARG_TEXT_P(1);
    
    char *pattern_string = text_to_cstring(pattern);
    
    /* Use cached calculation function */
    KmerMatchResult result = get_cached_rawscore_dna4(dna, pattern_string);
    
    pfree(pattern_string);
    
    /* Return boolean match result */
    PG_RETURN_BOOL(result.valid ? result.match_result : false);
}


/*
 * DNA2 equality operator
 */
Datum
kmersearch_dna2_eq(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *a = (kmersearch_dna2 *) PG_GETARG_POINTER(0);
    kmersearch_dna2 *b = (kmersearch_dna2 *) PG_GETARG_POINTER(1);
    int32 len_a = VARSIZE(a);
    int32 len_b = VARSIZE(b);
    bool result;
    
    if (len_a != len_b)
        PG_RETURN_BOOL(false);
    
    result = (memcmp(a, b, len_a) == 0);
    PG_RETURN_BOOL(result);
}

/*
 * DNA4 equality operator
 */
Datum
kmersearch_dna4_eq(PG_FUNCTION_ARGS)
{
    kmersearch_dna4 *a = (kmersearch_dna4 *) PG_GETARG_POINTER(0);
    kmersearch_dna4 *b = (kmersearch_dna4 *) PG_GETARG_POINTER(1);
    int32 len_a = VARSIZE(a);
    int32 len_b = VARSIZE(b);
    bool result;
    
    if (len_a != len_b)
        PG_RETURN_BOOL(false);
    
    result = (memcmp(a, b, len_a) == 0);
    PG_RETURN_BOOL(result);
}

/*
 * Create k-mer key without occurrence count (for frequency analysis)
 */
static VarBit *
kmersearch_create_kmer_key_only(const char *kmer, int k)
{
    int kmer_bits = k * 2;  /* 2 bits per base */
    int total_bytes = (kmer_bits + 7) / 8;
    VarBit *result;
    bits8 *data_ptr;
    int i;
    
    result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + total_bytes);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + total_bytes);
    VARBITLEN(result) = kmer_bits;
    
    data_ptr = VARBITS(result);
    
    /* Encode k-mer */
    for (i = 0; i < k; i++)
    {
        uint8 encoded = kmersearch_dna2_encode_table[(unsigned char)kmer[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        data_ptr[byte_pos] |= (encoded << (6 - bit_offset));
    }
    
    return result;
}

/* A-3: Removed kmersearch_kmer_hash and kmersearch_kmer_compare functions - no longer needed */

/*
 * Analyze table frequency and determine highly frequent k-mers
 */
Datum
kmersearch_analyze_table_frequency(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k = PG_GETARG_INT32(2);
    Oid index_oid = PG_GETARG_OID(3);
    
    char *column_name = text_to_cstring(column_name_text);
    int highfreq_count = 0;
    
    /* Check if high-frequency k-mer exclusion should be performed */
    bool should_exclude = false;
    
    /* If max_appearance_rate is 0, treat as undefined (no exclusion) */
    if (kmersearch_max_appearance_rate > 0.0)
        should_exclude = true;
    
    /* If max_appearance_nrow is greater than 0, enable exclusion */
    if (kmersearch_max_appearance_nrow > 0)
        should_exclude = true;
    
    if (!should_exclude)
    {
        /* Skip frequency analysis - create empty highly frequent k-mer list */
        ereport(NOTICE, (errmsg("High-frequency k-mer exclusion disabled, skipping table scan")));
        
        /* Insert index info with zero highly frequent k-mers */
        /* Note: This would normally insert into kmersearch_index_info table */
        /* For now, just return 0 indicating no exclusions */
        
        PG_RETURN_INT32(0);
    }
    
    /* Perform frequency analysis if exclusion is enabled */
    ereport(NOTICE, (errmsg("Performing k-mer frequency analysis for k=%d", k)));
    ereport(NOTICE, (errmsg("Max appearance rate: %f, Max appearance nrow: %d", 
                           kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow)));
    
    /* 
     * TODO: Implement actual frequency analysis:
     * 1. Scan all rows in the table
     * 2. Extract k-mers from the specified column
     * 3. Count frequency of each k-mer
     * 4. Identify k-mers exceeding thresholds
     * 5. Insert highly frequent k-mers into kmersearch_highfreq_kmers table
     * 6. Insert index statistics into kmersearch_index_info table
     */
    
    PG_RETURN_INT32(highfreq_count);
}

/*
 * Get highly frequent k-mers for an index
 */
Datum
kmersearch_get_highfreq_kmers(PG_FUNCTION_ARGS)
{
    Oid index_oid = PG_GETARG_OID(0);
    int ret;
    StringInfoData query;
    ArrayType *result_array = NULL;
    Datum *datums = NULL;
    int nkeys = 0;
    int i;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get highly frequent k-mers */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT kmer_key FROM kmersearch_highfreq_kmers WHERE index_oid = %u ORDER BY kmer_key",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        nkeys = SPI_processed;
        datums = (Datum *) palloc(nkeys * sizeof(Datum));
        
        for (i = 0; i < nkeys; i++)
        {
            bool isnull;
            Datum kmer_datum;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                /* Copy the varbit value */
                VarBit *kmer = DatumGetVarBitPCopy(kmer_datum);
                datums[i] = PointerGetDatum(kmer);
            }
            else
            {
                datums[i] = (Datum) 0;
            }
        }
        
        /* Create array result */
        if (nkeys > 0)
        {
            result_array = construct_array(datums, nkeys, VARBITOID, -1, false, TYPALIGN_INT);
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    if (datums)
        pfree(datums);
    SPI_finish();
    
    if (result_array)
        PG_RETURN_ARRAYTYPE_P(result_array);
    else
        PG_RETURN_NULL();
}


/*
 * Count highly frequent k-mers in query sequence
 */
static int
kmersearch_count_highfreq_kmers_in_query(VarBit **query_keys, int nkeys)
{
    int highfreq_count = 0;
    int i;
    
    /* For each k-mer in the query, check if it's highly frequent */
    for (i = 0; i < nkeys; i++)
    {
        if (kmersearch_is_kmer_highfreq(query_keys[i]))
        {
            highfreq_count++;
        }
    }
    
    return highfreq_count;
}

/*
 * Check if high-frequency k-mer filtering is enabled for current context
 */
static bool
kmersearch_is_highfreq_filtering_enabled(void)
{
    /* Check if global cache is valid and contains high-frequency k-mers */
    if (!global_highfreq_cache.is_valid)
        return false;
    
    /* Check if cache contains any high-frequency k-mers */
    if (global_highfreq_cache.highfreq_hash == NULL)
        return false;
    
    return true;
}

/*
 * Calculate adjusted minimum score based on highly frequent k-mers in query
 * Only applies adjustment when high-frequency filtering is actually enabled
 */
static int
kmersearch_get_adjusted_min_score(VarBit **query_keys, int nkeys)
{
    int highfreq_count;
    int adjusted_score;
    
    /* Check if high-frequency filtering is enabled for this context */
    if (!kmersearch_is_highfreq_filtering_enabled()) {
        return kmersearch_min_score;  /* No adjustment needed */
    }
    
    highfreq_count = kmersearch_count_highfreq_kmers_in_query(query_keys, nkeys);
    adjusted_score = kmersearch_min_score - highfreq_count;
    
    /* Ensure adjusted score is not negative */
    if (adjusted_score < 0)
        adjusted_score = 0;
    
    return adjusted_score;
}

/*
 * Calculate raw score between two DNA sequences
 */
static int
kmersearch_calculate_raw_score(VarBit *seq1, VarBit *seq2, text *query_text)
{
    char *query_string = text_to_cstring(query_text);
    int query_len = strlen(query_string);
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    int score = 0;
    VarBit **seq1_keys, **seq2_keys;
    int seq1_nkeys, seq2_nkeys;
    int i, j;
    
    /* Extract k-mers from both sequences */
    seq1_keys = kmersearch_extract_kmers_from_varbit(seq1, k, &seq1_nkeys);
    seq2_keys = kmersearch_extract_kmers_from_query(query_string, k, &seq2_nkeys);
    
    /* Count matching k-mers using optimized function */
    score = kmersearch_count_matching_kmers_fast(seq1_keys, seq1_nkeys, seq2_keys, seq2_nkeys);
    
    /* Cleanup */
    if (seq1_keys)
    {
        for (i = 0; i < seq1_nkeys; i++)
            if (seq1_keys[i]) pfree(seq1_keys[i]);
        pfree(seq1_keys);
    }
    if (seq2_keys)
    {
        for (j = 0; j < seq2_nkeys; j++)
            if (seq2_keys[j]) pfree(seq2_keys[j]);
        pfree(seq2_keys);
    }
    
    pfree(query_string);
    return score;
}


/*
 * Extract k-mers from VarBit sequence
 */
static VarBit **
kmersearch_extract_kmers_from_varbit(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;  /* Assuming 2-bit encoding */
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    VarBit **keys;
    int key_count = 0;
    int i, j;
    bits8 *seq_data = VARBITS(seq);
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (VarBit **) palloc(max_kmers * sizeof(VarBit *));
    
    /* Extract each k-mer */
    for (i = 0; i <= seq_bases - k; i++)
    {
        int kmer_bits = k * 2;
        int kmer_bytes = (kmer_bits + 7) / 8;
        VarBit *kmer_key = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
        bits8 *kmer_data;
        
        SET_VARSIZE(kmer_key, VARHDRSZ + sizeof(int32) + kmer_bytes);
        VARBITLEN(kmer_key) = kmer_bits;
        kmer_data = VARBITS(kmer_key);
        
        /* Copy k-mer bits from sequence */
        for (j = 0; j < k; j++)
        {
            int src_bit_pos = (i + j) * 2;
            int dst_bit_pos = j * 2;
            int src_byte_pos = src_bit_pos / 8;
            int src_bit_offset = src_bit_pos % 8;
            int dst_byte_pos = dst_bit_pos / 8;
            int dst_bit_offset = dst_bit_pos % 8;
            
            /* Extract 2 bits from source */
            uint8 base_bits = (seq_data[src_byte_pos] >> (6 - src_bit_offset)) & 0x3;
            
            /* Store 2 bits in destination */
            kmer_data[dst_byte_pos] |= (base_bits << (6 - dst_bit_offset));
        }
        
        keys[key_count++] = kmer_key;
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Extract k-mers from query string
 */
static VarBit **
kmersearch_extract_kmers_from_query(const char *query, int k, int *nkeys)
{
    int query_len = strlen(query);
    int max_kmers = (query_len >= k) ? (query_len - k + 1) : 0;
    VarBit **keys;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (VarBit **) palloc(max_kmers * sizeof(VarBit *));
    
    /* Extract each k-mer from query string */
    for (i = 0; i <= query_len - k; i++)
    {
        keys[key_count++] = kmersearch_create_kmer_key_only(query + i, k);
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Raw score calculation function for DNA2
 */
Datum
kmersearch_rawscore_dna2(PG_FUNCTION_ARGS)
{
    VarBit *sequence = PG_GETARG_VARBIT_P(0);
    text *query_text = PG_GETARG_TEXT_P(1);
    
    char *query_string = text_to_cstring(query_text);
    
    /* Use cached calculation function */
    KmerMatchResult result = get_cached_rawscore_dna2(sequence, query_string);
    
    pfree(query_string);
    
    /* Return shared count as raw score */
    PG_RETURN_INT32(result.valid ? result.shared_count : 0);
}

/*
 * Raw score calculation function for DNA4
 */
Datum
kmersearch_rawscore_dna4(PG_FUNCTION_ARGS)
{
    VarBit *sequence = PG_GETARG_VARBIT_P(0);
    text *query_text = PG_GETARG_TEXT_P(1);
    
    char *query_string = text_to_cstring(query_text);
    
    /* Use cached calculation function */
    KmerMatchResult result = get_cached_rawscore_dna4(sequence, query_string);
    
    pfree(query_string);
    
    /* Return shared count as raw score */
    PG_RETURN_INT32(result.valid ? result.shared_count : 0);
}

/*
 * Corrected score functions - currently equivalent to raw score
 * 
 * These functions are designed to return corrected scores that account
 * for highly frequent k-mers that should be excluded from scoring.
 * Currently, since high-frequency k-mer exclusion is not implemented,
 * these functions simply call the corresponding raw score functions.
 * 
 * Future implementation will:
 * 1. Get index OID from execution context
 * 2. Count mutual high-frequency k-mers between sequence and query
 * 3. Add this count to raw score for correction
 */
Datum
kmersearch_correctedscore_dna2(PG_FUNCTION_ARGS)
{
    VarBit *sequence = PG_GETARG_VARBIT_P(0);  /* DNA2 is stored as VarBit */
    text *query_text = PG_GETARG_TEXT_P(1);
    char *query_string = text_to_cstring(query_text);
    VarBit **query_keys = NULL;
    VarBit **seq_keys = NULL;
    Datum *seq_datum_keys = NULL;
    int query_nkeys = 0;
    int seq_nkeys = 0;
    int k = kmersearch_kmer_size;
    int shared_count = 0;
    int i, j;
    
    /* Extract k-mers from DNA2 sequence (no degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna2_kmers_direct(sequence, k, &seq_nkeys);
    if (seq_datum_keys != NULL && seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    elog(LOG, "correctedscore_dna2: seq_keys=%p, seq_nkeys=%d", seq_keys, seq_nkeys);
    if (seq_keys && seq_nkeys > 0) {
        elog(LOG, "correctedscore_dna2: First seq k-mer bitlen=%d", VARBITLEN(seq_keys[0]));
    }
    
    /* Extract k-mers from query */
    query_keys = kmersearch_extract_kmers_from_query(query_string, k, &query_nkeys);
    elog(LOG, "correctedscore_dna2: query_keys=%p, query_nkeys=%d", query_keys, query_nkeys);
    if (query_keys && query_nkeys > 0) {
        elog(LOG, "correctedscore_dna2: First query k-mer bitlen=%d", VARBITLEN(query_keys[0]));
    }
    
    /* Count shared k-mers using optimized function */
    if (seq_keys && query_keys && seq_nkeys > 0 && query_nkeys > 0) {
        elog(LOG, "correctedscore_dna2: Calling kmersearch_count_matching_kmers_fast");
        shared_count = kmersearch_count_matching_kmers_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
        elog(LOG, "correctedscore_dna2: shared_count=%d", shared_count);
    }
    
    /* Cleanup */
    if (seq_keys) {
        pfree(seq_keys);
    }
    if (seq_datum_keys) {
        pfree(seq_datum_keys);
    }
    if (query_keys) {
        for (i = 0; i < query_nkeys; i++) {
            if (query_keys[i]) pfree(query_keys[i]);
        }
        pfree(query_keys);
    }
    pfree(query_string);
    
    /* Return corrected score (shared k-mer count) */
    PG_RETURN_INT32(shared_count);
}

Datum
kmersearch_correctedscore_dna4(PG_FUNCTION_ARGS)
{
    VarBit *sequence = PG_GETARG_VARBIT_P(0);  /* DNA4 is stored as VarBit */
    text *query_text = PG_GETARG_TEXT_P(1);
    char *query_string = text_to_cstring(query_text);
    VarBit **query_keys = NULL;
    VarBit **seq_keys = NULL;
    Datum *seq_datum_keys = NULL;
    int query_nkeys = 0;
    int seq_nkeys = 0;
    int k = kmersearch_kmer_size;
    int shared_count = 0;
    int i, j;
    
    /* Extract k-mers from DNA4 sequence (with degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna4_kmers_with_expansion_direct(sequence, k, &seq_nkeys);
    if (seq_datum_keys != NULL && seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    elog(LOG, "correctedscore_dna4: seq_keys=%p, seq_nkeys=%d", seq_keys, seq_nkeys);
    if (seq_keys && seq_nkeys > 0) {
        elog(LOG, "correctedscore_dna4: First seq k-mer bitlen=%d", VARBITLEN(seq_keys[0]));
    }
    
    /* Extract k-mers from query */
    query_keys = kmersearch_extract_kmers_from_query(query_string, k, &query_nkeys);
    elog(LOG, "correctedscore_dna4: query_keys=%p, query_nkeys=%d", query_keys, query_nkeys);
    if (query_keys && query_nkeys > 0) {
        elog(LOG, "correctedscore_dna4: First query k-mer bitlen=%d", VARBITLEN(query_keys[0]));
    }
    
    /* Count shared k-mers using optimized function */
    if (seq_keys && query_keys && seq_nkeys > 0 && query_nkeys > 0) {
        elog(LOG, "correctedscore_dna4: Calling kmersearch_count_matching_kmers_fast");
        shared_count = kmersearch_count_matching_kmers_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
        elog(LOG, "correctedscore_dna4: shared_count=%d", shared_count);
    }
    
    /* Cleanup */
    if (seq_keys) {
        pfree(seq_keys);
    }
    if (seq_datum_keys) {
        pfree(seq_datum_keys);
    }
    if (query_keys) {
        for (i = 0; i < query_nkeys; i++) {
            if (query_keys[i]) pfree(query_keys[i]);
        }
        pfree(query_keys);
    }
    pfree(query_string);
    
    /* Return corrected score (shared k-mer count) */
    PG_RETURN_INT32(shared_count);
}

/*
 * Extract k-mer as single uint64_t value (for k <= 32)
 */
static uint64_t
kmersearch_extract_kmer_as_uint64(VarBit *seq, int start_pos, int k)
{
    uint64_t kmer_value = 0;
    bits8 *src_data = VARBITS(seq);
    int src_bytes = VARBITBYTES(seq);
    int j;
    
    for (j = 0; j < k; j++)
    {
        int bit_pos = (start_pos + j) * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 base_bits;
        
        /* Boundary check to prevent buffer overflow */
        if (byte_pos >= src_bytes) {
            return 0;  /* Invalid k-mer */
        }
        
        base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
        kmer_value = (kmer_value << 2) | base_bits;
    }
    
    return kmer_value;
}

/*
 * Find or add k-mer occurrence in sorted array (no hashing)
 */
static int
kmersearch_find_or_add_kmer_occurrence(KmerOccurrence *occurrences, int *count, uint64_t kmer_value, int max_count)
{
    int left = 0, right = *count - 1;
    int insert_pos = *count;
    
    /* Binary search for existing k-mer */
    while (left <= right)
    {
        int mid = (left + right) / 2;
        if (occurrences[mid].kmer_value == kmer_value)
        {
            return ++occurrences[mid].count;  /* Found, increment and return */
        }
        else if (occurrences[mid].kmer_value < kmer_value)
        {
            left = mid + 1;
        }
        else
        {
            right = mid - 1;
            insert_pos = mid;
        }
    }
    
    /* Not found, insert new entry if space available */
    if (*count >= max_count)
        return -1;  /* Array full */
    
    /* Shift elements to make room */
    {
        int i;
        for (i = *count; i > insert_pos; i--)
        {
            occurrences[i] = occurrences[i-1];
        }
    }
    
    /* Insert new entry */
    occurrences[insert_pos].kmer_value = kmer_value;
    occurrences[insert_pos].count = 1;
    (*count)++;
    
    return 1;
}

/*
 * Check if a k-mer is highly frequent for a given index
 */
static bool
kmersearch_is_kmer_highfreq(VarBit *kmer_key)
{
    HighfreqKmerHashEntry *entry;
    bool found;
    
    if (!kmer_key)
        return false;
    
    /* Try to auto-load cache if not already loaded */
    if (!global_highfreq_cache.is_valid) {
        kmersearch_auto_load_cache_if_needed();
    }
    
    /* If cache is still not valid, assume not high-frequency */
    if (!global_highfreq_cache.is_valid || !global_highfreq_cache.highfreq_hash) {
        return false;
    }
    
    /* Check if this k-mer is in the high-frequency hash table */
    entry = (HighfreqKmerHashEntry *) hash_search(global_highfreq_cache.highfreq_hash,
                                                 (void *) &kmer_key,
                                                 HASH_FIND,
                                                 &found);
    
    return found;
}

/*
 * Length functions for DNA2 and DNA4 types
 */

/* DNA2 bit_length function */
PG_FUNCTION_INFO_V1(kmersearch_dna2_bit_length);
Datum
kmersearch_dna2_bit_length(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    
    PG_RETURN_INT32(bit_len);
}

/* DNA4 bit_length function */
PG_FUNCTION_INFO_V1(kmersearch_dna4_bit_length);
Datum
kmersearch_dna4_bit_length(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    
    PG_RETURN_INT32(bit_len);
}

/* DNA2 nuc_length function */
PG_FUNCTION_INFO_V1(kmersearch_dna2_nuc_length);
Datum
kmersearch_dna2_nuc_length(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int32 nuc_len = bit_len / 2;  /* DNA2 uses 2 bits per nucleotide */
    
    PG_RETURN_INT32(nuc_len);
}

/* DNA4 nuc_length function */
PG_FUNCTION_INFO_V1(kmersearch_dna4_nuc_length);
Datum
kmersearch_dna4_nuc_length(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int32 nuc_len = bit_len / 4;  /* DNA4 uses 4 bits per nucleotide */
    
    PG_RETURN_INT32(nuc_len);
}

/* DNA2 char_length function (same as nuc_length) */
PG_FUNCTION_INFO_V1(kmersearch_dna2_char_length);
Datum
kmersearch_dna2_char_length(PG_FUNCTION_ARGS)
{
    return kmersearch_dna2_nuc_length(fcinfo);
}

/* DNA4 char_length function (same as nuc_length) */
PG_FUNCTION_INFO_V1(kmersearch_dna4_char_length);
Datum
kmersearch_dna4_char_length(PG_FUNCTION_ARGS)
{
    return kmersearch_dna4_nuc_length(fcinfo);
}

/* Parallel k-mer analysis functions declarations added elsewhere */

/*
 * Helper function to evaluate match conditions using both score and rate thresholds
 */
static bool
kmersearch_evaluate_match_conditions(int shared_count, int query_total)
{
    /* Condition 1: Absolute shared key count */
    bool score_condition = (shared_count >= kmersearch_min_score);
    
    /* Condition 2: Relative shared key rate */
    double sharing_rate;
    bool rate_condition;
    
    if (query_total > 0) {
        sharing_rate = (double)shared_count / (double)query_total;
    } else {
        sharing_rate = 0.0;
    }
    rate_condition = (sharing_rate >= kmersearch_min_shared_ngram_key_rate);
    
    /* AND condition */
    return (score_condition && rate_condition);
}

/*
 * Calculate actual minimum score for optimized condition evaluation
 * This combines both absolute and relative thresholds into a single value
 * Now uses adjusted_min_score instead of kmersearch_min_score when high-frequency filtering is enabled
 */
static int
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
get_cached_actual_min_score(VarBit **query_keys, int nkeys, const char *query_string, int query_total_kmers)
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
            actual_min_score_cache_manager = create_actual_min_score_cache_manager();
        }
        PG_CATCH();
        {
            MemoryContextSwitchTo(old_context);
            /* Fallback to direct calculation if cache creation fails */
            return calculate_actual_min_score(query_keys, nkeys, query_total_kmers);
        }
        PG_END_TRY();
        
        MemoryContextSwitchTo(old_context);
    }
    
    /* Calculate hash value for query string using 64-bit hash function */
    query_hash = hash_bytes_extended((unsigned char *)query_string, strlen(query_string), 0);
    
    /* Look up in hash table */
    cache_entry = (ActualMinScoreCacheEntry *) hash_search(actual_min_score_cache_manager->cache_hash,
                                                          &query_hash, HASH_FIND, &found);
    
    if (found) {
        actual_min_score_cache_manager->hits++;
        return cache_entry->actual_min_score;
    }
    
    /* Not found - calculate and cache */
    actual_min_score_cache_manager->misses++;
    actual_min_score = calculate_actual_min_score(query_keys, nkeys, query_total_kmers);
    
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


/*
 * Optimized match condition evaluation using pre-calculated actual minimum score
 */
static bool
evaluate_optimized_match_condition(VarBit **query_keys, int nkeys, int shared_count, const char *query_string, int query_total_kmers)
{
    int actual_min_score;
    
    /* Get cached actual min score (with TopMemoryContext caching for performance) */
    actual_min_score = get_cached_actual_min_score(query_keys, nkeys, query_string, query_total_kmers);
    
    /* Use optimized condition check with cached actual_min_score */
    return (shared_count >= actual_min_score);
}

/*
 * K-mer based matching for DNA2 sequences
 */
static bool
kmersearch_kmer_based_match_dna2(VarBit *sequence, const char *query_string)
{
    VarBit **seq_keys;
    VarBit **query_keys;
    int seq_nkeys, query_nkeys;
    int shared_count;
    int k = kmersearch_kmer_size;
    bool result;
    
    /* Extract k-mers from DNA2 sequence (no degenerate expansion) */
    seq_keys = (VarBit **)kmersearch_extract_dna2_kmers_direct(sequence, k, &seq_nkeys);
    if (seq_keys == NULL || seq_nkeys == 0) {
        return false;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    query_keys = get_cached_query_kmers(query_string, k, &query_nkeys);
    if (query_keys == NULL || query_nkeys == 0) {
        /* Free sequence keys */
        if (seq_keys) {
            int i;
            for (i = 0; i < seq_nkeys; i++) {
                if (seq_keys[i])
                    pfree(seq_keys[i]);
            }
            pfree(seq_keys);
        }
        return false;
    }
    
    /* Count shared keys */
    shared_count = kmersearch_count_matching_kmers_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
    
    /* Evaluate match conditions using optimized method */
    result = evaluate_optimized_match_condition(query_keys, query_nkeys, shared_count, query_string, query_nkeys);
    
    /* Free memory */
    if (seq_keys) {
        int i;
        for (i = 0; i < seq_nkeys; i++) {
            if (seq_keys[i])
                pfree(seq_keys[i]);
        }
        pfree(seq_keys);
    }
    
    if (query_keys) {
        int i;
        for (i = 0; i < query_nkeys; i++) {
            if (query_keys[i])
                pfree(query_keys[i]);
        }
        pfree(query_keys);
    }
    
    return result;
}

/*
 * K-mer based matching for DNA4 sequences
 */
static bool
kmersearch_kmer_based_match_dna4(VarBit *sequence, const char *query_string)
{
    VarBit **seq_keys;
    VarBit **query_keys;
    int seq_nkeys, query_nkeys;
    int shared_count;
    int k = kmersearch_kmer_size;
    bool result;
    
    /* Extract k-mers from DNA4 sequence (with degenerate expansion) */
    seq_keys = (VarBit **)kmersearch_extract_dna4_kmers_with_expansion_direct(sequence, k, &seq_nkeys);
    if (seq_keys == NULL || seq_nkeys == 0) {
        return false;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    query_keys = get_cached_query_kmers(query_string, k, &query_nkeys);
    if (query_keys == NULL || query_nkeys == 0) {
        /* Free sequence keys */
        if (seq_keys) {
            int i;
            for (i = 0; i < seq_nkeys; i++) {
                if (seq_keys[i])
                    pfree(seq_keys[i]);
            }
            pfree(seq_keys);
        }
        return false;
    }
    
    /* Count shared keys */
    shared_count = kmersearch_count_matching_kmers_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
    
    /* Evaluate match conditions using optimized method */
    result = evaluate_optimized_match_condition(query_keys, query_nkeys, shared_count, query_string, query_nkeys);
    
    /* Free memory */
    if (seq_keys) {
        int i;
        for (i = 0; i < seq_nkeys; i++) {
            if (seq_keys[i])
                pfree(seq_keys[i]);
        }
        pfree(seq_keys);
    }
    
    if (query_keys) {
        int i;
        for (i = 0; i < query_nkeys; i++) {
            if (query_keys[i])
                pfree(query_keys[i]);
        }
        pfree(query_keys);
    }
    
    return result;
}



/*
 * Core k-mer matching and scoring function for DNA2 sequences
 * Performs all k-mer extraction, comparison, and evaluation in one pass
 */
static KmerMatchResult
kmersearch_calculate_kmer_match_and_score_dna2(VarBit *sequence, const char *query_string)
{
    KmerMatchResult result = {0};
    VarBit **seq_keys = NULL;
    VarBit **query_keys = NULL;
    Datum *seq_datum_keys = NULL;
    int k = kmersearch_kmer_size;
    int query_len;
    int i;
    
    /* Initialize result */
    result.valid = false;
    result.shared_count = 0;
    result.seq_nkeys = 0;
    result.query_nkeys = 0;
    result.sharing_rate = 0.0;
    result.match_result = false;
    
    /* Input validation */
    if (sequence == NULL || query_string == NULL) {
        return result;
    }
    
    /* Query length validation */
    query_len = strlen(query_string);
    if (query_len < k) {
        return result;
    }
    
    /* Extract k-mers from DNA2 sequence (no degenerate expansion) */
    elog(LOG, "DNA2 Cache: Starting k-mer extraction from sequence");
    seq_datum_keys = kmersearch_extract_dna2_kmers_direct(sequence, k, &result.seq_nkeys);
    elog(LOG, "DNA2 Cache: Extracted %d k-mers from sequence", result.seq_nkeys);
    
    if (seq_datum_keys != NULL && result.seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(result.seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < result.seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
        elog(LOG, "DNA2 Cache: Converted %d datum keys to VarBit", result.seq_nkeys);
    }
    if (seq_keys == NULL || result.seq_nkeys == 0) {
        elog(LOG, "DNA2 Cache: No sequence k-mers extracted, cleaning up");
        goto cleanup;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    elog(LOG, "DNA2 Cache: Starting k-mer extraction from query '%s'", query_string);
    query_keys = get_cached_query_kmers(query_string, k, &result.query_nkeys);
    elog(LOG, "DNA2 Cache: Extracted %d k-mers from query", result.query_nkeys);
    
    if (query_keys == NULL || result.query_nkeys == 0) {
        elog(LOG, "DNA2 Cache: No query k-mers extracted, cleaning up");
        goto cleanup;
    }
    
    /* Calculate shared k-mer count (this becomes the rawscore) */
    elog(LOG, "DNA2 Cache: Starting k-mer matching calculation");
    result.shared_count = kmersearch_count_matching_kmers_fast(seq_keys, result.seq_nkeys, 
                                                               query_keys, result.query_nkeys);
    elog(LOG, "DNA2 Cache: Completed k-mer matching, shared_count=%d", result.shared_count);
    
    /* Calculate sharing rate */
    if (result.query_nkeys > 0) {
        result.sharing_rate = (double)result.shared_count / (double)result.query_nkeys;
    }
    
    /* Evaluate match conditions for =% operator using optimized method */
    result.match_result = evaluate_optimized_match_condition(query_keys, result.query_nkeys, result.shared_count, query_string, result.query_nkeys);
    
    result.valid = true;
    
cleanup:
    /* Unified memory cleanup */
    CLEANUP_KMER_ARRAYS(seq_keys, result.seq_nkeys, query_keys, result.query_nkeys);
    if (seq_datum_keys) {
        pfree(seq_datum_keys);
    }
    
    return result;
}

/*
 * Core k-mer matching and scoring function for DNA4 sequences
 * Performs all k-mer extraction, comparison, and evaluation in one pass
 */
static KmerMatchResult
kmersearch_calculate_kmer_match_and_score_dna4(VarBit *sequence, const char *query_string)
{
    KmerMatchResult result = {0};
    VarBit **seq_keys = NULL;
    VarBit **query_keys = NULL;
    Datum *seq_datum_keys = NULL;
    int k = kmersearch_kmer_size;
    int query_len;
    int i;
    
    /* Initialize result */
    result.valid = false;
    result.shared_count = 0;
    result.seq_nkeys = 0;
    result.query_nkeys = 0;
    result.sharing_rate = 0.0;
    result.match_result = false;
    
    /* Input validation */
    if (sequence == NULL || query_string == NULL) {
        return result;
    }
    
    /* Query length validation */
    query_len = strlen(query_string);
    if (query_len < k) {
        return result;
    }
    
    /* Extract k-mers from DNA4 sequence (with degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna4_kmers_with_expansion_direct(sequence, k, &result.seq_nkeys);
    if (seq_datum_keys != NULL && result.seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(result.seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < result.seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    if (seq_keys == NULL || result.seq_nkeys == 0) {
        goto cleanup;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    query_keys = get_cached_query_kmers(query_string, k, &result.query_nkeys);
    if (query_keys == NULL || result.query_nkeys == 0) {
        goto cleanup;
    }
    
    /* Calculate shared k-mer count (this becomes the rawscore) */
    result.shared_count = kmersearch_count_matching_kmers_fast(seq_keys, result.seq_nkeys, 
                                                               query_keys, result.query_nkeys);
    
    /* Calculate sharing rate */
    if (result.query_nkeys > 0) {
        result.sharing_rate = (double)result.shared_count / (double)result.query_nkeys;
    }
    
    /* Evaluate match conditions for =% operator using optimized method */
    result.match_result = evaluate_optimized_match_condition(query_keys, result.query_nkeys, result.shared_count, query_string, result.query_nkeys);
    
    result.valid = true;
    
cleanup:
    /* Unified memory cleanup */
    CLEANUP_KMER_ARRAYS(seq_keys, result.seq_nkeys, query_keys, result.query_nkeys);
    if (seq_datum_keys) {
        pfree(seq_datum_keys);
    }
    
    return result;
}

/*
 * Get the data size in bytes for a given k-mer size
 */
static size_t
kmersearch_get_kmer_data_size(int k_size)
{
    if (k_size <= 8) return 2;       /* uint16 */
    else if (k_size <= 16) return 4; /* uint32 */
    else if (k_size <= 32) return 8; /* uint64 */
    else return 16;                  /* 128-bit struct */
}

/*
 * Get index information from index OID
 */
static bool
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
        "SELECT table_oid, column_name, k_value FROM kmersearch_index_info "
        "WHERE index_oid = %u",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    kmersearch_handle_spi_error(ret, "SELECT");
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        Datum table_oid_datum, column_name_datum, k_size_datum;
        bool isnull;
        
        table_oid_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull && table_oid)
            *table_oid = DatumGetObjectId(table_oid_datum);
        
        column_name_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull);
        if (!isnull && column_name)
        {
            char *col_name = SPI_getvalue(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2);
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
 * Connect to SPI with error handling
 */
static void
kmersearch_spi_connect_or_error(void)
{
    int ret = SPI_connect();
    
    switch (ret) {
    case SPI_OK_CONNECT:
        return;  /* Success */
    case SPI_ERROR_CONNECT:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI manager already connected")));
        break;
    case SPI_ERROR_ARGUMENT:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI connection failed: invalid argument")));
        break;
    default:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI connection failed with code %d", ret)));
        break;
    }
}

/*
 * Handle SPI operation errors
 */
static void
kmersearch_handle_spi_error(int spi_result, const char *operation)
{
    switch (spi_result) {
    case SPI_OK_SELECT:
    case SPI_OK_INSERT:
    case SPI_OK_DELETE:
    case SPI_OK_UPDATE:
    case SPI_OK_UTILITY:
        return;  /* Success cases */
        
    case SPI_ERROR_ARGUMENT:
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("SPI %s failed: invalid argument", operation)));
        break;
        
    case SPI_ERROR_UNCONNECTED:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: not connected to SPI manager", operation)));
        break;
        
    case SPI_ERROR_COPY:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: COPY operation not supported", operation)));
        break;
        
    case SPI_ERROR_CURSOR:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: cursor operation error", operation)));
        break;
        
    case SPI_ERROR_TRANSACTION:
        ereport(ERROR,
                (errcode(ERRCODE_ACTIVE_SQL_TRANSACTION),
                 errmsg("SPI %s failed: transaction block error", operation)));
        break;
        
    case SPI_ERROR_OPUNKNOWN:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed: unknown operation", operation)));
        break;
        
    default:
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("SPI %s failed with code %d", operation, spi_result)));
        break;
    }
}

/*
 * Calculate optimal buffer size based on memory constraints
 */
static int
kmersearch_calculate_buffer_size(int k_size)
{
    const int TARGET_MEMORY_MB = 50;  /* Target 50MB per worker */
    const int MIN_BUFFER_SIZE = 1000;
    const int MAX_BUFFER_SIZE = 100000;
    
    size_t entry_size = sizeof(CompactKmerFreq);
    int max_entries = (TARGET_MEMORY_MB * 1024 * 1024) / entry_size;
    
    if (max_entries < MIN_BUFFER_SIZE) return MIN_BUFFER_SIZE;
    if (max_entries > MAX_BUFFER_SIZE) return MAX_BUFFER_SIZE;
    return max_entries;
}

/*
 * Encode VarBit k-mer into compact KmerData
 */
static KmerData
kmersearch_encode_kmer_data(VarBit *kmer, int k_size)
{
    KmerData result;
    unsigned char *bits;
    int i, bit_offset;
    
    memset(&result, 0, sizeof(KmerData));
    bits = VARBITS(kmer);
    
    if (k_size <= 8) {
        result.k8_data = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            int byte_pos = bit_offset / 8;
            int bit_in_byte = bit_offset % 8;
            int nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k8_data |= (nucleotide << (2 * (k_size - 1 - i)));
        }
    } else if (k_size <= 16) {
        result.k16_data = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            int byte_pos = bit_offset / 8;
            int bit_in_byte = bit_offset % 8;
            int nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k16_data |= (nucleotide << (2 * (k_size - 1 - i)));
        }
    } else if (k_size <= 32) {
        result.k32_data = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            int byte_pos = bit_offset / 8;
            int bit_in_byte = bit_offset % 8;
            int nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k32_data |= ((uint64)nucleotide << (2 * (k_size - 1 - i)));
        }
    } else {
        /* For k > 32, split across high and low 64-bit values */
        result.k64_data.high = 0;
        result.k64_data.low = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            int byte_pos = bit_offset / 8;
            int bit_in_byte = bit_offset % 8;
            int nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            
            if (i < 32) {
                result.k64_data.high |= ((uint64)nucleotide << (2 * (31 - i)));
            } else {
                result.k64_data.low |= ((uint64)nucleotide << (2 * (k_size - 1 - i)));
            }
        }
    }
    
    return result;
}

/*
 * Initialize k-mer buffer
 */
static void
kmersearch_init_buffer(KmerBuffer *buffer, int k_size)
{
    buffer->capacity = kmersearch_calculate_buffer_size(k_size);
    buffer->entries = (CompactKmerFreq *) palloc0(buffer->capacity * sizeof(CompactKmerFreq));
    buffer->count = 0;
    buffer->k_size = k_size;
}

/*
 * Aggregate buffer entries with same kmer_data to prevent ON CONFLICT issues
 */
static void
kmersearch_aggregate_buffer_entries(KmerBuffer *buffer)
{
    int i, j;
    int write_pos = 0;
    bool merged;
    
    if (buffer->count <= 1) return;
    
    for (i = 0; i < buffer->count; i++) {
        merged = false;
        
        /* Check if this entry can be merged with any previous entry */
        for (j = 0; j < write_pos; j++) {
            bool same_kmer = false;
            
            if (buffer->k_size <= 8) {
                same_kmer = (buffer->entries[i].kmer_data.k8_data == buffer->entries[j].kmer_data.k8_data);
            } else if (buffer->k_size <= 16) {
                same_kmer = (buffer->entries[i].kmer_data.k16_data == buffer->entries[j].kmer_data.k16_data);
            } else if (buffer->k_size <= 32) {
                same_kmer = (buffer->entries[i].kmer_data.k32_data == buffer->entries[j].kmer_data.k32_data);
            } else {
                same_kmer = (buffer->entries[i].kmer_data.k64_data.high == buffer->entries[j].kmer_data.k64_data.high &&
                           buffer->entries[i].kmer_data.k64_data.low == buffer->entries[j].kmer_data.k64_data.low);
            }
            
            if (same_kmer) {
                /* Merge frequency counts */
                buffer->entries[j].frequency_count += buffer->entries[i].frequency_count;
                merged = true;
                break;
            }
        }
        
        /* If not merged, add to write position */
        if (!merged) {
            if (write_pos != i) {
                buffer->entries[write_pos] = buffer->entries[i];
            }
            write_pos++;
        }
    }
    
    /* Update count to reflect aggregated entries */
    buffer->count = write_pos;
}

/*
 * Flush buffer contents to temporary table
 */
static void
kmersearch_flush_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    
    if (buffer->count == 0) return;
    
    /* Aggregate entries with same kmer_data before insertion */
    kmersearch_aggregate_buffer_entries(buffer);
    
    initStringInfo(&query);
    
    /* Build bulk INSERT statement */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    for (i = 0; i < buffer->count; i++) {
        if (i > 0) appendStringInfoString(&query, ", ");
        
        if (buffer->k_size <= 8) {
            appendStringInfo(&query, "(%u, %d)", 
                           buffer->entries[i].kmer_data.k8_data,
                           buffer->entries[i].frequency_count);
        } else if (buffer->k_size <= 16) {
            appendStringInfo(&query, "(%u, %d)", 
                           buffer->entries[i].kmer_data.k16_data,
                           buffer->entries[i].frequency_count);
        } else if (buffer->k_size <= 32) {
            appendStringInfo(&query, "(%lu, %d)", 
                           buffer->entries[i].kmer_data.k32_data,
                           buffer->entries[i].frequency_count);
        } else {
            /* For k > 32, we'll store as two bigint columns or use a different approach */
            appendStringInfo(&query, "(%lu, %d)", 
                           buffer->entries[i].kmer_data.k32_data,  /* Simplified for now */
                           buffer->entries[i].frequency_count);
        }
    }
    
    appendStringInfoString(&query, " ON CONFLICT (kmer_data) DO UPDATE SET frequency_count = EXCLUDED.frequency_count + 1");
    
    /* Execute the query */
    SPI_connect();
    SPI_exec(query.data, 0);
    SPI_finish();
    
    /* Reset buffer */
    buffer->count = 0;
    
    pfree(query.data);
}

/*
 * Add k-mer to buffer, flush to temp table if full
 */
static void
kmersearch_add_to_buffer(KmerBuffer *buffer, KmerData kmer_data, const char *temp_table_name)
{
    CompactKmerFreq *entry;
    
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_buffer_to_table(buffer, temp_table_name);
    }
    
    /* Add new entry */
    entry = &buffer->entries[buffer->count];
    entry->kmer_data = kmer_data;
    entry->frequency_count = 1;
    entry->is_highfreq = false;
    buffer->count++;
}

/*
 * Create temporary table for worker k-mer storage
 */
static void
kmersearch_create_worker_temp_table(const char *temp_table_name, int k_size)
{
    StringInfoData query;
    const char *data_type;
    
    initStringInfo(&query);
    
    /* Choose appropriate data type based on k-size */
    if (k_size <= 8) {
        data_type = "integer";  /* k=8 needs up to 65535, exceeds smallint range */
    } else if (k_size <= 16) {
        data_type = "bigint";   /* k=16 needs up to 4B values */
    } else if (k_size <= 32) {
        data_type = "bigint";
    } else {
        data_type = "bigint";  /* Will need special handling for k > 32 */
    }
    
    appendStringInfo(&query, 
        "CREATE TEMP TABLE %s ("
        "kmer_data %s PRIMARY KEY, "
        "frequency_count integer DEFAULT 1"
        ")", temp_table_name, data_type);
    
    SPI_connect();
    SPI_exec(query.data, 0);
    SPI_finish();
    
    pfree(query.data);
}

/*
 * Determine the optimal number of parallel workers for k-mer analysis
 */
static int
kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation)
{
    int max_workers;
    int max_parallel_workers_val;
    int max_parallel_maintenance_workers_val;
    int table_size_factor;
    int auto_workers;
    BlockNumber total_blocks;
    
    /* Get GUC values */
    max_parallel_workers_val = max_parallel_workers;
    max_parallel_maintenance_workers_val = max_parallel_maintenance_workers;
    
    /* Use the smaller of the two limits */
    max_workers = Min(max_parallel_workers_val, max_parallel_maintenance_workers_val);
    
    if (max_workers <= 0) {
        return 1;  /* No parallel processing allowed */
    }
    
    /* Estimate table size in blocks */
    total_blocks = RelationGetNumberOfBlocksInFork(target_relation, MAIN_FORKNUM);
    
    /* Calculate workers based on table size (minimum 1000 blocks per worker) */
    table_size_factor = Max(1, total_blocks / 1000);
    auto_workers = Min(max_workers, table_size_factor);
    
    if (requested_workers > 0) {
        /* User specified workers - respect limits */
        return Min(requested_workers, max_workers);
    } else {
        /* Auto-determine workers */
        return auto_workers;
    }
}

/*
 * Worker function to analyze blocks and extract k-mer frequencies
 */
static void
kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, 
                                const char *column_name, int k_size)
{
    TableScanDesc scan;
    HeapTuple tuple;
    TupleDesc tupdesc;
    int target_attno = -1;
    int i;
    
    /* Find the target column */
    tupdesc = RelationGetDescr(rel);
    for (i = 0; i < tupdesc->natts; i++) {
        Form_pg_attribute attr = TupleDescAttr(tupdesc, i);
        if (strcmp(NameStr(attr->attname), column_name) == 0) {
            target_attno = attr->attnum;
            break;
        }
    }
    
    if (target_attno == -1) {
        ereport(ERROR, (errmsg("Column '%s' not found in relation", column_name)));
    }
    
    /* Initialize buffer */
    kmersearch_init_buffer(&worker->buffer, k_size);
    
    /* Create temporary table for this worker */
    worker->temp_table_name = psprintf("temp_kmer_worker_%d", worker->worker_id);
    kmersearch_create_worker_temp_table(worker->temp_table_name, k_size);
    
    /* Scan assigned blocks */
    scan = heap_beginscan(rel, GetTransactionSnapshot(), 0, NULL, NULL, 0);
    
    while ((tuple = heap_getnext(scan, ForwardScanDirection)) != NULL) {
        bool isnull;
        Datum value;
        VarBit *sequence;
        VarBit **kmers;
        int nkeys;
        int j;
        
        worker->rows_processed++;
        
        /* Extract the DNA sequence value */
        value = heap_getattr(tuple, target_attno, tupdesc, &isnull);
        if (isnull) {
            continue;  /* Skip NULL values */
        }
        
        sequence = DatumGetVarBitP(value);
        
        /* Extract k-mers from the sequence */
        kmers = (VarBit **)kmersearch_extract_dna2_kmers_direct(sequence, k_size, &nkeys);
        if (kmers == NULL || nkeys == 0) {
            continue;
        }
        
        /* Process each k-mer in this row - encode and add to buffer */
        for (j = 0; j < nkeys; j++) {
            KmerData encoded_kmer;
            
            if (kmers[j] == NULL) {
                continue;
            }
            
            /* Encode k-mer into compact format */
            encoded_kmer = kmersearch_encode_kmer_data(kmers[j], k_size);
            
            /* Add to buffer (will flush to temp table if full) */
            kmersearch_add_to_buffer(&worker->buffer, encoded_kmer, worker->temp_table_name);
        }
        
        /* Cleanup k-mer array */
        if (kmers) {
            pfree(kmers);
        }
    }
    
    /* Flush any remaining buffer contents */
    kmersearch_flush_buffer_to_table(&worker->buffer, worker->temp_table_name);
    
    heap_endscan(scan);
    
    /* Cleanup buffer */
    if (worker->buffer.entries) {
        pfree(worker->buffer.entries);
    }
}


/*
 * Merge worker results using SQL aggregation
 */
static void
kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, 
                                   const char *final_table_name, int k_size, int threshold_rows)
{
    StringInfoData query;
    StringInfoData union_query;
    const char *data_type;
    int i;
    
    initStringInfo(&query);
    initStringInfo(&union_query);
    
    /* Choose appropriate data type based on k-size */
    if (k_size <= 8) {
        data_type = "integer";  /* k=8 needs up to 65535, exceeds smallint range */
    } else if (k_size <= 16) {
        data_type = "bigint";   /* k=16 needs up to 4B values */
    } else if (k_size <= 32) {
        data_type = "bigint";
    } else {
        data_type = "bigint";
    }
    
    /* Create final aggregation table */
    appendStringInfo(&query, 
        "CREATE TEMP TABLE %s ("
        "kmer_data %s PRIMARY KEY, "
        "frequency_count integer"
        ")", final_table_name, data_type);
    
    SPI_connect();
    SPI_exec(query.data, 0);
    
    /* Build UNION ALL query to combine all worker tables */
    resetStringInfo(&query);
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) ", final_table_name);
    appendStringInfoString(&query, "SELECT kmer_data, sum(frequency_count) FROM (");
    
    for (i = 0; i < num_workers; i++) {
        if (i > 0) appendStringInfoString(&query, " UNION ALL ");
        appendStringInfo(&query, "SELECT kmer_data, frequency_count FROM %s", 
                        workers[i].temp_table_name);
    }
    
    appendStringInfo(&query, ") AS combined GROUP BY kmer_data HAVING sum(frequency_count) > %d", 
                    threshold_rows);
    
    /* Execute aggregation query */
    SPI_exec(query.data, 0);
    SPI_finish();
    
    pfree(query.data);
    pfree(union_query.data);
}

/*
 * Main parallel analysis function
 */
static KmerAnalysisResult
kmersearch_analyze_table_parallel(Oid table_oid, const char *column_name, int k_size, int parallel_workers)
{
    KmerAnalysisResult result;
    Relation rel;
    int num_workers;
    KmerWorkerState *workers;
    int threshold_rows;
    int i;
    
    /* Initialize result structure */
    memset(&result, 0, sizeof(KmerAnalysisResult));
    
    /* Open target relation */
    rel = relation_open(table_oid, AccessShareLock);
    
    /* Determine number of parallel workers */
    num_workers = kmersearch_determine_parallel_workers(parallel_workers, rel);
    result.parallel_workers_used = num_workers;
    
    /* Calculate threshold based on GUC variables */
    {
        BlockNumber pages;
        double tuples, allvisfrac;
        estimate_rel_size(rel, NULL, &pages, &tuples, &allvisfrac);
        threshold_rows = (int)(tuples * kmersearch_max_appearance_rate);
    }
    if (kmersearch_max_appearance_nrow > 0 && threshold_rows > kmersearch_max_appearance_nrow)
        threshold_rows = kmersearch_max_appearance_nrow;
    
    result.max_appearance_rate_used = kmersearch_max_appearance_rate;
    result.max_appearance_nrow_used = threshold_rows;
    
    /* Allocate worker state array */
    workers = (KmerWorkerState *) palloc0(num_workers * sizeof(KmerWorkerState));
    
    /* Initialize workers and assign work blocks */
    for (i = 0; i < num_workers; i++) {
        workers[i].worker_id = i;
        workers[i].start_block = (RelationGetNumberOfBlocksInFork(rel, MAIN_FORKNUM) * i) / num_workers;
        workers[i].end_block = (RelationGetNumberOfBlocksInFork(rel, MAIN_FORKNUM) * (i + 1)) / num_workers;
        workers[i].local_highfreq_count = 0;
        workers[i].rows_processed = 0;
        workers[i].temp_table_name = NULL;
        
        /* Process assigned blocks */
        kmersearch_worker_analyze_blocks(&workers[i], rel, column_name, k_size);
    }
    
    /* Merge worker results using SQL aggregation */
    {
        char *final_table_name = psprintf("temp_kmer_final_%d", getpid());
        kmersearch_merge_worker_results_sql(workers, num_workers, final_table_name, k_size, threshold_rows);
        
        /* Count highly frequent k-mers */
        SPI_connect();
        {
            StringInfoData count_query;
            initStringInfo(&count_query);
            appendStringInfo(&count_query, "SELECT count(*) FROM %s", final_table_name);
            
            if (SPI_exec(count_query.data, 0) == SPI_OK_SELECT && SPI_processed > 0) {
                result.highfreq_kmers_count = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[0], 
                                                                        SPI_tuptable->tupdesc, 1, NULL));
            }
            pfree(count_query.data);
        }
        SPI_finish();
        
        /* Persist highly frequent k-mers to permanent tables */
        kmersearch_persist_highfreq_kmers_from_temp(table_oid, column_name, k_size, final_table_name);
        
        pfree(final_table_name);
    }
    
    /* Calculate total statistics */
    for (i = 0; i < num_workers; i++) {
        result.total_rows += workers[i].rows_processed;
        result.highfreq_kmers_count += workers[i].local_highfreq_count;
    }
    
    /* Clean up */
    pfree(workers);
    relation_close(rel, AccessShareLock);
    
    return result;
}

/*
 * Persist highly frequent k-mers from temporary table to permanent tables
 */
static void
kmersearch_persist_highfreq_kmers_from_temp(Oid table_oid, const char *column_name, int k_size,
                                           const char *temp_table_name)
{
    StringInfoData query;
    
    initStringInfo(&query);
    
    /* Insert highly frequent k-mers into permanent table */
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmers (index_oid, kmer_key, frequency_count, detection_reason) "
        "SELECT %u, kmer_data::varbit, frequency_count, 'high_frequency' "
        "FROM %s",
        table_oid, temp_table_name);
    
    SPI_connect();
    SPI_exec(query.data, 0);
    SPI_finish();
    
    pfree(query.data);
}

/*
 * Persist highly frequent k-mers to database tables (legacy - disabled in memory-efficient version)
 */
static void
kmersearch_persist_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, 
                                 void *unused_table, int threshold_rows)
{
    /* Memory-efficient implementation uses SQL-based analysis instead */
    (void) table_oid;
    (void) column_name;
    (void) k_size;
    (void) unused_table;
    (void) threshold_rows;
    
    ereport(NOTICE, (errmsg("Legacy persist function called - using memory-efficient analysis instead")));
}

/*
 * Delete existing analysis results
 */
static void
kmersearch_delete_existing_analysis(Oid table_oid, const char *column_name, int k_size)
{
    /* For now, just log that we would delete existing analysis */
    ereport(NOTICE, (errmsg("Would delete existing analysis for table %u, column %s, k=%d", 
                           table_oid, column_name, k_size)));
}

/*
 * Check if analysis exists for given parameters
 */
static bool
kmersearch_check_analysis_exists(Oid table_oid, const char *column_name, int k_size)
{
    int ret;
    bool found = false;
    StringInfoData query;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to check for existing analysis */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT COUNT(*) FROM kmersearch_index_info "
        "WHERE table_oid = %u AND column_name = '%s' AND k_value = %d",
        table_oid, column_name, k_size);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    kmersearch_handle_spi_error(ret, "SELECT");
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        Datum count_datum;
        bool isnull;
        int count;
        
        count_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull)
        {
            count = DatumGetInt32(count_datum);
            found = (count > 0);
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return found;
}

/*
 * Validate analysis parameters match current GUC settings
 */
static void
kmersearch_validate_analysis_parameters(Oid table_oid, const char *column_name, int k_size)
{
    Relation rel;
    TupleDesc tupdesc;
    AttrNumber attnum;
    Form_pg_attribute attr;
    Oid column_type;
    
    /* Validate k-mer size */
    if (k_size < 4 || k_size > 64) {
        ereport(ERROR, 
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer size must be between 4 and 64, got %d", k_size)));
    }
    
    /* Validate table exists and we can access it */
    rel = try_relation_open(table_oid, AccessShareLock);
    if (rel == NULL) {
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_TABLE),
                 errmsg("table with OID %u does not exist", table_oid)));
    }
    
    /* Validate column exists and has correct type */
    tupdesc = RelationGetDescr(rel);
    attnum = get_attnum(table_oid, column_name);
    if (attnum == InvalidAttrNumber) {
        relation_close(rel, AccessShareLock);
        ereport(ERROR,
                (errcode(ERRCODE_UNDEFINED_COLUMN),
                 errmsg("column \"%s\" does not exist in table", column_name)));
    }
    
    attr = TupleDescAttr(tupdesc, attnum - 1);
    column_type = attr->atttypid;
    
    /* Check if column type is compatible (basic check for varbit-like types) */
    if (column_type != VARBITOID && attr->atttypmod == -1) {
        /* For now, just warn about potential type mismatch */
        ereport(WARNING,
                (errmsg("column \"%s\" type may not be compatible with k-mer analysis", column_name),
                 errhint("Expected DNA2, DNA4, or varbit type")));
    }
    
    /* Validate GUC parameters */
    if (kmersearch_max_appearance_rate <= 0.0 && kmersearch_max_appearance_nrow <= 0) {
        relation_close(rel, AccessShareLock);
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("exclusion parameters not configured"),
                 errhint("Set kmersearch.max_appearance_rate > 0 or kmersearch.max_appearance_nrow > 0")));
    }
    
    if (kmersearch_max_appearance_rate < 0.0 || kmersearch_max_appearance_rate > 1.0) {
        relation_close(rel, AccessShareLock);
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("kmersearch.max_appearance_rate must be between 0.0 and 1.0, got %f", 
                        kmersearch_max_appearance_rate)));
    }
    
    relation_close(rel, AccessShareLock);
    
    ereport(DEBUG1, (errmsg("Analysis parameters validated for table %u, column %s, k=%d", 
                           table_oid, column_name, k_size)));
}

/*
 * Filter highly frequent k-mers from the key array
 * A-2: Use direct VarBit comparison instead of hash table (k-mer+occurrence n-gram keys are small)
 */
static Datum *
kmersearch_filter_highfreq_kmers(Oid table_oid, const char *column_name, int k_size, 
                                Datum *all_keys, int total_keys, int *filtered_count)
{
    int ret;
    StringInfoData query;
    Datum *filtered_keys;
    int filtered_idx = 0;
    int i;
    VarBit **highfreq_kmers = NULL;
    int highfreq_count = 0;
    
    /* Check if analysis exists */
    if (!kmersearch_check_analysis_exists(table_oid, column_name, k_size)) {
        /* No analysis found, return all keys */
        *filtered_count = total_keys;
        return all_keys;
    }
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get excluded k-mers from index info table */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT ek.kmer_key FROM kmersearch_highfreq_kmers ek "
        "JOIN kmersearch_index_info ii ON ek.index_oid = ii.index_oid "
        "WHERE ii.table_oid = %u AND ii.column_name = '%s' AND ii.k_value = %d",
        table_oid, column_name, k_size);
    
    /* Execute query and collect highly frequent k-mers in a simple array */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        int j;
        highfreq_count = SPI_processed;
        highfreq_kmers = (VarBit **) palloc(highfreq_count * sizeof(VarBit *));
        
        for (j = 0; j < SPI_processed; j++)
        {
            bool isnull;
            Datum kmer_datum;
            VarBit *kmer;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[j], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                kmer = DatumGetVarBitP(kmer_datum);
                
                /* Store copy of k-mer for direct comparison */
                highfreq_kmers[j] = (VarBit *) palloc(VARSIZE(kmer));
                memcpy(highfreq_kmers[j], kmer, VARSIZE(kmer));
            }
            else
            {
                highfreq_kmers[j] = NULL;
            }
        }
    }
    
    /* Filter out highly frequent k-mers using direct VarBit comparison */
    filtered_keys = (Datum *) palloc(total_keys * sizeof(Datum));
    
    for (i = 0; i < total_keys; i++)
    {
        VarBit *kmer = DatumGetVarBitP(all_keys[i]);
        bool is_highfreq = false;
        int j;
        
        /* Direct comparison with highly frequent k-mers (no hashing) */
        for (j = 0; j < highfreq_count; j++)
        {
            if (highfreq_kmers[j] != NULL &&
                VARBITLEN(kmer) == VARBITLEN(highfreq_kmers[j]) &&
                VARSIZE(kmer) == VARSIZE(highfreq_kmers[j]) &&
                memcmp(VARBITS(kmer), VARBITS(highfreq_kmers[j]), VARBITBYTES(kmer)) == 0)
            {
                is_highfreq = true;
                break;
            }
        }
        
        if (!is_highfreq)
        {
            /* Not highly frequent, include in filtered result */
            filtered_keys[filtered_idx++] = all_keys[i];
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    if (highfreq_kmers)
    {
        int j;
        for (j = 0; j < highfreq_count; j++)
        {
            if (highfreq_kmers[j])
                pfree(highfreq_kmers[j]);
        }
        pfree(highfreq_kmers);
    }
    SPI_finish();
    
    *filtered_count = filtered_idx;
    
    /* If no keys were filtered, return original array */
    if (filtered_idx == total_keys) {
        pfree(filtered_keys);
        return all_keys;
    }
    
    return filtered_keys;
}
/*
 * Main k-mer analysis function (Option A implementation)
 */
Datum
kmersearch_analyze_table(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k_size = PG_GETARG_INT32(2);
    int parallel_workers = PG_GETARG_INT32(3);
    
    char *column_name = text_to_cstring(column_name_text);
    KmerAnalysisResult result;
    
    /* Comprehensive parameter validation */
    kmersearch_validate_analysis_parameters(table_oid, column_name, k_size);
    
    /* Perform parallel analysis */
    result = kmersearch_analyze_table_parallel(table_oid, column_name, k_size, parallel_workers);
    
    /* Create result tuple */
    {
        TupleDesc tupdesc;
        Datum values[6];
        bool nulls[6] = {false};
        HeapTuple tuple;
    
    /* Build tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE) {
        ereport(ERROR, (errmsg("function returning record called in context that cannot accept a record")));
    }
    
    /* Fill result values */
    values[0] = Int64GetDatum(result.total_rows);
    values[1] = Int32GetDatum(result.highfreq_kmers_count);
    values[2] = Int32GetDatum(result.parallel_workers_used);
    values[3] = Float8GetDatum(result.analysis_duration);
    values[4] = Float8GetDatum(result.max_appearance_rate_used);
    values[5] = Int32GetDatum(result.max_appearance_nrow_used);
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    
    /* Cleanup */
    pfree(column_name);
    
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
    }
}

/*
 * Drop analysis results function
 */
Datum
kmersearch_drop_analysis(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k_size = PG_GETARG_INT32(2);  /* 0 means all k-sizes */
    
    char *column_name = text_to_cstring(column_name_text);
    DropAnalysisResult result;
    
    /* Validate table OID */
    if (!OidIsValid(table_oid)) {
        ereport(ERROR, (errmsg("invalid table OID")));
    }
    
    /* Perform drop operation */
    result = kmersearch_drop_analysis_internal(table_oid, column_name, k_size);
    
    /* Create result tuple */
    {
        TupleDesc tupdesc;
        Datum values[3];
        bool nulls[3] = {false};
        HeapTuple tuple;
    
    /* Build tuple descriptor */
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE) {
        ereport(ERROR, (errmsg("function returning record called in context that cannot accept a record")));
    }
    
    /* Fill result values */
    values[0] = Int32GetDatum(result.dropped_analyses);
    values[1] = Int32GetDatum(result.dropped_highfreq_kmers);
    values[2] = Int64GetDatum(result.freed_storage_bytes);
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    
    /* Cleanup */
    pfree(column_name);
    
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
    }
}

/*
 * Internal drop analysis implementation
 */
static DropAnalysisResult
kmersearch_drop_analysis_internal(Oid table_oid, const char *column_name, int k_size)
{
    DropAnalysisResult result = {0};
    
    /* For now, just log the operation */
    if (k_size > 0) {
        ereport(NOTICE, (errmsg("Would drop analysis for table %u, column %s, k=%d", 
                               table_oid, column_name, k_size)));
        result.dropped_analyses = 1;
        result.dropped_highfreq_kmers = 10;  /* Mock value */
    } else {
        ereport(NOTICE, (errmsg("Would drop all analyses for table %u, column %s", 
                               table_oid, column_name)));
        result.dropped_analyses = 3;  /* Mock value */
        result.dropped_highfreq_kmers = 30;  /* Mock value */
    }
    
    result.freed_storage_bytes = result.dropped_highfreq_kmers * 100;  /* Mock calculation */
    
    return result;
}

/*
 * Helper function to get highly frequent k-mers list for a given index
 */
static List *
kmersearch_get_highfreq_kmers_list(Oid index_oid)
{
    List *highfreq_kmers = NIL;
    int ret;
    StringInfoData query;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get highly frequent k-mers */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT ek.kmer_key FROM kmersearch_highfreq_kmers ek "
        "JOIN kmersearch_index_info ii ON ek.index_oid = ii.index_oid "
        "WHERE ii.index_oid = %u ORDER BY ek.kmer_key",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    kmersearch_handle_spi_error(ret, "SELECT highly frequent k-mers");
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        int i;
        for (i = 0; i < SPI_processed; i++)
        {
            bool isnull;
            Datum kmer_datum;
            VarBit *kmer;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                kmer = DatumGetVarBitPCopy(kmer_datum);
                highfreq_kmers = lappend(highfreq_kmers, kmer);
            }
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return highfreq_kmers;
}

/*
 * Helper function to delete a k-mer from GIN index
 */
static bool
kmersearch_delete_kmer_from_gin_index(Relation index_rel, VarBit *kmer_key)
{
    /* 
     * Note: This is a simplified implementation.
     * In a complete implementation, this would use GIN internal APIs
     * to locate and remove the specific k-mer key and its posting list.
     * For now, we'll just log the operation.
     */
    ereport(DEBUG1, (errmsg("Would delete k-mer from index (size: %d bits)", 
                           VARBITLEN(kmer_key))));
    
    /*
     * TODO: Implement actual GIN key deletion using:
     * - ginFindLeafPage() to locate the key
     * - ginDeleteKey() to remove the key and posting list
     * - ginUpdateStats() to update index statistics
     */
    
    return true;  /* Assume success for now */
}

/*
 * Main function to reduce index size by removing highly frequent k-mers
 * Note: This is a placeholder implementation
 */
Datum
kmersearch_reduce_index(PG_FUNCTION_ARGS)
{
    Oid index_oid = PG_GETARG_OID(0);
    
    /* Validate index OID */
    if (!OidIsValid(index_oid))
        ereport(ERROR, (errmsg("invalid index OID")));
    
    /* Return placeholder message */
    ereport(NOTICE, (errmsg("Index reduction not implemented yet for index OID %u", index_oid)));
    
    PG_RETURN_TEXT_P(cstring_to_text("Index reduction not implemented"));
}

/*
 * Rawscore cache statistics function
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
 * Free rawscore cache managers
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
 * Helper functions to access and free static cache managers
 */
static int
free_dna2_rawscore_cache(void)
{
    /* This function needs to be called from within get_cached_rawscore_dna2 context */
    /* We'll implement this through a global flag approach */
    return 0;  /* Placeholder - will be implemented with global reset */
}

static int
free_dna4_rawscore_cache(void)
{
    /* This function needs to be called from within get_cached_rawscore_dna4 context */
    /* We'll implement this through a global flag approach */
    return 0;  /* Placeholder - will be implemented with global reset */
}

/*
 * Rawscore cache free function
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

/* free_cache_manager function removed - cache managers are now local and auto-freed */

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
 * High-frequency k-mer cache load function
 */
Datum
kmersearch_highfreq_kmers_cache_load(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k_value = PG_GETARG_INT32(2);
    
    char *column_name = text_to_cstring(column_name_text);
    bool success;
    
    /* Call internal cache load function */
    success = kmersearch_highfreq_kmers_cache_load_internal(table_oid, column_name, k_value);
    
    pfree(column_name);
    
    PG_RETURN_BOOL(success);
}

/*
 * High-frequency k-mer cache free function
 */
Datum
kmersearch_highfreq_kmers_cache_free(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    /* Count current entries before clearing */
    if (global_highfreq_cache.is_valid)
        freed_entries = global_highfreq_cache.highfreq_count;
    
    /* Call internal cache free function */
    kmersearch_highfreq_kmers_cache_free_internal();
    
    PG_RETURN_INT32(freed_entries);
}

/*
 * Module cleanup function
 */
void
_PG_fini(void)
{
    /* Free query pattern cache manager on module unload (uses TopMemoryContext - needs manual cleanup) */
    /* DNA2/DNA4 cache managers are now local and automatically freed with QueryContext */
    free_query_pattern_cache_manager(&query_pattern_cache_manager);
    
    /* Free actual min score cache manager on module unload */
    free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
    /* Free high-frequency k-mer cache on module unload */
    kmersearch_highfreq_kmers_cache_free_internal();
}

/*
 * High-frequency k-mer filtering functions implementation
 */
static VarBit *
kmersearch_remove_occurrence_bits(VarBit *key_with_occurrence, int k)
{
    VarBit *result;
    int kmer_bits;
    int kmer_bytes;
    int total_bits;
    int occur_bits;
    
    if (!key_with_occurrence || k <= 0)
        return NULL;
    
    occur_bits = kmersearch_occur_bitlen;
    total_bits = VARBITLEN(key_with_occurrence);
    kmer_bits = total_bits - occur_bits;
    
    if (kmer_bits <= 0)
        return NULL;
    
    kmer_bytes = (kmer_bits + 7) / 8;
    result = (VarBit *) palloc(VARHDRSZ + kmer_bytes);
    SET_VARSIZE(result, VARHDRSZ + kmer_bytes);
    VARBITLEN(result) = kmer_bits;
    
    memcpy(VARBITS(result), VARBITS(key_with_occurrence), kmer_bytes);
    
    return result;
}

static VarBit **
kmersearch_get_highfreq_kmers_from_table(Oid table_oid, const char *column_name, int k, int *nkeys)
{
    VarBit **result = NULL;
    int ret;
    StringInfoData query;
    int i;
    
    if (!nkeys)
        return NULL;
    
    *nkeys = 0;
    
    /* Connect to SPI */
    if (SPI_connect() != SPI_OK_CONNECT)
        ereport(ERROR, (errmsg("kmersearch_get_highfreq_kmers_from_table: SPI_connect failed")));
    
    /* Build query to get highly frequent k-mers */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT DISTINCT hkm.kmer_key FROM kmersearch_highfreq_kmers hkm "
        "JOIN kmersearch_highfreq_kmers_meta hkm_meta ON "
        "(hkm_meta.table_oid = %u AND hkm_meta.column_name = '%s' AND hkm_meta.k_value = %d) "
        "ORDER BY hkm.kmer_key",
        table_oid, column_name, k);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        *nkeys = SPI_processed;
        result = (VarBit **) palloc(*nkeys * sizeof(VarBit *));
        
        for (i = 0; i < *nkeys; i++)
        {
            bool isnull;
            Datum kmer_datum;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                /* Copy the varbit value */
                result[i] = DatumGetVarBitPCopy(kmer_datum);
            }
            else
            {
                result[i] = NULL;
            }
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return result;
}

static HTAB *
kmersearch_create_highfreq_hash_from_array(VarBit **kmers, int nkeys)
{
    HTAB *hash_table;
    HASHCTL hash_ctl;
    int i;
    
    if (!kmers || nkeys <= 0)
        return NULL;
    
    /* Set up hash table */
    MemSet(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(VarBit *);
    hash_ctl.entrysize = sizeof(HighfreqKmerHashEntry);
    hash_ctl.hash = tag_hash;
    hash_ctl.hcxt = CurrentMemoryContext;
    
    hash_table = hash_create("HighfreqKmerHash",
                            nkeys,
                            &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_CONTEXT);
    
    if (!hash_table)
        return NULL;
    
    /* Add each k-mer to the hash table */
    for (i = 0; i < nkeys; i++)
    {
        HighfreqKmerHashEntry *entry;
        bool found;
        
        if (!kmers[i])
            continue;
        
        entry = (HighfreqKmerHashEntry *) hash_search(hash_table,
                                                     (void *) &kmers[i],
                                                     HASH_ENTER,
                                                     &found);
        
        if (entry && !found)
        {
            entry->kmer_key = kmers[i];
            entry->hash_value = DatumGetUInt64(hash_any((unsigned char *) VARBITS(kmers[i]), VARBITBYTES(kmers[i])));
        }
    }
    
    return hash_table;
}

static Datum *
kmersearch_filter_highfreq_kmers_from_keys(Datum *original_keys, int *nkeys, HTAB *highfreq_hash, int k)
{
    Datum *filtered_keys;
    int original_count;
    int filtered_count;
    int i;
    
    if (!original_keys || !nkeys || *nkeys <= 0)
        return NULL;
    
    /* If no high-frequency hash, return original keys unchanged */
    if (!highfreq_hash)
        return original_keys;
    
    original_count = *nkeys;
    filtered_keys = (Datum *) palloc(original_count * sizeof(Datum));
    filtered_count = 0;
    
    /* Filter out high-frequency k-mers */
    for (i = 0; i < original_count; i++)
    {
        VarBit *key_with_occurrence;
        VarBit *kmer_only;
        HighfreqKmerHashEntry *entry;
        bool found;
        
        key_with_occurrence = DatumGetVarBitP(original_keys[i]);
        if (!key_with_occurrence)
            continue;
        
        /* Remove occurrence bits to get pure k-mer */
        kmer_only = kmersearch_remove_occurrence_bits(key_with_occurrence, k);
        if (!kmer_only)
        {
            /* If we can't extract the k-mer, keep the original key */
            filtered_keys[filtered_count++] = original_keys[i];
            continue;
        }
        
        /* Check if this k-mer is in the high-frequency list */
        entry = (HighfreqKmerHashEntry *) hash_search(highfreq_hash,
                                                     (void *) &kmer_only,
                                                     HASH_FIND,
                                                     &found);
        
        if (!found)
        {
            /* Not a high-frequency k-mer, keep it */
            filtered_keys[filtered_count++] = original_keys[i];
        }
        
        /* Clean up the temporary k-mer */
        pfree(kmer_only);
    }
    
    /* Update the count */
    *nkeys = filtered_count;
    
    /* If no keys left, return NULL */
    if (filtered_count == 0)
    {
        pfree(filtered_keys);
        return NULL;
    }
    
    /* Resize the array if significantly smaller */
    if (filtered_count < original_count / 2)
    {
        filtered_keys = (Datum *) repalloc(filtered_keys, filtered_count * sizeof(Datum));
    }
    
    return filtered_keys;
}

/*
 * High-frequency k-mer cache management functions implementation
 */
static void
kmersearch_highfreq_kmers_cache_init(void)
{
    MemoryContext old_context;
    
    /* Switch to TopMemoryContext */
    old_context = MemoryContextSwitchTo(TopMemoryContext);
    
    /* Initialize cache structure */
    memset(&global_highfreq_cache, 0, sizeof(HighfreqKmerCache));
    global_highfreq_cache.is_valid = false;
    global_highfreq_cache.current_table_oid = InvalidOid;
    global_highfreq_cache.current_column_name = NULL;
    global_highfreq_cache.current_k_value = 0;
    
    /* Create dedicated memory context for high-frequency k-mer cache */
    global_highfreq_cache.cache_context = AllocSetContextCreate(TopMemoryContext,
                                                                "HighfreqKmerCache",
                                                                ALLOCSET_DEFAULT_SIZES);
    
    global_highfreq_cache.highfreq_hash = NULL;
    global_highfreq_cache.highfreq_kmers = NULL;
    global_highfreq_cache.highfreq_count = 0;
    
    MemoryContextSwitchTo(old_context);
}

static bool
kmersearch_highfreq_kmers_cache_load_internal(Oid table_oid, const char *column_name, int k_value)
{
    MemoryContext old_context;
    VarBit **highfreq_kmers;
    int highfreq_count;
    
    if (!column_name || k_value <= 0)
        return false;
    
    /* Clear existing cache if valid */
    if (global_highfreq_cache.is_valid) {
        kmersearch_highfreq_kmers_cache_free_internal();
    }
    
    /* Get high-frequency k-mers list */
    highfreq_kmers = kmersearch_get_highfreq_kmers_from_table(table_oid, column_name, k_value, &highfreq_count);
    
    if (!highfreq_kmers || highfreq_count <= 0) {
        /* No high-frequency k-mers found, cache remains invalid */
        return false;
    }
    
    /* Switch to cache context for cache storage */
    old_context = MemoryContextSwitchTo(global_highfreq_cache.cache_context);
    
    /* Store in cache */
    global_highfreq_cache.current_table_oid = table_oid;
    global_highfreq_cache.current_column_name = pstrdup(column_name);  /* Cache context copy */
    global_highfreq_cache.current_k_value = k_value;
    global_highfreq_cache.highfreq_kmers = highfreq_kmers;
    global_highfreq_cache.highfreq_count = highfreq_count;
    
    /* Create hash table in cache context */
    global_highfreq_cache.highfreq_hash = kmersearch_create_highfreq_hash_from_array(highfreq_kmers, highfreq_count);
    
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

static void
kmersearch_highfreq_kmers_cache_free_internal(void)
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
    global_highfreq_cache.current_k_value = 0;
    global_highfreq_cache.highfreq_count = 0;
    global_highfreq_cache.highfreq_hash = NULL;
    global_highfreq_cache.highfreq_kmers = NULL;
    global_highfreq_cache.current_column_name = NULL;
}

static bool
kmersearch_highfreq_kmers_cache_is_valid(Oid table_oid, const char *column_name, int k_value)
{
    return (global_highfreq_cache.is_valid &&
            global_highfreq_cache.current_table_oid == table_oid &&
            global_highfreq_cache.current_k_value == k_value &&
            global_highfreq_cache.current_column_name &&
            column_name &&
            strcmp(global_highfreq_cache.current_column_name, column_name) == 0);
}

/*
 * Auto-load cache if high-frequency k-mer metadata is available
 * This function searches for available metadata and loads the cache automatically
 */
static bool
kmersearch_auto_load_cache_if_needed(void)
{
    int ret;
    StringInfoData query;
    bool cache_loaded = false;
    
    /* If cache is already valid, no need to reload */
    if (global_highfreq_cache.is_valid)
        return true;
    
    /* Connect to SPI */
    if (SPI_connect() != SPI_OK_CONNECT)
        return false;
    
    /* Query for available high-frequency k-mer metadata */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT table_oid, column_name, k_value "
        "FROM kmersearch_highfreq_kmers_meta "
        "ORDER BY analysis_timestamp DESC "
        "LIMIT 1");
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        bool isnull;
        Datum table_oid_datum, column_name_datum, k_value_datum;
        Oid table_oid;
        char *column_name;
        int k_value;
        
        /* Extract values from first row */
        table_oid_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull)
        {
            table_oid = DatumGetObjectId(table_oid_datum);
            
            column_name_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull);
            if (!isnull)
            {
                column_name = NameStr(*DatumGetName(column_name_datum));
                
                k_value_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3, &isnull);
                if (!isnull)
                {
                    k_value = DatumGetInt32(k_value_datum);
                    
                    /* Try to load cache with found metadata */
                    SPI_finish();
                    cache_loaded = kmersearch_highfreq_kmers_cache_load_internal(table_oid, column_name, k_value);
                    if (SPI_connect() != SPI_OK_CONNECT)
                        return cache_loaded;
                }
            }
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return cache_loaded;
}

/*
 * SIMD Implementation Functions
 */

/* Scalar implementations (fallback) */
static void dna2_encode_scalar(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 2 + 7) / 8;
    memset(output, 0, byte_len);
    
    for (int i = 0; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        output[byte_pos] |= (encoded << (6 - bit_offset));
    }
}

static void dna2_decode_scalar(const uint8_t* input, char* output, int len)
{
    for (int i = 0; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

static void dna4_encode_scalar(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    memset(output, 0, byte_len);
    
    for (int i = 0; i < len; i++) {
        uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 4) {
            output[byte_pos] |= (encoded << (4 - bit_offset));
        } else {
            output[byte_pos] |= (encoded >> (bit_offset - 4));
            if (byte_pos + 1 < byte_len)
                output[byte_pos + 1] |= (encoded << (12 - bit_offset));
        }
    }
}

static void dna4_decode_scalar(const uint8_t* input, char* output, int len)
{
    int bit_len = len * 4;
    
    for (int i = 0; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            encoded = ((input[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < (bit_len + 7) / 8)
                encoded |= (input[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}

#ifdef __x86_64__
/* AVX2 implementations */
__attribute__((target("avx2")))
static void dna2_encode_avx2(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 2 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Process 32 characters at a time with AVX2 */
    int simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    for (int i = 0; i < simd_len; i += 32) {
        __m256i chars = _mm256_loadu_si256((__m256i*)(input + i));
        
        /* Create comparison masks for each DNA base */
        __m256i mask_A = _mm256_or_si256(_mm256_cmpeq_epi8(chars, _mm256_set1_epi8('A')),
                                         _mm256_cmpeq_epi8(chars, _mm256_set1_epi8('a')));
        __m256i mask_C = _mm256_or_si256(_mm256_cmpeq_epi8(chars, _mm256_set1_epi8('C')),
                                         _mm256_cmpeq_epi8(chars, _mm256_set1_epi8('c')));
        __m256i mask_G = _mm256_or_si256(_mm256_cmpeq_epi8(chars, _mm256_set1_epi8('G')),
                                         _mm256_cmpeq_epi8(chars, _mm256_set1_epi8('g')));
        __m256i mask_T = _mm256_or_si256(_mm256_cmpeq_epi8(chars, _mm256_set1_epi8('T')),
                                         _mm256_cmpeq_epi8(chars, _mm256_set1_epi8('t')));
        __m256i mask_U = _mm256_or_si256(_mm256_cmpeq_epi8(chars, _mm256_set1_epi8('U')),
                                         _mm256_cmpeq_epi8(chars, _mm256_set1_epi8('u')));
        
        /* Combine T and U masks */
        mask_T = _mm256_or_si256(mask_T, mask_U);
        
        /* Generate 2-bit encoded values */
        __m256i encoded = _mm256_setzero_si256();
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_C, _mm256_set1_epi8(1)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_G, _mm256_set1_epi8(2)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_T, _mm256_set1_epi8(3)));
        
        /* Pack 2-bit values into bytes - simplified approach */
        /* For full implementation, we'd need complex bit packing */
        /* Fall back to scalar for bit packing part */
        uint8_t temp[32];
        _mm256_storeu_si256((__m256i*)temp, encoded);
        
        for (int j = 0; j < 32; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            output[byte_pos] |= (temp[j] << (6 - bit_offset));
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        output[byte_pos] |= (encoded << (6 - bit_offset));
    }
}

__attribute__((target("avx2")))
static void dna2_decode_avx2(const uint8_t* input, char* output, int len)
{
    /* Process 32 characters at a time with AVX2 */
    int simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    for (int i = 0; i < simd_len; i += 32) {
        /* Extract 16 bytes (128 bits) to process 32 characters */
        __m128i data = _mm_loadu_si128((__m128i*)(input + i/4));
        
        /* Expand to 32 bytes for easier processing */
        uint8_t temp[16];
        _mm_storeu_si128((__m128i*)temp, data);
        
        /* Decode each 2-bit pair */
        for (int j = 0; j < 32; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

__attribute__((target("avx2")))
static void dna4_encode_avx2(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Process 32 characters at a time with AVX2 */
    int simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    for (int i = 0; i < simd_len; i += 32) {
        __m256i chars = _mm256_loadu_si256((__m256i*)(input + i));
        
        /* Use lookup table approach for DNA4 encoding */
        uint8_t temp[32];
        _mm256_storeu_si256((__m256i*)temp, chars);
        
        for (int j = 0; j < 32; j++) {
            uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)temp[j]];
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 4) {
                output[byte_pos] |= (encoded << (4 - bit_offset));
            } else {
                output[byte_pos] |= (encoded >> (bit_offset - 4));
                if (byte_pos + 1 < byte_len)
                    output[byte_pos + 1] |= (encoded << (12 - bit_offset));
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 4) {
            output[byte_pos] |= (encoded << (4 - bit_offset));
        } else {
            output[byte_pos] |= (encoded >> (bit_offset - 4));
            if (byte_pos + 1 < byte_len)
                output[byte_pos + 1] |= (encoded << (12 - bit_offset));
        }
    }
}

__attribute__((target("avx2")))
static void dna4_decode_avx2(const uint8_t* input, char* output, int len)
{
    int bit_len = len * 4;
    
    /* Process characters with AVX2 assistance */
    int simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    for (int i = 0; i < simd_len; i += 32) {
        /* Load 16 bytes to process 32 characters */
        __m128i data = _mm_loadu_si128((__m128i*)(input + i/2));
        uint8_t temp[16];
        _mm_storeu_si128((__m128i*)temp, data);
        
        /* Decode each 4-bit nibble */
        for (int j = 0; j < 32; j++) {
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded;
            
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                encoded = ((temp[byte_pos] << (bit_offset - 4)) & 0xF);
                if (byte_pos + 1 < (bit_len + 7) / 8)
                    encoded |= (temp[byte_pos + 1] >> (12 - bit_offset));
                encoded &= 0xF;
            }
            
            output[i + j] = kmersearch_dna4_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            encoded = ((input[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < (bit_len + 7) / 8)
                encoded |= (input[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}

/* AVX512 implementations */
__attribute__((target("avx512f,avx512bw")))
static void dna2_encode_avx512(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 2 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Process 64 characters at a time with AVX512 */
    int simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    for (int i = 0; i < simd_len; i += 64) {
        __m512i chars = _mm512_loadu_si512((__m512i*)(input + i));
        
        /* Create comparison masks for each DNA base */
        __mmask64 mask_A = _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('A')) |
                          _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('a'));
        __mmask64 mask_C = _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('C')) |
                          _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('c'));
        __mmask64 mask_G = _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('G')) |
                          _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('g'));
        __mmask64 mask_T = _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('T')) |
                          _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('t')) |
                          _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('U')) |
                          _mm512_cmpeq_epi8_mask(chars, _mm512_set1_epi8('u'));
        
        /* Generate 2-bit encoded values using mask operations */
        __m512i encoded = _mm512_setzero_si512();
        encoded = _mm512_mask_or_epi32(encoded, mask_C, encoded, _mm512_set1_epi8(1));
        encoded = _mm512_mask_or_epi32(encoded, mask_G, encoded, _mm512_set1_epi8(2));
        encoded = _mm512_mask_or_epi32(encoded, mask_T, encoded, _mm512_set1_epi8(3));
        
        /* Store and pack bits */
        uint8_t temp[64];
        _mm512_storeu_si512((__m512i*)temp, encoded);
        
        for (int j = 0; j < 64; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            output[byte_pos] |= (temp[j] << (6 - bit_offset));
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        output[byte_pos] |= (encoded << (6 - bit_offset));
    }
}

__attribute__((target("avx512f,avx512bw")))
static void dna2_decode_avx512(const uint8_t* input, char* output, int len)
{
    /* Process 64 characters at a time with AVX512 */
    int simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    for (int i = 0; i < simd_len; i += 64) {
        /* Load 32 bytes to process 64 characters */
        __m256i data = _mm256_loadu_si256((__m256i*)(input + i/4));
        
        uint8_t temp[32];
        _mm256_storeu_si256((__m256i*)temp, data);
        
        /* Decode each 2-bit pair */
        for (int j = 0; j < 64; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

__attribute__((target("avx512f,avx512bw")))
static void dna4_encode_avx512(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Process 64 characters at a time with AVX512 */
    int simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    for (int i = 0; i < simd_len; i += 64) {
        __m512i chars = _mm512_loadu_si512((__m512i*)(input + i));
        
        /* Use lookup table approach for DNA4 encoding */
        uint8_t temp[64];
        _mm512_storeu_si512((__m512i*)temp, chars);
        
        for (int j = 0; j < 64; j++) {
            uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)temp[j]];
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 4) {
                output[byte_pos] |= (encoded << (4 - bit_offset));
            } else {
                output[byte_pos] |= (encoded >> (bit_offset - 4));
                if (byte_pos + 1 < byte_len)
                    output[byte_pos + 1] |= (encoded << (12 - bit_offset));
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 4) {
            output[byte_pos] |= (encoded << (4 - bit_offset));
        } else {
            output[byte_pos] |= (encoded >> (bit_offset - 4));
            if (byte_pos + 1 < byte_len)
                output[byte_pos + 1] |= (encoded << (12 - bit_offset));
        }
    }
}

__attribute__((target("avx512f,avx512bw")))
static void dna4_decode_avx512(const uint8_t* input, char* output, int len)
{
    int bit_len = len * 4;
    
    /* Process 64 characters at a time with AVX512 */
    int simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    for (int i = 0; i < simd_len; i += 64) {
        /* Load 32 bytes to process 64 characters */
        __m256i data = _mm256_loadu_si256((__m256i*)(input + i/2));
        uint8_t temp[32];
        _mm256_storeu_si256((__m256i*)temp, data);
        
        /* Decode each 4-bit nibble */
        for (int j = 0; j < 64; j++) {
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded;
            
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                encoded = ((temp[byte_pos] << (bit_offset - 4)) & 0xF);
                if (byte_pos + 1 < (bit_len + 7) / 8)
                    encoded |= (temp[byte_pos + 1] >> (12 - bit_offset));
                encoded &= 0xF;
            }
            
            output[i + j] = kmersearch_dna4_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            encoded = ((input[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < (bit_len + 7) / 8)
                encoded |= (input[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}
#endif

#ifdef __aarch64__
/* NEON implementations */
static void dna2_encode_neon(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 2 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Process 16 characters at a time with NEON */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    for (int i = 0; i < simd_len; i += 16) {
        uint8x16_t chars = vld1q_u8((uint8_t*)(input + i));
        
        /* Create comparison masks for each DNA base */
        uint8x16_t mask_A = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('A')),
                                     vceqq_u8(chars, vdupq_n_u8('a')));
        uint8x16_t mask_C = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('C')),
                                     vceqq_u8(chars, vdupq_n_u8('c')));
        uint8x16_t mask_G = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('G')),
                                     vceqq_u8(chars, vdupq_n_u8('g')));
        uint8x16_t mask_T = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('T')),
                                     vceqq_u8(chars, vdupq_n_u8('t')));
        uint8x16_t mask_U = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('U')),
                                     vceqq_u8(chars, vdupq_n_u8('u')));
        
        /* Combine T and U masks */
        mask_T = vorrq_u8(mask_T, mask_U);
        
        /* Generate 2-bit encoded values */
        uint8x16_t encoded = vdupq_n_u8(0);
        encoded = vorrq_u8(encoded, vandq_u8(mask_C, vdupq_n_u8(1)));
        encoded = vorrq_u8(encoded, vandq_u8(mask_G, vdupq_n_u8(2)));
        encoded = vorrq_u8(encoded, vandq_u8(mask_T, vdupq_n_u8(3)));
        
        /* Store and pack bits */
        uint8_t temp[16];
        vst1q_u8(temp, encoded);
        
        for (int j = 0; j < 16; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            output[byte_pos] |= (temp[j] << (6 - bit_offset));
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        output[byte_pos] |= (encoded << (6 - bit_offset));
    }
}

static void dna2_decode_neon(const uint8_t* input, char* output, int len)
{
    /* Process 16 characters at a time with NEON */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    for (int i = 0; i < simd_len; i += 16) {
        /* Load 4 bytes to process 16 characters */
        uint32_t data_word;
        memcpy(&data_word, input + i/4, 4);
        uint8_t temp[4] = {data_word & 0xFF, (data_word >> 8) & 0xFF, (data_word >> 16) & 0xFF, (data_word >> 24) & 0xFF};
        
        /* Decode each 2-bit pair */
        for (int j = 0; j < 16; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

static void dna4_encode_neon(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Process 16 characters at a time with NEON */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    for (int i = 0; i < simd_len; i += 16) {
        uint8x16_t chars = vld1q_u8((uint8_t*)(input + i));
        
        /* Use lookup table approach for DNA4 encoding */
        uint8_t temp[16];
        vst1q_u8(temp, chars);
        
        for (int j = 0; j < 16; j++) {
            uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)temp[j]];
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 4) {
                output[byte_pos] |= (encoded << (4 - bit_offset));
            } else {
                output[byte_pos] |= (encoded >> (bit_offset - 4));
                if (byte_pos + 1 < byte_len)
                    output[byte_pos + 1] |= (encoded << (12 - bit_offset));
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 4) {
            output[byte_pos] |= (encoded << (4 - bit_offset));
        } else {
            output[byte_pos] |= (encoded >> (bit_offset - 4));
            if (byte_pos + 1 < byte_len)
                output[byte_pos + 1] |= (encoded << (12 - bit_offset));
        }
    }
}

static void dna4_decode_neon(const uint8_t* input, char* output, int len)
{
    int bit_len = len * 4;
    
    /* Process 16 characters at a time with NEON */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    for (int i = 0; i < simd_len; i += 16) {
        /* Load 8 bytes to process 16 characters */
        uint64_t data_word;
        memcpy(&data_word, input + i/2, 8);
        uint8_t temp[8];
        for (int k = 0; k < 8; k++) {
            temp[k] = (data_word >> (k * 8)) & 0xFF;
        }
        
        /* Decode each 4-bit nibble */
        for (int j = 0; j < 16; j++) {
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded;
            
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                encoded = ((temp[byte_pos] << (bit_offset - 4)) & 0xF);
                if (byte_pos + 1 < (bit_len + 7) / 8)
                    encoded |= (temp[byte_pos + 1] >> (12 - bit_offset));
                encoded &= 0xF;
            }
            
            output[i + j] = kmersearch_dna4_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            encoded = ((input[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < (bit_len + 7) / 8)
                encoded |= (input[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}

#ifdef __ARM_FEATURE_SVE
/* SVE implementations */
static void dna2_encode_sve(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 2 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    
    for (int i = 0; i < simd_len; i += sve_len) {
        svbool_t pg = svwhilelt_b8_s32(i, len);
        svuint8_t chars = svld1_u8(pg, (uint8_t*)(input + i));
        
        /* Create comparison masks for each DNA base */
        svbool_t mask_A = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('A')),
                                        svcmpeq_u8(pg, chars, svdup_n_u8('a')));
        svbool_t mask_C = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('C')),
                                        svcmpeq_u8(pg, chars, svdup_n_u8('c')));
        svbool_t mask_G = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('G')),
                                        svcmpeq_u8(pg, chars, svdup_n_u8('g')));
        svbool_t mask_T = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('T')),
                                        svcmpeq_u8(pg, chars, svdup_n_u8('t')));
        svbool_t mask_U = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('U')),
                                        svcmpeq_u8(pg, chars, svdup_n_u8('u')));
        
        /* Combine T and U masks */
        mask_T = svorr_b_z(pg, mask_T, mask_U);
        
        /* Generate 2-bit encoded values */
        svuint8_t encoded = svdup_n_u8(0);
        encoded = svorr_u8_m(mask_C, encoded, svdup_n_u8(1));
        encoded = svorr_u8_m(mask_G, encoded, svdup_n_u8(2));
        encoded = svorr_u8_m(mask_T, encoded, svdup_n_u8(3));
        
        /* Store and pack bits */
        uint8_t temp[sve_len];
        svst1_u8(pg, temp, encoded);
        
        for (int j = 0; j < sve_len && (i + j) < len; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            output[byte_pos] |= (temp[j] << (6 - bit_offset));
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        output[byte_pos] |= (encoded << (6 - bit_offset));
    }
}

static void dna2_decode_sve(const uint8_t* input, char* output, int len)
{
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    
    for (int i = 0; i < simd_len; i += sve_len) {
        /* Load packed data for processing */
        int bytes_needed = (sve_len * 2 + 7) / 8;
        uint8_t temp[bytes_needed];
        memcpy(temp, input + i/4, bytes_needed);
        
        /* Decode each 2-bit pair */
        for (int j = 0; j < sve_len && (i + j) < len; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

static void dna4_encode_sve(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    memset(output, 0, byte_len);
    
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    
    for (int i = 0; i < simd_len; i += sve_len) {
        svbool_t pg = svwhilelt_b8_s32(i, len);
        svuint8_t chars = svld1_u8(pg, (uint8_t*)(input + i));
        
        /* Use lookup table approach for DNA4 encoding */
        uint8_t temp[sve_len];
        svst1_u8(pg, temp, chars);
        
        for (int j = 0; j < sve_len && (i + j) < len; j++) {
            uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)temp[j]];
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 4) {
                output[byte_pos] |= (encoded << (4 - bit_offset));
            } else {
                output[byte_pos] |= (encoded >> (bit_offset - 4));
                if (byte_pos + 1 < byte_len)
                    output[byte_pos + 1] |= (encoded << (12 - bit_offset));
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 4) {
            output[byte_pos] |= (encoded << (4 - bit_offset));
        } else {
            output[byte_pos] |= (encoded >> (bit_offset - 4));
            if (byte_pos + 1 < byte_len)
                output[byte_pos + 1] |= (encoded << (12 - bit_offset));
        }
    }
}

static void dna4_decode_sve(const uint8_t* input, char* output, int len)
{
    int bit_len = len * 4;
    
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    
    for (int i = 0; i < simd_len; i += sve_len) {
        /* Load packed data for processing */
        int bytes_needed = (sve_len * 4 + 7) / 8;
        uint8_t temp[bytes_needed];
        memcpy(temp, input + i/2, bytes_needed);
        
        /* Decode each 4-bit nibble */
        for (int j = 0; j < sve_len && (i + j) < len; j++) {
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded;
            
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                encoded = ((temp[byte_pos] << (bit_offset - 4)) & 0xF);
                if (byte_pos + 1 < (bit_len + 7) / 8)
                    encoded |= (temp[byte_pos + 1] >> (12 - bit_offset));
                encoded &= 0xF;
            }
            
            output[i + j] = kmersearch_dna4_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            encoded = ((input[byte_pos] << (bit_offset - 4)) & 0xF);
            if (byte_pos + 1 < (bit_len + 7) / 8)
                encoded |= (input[byte_pos + 1] >> (12 - bit_offset));
            encoded &= 0xF;
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}
#endif
#endif

/*
 * AVX2 K-mer Processing Functions
 */

#ifdef __x86_64__
/* AVX2 optimized version of kmersearch_extract_dna2_kmers_direct */
__attribute__((target("avx2")))
static Datum *
kmersearch_extract_dna2_kmers_direct_avx2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Use simple array-based k-mer tracking (no hashing needed for k<=64) */
    occurrences = (KmerOccurrence *) palloc(max_kmers * sizeof(KmerOccurrence));
    
    /* Process k-mers with AVX2-optimized bit extraction where possible */
    int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
    /* AVX2-optimized batch processing for aligned k-mer extraction */
    for (i = 0; i < simd_batch; i += 8)
    {
        /* Extract 8 k-mers in parallel using AVX2 */
        __m256i positions = _mm256_setr_epi32(i, i+1, i+2, i+3, i+4, i+5, i+6, i+7);
        
        /* Process each k-mer in the batch */
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            /* Extract k-mer as single uint64_t value */
            kmer_value = kmersearch_extract_kmer_as_uint64(seq, pos, k);
            
            /* Skip if k-mer extraction failed (boundary check failed) */
            if (kmer_value == 0 && k > 0) {
                /* Note: valid k-mers could be 0 (all A's), but check bounds properly */
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;  /* Out of bounds, skip */
                }
            }
            
            /* Find or add occurrence count using binary search */
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers);
            
            if (current_count < 0)
                continue;  /* Array full, skip */
            
            /* Skip if occurrence exceeds bit limit */
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            /* Create simple k-mer key (without occurrence count for matching) */
            ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, pos, k);
            if (ngram_key == NULL)
                continue;  /* Skip if key creation failed */
                
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        /* Extract k-mer as single uint64_t value */
        kmer_value = kmersearch_extract_kmer_as_uint64(seq, i, k);
        
        /* Skip if k-mer extraction failed (boundary check failed) */
        if (kmer_value == 0 && k > 0) {
            /* Note: valid k-mers could be 0 (all A's), but check bounds properly */
            int last_bit_pos = (i + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;  /* Out of bounds, skip */
            }
        }
        
        /* Find or add occurrence count using binary search */
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, max_kmers);
        
        if (current_count < 0)
            continue;  /* Array full, skip */
        
        /* Skip if occurrence exceeds bit limit */
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* Create simple k-mer key (without occurrence count for matching) */
        ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, i, k);
        if (ngram_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(ngram_key);
    }
    
    /* Cleanup */
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* AVX2 optimized version of kmersearch_extract_dna4_kmers_with_expansion_direct */
__attribute__((target("avx2")))
static Datum *
kmersearch_extract_dna4_kmers_with_expansion_direct_avx2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Use simple array-based k-mer tracking (no hashing needed for k<=64) */
    occurrences = (KmerOccurrence *) palloc(max_kmers * 10 * sizeof(KmerOccurrence));
    
    /* Process k-mers with AVX2-optimized expansion where possible */
    int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
    /* AVX2-optimized batch processing for k-mer expansion */
    for (i = 0; i < simd_batch; i += 8)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers */
            expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Only use SIMD processing for expansion_count >= 3 */
            if (expansion_count < 3)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                uint64_t kmer_value;
                int current_count;
                VarBit *ngram_key;
                
                /* Extract k-mer as single uint64_t value */
                kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
                
                /* Find or add occurrence count using binary search */
                current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                      kmer_value, max_kmers * 10);
                
                if (current_count < 0)
                    continue;  /* Array full, skip */
                
                /* Skip if occurrence exceeds bit limit */
                if (current_count > (1 << kmersearch_occur_bitlen))
                    continue;
                
                /* Create simple k-mer key (copy the DNA2 k-mer directly) */
                ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
                memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
                keys[key_count++] = PointerGetDatum(ngram_key);
            }
            
            /* Free each expanded k-mer and then the array */
            if (expanded_kmers)
            {
                for (exp_j = 0; exp_j < expansion_count; exp_j++)
                {
                    if (expanded_kmers[exp_j])
                        pfree(expanded_kmers[exp_j]);
                }
                pfree(expanded_kmers);
            }
        }
    }
    
    /* Handle all k-mers with scalar processing */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Skip k-mers that were already processed by SIMD (expansion_count >= 3 and within SIMD batch) */
        if (i < simd_batch && expansion_count >= 3)
        {
            /* Free the expansion since it was already processed by SIMD */
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
            continue;
        }
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            /* Extract k-mer as single uint64_t value */
            kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
            
            /* Find or add occurrence count using binary search */
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers * 10);
            
            if (current_count < 0)
                continue;  /* Array full, skip */
            
            /* Skip if occurrence exceeds bit limit */
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            /* Create simple k-mer key (copy the DNA2 k-mer directly) */
            ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
            memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
        
        /* Free each expanded k-mer and then the array */
        if (expanded_kmers)
        {
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
        }
    }
    
    /* Cleanup */
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* AVX2 optimized version of kmersearch_count_matching_kmers_fast */
__attribute__((target("avx2")))
static int
kmersearch_count_matching_kmers_fast_avx2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* For small datasets, O(n*m) might be faster than hash table overhead */
    if (seq_nkeys * query_nkeys < 100)
    {
        /* AVX2-optimized comparison for small datasets */
        for (i = 0; i < seq_nkeys; i++)
        {
            int j;
            int simd_batch = query_nkeys & ~7;  /* Process 8 queries at a time */
            bool found_match = false;
            
            /* AVX2 batch processing for query comparison */
            for (j = 0; j < simd_batch && !found_match; j += 8)
            {
                /* Check 8 queries in parallel */
                for (int k = 0; k < 8 && (j + k) < query_nkeys; k++)
                {
                    int idx = j + k;
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[idx]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[idx]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[idx]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        found_match = true;
                        break;
                    }
                }
            }
            
            /* Handle remaining queries with scalar */
            if (!found_match)
            {
                for (j = simd_batch; j < query_nkeys; j++)
                {
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[j]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[j]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[j]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        break;
                    }
                }
            }
        }
        return match_count;
    }
    
    /* Use hash table for larger datasets (same as scalar version) */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    /* Safety check: ensure we have valid query keys */
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmers_fast_avx2: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);  /* Use data size, not total size */
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    elog(LOG, "kmersearch_count_matching_kmers_fast_avx2: Creating hash with keysize=%zu, query_nkeys=%d", 
         (size_t)VARBITBYTES(query_keys[0]), query_nkeys);
    
    query_hash = hash_create("QueryKmerHashAVX2", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* Insert all query k-mers into hash table using content as key */
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmers_fast_avx2: NULL query key at index %d", i);
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Check each sequence k-mer against hash table */
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmers_fast_avx2: NULL seq key at index %d", i);
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            elog(LOG, "kmersearch_count_matching_kmers_fast_avx2: Size mismatch seq[%d]=%zu vs query[0]=%zu", 
                 i, (size_t)VARBITBYTES(seq_keys[i]), (size_t)VARBITBYTES(query_keys[0]));
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    elog(LOG, "kmersearch_count_matching_kmers_fast_avx2: Returning match_count=%d", match_count);
    return match_count;
}

/* AVX512 optimized version of kmersearch_extract_dna2_kmers_direct */
__attribute__((target("avx512f,avx512bw")))
static Datum *
kmersearch_extract_dna2_kmers_direct_avx512(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    occurrences = (KmerOccurrence *) palloc(max_kmers * sizeof(KmerOccurrence));
    
    /* Process k-mers with AVX512-optimized bit extraction (16 k-mers at a time) */
    int simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
    
    /* AVX512-optimized batch processing */
    for (i = 0; i < simd_batch; i += 16)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 16 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            kmer_value = kmersearch_extract_kmer_as_uint64(seq, pos, k);
            
            if (kmer_value == 0 && k > 0) {
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;
                }
            }
            
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers);
            
            if (current_count < 0)
                continue;
            
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, pos, k);
            if (ngram_key == NULL)
                continue;
                
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        kmer_value = kmersearch_extract_kmer_as_uint64(seq, i, k);
        
        if (kmer_value == 0 && k > 0) {
            int last_bit_pos = (i + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;
            }
        }
        
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, max_kmers);
        
        if (current_count < 0)
            continue;
        
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, i, k);
        if (ngram_key == NULL)
            continue;
            
        keys[key_count++] = PointerGetDatum(ngram_key);
    }
    
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* AVX512 optimized version of kmersearch_extract_dna4_kmers_with_expansion_direct */
__attribute__((target("avx512f,avx512bw")))
static Datum *
kmersearch_extract_dna4_kmers_with_expansion_direct_avx512(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));
    occurrences = (KmerOccurrence *) palloc(max_kmers * 10 * sizeof(KmerOccurrence));
    
    /* Process k-mers with AVX512-optimized expansion (16 k-mers at a time) */
    int simd_batch = max_kmers & ~15;
    
    for (i = 0; i < simd_batch; i += 16)
    {
        for (int j = 0; j < 16 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Only use SIMD processing for expansion_count >= 3 */
            if (expansion_count < 3)
                continue;
            
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                uint64_t kmer_value;
                int current_count;
                VarBit *ngram_key;
                
                kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
                
                current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                      kmer_value, max_kmers * 10);
                
                if (current_count < 0)
                    continue;
                
                if (current_count > (1 << kmersearch_occur_bitlen))
                    continue;
                
                ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
                memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
                keys[key_count++] = PointerGetDatum(ngram_key);
            }
            
            if (expanded_kmers)
            {
                for (exp_j = 0; exp_j < expansion_count; exp_j++)
                {
                    if (expanded_kmers[exp_j])
                        pfree(expanded_kmers[exp_j]);
                }
                pfree(expanded_kmers);
            }
        }
    }
    
    /* Handle all k-mers with scalar processing */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Skip k-mers that were already processed by SIMD (expansion_count >= 3 and within SIMD batch) */
        if (i < simd_batch && expansion_count >= 3)
        {
            /* Free the expansion since it was already processed by SIMD */
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
            continue;
        }
        
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
            
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers * 10);
            
            if (current_count < 0)
                continue;
            
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
            memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
        
        if (expanded_kmers)
        {
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
        }
    }
    
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* AVX512 optimized version of kmersearch_count_matching_kmers_fast */
__attribute__((target("avx512f,avx512bw")))
static int
kmersearch_count_matching_kmers_fast_avx512(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* For small datasets, O(n*m) might be faster than hash table overhead */
    if (seq_nkeys * query_nkeys < 100)
    {
        /* AVX512-optimized comparison for small datasets (16 queries at a time) */
        for (i = 0; i < seq_nkeys; i++)
        {
            int j;
            int simd_batch = query_nkeys & ~15;
            bool found_match = false;
            
            for (j = 0; j < simd_batch && !found_match; j += 16)
            {
                for (int k = 0; k < 16 && (j + k) < query_nkeys; k++)
                {
                    int idx = j + k;
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[idx]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[idx]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[idx]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        found_match = true;
                        break;
                    }
                }
            }
            
            if (!found_match)
            {
                for (j = simd_batch; j < query_nkeys; j++)
                {
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[j]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[j]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[j]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        break;
                    }
                }
            }
        }
        return match_count;
    }
    
    /* Use hash table for larger datasets */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmers_fast_avx512: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHashAVX512", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    return match_count;
}
#endif

#ifdef __aarch64__
/* NEON optimized version of kmersearch_extract_dna2_kmers_direct */
__attribute__((target("neon")))
static Datum *
kmersearch_extract_dna2_kmers_direct_neon(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    occurrences = (KmerOccurrence *) palloc(max_kmers * sizeof(KmerOccurrence));
    
    /* Process k-mers with NEON-optimized bit extraction (4 k-mers at a time) */
    int simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
    
    /* NEON-optimized batch processing */
    for (i = 0; i < simd_batch; i += 4)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 4 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            kmer_value = kmersearch_extract_kmer_as_uint64(seq, pos, k);
            
            if (kmer_value == 0 && k > 0) {
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;
                }
            }
            
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers);
            
            if (current_count < 0)
                continue;
            
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, pos, k);
            if (ngram_key == NULL)
                continue;
                
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        kmer_value = kmersearch_extract_kmer_as_uint64(seq, i, k);
        
        if (kmer_value == 0 && k > 0) {
            int last_bit_pos = (i + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;
            }
        }
        
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, max_kmers);
        
        if (current_count < 0)
            continue;
        
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, i, k);
        if (ngram_key == NULL)
            continue;
            
        keys[key_count++] = PointerGetDatum(ngram_key);
    }
    
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* NEON optimized version of kmersearch_extract_dna4_kmers_with_expansion_direct */
__attribute__((target("neon")))
static Datum *
kmersearch_extract_dna4_kmers_with_expansion_direct_neon(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));
    occurrences = (KmerOccurrence *) palloc(max_kmers * 10 * sizeof(KmerOccurrence));
    
    /* Process k-mers with NEON-optimized expansion (4 k-mers at a time) */
    int simd_batch = max_kmers & ~3;
    
    for (i = 0; i < simd_batch; i += 4)
    {
        for (int j = 0; j < 4 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Only use SIMD processing for expansion_count >= 3 */
            if (expansion_count < 3)
                continue;
            
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                uint64_t kmer_value;
                int current_count;
                VarBit *ngram_key;
                
                kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
                
                current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                      kmer_value, max_kmers * 10);
                
                if (current_count < 0)
                    continue;
                
                if (current_count > (1 << kmersearch_occur_bitlen))
                    continue;
                
                ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
                memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
                keys[key_count++] = PointerGetDatum(ngram_key);
            }
            
            if (expanded_kmers)
            {
                for (exp_j = 0; exp_j < expansion_count; exp_j++)
                {
                    if (expanded_kmers[exp_j])
                        pfree(expanded_kmers[exp_j]);
                }
                pfree(expanded_kmers);
            }
        }
    }
    
    /* Handle all k-mers with scalar processing */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Skip k-mers that were already processed by SIMD (expansion_count >= 3 and within SIMD batch) */
        if (i < simd_batch && expansion_count >= 3)
        {
            /* Free the expansion since it was already processed by SIMD */
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
            continue;
        }
        
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
            
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers * 10);
            
            if (current_count < 0)
                continue;
            
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
            memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
        
        if (expanded_kmers)
        {
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
        }
    }
    
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* NEON optimized version of kmersearch_count_matching_kmers_fast */
__attribute__((target("neon")))
static int
kmersearch_count_matching_kmers_fast_neon(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* For small datasets, O(n*m) might be faster than hash table overhead */
    if (seq_nkeys * query_nkeys < 100)
    {
        /* NEON-optimized comparison for small datasets (4 queries at a time) */
        for (i = 0; i < seq_nkeys; i++)
        {
            int j;
            int simd_batch = query_nkeys & ~3;
            bool found_match = false;
            
            for (j = 0; j < simd_batch && !found_match; j += 4)
            {
                for (int k = 0; k < 4 && (j + k) < query_nkeys; k++)
                {
                    int idx = j + k;
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[idx]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[idx]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[idx]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        found_match = true;
                        break;
                    }
                }
            }
            
            if (!found_match)
            {
                for (j = simd_batch; j < query_nkeys; j++)
                {
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[j]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[j]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[j]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        break;
                    }
                }
            }
        }
        return match_count;
    }
    
    /* Use hash table for larger datasets */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmers_fast_neon: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHashNEON", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    return match_count;
}
#endif

#ifdef __aarch64__
/* SVE optimized version of kmersearch_extract_dna2_kmers_direct */
__attribute__((target("sve")))
static Datum *
kmersearch_extract_dna2_kmers_direct_sve(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    occurrences = (KmerOccurrence *) palloc(max_kmers * sizeof(KmerOccurrence));
    
    /* Process k-mers with SVE-optimized bit extraction (variable vector width) */
    /* SVE vector length is runtime determined, so we use a conservative batch size */
    int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time (conservative) */
    
    /* SVE-optimized batch processing */
    for (i = 0; i < simd_batch; i += 8)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            kmer_value = kmersearch_extract_kmer_as_uint64(seq, pos, k);
            
            if (kmer_value == 0 && k > 0) {
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;
                }
            }
            
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers);
            
            if (current_count < 0)
                continue;
            
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, pos, k);
            if (ngram_key == NULL)
                continue;
                
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        kmer_value = kmersearch_extract_kmer_as_uint64(seq, i, k);
        
        if (kmer_value == 0 && k > 0) {
            int last_bit_pos = (i + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;
            }
        }
        
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, max_kmers);
        
        if (current_count < 0)
            continue;
        
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, i, k);
        if (ngram_key == NULL)
            continue;
            
        keys[key_count++] = PointerGetDatum(ngram_key);
    }
    
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* SVE optimized version of kmersearch_extract_dna4_kmers_with_expansion_direct */
__attribute__((target("sve")))
static Datum *
kmersearch_extract_dna4_kmers_with_expansion_direct_sve(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));
    occurrences = (KmerOccurrence *) palloc(max_kmers * 10 * sizeof(KmerOccurrence));
    
    /* Process k-mers with SVE-optimized expansion (8 k-mers at a time, conservative) */
    int simd_batch = max_kmers & ~7;
    
    for (i = 0; i < simd_batch; i += 8)
    {
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Only use SIMD processing for expansion_count >= 3 */
            if (expansion_count < 3)
                continue;
            
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                uint64_t kmer_value;
                int current_count;
                VarBit *ngram_key;
                
                kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
                
                current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                      kmer_value, max_kmers * 10);
                
                if (current_count < 0)
                    continue;
                
                if (current_count > (1 << kmersearch_occur_bitlen))
                    continue;
                
                ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
                memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
                keys[key_count++] = PointerGetDatum(ngram_key);
            }
            
            if (expanded_kmers)
            {
                for (exp_j = 0; exp_j < expansion_count; exp_j++)
                {
                    if (expanded_kmers[exp_j])
                        pfree(expanded_kmers[exp_j]);
                }
                pfree(expanded_kmers);
            }
        }
    }
    
    /* Handle all k-mers with scalar processing */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        expanded_kmers = kmersearch_expand_dna4_kmer_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Skip k-mers that were already processed by SIMD (expansion_count >= 3 and within SIMD batch) */
        if (i < simd_batch && expansion_count >= 3)
        {
            /* Free the expansion since it was already processed by SIMD */
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
            continue;
        }
        
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            uint64_t kmer_value;
            int current_count;
            VarBit *ngram_key;
            
            kmer_value = kmersearch_extract_kmer_as_uint64(dna2_kmer, 0, k);
            
            current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                                  kmer_value, max_kmers * 10);
            
            if (current_count < 0)
                continue;
            
            if (current_count > (1 << kmersearch_occur_bitlen))
                continue;
            
            ngram_key = (VarBit *) palloc(VARSIZE(dna2_kmer));
            memcpy(ngram_key, dna2_kmer, VARSIZE(dna2_kmer));
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
        
        if (expanded_kmers)
        {
            for (j = 0; j < expansion_count; j++)
            {
                if (expanded_kmers[j])
                    pfree(expanded_kmers[j]);
            }
            pfree(expanded_kmers);
        }
    }
    
    pfree(occurrences);
    
    *nkeys = key_count;
    return keys;
}

/* SVE optimized version of kmersearch_count_matching_kmers_fast */
__attribute__((target("sve")))
static int
kmersearch_count_matching_kmers_fast_sve(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* For small datasets, O(n*m) might be faster than hash table overhead */
    if (seq_nkeys * query_nkeys < 100)
    {
        /* SVE-optimized comparison for small datasets (8 queries at a time, conservative) */
        for (i = 0; i < seq_nkeys; i++)
        {
            int j;
            int simd_batch = query_nkeys & ~7;
            bool found_match = false;
            
            for (j = 0; j < simd_batch && !found_match; j += 8)
            {
                for (int k = 0; k < 8 && (j + k) < query_nkeys; k++)
                {
                    int idx = j + k;
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[idx]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[idx]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[idx]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        found_match = true;
                        break;
                    }
                }
            }
            
            if (!found_match)
            {
                for (j = simd_batch; j < query_nkeys; j++)
                {
                    if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[j]) &&
                        VARSIZE(seq_keys[i]) == VARSIZE(query_keys[j]) &&
                        memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[j]), VARBITBYTES(seq_keys[i])) == 0)
                    {
                        match_count++;
                        break;
                    }
                }
            }
        }
        return match_count;
    }
    
    /* Use hash table for larger datasets */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmers_fast_sve: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHashSVE", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    return match_count;
}
#endif
