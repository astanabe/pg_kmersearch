/*-------------------------------------------------------------------------
 *
 * kmersearch.h
 *		PostgreSQL extension for DNA sequence k-mer search
 *
 * This header contains shared type definitions, constants, and function
 * declarations for the pg_kmersearch extension.
 *
 * Copyright (c) 2024, pg_kmersearch contributors
 *
 *-------------------------------------------------------------------------
 */
#ifndef KMERSEARCH_H
#define KMERSEARCH_H

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
#include "access/parallel.h"
#include "storage/dsm.h"
#include "utils/dsa.h"
#include "lib/dshash.h"
#include "storage/ipc.h"
#include "postmaster/bgworker.h"
#include "utils/backend_status.h"
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

/*
 * SIMD capability detection
 */
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

/*
 * K-mer frequency entry for analysis
 */
typedef struct KmerFreqEntry
{
    VarBit      *kmer_key;      /* K-mer binary key (without occurrence count) */
    int         row_count;      /* Number of rows where this k-mer appears */
    bool        highfreq;       /* Whether this k-mer is highly frequent */
} KmerFreqEntry;

/*
 * K-mer match result structure
 */
typedef struct KmerMatchResult
{
    int         shared_count;       /* Number of shared k-mers (rawscore) */
    int         seq_nkeys;          /* Number of k-mers in sequence */
    int         query_nkeys;        /* Number of k-mers in query */
    double      sharing_rate;       /* Calculated sharing rate */
    bool        match_result;       /* Boolean match result for =% operator */
    bool        valid;              /* Whether the result is valid */
} KmerMatchResult;

/*
 * Actual min score cache entry
 */
typedef struct ActualMinScoreCacheEntry
{
    uint64      query_hash;                /* Hash value of query string */
    int         actual_min_score;          /* Actual minimum score for this query */
} ActualMinScoreCacheEntry;

/*
 * Actual min score cache manager
 */
typedef struct ActualMinScoreCacheManager
{
    HTAB        *cache_hash;               /* Hash table for actual min score cache */
    MemoryContext cache_context;           /* Memory context for cache */
    int         hits;                      /* Cache hit count */
    int         misses;                    /* Cache miss count */
    int         max_entries;               /* Maximum number of entries */
    int         current_entries;           /* Current number of entries */
} ActualMinScoreCacheManager;

/*
 * K-mer occurrence tracking
 */
typedef struct KmerOccurrence
{
    uint64_t    kmer_value;         /* K-mer as single 64-bit value */
    int         count;              /* Occurrence count */
} KmerOccurrence;

/*
 * K-mer analysis result
 */
typedef struct KmerAnalysisResult
{
    int64       total_rows;                /* Total number of rows analyzed */
    int         highfreq_kmers_count;      /* Number of highly frequent k-mers */
    int         parallel_workers_used;     /* Number of parallel workers used */
    double      analysis_duration;        /* Analysis duration in seconds */
    double      max_appearance_rate_used; /* Max appearance rate used */
    int         max_appearance_nrow_used;  /* Max appearance nrow used */
} KmerAnalysisResult;

/*
 * Drop analysis result
 */
typedef struct DropAnalysisResult
{
    int         dropped_analyses;          /* Number of dropped analyses */
    int         dropped_highfreq_kmers;    /* Number of dropped highly frequent k-mers */
    int64       freed_storage_bytes;       /* Freed storage in bytes */
} DropAnalysisResult;

/*
 * Rawscore cache entry
 */
typedef struct RawscoreCacheEntry
{
    uint64      hash_key;                  /* Hash key for this entry */
    VarBit      *sequence_copy;            /* Copy of sequence for exact match */
    char        *query_string_copy;        /* Copy of query string for exact match */
    KmerMatchResult result;                /* Cached rawscore calculation result */
    int         heap_index;                /* Index in min-heap array (-1 if not in heap) */
} RawscoreCacheEntry;

/*
 * Rawscore cache manager
 */
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

/*
 * Query pattern cache entry
 */
typedef struct QueryPatternCacheEntry
{
    uint64      hash_key;                  /* Hash key for this entry */
    char        *query_string_copy;        /* Copy of query string */
    int         kmer_size;                 /* K-mer size for this pattern */
    VarBit      **extracted_kmers;         /* Cached extracted k-mers */
    int         kmer_count;                /* Number of extracted k-mers */
    struct QueryPatternCacheEntry *next;   /* For LRU chain */
    struct QueryPatternCacheEntry *prev;   /* For LRU chain */
} QueryPatternCacheEntry;

/*
 * Query pattern cache manager
 */
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

/*
 * K-mer data union for different k-values
 */
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

/*
 * Compact k-mer frequency entry
 */
typedef struct CompactKmerFreq
{
    KmerData    kmer_data;               /* Raw k-mer data */
    int         frequency_count;          /* Frequency count */
    bool        is_highfreq;              /* Whether this k-mer is highly frequent */
} CompactKmerFreq;

/*
 * K-mer buffer for batch processing
 */
typedef struct KmerBuffer
{
    CompactKmerFreq *entries;            /* Buffer entries */
    int             count;               /* Current count */
    int             capacity;            /* Buffer capacity */
    int             kmer_size;           /* K-mer size */
} KmerBuffer;

/*
 * High-frequency k-mer hash entry
 */
typedef struct HighfreqKmerHashEntry
{
    VarBit         *kmer_key;             /* K-mer key for hashing */
    uint64          hash_value;           /* Precomputed hash value */
} HighfreqKmerHashEntry;

/*
 * High-frequency k-mer cache
 */
typedef struct HighfreqKmerCache
{
    Oid         current_table_oid;        /* Current table OID */
    char       *current_column_name;      /* Current column name */
    int         current_kmer_size;        /* Current k-mer size */
    MemoryContext cache_context;          /* Memory context for cache data */
    HTAB       *highfreq_hash;           /* Hash table for fast lookup */
    VarBit    **highfreq_kmers;          /* Array of high-frequency k-mers */
    int         highfreq_count;          /* Number of high-frequency k-mers */
    bool        is_valid;                /* Cache validity flag */
} HighfreqKmerCache;

/*
 * Parallel high-frequency k-mer cache entry
 */
typedef struct ParallelHighfreqKmerCacheEntry
{
    uint64      kmer_hash;               /* k-mer hash value (key) */
    int32       frequency_count;         /* frequency count */
    Oid         table_oid;               /* table OID */
    int32       kmer_size;               /* k-mer size */
} ParallelHighfreqKmerCacheEntry;

/*
 * Parallel high-frequency k-mer cache
 */
typedef struct ParallelHighfreqKmerCache
{
    dshash_table_handle hash_handle;     /* dshash table handle */
    int32               num_entries;     /* number of entries */
    Size                segment_size;    /* DSM segment size */
    bool                is_initialized;  /* initialization flag */
    dsm_handle          dsm_handle;      /* DSM segment handle */
    Oid                 table_oid;       /* table OID */
    int32               kmer_size;       /* k-mer size */
} ParallelHighfreqKmerCache;

/*
 * Worker state for parallel k-mer analysis
 */
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

/* DNA type definitions */
typedef struct
{
    int32 vl_len_;
    bits8 data[FLEXIBLE_ARRAY_MEMBER];
} kmersearch_dna2;

typedef struct
{
    int32 vl_len_;
    bits8 data[FLEXIBLE_ARRAY_MEMBER];
} kmersearch_dna4;

/*
 * External variable declarations
 */
extern simd_dispatch_table_t simd_dispatch;
extern simd_capability_t simd_capability;

/* Global configuration variables */
extern int kmersearch_occur_bitlen;
extern int kmersearch_kmer_size;
extern double kmersearch_max_appearance_rate;
extern int kmersearch_max_appearance_nrow;
extern int kmersearch_min_score;
extern double kmersearch_min_shared_ngram_key_rate;
extern bool kmersearch_preclude_highfreq_kmer;

/* Cache configuration variables */
extern int kmersearch_rawscore_cache_max_entries;
extern int kmersearch_query_pattern_cache_max_entries;
extern int kmersearch_actual_min_score_cache_max_entries;

/* Global cache managers */
extern ActualMinScoreCacheManager *actual_min_score_cache_manager;
extern QueryPatternCacheManager *query_pattern_cache_manager;
extern RawscoreCacheManager *rawscore_cache_manager;

/* Global high-frequency k-mer cache */
extern HighfreqKmerCache global_highfreq_cache;

/* Global testing variable for dshash usage */
extern bool kmersearch_force_use_dshash;

/* Global parallel cache state */
extern ParallelHighfreqKmerCache *parallel_highfreq_cache;
extern dsm_segment *parallel_cache_segment;
extern dsa_area *parallel_cache_dsa;
extern dshash_table *parallel_cache_hash;
extern bool parallel_cache_exit_callback_registered;

/* DNA encoding/decoding tables */
extern const uint8 kmersearch_dna2_encode_table[256];
extern const char kmersearch_dna2_decode_table[4];
extern const uint8 kmersearch_dna4_encode_table[256];
extern const char kmersearch_dna4_decode_table[16];

/*
 * Function declarations
 */

/* Module initialization */
void _PG_init(void);
void _PG_fini(void);

/* DNA datatype functions (implemented in kmersearch_datatype.c) */
Datum kmersearch_dna2_in(PG_FUNCTION_ARGS);
Datum kmersearch_dna2_out(PG_FUNCTION_ARGS);
Datum kmersearch_dna2_recv(PG_FUNCTION_ARGS);
Datum kmersearch_dna2_send(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_in(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_out(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_recv(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_send(PG_FUNCTION_ARGS);
Datum kmersearch_dna2_eq(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_eq(PG_FUNCTION_ARGS);

/* DNA datatype length functions */
Datum kmersearch_dna2_char_length(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_char_length(PG_FUNCTION_ARGS);

/* DNA datatype utility functions */
char *kmersearch_dna2_to_string(VarBit *dna);
char *kmersearch_dna4_to_string(VarBit *dna);

/* GIN operator class functions */
Datum kmersearch_extract_value_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna4(PG_FUNCTION_ARGS);
Datum kmersearch_extract_query_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_extract_query_dna4(PG_FUNCTION_ARGS);
Datum kmersearch_consistent_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_consistent_dna4(PG_FUNCTION_ARGS);
Datum kmersearch_compare_partial_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_compare_partial_dna4(PG_FUNCTION_ARGS);

/* Search operator functions */
Datum kmersearch_match_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_match_dna4(PG_FUNCTION_ARGS);

/* Scoring functions */
Datum kmersearch_rawscore_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_rawscore_dna4(PG_FUNCTION_ARGS);
Datum kmersearch_correctedscore_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_correctedscore_dna4(PG_FUNCTION_ARGS);

/* Cache management functions */
Datum kmersearch_actual_min_score_cache_stats(PG_FUNCTION_ARGS);
Datum kmersearch_actual_min_score_cache_free(PG_FUNCTION_ARGS);
Datum kmersearch_rawscore_cache_stats(PG_FUNCTION_ARGS);
Datum kmersearch_rawscore_cache_free(PG_FUNCTION_ARGS);
Datum kmersearch_query_pattern_cache_stats(PG_FUNCTION_ARGS);
Datum kmersearch_query_pattern_cache_free(PG_FUNCTION_ARGS);

/* High-frequency k-mer cache functions */
Datum kmersearch_highfreq_kmer_cache_load(PG_FUNCTION_ARGS);
Datum kmersearch_highfreq_kmer_cache_free(PG_FUNCTION_ARGS);

/* Parallel cache functions */
Datum kmersearch_parallel_highfreq_kmer_cache_load(PG_FUNCTION_ARGS);
Datum kmersearch_parallel_highfreq_kmer_cache_free(PG_FUNCTION_ARGS);

/* Analysis functions */
Datum kmersearch_analyze_table(PG_FUNCTION_ARGS);
Datum kmersearch_drop_highfreq_analysis(PG_FUNCTION_ARGS);

/* K-mer utility functions */
void kmersearch_expand_degenerate_sequence(const char *kmer, int k, char **expanded, int *expand_count);
VarBit *kmersearch_create_ngram_key2(const char *kmer, int k, int occurrence);
Datum *kmersearch_extract_dna2_kmer2_direct(VarBit *dna, int k, int *nkeys);
Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct(VarBit *dna, int k, int *nkeys);
Datum *kmersearch_filter_highfreq_kmers_from_keys(Datum *keys, int *nkeys, HTAB *highfreq_hash, int k);
Datum *kmersearch_filter_highfreq_kmers_from_keys_parallel(Datum *keys, int *nkeys, int k);
int calculate_actual_min_score(VarBit **query_keys, int nkeys, int original_nkeys);
int kmersearch_count_degenerate_combinations(const char *kmer, int k);
void kmersearch_set_bit_at(bits8 *data, int bit_pos, int value);
bool kmersearch_will_exceed_degenerate_limit_dna4_bits(VarBit *seq, int start_pos, int k);
VarBit *kmersearch_create_ngram_key2(const char *kmer, int k, int occurrence);
VarBit *kmersearch_create_kmer2_key_only(const char *kmer, int k);
VarBit *kmersearch_create_kmer2_key_from_dna2_bits(VarBit *seq, int start_pos, int k);
VarBit *kmersearch_create_ngram_key2_from_dna2_bits(VarBit *seq, int start_pos, int k, int occurrence_count);
VarBit *kmersearch_create_ngram_key2_from_dna4_bits(VarBit *seq, int start_pos, int k, int occurrence_count);
VarBit *kmersearch_create_ngram_key2_with_occurrence_from_dna2(VarBit *dna2_kmer, int k, int occurrence);
VarBit **kmersearch_expand_dna4_kmer2_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count);

/* GIN index support functions (implemented in kmersearch_gin.c) */
Datum kmersearch_extract_value_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna4(PG_FUNCTION_ARGS);
Datum kmersearch_extract_query(PG_FUNCTION_ARGS);
Datum kmersearch_consistent(PG_FUNCTION_ARGS);
Datum kmersearch_compare_partial(PG_FUNCTION_ARGS);

#endif   /* KMERSEARCH_H */