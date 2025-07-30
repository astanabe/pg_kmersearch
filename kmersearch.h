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
#include "storage/lmgr.h"
#include "utils/dsa.h"
#include "lib/dshash.h"
#include "storage/ipc.h"
#include "postmaster/bgworker.h"
#include "utils/backend_status.h"
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

/* SIMD intrinsics headers */
#ifdef __x86_64__
#include <cpuid.h>
#include <immintrin.h>
#elif defined(__aarch64__)
#include <stdio.h>
#include <signal.h>
#include <setjmp.h>
#include <arm_neon.h>
#include <arm_sve.h>
#endif
#include <unistd.h>

/*
 * LWLock tranche for parallel cache and analysis
 */
#define LWTRANCHE_KMERSEARCH_CACHE    (LWTRANCHE_FIRST_USER_DEFINED + 100)
#define LWTRANCHE_KMERSEARCH_ANALYSIS (LWTRANCHE_FIRST_USER_DEFINED + 101)

/*
 * SIMD capability detection
 */
typedef enum {
    SIMD_NONE = 0,              /* No SIMD support */
    
    /* x86/x64 SIMD capabilities */
    SIMD_AVX2 = 1,              /* AVX2 support */
    SIMD_BMI2 = 2,              /* AVX2 + BMI2 support */
    SIMD_AVX512F = 3,           /* AVX512F support */
    SIMD_AVX512BW = 4,          /* AVX512BW support */
    SIMD_AVX512VBMI = 5,        /* AVX512 + VBMI support */
    SIMD_AVX512VBMI2 = 6,       /* AVX512 + VBMI + VBMI2 support */
    
    /* ARM SIMD capabilities */
    SIMD_NEON = 21,             /* ARM NEON support */
    SIMD_SVE = 22,              /* ARM SVE support */
    SIMD_SVE2 = 23              /* ARM SVE2 support */
} simd_capability_t;

/*
 * SIMD comparison thresholds (bit length)
 * These thresholds determine when to use SIMD over memcmp
 */
#define SIMD_COMPARE_AVX2_THRESHOLD    128     /* 128 bits: Use AVX2 */
#define SIMD_COMPARE_AVX512_THRESHOLD  256     /* 256 bits: Use AVX512 */
#define SIMD_COMPARE_NEON_THRESHOLD    64      /* 64 bits: Use NEON */
#define SIMD_COMPARE_SVE_THRESHOLD     128     /* 128 bits: Use SVE */
#define SIMD_COMPARE_SVE2_THRESHOLD    128     /* 128 bits: Use SVE2 */

/*
 * SIMD k-mer extraction thresholds (sequence bit length)
 * These thresholds determine when to use SIMD for k-mer extraction
 */
#define SIMD_EXTRACT_AVX2_THRESHOLD    512     /* 512 bits: Use AVX2 for extraction */
#define SIMD_EXTRACT_AVX512_THRESHOLD  1024    /* 1024 bits: Use AVX512 for extraction */
#define SIMD_EXTRACT_NEON_THRESHOLD    256     /* 256 bits: Use NEON for extraction */
#define SIMD_EXTRACT_SVE_THRESHOLD     512     /* 512 bits: Use SVE for extraction */
#define SIMD_EXTRACT_SVE2_THRESHOLD    512     /* 512 bits: Use SVE2 for extraction */

/*
 * SIMD k-mer matching thresholds (key combination count)
 * These thresholds determine when to use SIMD for k-mer matching
 */
#define SIMD_KEYCOMB_AVX2_THRESHOLD    128     /* 128 combinations: Use AVX2 for matching */
#define SIMD_KEYCOMB_AVX512_THRESHOLD  256     /* 256 combinations: Use AVX512 for matching */
#define SIMD_KEYCOMB_NEON_THRESHOLD    64      /* 64 combinations: Use NEON for matching */
#define SIMD_KEYCOMB_SVE_THRESHOLD     128     /* 128 combinations: Use SVE for matching */
#define SIMD_KEYCOMB_SVE2_THRESHOLD    128     /* 128 combinations: Use SVE2 for matching */

/*
 * SIMD encoding thresholds (input character length)
 * Initially set to same values as SIMD_EXTRACT thresholds
 * TODO: These values need performance testing for optimal settings
 */
#define SIMD_ENCODE_AVX2_THRESHOLD     512     /* 512 chars: Use AVX2 for encoding */
#define SIMD_ENCODE_AVX512_THRESHOLD   1024    /* 1024 chars: Use AVX512 for encoding */
#define SIMD_ENCODE_NEON_THRESHOLD     256     /* 256 chars: Use NEON for encoding */
#define SIMD_ENCODE_SVE_THRESHOLD      512     /* 512 chars: Use SVE for encoding */
#define SIMD_ENCODE_SVE2_THRESHOLD     512     /* 512 chars: Use SVE2 for encoding */

/*
 * SIMD decoding thresholds (bit length)
 * Initially set to same values as SIMD_EXTRACT thresholds
 * TODO: These values need performance testing for optimal settings
 */
#define SIMD_DECODE_AVX2_THRESHOLD     512     /* 512 bits: Use AVX2 for decoding */
#define SIMD_DECODE_AVX512_THRESHOLD   1024    /* 1024 bits: Use AVX512 for decoding */
#define SIMD_DECODE_NEON_THRESHOLD     256     /* 256 bits: Use NEON for decoding */
#define SIMD_DECODE_SVE_THRESHOLD      512     /* 512 bits: Use SVE for decoding */
#define SIMD_DECODE_SVE2_THRESHOLD     512     /* 512 bits: Use SVE2 for decoding */

/*
 * SIMD DNA4 k-mer expansion thresholds (bit length)
 * These thresholds determine when to use SIMD for DNA4 k-mer expansion
 * 
 * NOTE: These values are initially set to match SIMD_EXTRACT_*_THRESHOLD values
 * for consistency, but should be re-evaluated through performance testing
 * as DNA4 expansion has different characteristics than simple extraction.
 * TODO: Conduct performance benchmarks to determine optimal threshold values
 */
#define SIMD_DNA4KMER2_AVX2_THRESHOLD    512     /* 512 bits: Use AVX2 for expansion */
#define SIMD_DNA4KMER2_AVX512_THRESHOLD  1024    /* 1024 bits: Use AVX512 for expansion */
#define SIMD_DNA4KMER2_NEON_THRESHOLD    256     /* 256 bits: Use NEON for expansion */
#define SIMD_DNA4KMER2_SVE_THRESHOLD     512     /* 512 bits: Use SVE for expansion */
#define SIMD_DNA4KMER2_SVE2_THRESHOLD    512     /* 512 bits: Use SVE2 for expansion */

/* Removed SIMD dispatch table - now using direct threshold-based dispatch */

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
 * Shared state for parallel k-mer analysis
 */
typedef struct KmerAnalysisSharedState
{
    LWLockPadded mutex;                   /* Exclusive control mutex */
    int         num_workers;              /* Number of parallel workers */
    Oid         table_oid;                /* Target table OID */
    AttrNumber  column_attnum;            /* Target column attnum */
    Oid         column_type_oid;          /* Target column data type OID */
    int         kmer_size;                /* K-mer size (kmersearch_kmer_size) */
    int         batch_size;               /* Batch size (kmersearch_highfreq_analysis_batch_size) */
    bool        all_processed;            /* All rows processed flag */
    BlockNumber next_block;               /* Next block number to process */
    BlockNumber total_blocks;             /* Total number of blocks in table */
    bool        worker_error_occurred;    /* Parallel worker error flag */
    char        error_message[256];       /* Error message buffer */
} KmerAnalysisSharedState;

/* K-mer entry structures for different sizes */
typedef struct KmerEntry16
{
    uint16      kmer;                     /* kmer2_as_uint (key) */
    int         count;                    /* Number of rows containing this k-mer */
} KmerEntry16;

typedef struct KmerEntry32
{
    uint32      kmer;                     /* kmer2_as_uint (key) */
    int         count;                    /* Number of rows containing this k-mer */
} KmerEntry32;

typedef struct KmerEntry64
{
    uint64      kmer;                     /* kmer2_as_uint (key) */
    int         count;                    /* Number of rows containing this k-mer */
} KmerEntry64;

/* Shared memory keys for parallel processing */
#define KMERSEARCH_KEY_SHARED_STATE  1
#define KMERSEARCH_KEY_HANDLES       2  /* Combined DSM and hash handles */

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
 * High-frequency k-mer cache key structure (24 bytes total)
 */
typedef struct HighfreqCacheKey
{
    Oid         table_oid;               /* Table OID (4 bytes) */
    uint32      column_name_hash;        /* Column name hash (4 bytes) */
    int         kmer_size;               /* K-mer size (4 bytes) */
    int         occur_bitlen;            /* Occurrence bit length (4 bytes) */
    float       max_appearance_rate;     /* Max appearance rate (4 bytes) */
    int         max_appearance_nrow;     /* Max appearance nrow (4 bytes) */
} HighfreqCacheKey;

/*
 * High-frequency k-mer cache
 */
typedef struct HighfreqKmerCache
{
    HighfreqCacheKey current_cache_key;   /* Current cache key */
    MemoryContext cache_context;          /* Memory context for cache data */
    HTAB       *highfreq_hash;           /* Hash table for fast lookup */
    uint64     *highfreq_kmers;          /* Array of high-frequency k-mers as kmer2_as_uint */
    int         highfreq_count;          /* Number of high-frequency k-mers */
    bool        is_valid;                /* Cache validity flag */
} HighfreqKmerCache;

/*
 * Parallel high-frequency k-mer cache entry
 */
typedef struct ParallelHighfreqKmerCacheEntry
{
    uint64      kmer2_as_uint;           /* Direct kmer2 as kmer2_as_uint (key) */
    int32       frequency_count;         /* frequency count */
    HighfreqCacheKey cache_key;          /* cache key for validation */
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
    HighfreqCacheKey    cache_key;       /* cache key for validation */
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
extern int kmersearch_query_pattern_cache_max_entries;
extern int kmersearch_actual_min_score_cache_max_entries;
extern int kmersearch_highfreq_kmer_cache_load_batch_size;
extern int kmersearch_highfreq_analysis_batch_size;

/* Global cache managers */
extern ActualMinScoreCacheManager *actual_min_score_cache_manager;
extern QueryPatternCacheManager *query_pattern_cache_manager;

/* Global high-frequency k-mer cache */
extern HighfreqKmerCache global_highfreq_cache;

/* Global testing variable for dshash usage */
extern bool kmersearch_force_use_parallel_highfreq_kmer_cache;

/* Global SIMD variables */
extern int kmersearch_force_simd_capability;
extern simd_capability_t simd_capability_auto;  /* Auto-detected capability */

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
void check_guc_initialization(void);

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

/* DNA datatype hash functions */
Datum kmersearch_dna2_hash(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_hash(PG_FUNCTION_ARGS);
Datum kmersearch_dna2_hash_extended(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_hash_extended(PG_FUNCTION_ARGS);

/* DNA datatype length functions */
Datum kmersearch_dna2_bit_length(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_bit_length(PG_FUNCTION_ARGS);
Datum kmersearch_dna2_nuc_length(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_nuc_length(PG_FUNCTION_ARGS);
Datum kmersearch_dna2_char_length(PG_FUNCTION_ARGS);
Datum kmersearch_dna4_char_length(PG_FUNCTION_ARGS);

/* DNA datatype utility functions */
char *kmersearch_dna2_to_string(VarBit *dna);
char *kmersearch_dna4_to_string(VarBit *dna);

/* Main dispatch functions with threshold-based SIMD selection */
void dna2_encode(const char* input, uint8_t* output, int len);
void dna2_decode(const uint8_t* input, char* output, int len);
void dna4_encode(const char* input, uint8_t* output, int len);
void dna4_decode(const uint8_t* input, char* output, int len);
int dna_compare(const uint8_t* a, const uint8_t* b, int bit_len);

/* DNA encoding/decoding functions */
void dna2_encode_scalar(const char* input, uint8_t* output, int len);
void dna2_decode_scalar(const uint8_t* input, char* output, int len);
void dna4_encode_scalar(const char* input, uint8_t* output, int len);
void dna4_decode_scalar(const uint8_t* input, char* output, int len);
#ifdef __x86_64__
void dna2_encode_avx2(const char* input, uint8_t* output, int len);
void dna2_encode_avx512(const char* input, uint8_t* output, int len);
void dna2_decode_avx2(const uint8_t* input, char* output, int len);
void dna2_decode_avx512(const uint8_t* input, char* output, int len);
void dna4_encode_avx2(const char* input, uint8_t* output, int len);
void dna4_decode_avx2(const uint8_t* input, char* output, int len);
void dna4_encode_avx512(const char* input, uint8_t* output, int len);
void dna4_decode_avx512(const uint8_t* input, char* output, int len);
#endif
#ifdef __aarch64__
void dna2_encode_neon(const char* input, uint8_t* output, int len);
void dna2_encode_sve(const char* input, uint8_t* output, int len);
void dna2_decode_neon(const uint8_t* input, char* output, int len);
void dna2_decode_sve(const uint8_t* input, char* output, int len);
void dna2_decode_sve2(const uint8_t* input, char* output, int len);
void dna4_encode_neon(const char* input, uint8_t* output, int len);
void dna4_decode_neon(const uint8_t* input, char* output, int len);
void dna4_encode_sve(const char* input, uint8_t* output, int len);
void dna4_decode_sve(const uint8_t* input, char* output, int len);
#endif

/* DNA comparison functions */
int dna_compare_scalar(const uint8_t* a, const uint8_t* b, int bit_len);
#ifdef __x86_64__
int dna_compare_avx2(const uint8_t* a, const uint8_t* b, int bit_len);
int dna_compare_avx512(const uint8_t* a, const uint8_t* b, int bit_len);
#endif
#ifdef __aarch64__
int dna_compare_neon(const uint8_t* a, const uint8_t* b, int bit_len);
int dna_compare_sve(const uint8_t* a, const uint8_t* b, int bit_len);
#endif

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
Datum kmersearch_correctedscore_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_correctedscore_dna4(PG_FUNCTION_ARGS);

/* Cache management functions */
Datum kmersearch_actual_min_score_cache_stats(PG_FUNCTION_ARGS);
Datum kmersearch_actual_min_score_cache_free(PG_FUNCTION_ARGS);
Datum kmersearch_query_pattern_cache_stats(PG_FUNCTION_ARGS);
Datum kmersearch_query_pattern_cache_free(PG_FUNCTION_ARGS);

/* Cache manager functions (defined in kmersearch_cache.c) */
void free_query_pattern_cache_manager(QueryPatternCacheManager **manager);
void free_actual_min_score_cache_manager(ActualMinScoreCacheManager **manager);

/* High-frequency k-mer cache global variables (defined in kmersearch_cache.c) */

/* High-frequency k-mer cache internal functions */
bool kmersearch_validate_guc_against_metadata(Oid table_oid, const char *column_name, int k_value);
HTAB *kmersearch_create_highfreq_hash_from_array(VarBit **kmers, int nkeys);
uint64 kmersearch_ngram_key_to_hash(VarBit *ngram_key);
bool kmersearch_is_global_highfreq_cache_loaded(void);
bool kmersearch_lookup_in_global_cache(VarBit *kmer_key);
bool kmersearch_lookup_kmer2_as_uint_in_global_cache(uint64 kmer2_as_uint, const char *table_name, const char *column_name);
bool kmersearch_lookup_kmer2_as_uint_in_parallel_cache(uint64 kmer2_as_uint, const char *table_name, const char *column_name);
void kmersearch_highfreq_kmer_cache_init(void);
bool kmersearch_highfreq_kmer_cache_load_internal(Oid table_oid, const char *column_name, int k_value);
void kmersearch_highfreq_kmer_cache_free_internal(void);
bool kmersearch_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value);

/* Parallel high-frequency k-mer cache internal functions */
void kmersearch_parallel_highfreq_kmer_cache_init(void);
bool kmersearch_parallel_highfreq_kmer_cache_load_internal(Oid table_oid, const char *column_name, int k_value);
void kmersearch_parallel_highfreq_kmer_cache_free_internal(void);
void kmersearch_parallel_cache_cleanup_on_exit(int code, Datum arg);
bool kmersearch_parallel_cache_lookup(uint64 kmer_hash);

/* High-frequency k-mer filtering functions */
Datum *kmersearch_filter_highfreq_ngram_key2(Datum *original_keys, int *nkeys, HTAB *highfreq_hash, int k);
Datum *kmersearch_filter_highfreq_ngram_key2_parallel(Datum *original_keys, int *nkeys, int k);

/* High-frequency k-mer cache functions */
Datum kmersearch_highfreq_kmer_cache_load(PG_FUNCTION_ARGS);
Datum kmersearch_highfreq_kmer_cache_free(PG_FUNCTION_ARGS);
bool kmersearch_is_kmer_highfreq(VarBit *kmer_key);

/* Parallel cache functions */
Datum kmersearch_parallel_highfreq_kmer_cache_load(PG_FUNCTION_ARGS);
Datum kmersearch_parallel_highfreq_kmer_cache_free(PG_FUNCTION_ARGS);
/* kmersearch_is_highfreq_kmer_parallel now static in kmersearch_gin.c */

/* Analysis functions */
Datum kmersearch_perform_highfreq_analysis(PG_FUNCTION_ARGS);
Datum kmersearch_undo_highfreq_analysis(PG_FUNCTION_ARGS);

/* Analysis dshash functions */
bool kmersearch_is_kmer_hash_in_analysis_dshash(uint64 kmer_hash);

/* K-mer utility functions */
void kmersearch_expand_degenerate_sequence(const char *kmer, int k, char **expanded, int *expand_count);
VarBit *kmersearch_kmer2_as_uint_to_kmer2(uint64 kmer2_as_uint, int kmer_size);
VarBit *kmersearch_create_ngram_key2_from_kmer2_as_uint(uint64 kmer2_as_uint, int kmer_size, int occurrence);
VarBit *kmersearch_remove_occurrence_from_ngram_key2(VarBit *ngram_key2);
void kmersearch_convert_kmer2_to_uint16(VarBit *kmer2, uint16 *result);
void kmersearch_convert_kmer2_to_uint32(VarBit *kmer2, uint32 *result);
void kmersearch_convert_kmer2_to_uint64(VarBit *kmer2, uint64 *result);
Datum *kmersearch_extract_dna2_kmer2_direct(VarBit *dna, int k, int *nkeys);
Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct(VarBit *dna, int k, int *nkeys);
Datum *kmersearch_extract_dna2_ngram_key2_direct(VarBit *dna, int k, int *nkeys);
Datum *kmersearch_extract_dna4_ngram_key2_with_expansion_direct(VarBit *dna, int k, int *nkeys);
int kmersearch_count_degenerate_combinations(const char *kmer, int k);
void kmersearch_set_bit_at(bits8 *data, int bit_pos, int value);
bool kmersearch_will_exceed_degenerate_limit_dna4_bits(VarBit *seq, int start_pos, int k);
VarBit *kmersearch_create_kmer2_key_only(const char *kmer, int k);
VarBit *kmersearch_create_kmer2_key_from_dna2_bits(VarBit *seq, int start_pos, int k);
VarBit *kmersearch_create_ngram_key2_from_dna2_bits(VarBit *seq, int start_pos, int k, int occurrence_count);
VarBit *kmersearch_create_ngram_key2_from_dna4_bits(VarBit *seq, int start_pos, int k, int occurrence_count);
VarBit *kmersearch_create_ngram_key2_with_occurrence_from_dna2(VarBit *dna2_kmer, int k, int occurrence);
VarBit *create_ngram_key2_from_kmer2_and_count(uint64_t kmer2_value, int k_size, int occurrence_count);
VarBit **kmersearch_expand_dna4_kmer2_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count);
uint64_t kmersearch_get_kmer_hash(VarBit *seq, int start_pos, int k);
int kmersearch_find_or_add_kmer_occurrence(KmerOccurrence *occurrences, int *count, uint64_t kmer_value, int max_count);
VarBit **kmersearch_extract_kmer_from_varbit(VarBit *seq, int k, int *nkeys);
VarBit **kmersearch_extract_kmer_from_query(const char *query, int k, int *nkeys);
VarBit **kmersearch_extract_query_ngram_key2(const char *query, int k, int *nkeys);
int kmersearch_count_matching_kmer_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
uint8 kmersearch_get_bit_at(bits8 *data, int bit_pos);
bool kmersearch_will_exceed_degenerate_limit(const char *seq, int len);


/* Parallel analysis functions (implemented in kmersearch.c) */
void kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size, int target_attno, bool is_dna4_type);
void kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, const char *final_table_name, int k_size, int threshold_rows);

/* Parallel worker function for k-mer analysis */
PGDLLEXPORT void kmersearch_analysis_worker(dsm_segment *seg, shm_toc *toc);

/* Frequency analysis persistence function (implemented in kmersearch_freq.c) */
void kmersearch_persist_highfreq_kmers_from_temp(Oid table_oid, const char *column_name, int k_size, const char *temp_table_name);
void create_worker_ngram_temp_table(const char *table_name);
void kmersearch_aggregate_buffer_entries(KmerBuffer *buffer);
void kmersearch_flush_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name);
void kmersearch_add_to_buffer(KmerBuffer *buffer, KmerData kmer_data, const char *temp_table_name);

/* Cache validity check function (implemented in kmersearch_cache.c) */
bool kmersearch_parallel_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value);

/* GIN index support function (implemented in kmersearch_gin.c) */
bool kmersearch_get_index_info(Oid index_oid, Oid *table_oid, char **column_name, int *k_size);
bool evaluate_optimized_match_condition(VarBit **query_keys, int nkeys, int shared_count, const char *query_string, int query_total_kmers);

/* Type OID helper functions */
Oid get_dna2_type_oid(void);
Oid get_dna4_type_oid(void);

/* Cache management functions (implemented in kmersearch_cache.c) */
int calculate_actual_min_score(VarBit **query_keys, int nkeys, int query_total_kmers);

/* Query pattern cache functions (implemented in kmersearch_cache.c) */
VarBit **get_cached_query_kmer(const char *query_string, int k_size, int *nkeys);

/* Actual min score cache functions (implemented in kmersearch_cache.c) */  
int get_cached_actual_min_score(VarBit **query_keys, int nkeys);
int get_cached_actual_min_score_or_error(VarBit **query_keys, int nkeys);

/* High-frequency k-mer filtering functions (implemented in kmersearch_gin.c) */
VarBit **filter_ngram_key2_and_set_actual_min_score(VarBit **query_keys, int *nkeys, const char *query_string);

/* Cache functions (implemented in kmersearch_cache.c) */
void kmersearch_query_pattern_cache_max_entries_assign_hook(int newval, void *extra);

/* Internal functions that should be declared (implemented in kmersearch_freq.c) */
bool kmersearch_is_highfreq_filtering_enabled(void);
int kmersearch_count_highfreq_kmer_in_query(VarBit **query_keys, int nkeys);

/* Cache key validation functions (implemented in kmersearch_cache.c) */
bool kmersearch_validate_cache_key_match(Oid table_oid, const char *column_name);
bool kmersearch_validate_parallel_cache_key_match(Oid table_oid, const char *column_name);
bool kmersearch_validate_guc_against_all_metadata(void);
bool kmersearch_is_parallel_highfreq_cache_loaded(void);
bool kmersearch_lookup_in_parallel_cache(VarBit *kmer_key);

/* GIN index support functions (implemented in kmersearch_gin.c) */
Datum kmersearch_extract_value_dna2(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna4(PG_FUNCTION_ARGS);
Datum kmersearch_extract_query(PG_FUNCTION_ARGS);
Datum kmersearch_consistent(PG_FUNCTION_ARGS);
Datum kmersearch_compare_partial(PG_FUNCTION_ARGS);

/* Frequency analysis functions (implemented in kmersearch_freq.c) */
Datum kmersearch_perform_highfreq_analysis(PG_FUNCTION_ARGS);
Datum kmersearch_undo_highfreq_analysis(PG_FUNCTION_ARGS);

/* Internal frequency analysis functions (implemented in kmersearch_freq.c) */
DropAnalysisResult kmersearch_undo_highfreq_analysis_internal(Oid table_oid, const char *column_name, int k_size);
KmerAnalysisResult kmersearch_perform_highfreq_analysis_parallel(Oid table_oid, const char *column_name, int k_size, int parallel_workers);
void kmersearch_validate_analysis_parameters(Oid table_oid, const char *column_name, int k_size);

/* Utility functions for unique temporary table name generation */
char *kmersearch_generate_unique_temp_table_name(const char *prefix, int additional_id);

/* Native uint k-mer extraction functions */
void kmersearch_extract_dna2_kmer2_as_uint_direct(VarBit *seq, int k, void **output, int *nkeys);
void kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(VarBit *seq, int k, void **output, int *nkeys);

/* Direct bit manipulation helper functions */
size_t kmersearch_get_kmer_uint_size(int k);
void kmersearch_extract_dna4_kmer_expansions_direct_bits(VarBit *seq, int start_pos, int k, uint64 *output, int *count);
void *kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count);

/* Utility functions */
void kmersearch_spi_connect_or_error(void);

#endif   /* KMERSEARCH_H */