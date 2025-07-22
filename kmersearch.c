#include "kmersearch.h"

/*
 * Generate a unique temporary table name with unified format
 * prefix: table name prefix (e.g., "temp_kmer_worker", "temp_kmer_final")
 * additional_id: additional identifier (e.g., worker_id, 0 for no additional id)
 */
char *
kmersearch_generate_unique_temp_table_name(const char *prefix, int additional_id)
{
    static int call_counter = 0;
    struct timeval tv;
    
    gettimeofday(&tv, NULL);
    
    if (additional_id >= 0)
    {
        return psprintf("%s_%d_%d_%ld_%ld_%d", 
                       prefix, getpid(), additional_id, tv.tv_sec, tv.tv_usec, ++call_counter);
    }
    else
    {
        return psprintf("%s_%d_%ld_%ld_%d", 
                       prefix, getpid(), tv.tv_sec, tv.tv_usec, ++call_counter);
    }
}

PG_MODULE_MAGIC;

/* Global SIMD dispatch table */
simd_dispatch_table_t simd_dispatch;
simd_capability_t simd_capability = SIMD_NONE;

/* Global variables for k-mer search configuration */
int kmersearch_occur_bitlen = 8;  /* Default 8 bits for occurrence count */
int kmersearch_kmer_size = 16;  /* Default k-mer size */
double kmersearch_max_appearance_rate = 0.5;  /* Default max appearance rate */
int kmersearch_max_appearance_nrow = 0;  /* Default max appearance nrow (0 = undefined) */
int kmersearch_min_score = 1;  /* Default minimum score for GIN search */
double kmersearch_min_shared_ngram_key_rate = 0.9;  /* Default minimum shared n-gram key rate for =% operator */
bool kmersearch_preclude_highfreq_kmer = false;  /* Default to not exclude high-frequency k-mers */

/* Cache configuration variables */
int kmersearch_rawscore_cache_max_entries = 50000;  /* Default max rawscore cache entries */
int kmersearch_query_pattern_cache_max_entries = 50000;  /* Default max query pattern cache entries */
int kmersearch_actual_min_score_cache_max_entries = 50000;  /* Default max actual min score cache entries */
int kmersearch_highfreq_kmer_cache_load_batch_size = 10000;  /* Default batch size for loading high-frequency k-mers */

/* Global cache managers */
ActualMinScoreCacheManager *actual_min_score_cache_manager = NULL;

/* Global query pattern cache manager for cross-query sharing */
QueryPatternCacheManager *query_pattern_cache_manager = NULL;

/* Global rawscore cache manager for cross-query sharing */
RawscoreCacheManager *rawscore_cache_manager = NULL;


/* Macro for safe memory cleanup */
#define CLEANUP_KMER_ARRAYS(seq_keys, seq_nkeys, query_keys, query_nkeys) \
    do { \
        if (seq_keys) { \
            int cleanup_i; \
            for (cleanup_i = 0; cleanup_i < (seq_nkeys); cleanup_i++) { \
                if (seq_keys[cleanup_i]) \
                    pfree(seq_keys[cleanup_i]); \
            } \
            pfree(seq_keys); \
            seq_keys = NULL; \
        } \
        if (query_keys) { \
            int cleanup_j; \
            for (cleanup_j = 0; cleanup_j < (query_nkeys); cleanup_j++) { \
                if (query_keys[cleanup_j]) \
                    pfree(query_keys[cleanup_j]); \
            } \
            pfree(query_keys); \
            query_keys = NULL; \
        } \
    } while(0)

/* Forward declarations */
/* Functions moved to other modules - declarations remain for compatibility */

/* Helper function declarations */
Oid get_dna2_type_oid(void);
Oid get_dna4_type_oid(void);
static bool tuple_in_worker_range(HeapTuple tuple, KmerWorkerState *worker);
static VarBit *extract_sequence_from_tuple(HeapTuple tuple, int attno, TupleDesc tupdesc);
static void create_worker_ngram_temp_table(const char *table_name);
static bool is_kmer2_in_highfreq_table(uint64_t kmer2_value, const char *highfreq_table_name);
static bool is_kmer2_in_analysis_dshash(uint64_t kmer2_value);
static void process_extracted_kmer2(uint64_t *kmer2_values, int nkeys, VarBit *sequence, int k_size, const char *worker_table, const char *highfreq_table);
static VarBit* create_ngram_key2_from_kmer2_and_count(uint64_t kmer2_value, int k_size, int occurrence_count);
static int get_processed_row_count(KmerWorkerState *worker);

/* Parallel high-frequency k-mer cache management functions */
static bool kmersearch_parallel_highfreq_kmer_cache_is_valid(Oid table_oid, const char *column_name, int k_value);

/* Analysis dshash functions (defined in kmersearch_freq.c) */
/* B-2: Other functions */
Datum *kmersearch_extract_dna2_kmer2_direct(VarBit *seq, int k, int *nkeys);
Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct(VarBit *seq, int k, int *nkeys);
Datum *kmersearch_extract_dna2_ngram_key2_direct(VarBit *seq, int k, int *nkeys);
Datum *kmersearch_extract_dna4_ngram_key2_with_expansion_direct(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmer_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
static bool kmersearch_kmer_based_match_dna2(VarBit *sequence, const char *query_string);
static bool kmersearch_kmer_based_match_dna4(VarBit *sequence, const char *query_string);
static bool kmersearch_evaluate_match_conditions(int shared_count, int query_total);

/* Actual min score cache functions */
static bool evaluate_optimized_match_condition(VarBit **query_keys, int nkeys, int shared_count, const char *query_string, int query_total_kmers);

/* New parallel analysis functions */
void kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size, int target_attno, bool is_dna4_type);
void kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, const char *final_table_name, int k_size, int threshold_rows);
static void kmersearch_persist_highfreq_kmers_from_temp(Oid table_oid, const char *column_name, int k_size, const char *temp_table_name);

/* New memory-efficient k-mer functions */
static size_t kmersearch_get_kmer_data_size(int k_size);
static bool kmersearch_get_index_info(Oid index_oid, Oid *table_oid, char **column_name, int *k_size);
static int kmersearch_calculate_buffer_size(int k_size);
static void kmersearch_init_buffer(KmerBuffer *buffer, int k_size);
static void kmersearch_add_to_buffer(KmerBuffer *buffer, KmerData kmer_data, const char *temp_table_name);
static void kmersearch_add_hash_to_buffer(KmerBuffer *buffer, uint64_t kmer_hash, const char *temp_table_name);
static void kmersearch_flush_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name);
static void kmersearch_flush_hash_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name);
static void kmersearch_aggregate_buffer_entries(KmerBuffer *buffer);
static void kmersearch_create_worker_temp_table(const char *temp_table_name, int k_size);

/* Custom GUC variables */
void _PG_init(void);

/*
 * DNA2 type: 2-bit encoding for ACGT
 * A=00, C=01, G=10, T=11
 * Definition moved to kmersearch.h
 */

/*
 * DNA4 type: 4-bit encoding with degenerate codes
 * Uses standard varbit structure
 * Definition moved to kmersearch.h
 */

/* DNA2 encoding table */
const uint8 kmersearch_dna2_encode_table[256] = {
    ['A'] = 0, ['a'] = 0,
    ['C'] = 1, ['c'] = 1,
    ['G'] = 2, ['g'] = 2,
    ['T'] = 3, ['t'] = 3,
    ['U'] = 3, ['u'] = 3,  /* U is treated as T */
};

/* DNA2 decoding table */
const char kmersearch_dna2_decode_table[4] = {'A', 'C', 'G', 'T'};

/* DNA4 encoding table */
const uint8 kmersearch_dna4_encode_table[256] = {
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
const char kmersearch_dna4_decode_table[16] = {
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
/* K-mer and GIN index functions */
PG_FUNCTION_INFO_V1(kmersearch_dna2_match);
PG_FUNCTION_INFO_V1(kmersearch_dna4_match);

/* K-mer frequency analysis functions */

/* Score calculation functions */
PG_FUNCTION_INFO_V1(kmersearch_rawscore_dna2);
PG_FUNCTION_INFO_V1(kmersearch_rawscore_dna4);
PG_FUNCTION_INFO_V1(kmersearch_correctedscore_dna2);
PG_FUNCTION_INFO_V1(kmersearch_correctedscore_dna4);
/* SIMD capability detection functions */
static simd_capability_t detect_cpu_capabilities(void);
static void init_simd_dispatch_table(void);

/* SIMD implementation functions */
static void dna2_encode_scalar(const char* input, uint8_t* output, int len);
static void dna2_decode_scalar(const uint8_t* input, char* output, int len);
static void dna4_encode_scalar(const char* input, uint8_t* output, int len);
static void dna4_decode_scalar(const uint8_t* input, char* output, int len);

/* Scalar versions */
static Datum *kmersearch_extract_dna2_kmer2_direct_scalar(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_scalar(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmer_fast_scalar_simple(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
static int kmersearch_count_matching_kmer_fast_scalar_hashtable(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

#ifdef __x86_64__
static void dna2_encode_avx2(const char* input, uint8_t* output, int len);
static void dna2_decode_avx2(const uint8_t* input, char* output, int len);
static void dna4_encode_avx2(const char* input, uint8_t* output, int len);
static void dna4_decode_avx2(const uint8_t* input, char* output, int len);

/* K-mer processing functions with SIMD optimization */
static Datum *kmersearch_extract_dna2_kmer2_direct_avx2(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmer_fast_avx2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

static Datum *kmersearch_extract_dna2_kmer2_direct_avx512(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmer_fast_avx512(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

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

static Datum *kmersearch_extract_dna2_kmer2_direct_neon(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_neon(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmer_fast_neon(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

static Datum *kmersearch_extract_dna2_kmer2_direct_sve(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_sve(VarBit *seq, int k, int *nkeys);
static int kmersearch_count_matching_kmer_fast_sve(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);

static void dna2_encode_sve(const char* input, uint8_t* output, int len);
static void dna2_decode_sve(const uint8_t* input, char* output, int len);
static void dna4_encode_sve(const char* input, uint8_t* output, int len);
static void dna4_decode_sve(const uint8_t* input, char* output, int len);
#endif

/*
 * GUC assign hook functions for cache invalidation
 */

/* High-frequency cache clear with warning */
static void 
clear_highfreq_cache_with_warning(void)
{
    bool had_valid_cache = global_highfreq_cache.is_valid;
    
    kmersearch_highfreq_kmer_cache_free_internal();
    
    /* Only show warning if cache was actually valid and cleared */
    if (had_valid_cache)
    {
        elog(WARNING, "High-frequency k-mer cache has been cleared. "
                      "You may need to manually execute kmersearch_highfreq_kmer_cache_load() "
                      "to reload the cache if needed.");
    }
}

/* K-mer size change affects all caches */
static void
kmersearch_kmer_size_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    
    /* Clear query pattern cache */
    if (query_pattern_cache_manager)
        kmersearch_free_query_pattern_cache_internal();
    
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        kmersearch_free_actual_min_score_cache_internal();
    
    /* Clear high-frequency k-mer cache with conditional warning */
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
        kmersearch_free_actual_min_score_cache_internal();
    
    /* Clear high-frequency k-mer cache with conditional warning */
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
        kmersearch_free_actual_min_score_cache_internal();
    
    /* Clear high-frequency k-mer cache with conditional warning */
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
        kmersearch_free_actual_min_score_cache_internal();
}

/* Min shared ngram key rate change affects actual min score cache */
static void
kmersearch_min_shared_ngram_key_rate_assign_hook(double newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        kmersearch_free_actual_min_score_cache_internal();
}
/* Query pattern cache max entries change requires cache recreation */
/* Function moved to kmersearch_cache.c */

/* Occurrence bit length change affects rawscore and high-freq caches */
static void
kmersearch_occur_bitlen_assign_hook(int newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    
    /* Clear high-frequency k-mer cache with conditional warning */
    clear_highfreq_cache_with_warning();
}

/*
 * Module initialization
 */
/* Global flag to prevent duplicate GUC initialization */
static bool guc_variables_initialized = false;

/*
 * Check if GUC variables are properly initialized
 * If not, report an error with instructions for shared_preload_libraries
 */
void
check_guc_initialization(void)
{
    if (!guc_variables_initialized) {
        ereport(ERROR,
                (errcode(ERRCODE_OBJECT_NOT_IN_PREREQUISITE_STATE),
                 errmsg("pg_kmersearch extension not properly initialized"),
                 errhint("Add 'pg_kmersearch' to shared_preload_libraries in postgresql.conf and restart PostgreSQL.")));
    }
}

void
_PG_init(void)
{
    /* Prevent duplicate initialization */
    if (guc_variables_initialized) {
        return;
    }
    
    /* Initialize SIMD capabilities */
    simd_capability = detect_cpu_capabilities();
    init_simd_dispatch_table();
    
    /* Define custom GUC variables */
    DefineCustomRealVariable("kmersearch.max_appearance_rate",
                            "Maximum appearance rate for k-mers to be included in index",
                            "K-mers appearing in more than this fraction of rows will be identified as highly frequent",
                            &kmersearch_max_appearance_rate,
                            0.5,
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
                           "Length of k-mer sequences for similarity matching (4-32)",
                           &kmersearch_kmer_size,
                           16,
                           4,
                           32,
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

    DefineCustomBoolVariable("kmersearch.preclude_highfreq_kmer",
                            "Enable high-frequency k-mer exclusion during GIN index construction",
                            "When enabled, high-frequency k-mers will be excluded from GIN index to improve performance",
                            &kmersearch_preclude_highfreq_kmer,
                            false,
                            PGC_USERSET,
                            0,
                            NULL,
                            NULL,
                            NULL);

    DefineCustomBoolVariable("kmersearch.force_use_parallel_highfreq_kmer_cache",
                            "Force use of dshash-based parallel cache (for testing)",
                            "When enabled, forces the use of parallel high-frequency k-mer cache even for main processes",
                            &kmersearch_force_use_parallel_highfreq_kmer_cache,
                            false,
                            PGC_USERSET,
                            0,
                            NULL,
                            NULL,
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

    DefineCustomIntVariable("kmersearch.highfreq_kmer_cache_load_batch_size",
                           "Batch size for loading high-frequency k-mers into cache",
                           "Controls the number of k-mers loaded in each batch to reduce memory usage",
                           &kmersearch_highfreq_kmer_cache_load_batch_size,
                           10000,
                           1000,
                           1000000,
                           PGC_USERSET,
                           0,
                           NULL,
                           NULL,
                           NULL);
    
    /* Initialize high-frequency k-mer cache */
    kmersearch_highfreq_kmer_cache_init();
    
    /* Mark GUC variables as initialized */
    guc_variables_initialized = true;
}

/*
 * CPU capability detection
 */
#ifdef __aarch64__
static sigjmp_buf jmpbuf;

static void sigill_handler(int sig) {
    siglongjmp(jmpbuf, 1);
}
#endif

static simd_capability_t detect_cpu_capabilities(void)
{
#ifdef __x86_64__
    unsigned int eax, ebx, ecx, edx;
    
    /* Check for AVX512 support */
    if (__get_cpuid_max(0, NULL) >= 7) {
        __cpuid_count(7, 0, eax, ebx, ecx, edx);
        if (ebx & (1 << 16)) { /* AVX512F */
            if (ebx & (1 << 30)) { /* AVX512BW */
                return SIMD_AVX512BW;
            }
            return SIMD_AVX512F;
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
    struct sigaction sa, old_sa;
    sa.sa_handler = sigill_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGILL, &sa, &old_sa);
    
     /* Check for SVE support */
    if (sigsetjmp(jmpbuf, 1) == 0) {
        size_t vl = svcntb();
        (void)vl;
        sigaction(SIGILL, &old_sa, NULL);
        return SIMD_SVE;
    }
    
    /* ARM64 always has NEON */
    if (sigsetjmp(jmpbuf, 1) == 0) {
        volatile uint8x8_t a = vdup_n_u8(1);
        (void)a;
        sigaction(SIGILL, &old_sa, NULL);
        return SIMD_NEON;
    }
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
    /* Enable only AVX2 implementations for testing */
    switch (simd_capability) {
#ifdef __x86_64__
        case SIMD_AVX512BW:
            simd_dispatch.dna2_encode = dna2_encode_avx512;
            simd_dispatch.dna2_decode = dna2_decode_avx512;
            simd_dispatch.dna4_encode = dna4_encode_avx512;
            simd_dispatch.dna4_decode = dna4_decode_avx512;
            break;
        case SIMD_AVX512F:
            /* AVX512F without BW - use AVX2 fallback for encode/decode */
            simd_dispatch.dna2_encode = dna2_encode_avx2;
            simd_dispatch.dna2_decode = dna2_decode_avx2;
            simd_dispatch.dna4_encode = dna4_encode_avx2;
            simd_dispatch.dna4_decode = dna4_decode_avx2;
            break;
        case SIMD_AVX2:
            simd_dispatch.dna2_encode = dna2_encode_avx2;
            simd_dispatch.dna2_decode = dna2_decode_avx2;
            simd_dispatch.dna4_encode = dna4_encode_avx2;
            simd_dispatch.dna4_decode = dna4_decode_avx2;
            break;
#elif defined(__aarch64__)
        case SIMD_SVE:
            /* ENABLED: Fixed ARM64 SVE implementations */
            simd_dispatch.dna2_encode = dna2_encode_sve;
            simd_dispatch.dna2_decode = dna2_decode_sve;
            simd_dispatch.dna4_encode = dna4_encode_sve;
            simd_dispatch.dna4_decode = dna4_decode_sve;
            break;
        case SIMD_NEON:
            /* ENABLED: Fixed ARM64 NEON implementations */
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
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using bit operations
 */
VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[32][4];  /* Max k=32, max 4 expansions per base */
    int base_counts[32];
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
 * Extract k-mers directly from DNA2 bit sequence (with SIMD dispatch)
 */
Datum *
kmersearch_extract_dna2_ngram_key2_direct(VarBit *seq, int k, int *nkeys)
{
    Datum *kmer2_keys;
    int kmer2_count;
    Datum *ngram_keys;
    int ngram_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    /* First, extract all kmer2 keys without occurrence count */
    kmer2_keys = kmersearch_extract_dna2_kmer2_direct(seq, k, &kmer2_count);
    if (!kmer2_keys || kmer2_count == 0) {
        *nkeys = 0;
        return NULL;
    }
    
    /* Allocate arrays for occurrence tracking and result */
    occurrences = (KmerOccurrence *) palloc(kmer2_count * sizeof(KmerOccurrence));
    ngram_keys = (Datum *) palloc(kmer2_count * sizeof(Datum));
    
    /* Process each kmer2 key and add occurrence count */
    for (i = 0; i < kmer2_count; i++) {
        VarBit *kmer2_key = (VarBit *) DatumGetPointer(kmer2_keys[i]);
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        /* Extract k-mer hash value from kmer2 key */
        kmer_value = kmersearch_get_kmer_hash(kmer2_key, 0, k);
        
        /* Find or add occurrence count using binary search */
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, kmer2_count);
        
        if (current_count < 0)
            continue;  /* Array full, skip */
        
        /* Skip if occurrence exceeds bit limit */
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* Create n-gram key (k-mer + occurrence count) from kmer2 key */
        ngram_key = kmersearch_create_ngram_key2_from_dna2_bits(kmer2_key, 0, k, current_count);
        if (ngram_key == NULL)
            continue;  /* Skip if key creation failed */
            
        ngram_keys[ngram_count++] = PointerGetDatum(ngram_key);
    }
    
    /* Cleanup */
    for (i = 0; i < kmer2_count; i++) {
        pfree(DatumGetPointer(kmer2_keys[i]));
    }
    pfree(kmer2_keys);
    pfree(occurrences);
    
    *nkeys = ngram_count;
    return ngram_keys;
}

/*
 * Extract k-mers directly from DNA2 bit sequence (kmer2 output without occurrence count)
 */
Datum *
kmersearch_extract_dna2_kmer2_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_avx512(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_AVX2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_avx2(seq, k, nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_sve(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_neon(seq, k, nkeys);
    }
#endif
    return kmersearch_extract_dna2_kmer2_direct_scalar(seq, k, nkeys);
}

/*
 * Scalar version: Extract k-mers directly from DNA2 bit sequence
 */
static Datum *
kmersearch_extract_dna2_kmer2_direct_scalar(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Extract k-mers without occurrence count (all k-mers, no deduplication) */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit *kmer_key;
        
        /* Check bounds before extraction to avoid invalid access */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;  /* Out of bounds, skip */
        }
        
        /* Create k-mer key (without occurrence count) */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Extract k-mers only (without occurrence counts) from DNA2 bit sequence
 * This function is used for frequency analysis phase
 */

/*
 * Encode k-mer-only VarBit into compact KmerData (ignoring occurrence count bits)
 */

/*
 * Extract k-mers directly from DNA4 bit sequence with degenerate expansion (with SIMD dispatch)
 */
Datum *
kmersearch_extract_dna4_ngram_key2_with_expansion_direct(VarBit *seq, int k, int *nkeys)
{
    Datum *kmer2_keys;
    int kmer2_count;
    Datum *ngram_keys;
    int ngram_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    /* First, extract all kmer2 keys without occurrence count */
    kmer2_keys = kmersearch_extract_dna4_kmer2_with_expansion_direct(seq, k, &kmer2_count);
    if (!kmer2_keys || kmer2_count == 0) {
        *nkeys = 0;
        return NULL;
    }
    
    /* Allocate arrays for occurrence tracking and result */
    occurrences = (KmerOccurrence *) palloc(kmer2_count * sizeof(KmerOccurrence));
    ngram_keys = (Datum *) palloc(kmer2_count * sizeof(Datum));
    
    /* Process each kmer2 key and add occurrence count */
    for (i = 0; i < kmer2_count; i++) {
        VarBit *kmer2_key = (VarBit *) DatumGetPointer(kmer2_keys[i]);
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        /* Extract k-mer hash value from kmer2 key */
        kmer_value = kmersearch_get_kmer_hash(kmer2_key, 0, k);
        
        /* Find or add occurrence count using binary search */
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, kmer2_count);
        
        if (current_count < 0)
            continue;  /* Array full, skip */
        
        /* Skip if occurrence exceeds bit limit */
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* Create n-gram key (k-mer + occurrence count) from kmer2 key */
        ngram_key = kmersearch_create_ngram_key2_from_dna2_bits(kmer2_key, 0, k, current_count);
        if (ngram_key == NULL)
            continue;  /* Skip if key creation failed */
            
        ngram_keys[ngram_count++] = PointerGetDatum(ngram_key);
    }
    
    /* Cleanup */
    for (i = 0; i < kmer2_count; i++) {
        pfree(DatumGetPointer(kmer2_keys[i]));
    }
    pfree(kmer2_keys);
    pfree(occurrences);
    
    *nkeys = ngram_count;
    return ngram_keys;
}

/*
 * Extract k-mers directly from DNA4 bit sequence with degenerate expansion (kmer2 output without occurrence count)
 */
Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_AVX2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2(seq, k, nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_sve(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_neon(seq, k, nkeys);
    }
#endif
    return kmersearch_extract_dna4_kmer2_with_expansion_direct_scalar(seq, k, nkeys);
}

/*
 * Scalar version: Extract k-mers directly from DNA4 bit sequence with degenerate expansion
 */
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_scalar(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Extract k-mers without occurrence count (all k-mers, no deduplication) */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Cache management functions
 */
/* Cache management functions moved to kmersearch_cache.c */

/*
 * Fast k-mer matching using hash table - optimized O(n+m) implementation (with SIMD dispatch)
 */
static int
kmersearch_count_matching_kmer_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int key_combinations;
    
    /* Validate input parameters */
    if (!seq_keys || !query_keys) {
        return 0;
    }
    
    if (seq_nkeys <= 0 || query_nkeys <= 0) {
        return 0;
    }
    
    key_combinations = seq_nkeys * query_nkeys;
    
    /* For small datasets, O(n*m) might be faster than hash table overhead */
    if (key_combinations < 100) {
        return kmersearch_count_matching_kmer_fast_scalar_simple(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    
    /* Use SIMD based on runtime capability and key combination thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && key_combinations >= SIMD_KEYCOMB_AVX512_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_avx512(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    if (simd_capability >= SIMD_AVX2 && key_combinations >= SIMD_KEYCOMB_AVX2_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_avx2(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && key_combinations >= SIMD_KEYCOMB_SVE_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_sve(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    if (simd_capability >= SIMD_NEON && key_combinations >= SIMD_KEYCOMB_NEON_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_neon(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
#endif
    return kmersearch_count_matching_kmer_fast_scalar_hashtable(seq_keys, seq_nkeys, query_keys, query_nkeys);
}

/*
 * Scalar version: Simple O(n*m) comparison for small datasets
 */
static int
kmersearch_count_matching_kmer_fast_scalar_simple(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i, j;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* Simple comparison for small datasets */
    for (i = 0; i < seq_nkeys; i++)
    {
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

/*
 * Scalar version: Fast k-mer matching using hash table - optimized O(n+m) implementation
 */
static int
kmersearch_count_matching_kmer_fast_scalar_hashtable(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* Create hash table using VarBit content as key */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    /* Safety check: ensure we have valid query keys */
    if (query_keys[0] == NULL) {
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);  /* Use data size, not total size */
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHash", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* Insert all query k-mers into hash table using content as key */
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Check each sequence k-mer against hash table */
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
    
    /* Cleanup */
    hash_destroy(query_hash);
    
    return match_count;
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
 * Corrected score functions - return unfiltered similarity scores
 * 
 * These functions calculate similarity scores without applying high-frequency k-mer exclusion.
 * Unlike rawscore functions which exclude high-frequency k-mers when analysis data is available,
 * correctedscore functions count all shared k-mers between sequence and query without filtering.
 * This provides the uncorrected baseline score for comparison with filtered rawscore results.
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
    seq_datum_keys = kmersearch_extract_dna2_ngram_key2_direct(sequence, k, &seq_nkeys);
    if (seq_datum_keys != NULL && seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    /* Extract k-mers from query as ngram_key2 format */
    query_keys = kmersearch_extract_query_ngram_key2(query_string, k, &query_nkeys);
    
    /* Count shared k-mers using optimized function */
    if (seq_keys && query_keys && seq_nkeys > 0 && query_nkeys > 0) {
        shared_count = kmersearch_count_matching_kmer_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
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
    seq_datum_keys = kmersearch_extract_dna4_ngram_key2_with_expansion_direct(sequence, k, &seq_nkeys);
    if (seq_datum_keys != NULL && seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    /* Extract k-mers from query as ngram_key2 format */
    query_keys = kmersearch_extract_query_ngram_key2(query_string, k, &query_nkeys);
    
    /* Count shared k-mers using optimized function */
    if (seq_keys && query_keys && seq_nkeys > 0 && query_nkeys > 0) {
        shared_count = kmersearch_count_matching_kmer_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
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
/* High-frequency k-mer validation and parallel cache functions moved to kmersearch_freq.c */

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
 * Optimized match condition evaluation using pre-calculated actual minimum score
 */
static bool
evaluate_optimized_match_condition(VarBit **query_keys, int nkeys, int shared_count, const char *query_string, int query_total_kmers)
{
    int actual_min_score;
    
    
    /* Get cached actual min score (with TopMemoryContext caching for performance) */
    actual_min_score = get_cached_actual_min_score(query_keys, nkeys);
    elog(LOG, "evaluate_optimized_match_condition: get_cached_actual_min_score returned %d", actual_min_score);
    
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
    seq_keys = (VarBit **)kmersearch_extract_dna2_ngram_key2_direct(sequence, k, &seq_nkeys);
    if (seq_keys == NULL || seq_nkeys == 0) {
        return false;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    query_keys = get_cached_query_kmer(query_string, k, &query_nkeys);
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
    shared_count = kmersearch_count_matching_kmer_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
    
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
    seq_keys = (VarBit **)kmersearch_extract_dna4_ngram_key2_with_expansion_direct(sequence, k, &seq_nkeys);
    if (seq_keys == NULL || seq_nkeys == 0) {
        return false;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    query_keys = get_cached_query_kmer(query_string, k, &query_nkeys);
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
    shared_count = kmersearch_count_matching_kmer_fast(seq_keys, seq_nkeys, query_keys, query_nkeys);
    
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
 * Helper function to clean up kmer match resources
 */
static void
cleanup_kmer_match_resources(VarBit **seq_keys, int seq_nkeys, 
                             VarBit **query_keys, int query_nkeys,
                             Datum *seq_datum_keys)
{
    CLEANUP_KMER_ARRAYS(seq_keys, seq_nkeys, query_keys, query_nkeys);
    if (seq_datum_keys) {
        pfree(seq_datum_keys);
    }
}

/*
 * Core k-mer matching and scoring function for DNA2 sequences
 * Performs all k-mer extraction, comparison, and evaluation in one pass
 */
KmerMatchResult
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
        ereport(ERROR, (errmsg("Query sequence must be at least %d bases long", k)));
    }
    
    /* Extract k-mers from DNA2 sequence (no degenerate expansion) */
    elog(LOG, "DNA2 Cache: Starting k-mer extraction from sequence");
    seq_datum_keys = kmersearch_extract_dna2_ngram_key2_direct(sequence, k, &result.seq_nkeys);
    elog(LOG, "DNA2 Cache: Extracted %d k-mers from sequence", result.seq_nkeys);
    
    if (seq_datum_keys != NULL && result.seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(result.seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < result.seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
        elog(LOG, "DNA2 Cache: Converted %d datum keys to VarBit", result.seq_nkeys);
    }
    
    if (seq_keys != NULL && result.seq_nkeys > 0) {
        /* Extract k-mers from query (with degenerate expansion) */
        elog(LOG, "DNA2 Cache: Starting k-mer extraction from query '%s'", query_string);
        query_keys = get_cached_query_kmer(query_string, k, &result.query_nkeys);
        elog(LOG, "DNA2 Cache: Extracted %d k-mers from query", result.query_nkeys);
        
        if (query_keys != NULL && result.query_nkeys > 0) {
            /* Calculate shared k-mer count (this becomes the rawscore) */
            elog(LOG, "DNA2 Cache: Starting k-mer matching calculation");
            result.shared_count = kmersearch_count_matching_kmer_fast(seq_keys, result.seq_nkeys, 
                                                                       query_keys, result.query_nkeys);
            elog(LOG, "DNA2 Cache: Completed k-mer matching, shared_count=%d", result.shared_count);
            
            /* Calculate sharing rate */
            if (result.query_nkeys > 0) {
                result.sharing_rate = (double)result.shared_count / (double)result.query_nkeys;
            }
            
            /* Evaluate match conditions for =% operator using optimized method */
            result.match_result = evaluate_optimized_match_condition(query_keys, result.query_nkeys, result.shared_count, query_string, result.query_nkeys);
            
            result.valid = true;
        } else {
            elog(LOG, "DNA2 Cache: No query k-mers extracted, cleaning up");
        }
    } else {
        elog(LOG, "DNA2 Cache: No sequence k-mers extracted, cleaning up");
    }
    
    /* Unified memory cleanup */
    cleanup_kmer_match_resources(seq_keys, result.seq_nkeys, query_keys, result.query_nkeys, seq_datum_keys);
    
    return result;
}

/*
 * Core k-mer matching and scoring function for DNA4 sequences
 * Performs all k-mer extraction, comparison, and evaluation in one pass
 */
KmerMatchResult
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
        ereport(ERROR, (errmsg("Query sequence must be at least %d bases long", k)));
    }
    
    /* Extract k-mers from DNA4 sequence (with degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna4_ngram_key2_with_expansion_direct(sequence, k, &result.seq_nkeys);
    if (seq_datum_keys != NULL && result.seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(result.seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < result.seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    
    if (seq_keys != NULL && result.seq_nkeys > 0) {
        /* Extract k-mers from query (with degenerate expansion) */
        query_keys = get_cached_query_kmer(query_string, k, &result.query_nkeys);
        
        if (query_keys != NULL && result.query_nkeys > 0) {
            /* Calculate shared k-mer count (this becomes the rawscore) */
            result.shared_count = kmersearch_count_matching_kmer_fast(seq_keys, result.seq_nkeys, 
                                                                       query_keys, result.query_nkeys);
            
            /* Calculate sharing rate */
            if (result.query_nkeys > 0) {
                result.sharing_rate = (double)result.shared_count / (double)result.query_nkeys;
            }
            
            /* Evaluate match conditions for =% operator using optimized method */
            result.match_result = evaluate_optimized_match_condition(query_keys, result.query_nkeys, result.shared_count, query_string, result.query_nkeys);
            
            result.valid = true;
        }
    }
    
    /* Unified memory cleanup */
    cleanup_kmer_match_resources(seq_keys, result.seq_nkeys, query_keys, result.query_nkeys, seq_datum_keys);
    
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
    else {
        /* k > 32 not supported */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k_size)));
    }
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
        "SELECT table_oid, column_name, kmer_size FROM kmersearch_index_info "
        "WHERE index_oid = %u",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 1);
    if (ret != SPI_OK_SELECT)
        ereport(ERROR, (errmsg("SPI_execute failed with code %d", ret)));
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        Datum table_oid_datum, column_name_datum, k_size_datum;
        bool isnull;
        
        table_oid_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull && table_oid)
            *table_oid = DatumGetObjectId(table_oid_datum);
        
        column_name_datum = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2, &isnull);
        (void) column_name_datum;  /* Suppress unused variable warning */
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
 * Initialize k-mer buffer
 */
static void
kmersearch_init_buffer(KmerBuffer *buffer, int k_size)
{
    buffer->capacity = kmersearch_calculate_buffer_size(k_size);
    buffer->entries = (CompactKmerFreq *) palloc0(buffer->capacity * sizeof(CompactKmerFreq));
    buffer->count = 0;
    buffer->kmer_size = k_size;
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
            
            if (buffer->kmer_size <= 8) {
                same_kmer = (buffer->entries[i].kmer_data.k8_data == buffer->entries[j].kmer_data.k8_data);
            } else if (buffer->kmer_size <= 16) {
                same_kmer = (buffer->entries[i].kmer_data.k16_data == buffer->entries[j].kmer_data.k16_data);
            } else if (buffer->kmer_size <= 32) {
                same_kmer = (buffer->entries[i].kmer_data.k32_data == buffer->entries[j].kmer_data.k32_data);
            } else {
                /* k > 32 not supported */
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                         errmsg("k-mer length must be between 4 and 32"),
                         errdetail("Provided k-mer length: %d", buffer->kmer_size)));
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
        
        if (buffer->kmer_size <= 8) {
            appendStringInfo(&query, "(%u, %d)", 
                           buffer->entries[i].kmer_data.k8_data,
                           buffer->entries[i].frequency_count);
        } else if (buffer->kmer_size <= 16) {
            appendStringInfo(&query, "(%u, %d)", 
                           buffer->entries[i].kmer_data.k16_data,
                           buffer->entries[i].frequency_count);
        } else if (buffer->kmer_size <= 32) {
            appendStringInfo(&query, "(%lu, %d)", 
                           buffer->entries[i].kmer_data.k32_data,
                           buffer->entries[i].frequency_count);
        } else {
            /* k > 32 not supported */
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                     errmsg("k-mer length must be between 4 and 32"),
                     errdetail("Provided k-mer length: %d", buffer->kmer_size)));
        }
    }
    
    /* Handle conflict resolution */
    appendStringInfo(&query, " ON CONFLICT (kmer_data) DO UPDATE SET frequency_count = %s.frequency_count + EXCLUDED.frequency_count", temp_table_name);
    
    /* Execute the query */
    SPI_connect();
    SPI_exec(query.data, 0);
    SPI_finish();
    
    /* Reset buffer */
    buffer->count = 0;
    
    pfree(query.data);
}

/*
 * Flush hash buffer to temporary table (simplified for Phase 1)
 */
static void
kmersearch_flush_hash_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    HTAB *aggregation_hash = NULL;
    HASHCTL hash_ctl;
    HASH_SEQ_STATUS status;
    bool first;
    typedef struct {
        uint64_t kmer_hash;
        int total_frequency;
    } AggregationEntry;
    AggregationEntry *entry;
    
    if (buffer->count == 0) {
        return;
    }
    
    initStringInfo(&query);
    
    /* Aggregate duplicates in buffer before inserting */
    
    /* Create hash table for aggregation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64_t);
    hash_ctl.entrysize = sizeof(AggregationEntry);
    hash_ctl.hcxt = CurrentMemoryContext;
    
    aggregation_hash = hash_create("buffer_aggregation", buffer->count, &hash_ctl, 
                                   HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
    
    /* Aggregate frequency counts for duplicate k-mers */
    for (i = 0; i < buffer->count; i++) {
        uint64_t kmer_hash = buffer->entries[i].kmer_data.k32_data;
        bool found;
        AggregationEntry *hash_entry = (AggregationEntry *) hash_search(aggregation_hash, 
                                                                        (void *) &kmer_hash, 
                                                                        HASH_ENTER, &found);
        if (!found) {
            hash_entry->kmer_hash = kmer_hash;
            hash_entry->total_frequency = buffer->entries[i].frequency_count;
        } else {
            hash_entry->total_frequency += buffer->entries[i].frequency_count;
        }
    }
    
    /* Build bulk INSERT statement with aggregated values */
    appendStringInfo(&query, "INSERT INTO %s (kmer_data, frequency_count) VALUES ", temp_table_name);
    
    /* Iterate through aggregated entries */
    first = true;
    
    hash_seq_init(&status, aggregation_hash);
    while ((entry = (AggregationEntry *) hash_seq_search(&status)) != NULL) {
        if (!first) appendStringInfoString(&query, ", ");
        
        appendStringInfo(&query, "(%lu, %d)", entry->kmer_hash, entry->total_frequency);
        first = false;
    }
    
    /* Clean up aggregation hash */
    hash_destroy(aggregation_hash);
    
    appendStringInfo(&query, " ON CONFLICT (kmer_data) DO UPDATE SET frequency_count = %s.frequency_count + EXCLUDED.frequency_count", temp_table_name);
    
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
 * Add hash value to buffer for Phase 1 processing
 */
static void
kmersearch_add_hash_to_buffer(KmerBuffer *buffer, uint64_t kmer_hash, const char *temp_table_name)
{
    CompactKmerFreq *entry;
    int i;
    
    /* No buffer-level deduplication - each row contribution should be counted */
    /* Row-level deduplication is already handled in the calling function */
    
    /* Check if buffer is full */
    if (buffer->count >= buffer->capacity) {
        kmersearch_flush_hash_buffer_to_table(buffer, temp_table_name);
        
        /* After flush, continue to add the new entry */
    }
    
    /* Add new entry using k32_data field to store uint64_t hash */
    entry = &buffer->entries[buffer->count];
    entry->kmer_data.k32_data = kmer_hash;  /* Store hash in k32_data field */
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
    
    initStringInfo(&query);
    
    if (k_size > 32)
    {
        /* k > 32 not supported */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k_size)));
    }
    
    /* For k <= 32, use single bigint column */
    appendStringInfo(&query, 
        "CREATE TEMP TABLE %s ("
        "kmer_data bigint PRIMARY KEY, "
        "frequency_count integer DEFAULT 1"
        ")", temp_table_name);
    
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
void
kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, 
                                const char *column_name, int k_size, int target_attno, bool is_dna4_type)
{
    TableScanDesc scan;
    HeapTuple tuple;
    TupleDesc tupdesc;
    int i;
    Datum *kmer_datums;
    BlockNumber current_block;
    bool isnull;
    Datum value;
    VarBit *sequence;
    VarBit **kmers;
    int nkeys;
    int j;
    
    /* Use passed parameters instead of determining them again */
    tupdesc = RelationGetDescr(rel);
    
    /* Initialize buffer */
    kmersearch_init_buffer(&worker->buffer, k_size);
    
    /* Create temporary table for this worker with unique name */
    worker->temp_table_name = kmersearch_generate_unique_temp_table_name("temp_kmer_worker", worker->worker_id);
    kmersearch_create_worker_temp_table(worker->temp_table_name, k_size);
    
    /* Scan assigned blocks only */
    
    for (current_block = worker->start_block; current_block < worker->end_block; current_block++) {
        Buffer buffer;
        Page page;
        OffsetNumber max_offset;
        OffsetNumber offset;
        
        /* Read the block */
        buffer = ReadBufferExtended(rel, MAIN_FORKNUM, current_block, RBM_NORMAL, NULL);
        LockBuffer(buffer, BUFFER_LOCK_SHARE);
        page = BufferGetPage(buffer);
        max_offset = PageGetMaxOffsetNumber(page);
        
        /* Process all tuples in this block */
        for (offset = FirstOffsetNumber; offset <= max_offset; offset = OffsetNumberNext(offset)) {
            ItemId item_id;
            HeapTupleData tuple_data;
            HeapTuple block_tuple = &tuple_data;
            
            /* Get the tuple */
            item_id = PageGetItemId(page, offset);
            if (!ItemIdIsNormal(item_id))
                continue;
                
            block_tuple->t_len = ItemIdGetLength(item_id);
            block_tuple->t_data = (HeapTupleHeader) PageGetItem(page, item_id);
            block_tuple->t_tableOid = RelationGetRelid(rel);
            block_tuple->t_self.ip_blkid.bi_hi = current_block >> 16;
            block_tuple->t_self.ip_blkid.bi_lo = current_block & 0xFFFF;
            block_tuple->t_self.ip_posid = offset;
            
            worker->rows_processed++;
            
            /* Extract the DNA sequence value */
            value = heap_getattr(block_tuple, target_attno, tupdesc, &isnull);
        if (isnull) {
            continue;  /* Skip NULL values */
        }
        
        /* Convert DNA data to VarBit representation */
        sequence = DatumGetVarBitP(value);
        
        /* Extract k-mers from the sequence using SIMD-optimized function based on DNA type */
        if (is_dna4_type) {
            kmer_datums = kmersearch_extract_dna4_kmer2_with_expansion_direct(sequence, k_size, &nkeys);
        } else {
            kmer_datums = kmersearch_extract_dna2_kmer2_direct(sequence, k_size, &nkeys);
        }
        
        if (kmer_datums == NULL || nkeys == 0) {
            continue;
        }
        
        /* Use hash set to track unique k-mers in this row to avoid counting duplicates */
        {
            HTAB *row_kmer_set = NULL;
            HASHCTL hash_ctl;
            
            /* Create hash table for unique k-mers in this row */
            memset(&hash_ctl, 0, sizeof(hash_ctl));
            hash_ctl.keysize = sizeof(uint64_t);
            hash_ctl.entrysize = sizeof(uint64_t);
            hash_ctl.hcxt = CurrentMemoryContext;
            
            row_kmer_set = hash_create("row_kmer_set", nkeys, &hash_ctl, 
                                       HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
            
            /* Process each k-mer in this row - convert to hash and deduplicate */
            for (j = 0; j < nkeys; j++) {
                uint64_t kmer_hash;
                VarBit *kmer;
                bool found;
                
                /* Convert Datum to VarBit */
                kmer = DatumGetVarBitP(kmer_datums[j]);
                if (kmer == NULL) {
                    continue;
                }
                
                /* Convert k-mer to consistent hash value */
                kmer_hash = kmersearch_get_kmer_hash(kmer, 0, k_size);
                
                /* Only add to buffer if not already seen in this row */
                hash_search(row_kmer_set, (void *) &kmer_hash, HASH_ENTER, &found);
                if (!found) {
                    /* Add to buffer (will flush to temp table if full) */
                    kmersearch_add_hash_to_buffer(&worker->buffer, kmer_hash, worker->temp_table_name);
                }
            }
            
            /* Clean up row hash table */
            hash_destroy(row_kmer_set);
        }
        
        /* Cleanup k-mer array */
        if (kmer_datums) {
            pfree(kmer_datums);
        }
        }
        
        /* Release buffer and lock */
        UnlockReleaseBuffer(buffer);
    }
    
    /* Flush any remaining buffer contents */
    kmersearch_flush_hash_buffer_to_table(&worker->buffer, worker->temp_table_name);
    
    /* Cleanup buffer */
    if (worker->buffer.entries) {
        pfree(worker->buffer.entries);
    }
}
/*
 * Merge worker results using SQL aggregation
 */
void
kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, 
                                   const char *final_table_name, int k_size, int threshold_rows)
{
    StringInfoData query;
    StringInfoData union_query;
    const char *data_type;
    int i;
    int ret;
    
    initStringInfo(&query);
    initStringInfo(&union_query);
    
    /* Create final aggregation table */
    if (k_size > 32) {
        /* k > 32 not supported */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k_size)));
    }
    
    appendStringInfo(&query, 
        "CREATE TEMP TABLE %s ("
        "kmer_data bigint PRIMARY KEY, "
        "frequency_count integer"
        ")", final_table_name);
    
    /* Use existing SPI connection from main function - don't call SPI_connect() again */
    SPI_exec(query.data, 0);
    
    /* Debug: Check worker table contents before merging */
    for (i = 0; i < num_workers; i++) {
        StringInfoData debug_query;
        initStringInfo(&debug_query);
        appendStringInfo(&debug_query, "SELECT count(*) FROM %s", workers[i].temp_table_name);
        
        ret = SPI_exec(debug_query.data, 0);
        if (ret == SPI_OK_SELECT && SPI_processed > 0) {
            bool isnull;
            int row_count = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull));
        }
        pfree(debug_query.data);
    }

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
    
    /* Execute aggregation query with detailed error handling */
    ret = SPI_exec(query.data, 0);
    
    if (ret != SPI_OK_INSERT) {
        ereport(ERROR,
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to execute k-mer aggregation query"),
                 errdetail("SPI_exec returned %d", ret),
                 errhint("Query was: %s", query.data)));
    }
    
    /* Don't call SPI_finish() - leave connection open for main function */
    
    pfree(query.data);
    pfree(union_query.data);
}
/*
 * Persist highly frequent k-mers from temporary table to permanent tables
 */
static void
kmersearch_persist_highfreq_kmers_from_temp(Oid table_oid, const char *column_name, int k_size,
                                           const char *temp_table_name)
{
    StringInfoData query;
    int ret;
    
    initStringInfo(&query);
    
    /* Insert highly frequent k-mers into permanent table */
    /* Note: Converting integer values to bit strings for ngram_key */
    ereport(DEBUG1, (errmsg("Building INSERT query for %s with occur_bitlen=%d", temp_table_name, kmersearch_occur_bitlen)));
    
    if (k_size > 32) {
        /* k > 32 not supported */
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k_size)));
    }
    
    /* For k <= 32, use single column approach */
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, ngram_key, detection_reason) "
        "SELECT %u, '%s', "
        "  CASE "
        "    WHEN %d <= 8 THEN (kmer_data::integer)::bit(32) || (frequency_count::integer)::bit(%d) "
        "    WHEN %d <= 16 THEN (kmer_data::bigint)::bit(64) || (frequency_count::integer)::bit(%d) "
        "    ELSE (kmer_data::bigint)::bit(64) || (frequency_count::integer)::bit(%d) "
        "  END AS ngram_key, "
        "  'high_frequency' "
        "FROM %s "
        "WHERE kmer_data IS NOT NULL AND frequency_count > 0",
        table_oid, column_name, 
        k_size, kmersearch_occur_bitlen,
        k_size, kmersearch_occur_bitlen,
        kmersearch_occur_bitlen,
        temp_table_name);
    
    ereport(DEBUG1, (errmsg("Generated INSERT query: %s", query.data)));
    
    SPI_connect();
    
    /* Validate source table has data before insertion */
    {
        StringInfoData count_query;
        int count_ret;
        
        initStringInfo(&count_query);
        appendStringInfo(&count_query, "SELECT COUNT(*) FROM %s", temp_table_name);
        
        count_ret = SPI_exec(count_query.data, 0);
        if (count_ret == SPI_OK_SELECT && SPI_processed == 1) {
            bool isnull;
            int64 source_count = DatumGetInt64(SPI_getbinval(SPI_tuptable->vals[0], 
                                                             SPI_tuptable->tupdesc, 1, &isnull));
            ereport(NOTICE, (errmsg("Source table %s contains %ld high-frequency k-mer records", 
                                    temp_table_name, source_count)));
            
            if (source_count == 0) {
                ereport(WARNING, (errmsg("No high-frequency k-mers found in source table %s", temp_table_name)));
            }
        }
        pfree(count_query.data);
    }
    
    /* Execute high-frequency k-mer insertion with error checking */
    ereport(DEBUG1, (errmsg("Executing high-frequency k-mer insertion query")));
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_INSERT && ret != SPI_OK_INSERT_RETURNING) {
        ereport(ERROR, 
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to insert high-frequency k-mers into kmersearch_highfreq_kmer"),
                 errdetail("SPI_exec returned %d for INSERT query", ret),
                 errhint("Check table structure and permissions")));
    }
    
    /* Validate insertion results */
    if (SPI_processed == 0) {
        ereport(WARNING, (errmsg("No records were inserted into kmersearch_highfreq_kmer. Check data compatibility and SQL query syntax.")));
    } else {
        ereport(NOTICE, (errmsg("Successfully inserted %lu high-frequency k-mer records", SPI_processed)));
    }
    
    /* Insert metadata record */
    pfree(query.data);
    initStringInfo(&query);
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer_meta "
        "(table_oid, column_name, kmer_size, occur_bitlen, max_appearance_rate, max_appearance_nrow) "
        "VALUES (%u, '%s', %d, %d, %f, %d) "
        "ON CONFLICT (table_oid, column_name, kmer_size) DO UPDATE SET "
        "occur_bitlen = EXCLUDED.occur_bitlen, "
        "max_appearance_rate = EXCLUDED.max_appearance_rate, "
        "max_appearance_nrow = EXCLUDED.max_appearance_nrow, "
        "analysis_timestamp = now()",
        table_oid, column_name, k_size, kmersearch_occur_bitlen, 
        kmersearch_max_appearance_rate, kmersearch_max_appearance_nrow);
    
    /* Execute metadata insertion with error checking */
    ereport(DEBUG1, (errmsg("Executing metadata insertion query")));
    ret = SPI_exec(query.data, 0);
    if (ret != SPI_OK_INSERT && ret != SPI_OK_UPDATE && ret != SPI_OK_INSERT_RETURNING) {
        ereport(ERROR, 
                (errcode(ERRCODE_INTERNAL_ERROR),
                 errmsg("Failed to insert/update metadata in kmersearch_highfreq_kmer_meta"),
                 errdetail("SPI_exec returned %d for metadata INSERT/UPDATE query", ret),
                 errhint("Check table structure and permissions")));
    }
    ereport(DEBUG1, (errmsg("Successfully inserted/updated %lu metadata records", SPI_processed)));
    
    SPI_finish();
    
    pfree(query.data);
}

/*
 * Helper functions implementation
 */

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

/*
 * Check if tuple is in worker's processing range
 */
static bool
tuple_in_worker_range(HeapTuple tuple, KmerWorkerState *worker)
{
    BlockNumber block_num = ItemPointerGetBlockNumber(&tuple->t_self);
    return (block_num >= worker->start_block && block_num < worker->end_block);
}

/*
 * Extract DNA sequence from tuple
 */
static VarBit *
extract_sequence_from_tuple(HeapTuple tuple, int attno, TupleDesc tupdesc)
{
    bool isnull;
    Datum value;
    
    value = heap_getattr(tuple, attno, tupdesc, &isnull);
    if (isnull)
        return NULL;
        
    return DatumGetVarBitP(value);
}

/*
 * Create worker temporary table for n-gram keys
 */
static void
create_worker_ngram_temp_table(const char *table_name)
{
    StringInfoData query;
    int ret;
    
    initStringInfo(&query);
    
    /* First try to drop table if it exists to avoid conflicts */
    appendStringInfo(&query, "DROP TABLE IF EXISTS %s", table_name);
    
    ret = SPI_exec(query.data, 0);
    if (ret < 0)
        ereport(ERROR, (errmsg("Failed to drop existing temp table %s", table_name)));
    
    /* Create the temp table */
    resetStringInfo(&query);
    appendStringInfo(&query,
        "CREATE TEMP TABLE %s ("
        "ngram_key varbit NOT NULL"
        ")", table_name);
    
    ret = SPI_exec(query.data, 0);
    if (ret < 0)
        ereport(ERROR, (errmsg("Failed to create temp table %s", table_name)));
    
    pfree(query.data);
}

/*
 * Check if kmer2 value corresponds to a high-frequency k-mer
 */
static bool
is_kmer2_in_highfreq_table(uint64_t kmer2_value, const char *highfreq_table_name)
{
    StringInfoData query;
    bool result = false;
    int ret;
    
    initStringInfo(&query);
    
    /* Check if kmer2_value exists in temp_kmer_final table */
    if (strstr(highfreq_table_name, "temp_kmer_final_") != NULL) {
        /* Final aggregation table uses kmer_data column (bigint) */
        appendStringInfo(&query,
            "SELECT 1 FROM %s WHERE kmer_data = %lu LIMIT 1",
            highfreq_table_name, kmer2_value);
    } else {
        /* 
         * Permanent table uses different structure - not supported yet
         * NOTE: This path is never reached in current implementation since
         * highfreq_table_name is always generated with "temp_kmer_final_" prefix
         * via kmersearch_generate_unique_temp_table_name()
         */
        ereport(ERROR, (errmsg("Permanent table lookup not supported in is_kmer2_in_highfreq_table")));
    }
    
    SPI_connect();
    ret = SPI_exec(query.data, 0);
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
        result = true;
        
    SPI_finish();
    pfree(query.data);
    
    return result;
}

/*
 * Check if kmer2 value corresponds to a high-frequency k-mer using analysis dshash
 * This function uses kmersearch_get_kmer_hash for consistent hash calculation
 */
static bool
is_kmer2_in_analysis_dshash(uint64_t kmer2_value)
{
    /* The kmer2_value is already a hash generated by kmersearch_get_kmer_hash(),
     * so we can directly use it for dshash lookup */
    return kmersearch_is_kmer_hash_in_analysis_dshash(kmer2_value);
}

/*
 * Create ngram_key2 from kmer2 value and occurrence count
 * Format: kmer2_bits + occurrence_count_bits
 */
static VarBit*
create_ngram_key2_from_kmer2_and_count(uint64_t kmer2_value, int k_size, int occurrence_count)
{
    int kmer2_bits = k_size * 2;  /* 2 bits per base */
    int occur_bits = kmersearch_occur_bitlen;
    int total_bits = kmer2_bits + occur_bits;
    int total_bytes = (total_bits + 7) / 8;
    VarBit *result;
    unsigned char *data;
    uint64_t combined_value;
    int i, bit_idx;
    
    /* Allocate VarBit structure */
    result = (VarBit *) palloc(VARHDRSZ + total_bytes);
    SET_VARSIZE(result, VARHDRSZ + total_bytes);
    VARBITLEN(result) = total_bits;
    
    /* Initialize data to zero */
    data = VARBITS(result);
    memset(data, 0, total_bytes);
    
    /* Combine kmer2 and occurrence count: kmer2 in high bits, occurrence in low bits */
    combined_value = (kmer2_value << occur_bits) | (occurrence_count & ((1ULL << occur_bits) - 1));
    
    /* Set bits in big-endian order */
    for (i = 0; i < total_bits; i++) {
        bit_idx = i / 8;
        if ((combined_value & (1ULL << (total_bits - 1 - i))) != 0) {
            data[bit_idx] |= (0x80 >> (i % 8));
        }
    }
    
    return result;
}


/*
 * Process extracted kmer2 values and generate ngram_key2 for high-frequency ones
 * This function counts occurrences of each kmer2 within the current sequence
 */
static void
process_extracted_kmer2(uint64_t *kmer2_values, int nkeys, VarBit *sequence, int k_size, const char *worker_table, const char *highfreq_table)
{
    StringInfoData query;
    int i, j;
    
    if (!kmer2_values || nkeys <= 0)
        return;
    
    initStringInfo(&query);
    
    /* Process each unique kmer2 */
    for (i = 0; i < nkeys; i++) {
        uint64_t current_kmer2 = kmer2_values[i];
        
        /* Check if this kmer2 is in the high-frequency list */
        bool is_highfreq = false;
        
        /* Use dshash lookup if highfreq_table is NULL (indicating dshash mode),
         * otherwise fall back to table lookup for backward compatibility */
        if (highfreq_table == NULL) {
            is_highfreq = is_kmer2_in_analysis_dshash(current_kmer2);
        } else {
            is_highfreq = is_kmer2_in_highfreq_table(current_kmer2, highfreq_table);
        }
        
        if (is_highfreq) {
            /* Count occurrences of this kmer2 in the sequence */
            int occurrence_count = 0;
            VarBit *ngram_key2;
            
            for (j = 0; j < nkeys; j++) {
                if (kmer2_values[j] == current_kmer2) {
                    occurrence_count++;
                }
            }
            
            /* Create ngram_key2 with occurrence count */
            ngram_key2 = create_ngram_key2_from_kmer2_and_count(current_kmer2, k_size, occurrence_count);
            
            /* Insert into worker table */
            resetStringInfo(&query);
            appendStringInfo(&query,
                "INSERT INTO %s (ngram_key) VALUES ('%s')",
                worker_table,
                DatumGetCString(DirectFunctionCall1(varbit_out, VarBitPGetDatum(ngram_key2))));
            
            SPI_exec(query.data, 0);
            
            /* Mark processed kmer2 values to avoid duplicates */
            for (j = i; j < nkeys; j++) {
                if (kmer2_values[j] == current_kmer2) {
                    kmer2_values[j] = UINT64_MAX; /* Mark as processed */
                }
            }
            
            pfree(ngram_key2);
        }
    }
    
    pfree(query.data);
}

/*
 * Get processed row count for worker
 */
static int
get_processed_row_count(KmerWorkerState *worker)
{
    StringInfoData query;
    int count = 0;
    int ret;
    
    initStringInfo(&query);
    appendStringInfo(&query, "SELECT COUNT(*) FROM %s", worker->temp_table_name);
    
    ret = SPI_exec(query.data, 0);
    
    if (ret == SPI_OK_SELECT && SPI_processed > 0) {
        bool isnull;
        Datum result = SPI_getbinval(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1, &isnull);
        if (!isnull)
            count = DatumGetInt32(result);
    }
    
    pfree(query.data);
    
    return count;
}




/*
 * Module cleanup function
 */
void
_PG_fini(void)
{
    /* Free query pattern cache manager on module unload (uses TopMemoryContext - needs manual cleanup) */
    /* DNA2/DNA4 cache managers are now local and automatically freed with QueryContext */
    kmersearch_free_query_pattern_cache_internal();
    
    /* Free actual min score cache manager on module unload */
    kmersearch_free_actual_min_score_cache_internal();
    
    /* Free high-frequency k-mer cache on module unload */
    kmersearch_highfreq_kmer_cache_free_internal();
}

/*
 * Parallel high-frequency k-mer cache internal functions
 */






/*
 * Check if parallel high-frequency k-mer cache is valid
 */
static bool
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
        
        if (bit_offset <= 6) {
            output[byte_pos] |= (encoded << (6 - bit_offset));
        } else {
            /* bit_offset == 7: 1st bit to current byte, 2nd bit to next byte */
            output[byte_pos] |= (encoded >> 1);
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded & 0x1) << 7;
            }
        }
    }
}

static void dna2_decode_scalar(const uint8_t* input, char* output, int len)
{
    for (int i = 0; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = 0;
        
        if (bit_offset <= 6) {
            encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        } else {
            /* bit_offset == 7: 1st bit from current byte, 2nd bit from next byte */
            encoded = (input[byte_pos] & 0x1) << 1;
            if (byte_pos + 1 < (len * 2 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> 7) & 0x1;
            }
        }
        
        /* Range check for safety */
        if (encoded >= 4) {
            encoded = 0; /* Default to 'A' */
        }
        
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
            /* bit_offset > 4: split across two bytes */
            int remaining_bits = 8 - bit_offset;
            output[byte_pos] |= (encoded >> (4 - remaining_bits));
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
            }
        }
    }
}

static void dna4_decode_scalar(const uint8_t* input, char* output, int len)
{
    for (int i = 0; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = 0;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            /* bit_offset > 4: get bits from two bytes */
            int remaining_bits = 8 - bit_offset;
            encoded = (input[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
            if (byte_pos + 1 < (len * 4 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
            }
        }
        
        /* Range check for safety */
        if (encoded >= 16) {
            encoded = 0; /* Default to '?' */
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
    int simd_len;
    
    memset(output, 0, byte_len);
    
    /* Process 32 characters at a time with AVX2 */
    simd_len = len & ~31;  /* Round down to multiple of 32 */
    
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
        __m256i encoded;
        uint8_t temp[32];
        
        mask_T = _mm256_or_si256(mask_T, mask_U);
        
        /* Generate 2-bit encoded values */
        encoded = _mm256_setzero_si256();
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_C, _mm256_set1_epi8(1)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_G, _mm256_set1_epi8(2)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_T, _mm256_set1_epi8(3)));
        
        /* Pack 2-bit values into bytes - simplified approach */
        /* For full implementation, we'd need complex bit packing */
        /* Fall back to scalar for bit packing part */
        _mm256_storeu_si256((__m256i*)temp, encoded);
        
        for (int j = 0; j < 32; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 6) {
                output[byte_pos] |= (temp[j] << (6 - bit_offset));
            } else {
                /* bit_offset == 7: 1st bit to current byte, 2nd bit to next byte */
                output[byte_pos] |= (temp[j] >> 1);
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (temp[j] & 0x1) << 7;
                }
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 6) {
            output[byte_pos] |= (encoded << (6 - bit_offset));
        } else {
            /* bit_offset == 7: 1st bit to current byte, 2nd bit to next byte */
            output[byte_pos] |= (encoded >> 1);
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded & 0x1) << 7;
            }
        }
    }
}

__attribute__((target("avx2")))
static void dna2_decode_avx2(const uint8_t* input, char* output, int len)
{
    /* Process 32 characters at a time with AVX2 */
    int simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    for (int i = 0; i < simd_len; i += 32) {
        /* Calculate how many bytes needed for 32 characters (32 * 2 bits = 64 bits = 8 bytes) */
        int bytes_needed = 8;
        uint8_t temp[8];
        
        /* Load 8 bytes for 32 characters */
        for (int b = 0; b < bytes_needed && (i * 2 / 8 + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[i * 2 / 8 + b];
        }
        
        /* Decode each 2-bit pair */
        for (int j = 0; j < 32; j++) {
            int rel_bit_pos = j * 2;  /* Relative bit position within this 32-char block */
            int byte_pos = rel_bit_pos / 8;
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded = 0;
            
            if (bit_offset <= 6) {
                encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            } else {
                /* bit_offset == 7: 1st bit from current byte, 2nd bit from next byte */
                encoded = (temp[byte_pos] & 0x1) << 1;
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> 7) & 0x1;
                }
            }
            
            /* Range check for safety */
            if (encoded >= 4) {
                encoded = 0; /* Default to 'A' */
            }
            
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = 0;
        
        if (bit_offset <= 6) {
            encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        } else {
            /* bit_offset == 7: 1st bit from current byte, 2nd bit from next byte */
            encoded = (input[byte_pos] & 0x1) << 1;
            if (byte_pos + 1 < (len * 2 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> 7) & 0x1;
            }
        }
        
        /* Range check for safety */
        if (encoded >= 4) {
            encoded = 0; /* Default to 'A' */
        }
        
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

__attribute__((target("avx2")))
static void dna4_encode_avx2(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    int simd_len;
    
    memset(output, 0, byte_len);
    
    /* Process 32 characters at a time with AVX2 */
    simd_len = len & ~31;  /* Round down to multiple of 32 */
    
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
                /* bit_offset > 4: split across two bytes */
                int remaining_bits = 8 - bit_offset;
                output[byte_pos] |= (encoded >> (4 - remaining_bits));
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
                }
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
            /* bit_offset > 4: split across two bytes */
            int remaining_bits = 8 - bit_offset;
            output[byte_pos] |= (encoded >> (4 - remaining_bits));
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
            }
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
        /* Calculate how many bytes needed for 32 characters (32 * 4 bits = 128 bits = 16 bytes) */
        int bytes_needed = 16;
        uint8_t temp[16];
        
        /* Load 16 bytes for 32 characters */
        for (int b = 0; b < bytes_needed && (i * 4 / 8 + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[i * 4 / 8 + b];
        }
        
        /* Decode each 4-bit nibble */
        for (int j = 0; j < 32; j++) {
            int rel_bit_pos = j * 4;  /* Relative bit position within this 32-char block */
            int byte_pos = rel_bit_pos / 8;
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded = 0;
            
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                /* bit_offset > 4: get bits from two bytes */
                int remaining_bits = 8 - bit_offset;
                encoded = (temp[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
                }
            }
            
            /* Range check for safety */
            if (encoded >= 16) {
                encoded = 0; /* Default to '?' */
            }
            
            output[i + j] = kmersearch_dna4_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = 0;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            /* bit_offset > 4: get bits from two bytes */
            int remaining_bits = 8 - bit_offset;
            encoded = (input[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
            if (byte_pos + 1 < (len * 4 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
            }
        }
        
        /* Range check for safety */
        if (encoded >= 16) {
            encoded = 0; /* Default to '?' */
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
    int simd_len;
    
    memset(output, 0, byte_len);
    
    /* Process 64 characters at a time with AVX512 */
    simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    for (int i = 0; i < simd_len; i += 64) {
        __m512i chars = _mm512_loadu_si512((__m512i*)(input + i));
        
        /* Generate 2-bit encoded values using lookup table approach */
        uint8_t temp[64];
        _mm512_storeu_si512((__m512i*)temp, chars);
        
        /* Encode each character using lookup table */
        for (int j = 0; j < 64; j++) {
            temp[j] = kmersearch_dna2_encode_table[(unsigned char)temp[j]];
        }
        
        for (int j = 0; j < 64; j++) {
            uint8_t encoded = temp[j];
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 6) {
                output[byte_pos] |= (encoded << (6 - bit_offset));
            } else {
                /* bit_offset == 7: handle byte boundary crossing */
                output[byte_pos] |= (encoded >> 1);
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (encoded & 0x1) << 7;
                }
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 6) {
            output[byte_pos] |= (encoded << (6 - bit_offset));
        } else {
            /* bit_offset == 7: handle byte boundary crossing */
            output[byte_pos] |= (encoded >> 1);
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded & 0x1) << 7;
            }
        }
    }
}

__attribute__((target("avx512f,avx512bw")))
static void dna2_decode_avx512(const uint8_t* input, char* output, int len)
{
    /* Process 64 characters at a time with AVX512 */
    int simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    for (int i = 0; i < simd_len; i += 64) {
        /* 64 characters need 128 bits = 16 bytes */
        int bytes_needed = 16;
        uint8_t temp[16];
        
        /* Clear temp array first */
        memset(temp, 0, bytes_needed);
        
        /* Load required bytes from correct position */
        for (int b = 0; b < bytes_needed && (i * 2 / 8 + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[i * 2 / 8 + b];
        }
        
        /* Decode each 2-bit pair using relative positioning */
        for (int j = 0; j < 64; j++) {
            int rel_bit_pos = j * 2;  /* Relative bit position in temp array */
            int byte_pos = rel_bit_pos / 8;
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded = 0;
            
            if (bit_offset <= 6) {
                encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            } else {
                /* bit_offset == 7: handle byte boundary crossing */
                encoded = (temp[byte_pos] & 0x1) << 1;
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> 7) & 0x1;
                }
            }
            
            if (encoded >= 4) {
                encoded = 0; /* Default to 'A' */
            }
            
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = 0;
        
        if (bit_offset <= 6) {
            encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        } else {
            /* bit_offset == 7: handle byte boundary crossing */
            encoded = (input[byte_pos] & 0x1) << 1;
            if (byte_pos + 1 < (len * 2 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> 7) & 0x1;
            }
        }
        
        if (encoded >= 4) {
            encoded = 0; /* Default to 'A' */
        }
        
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

__attribute__((target("avx512f,avx512bw")))
static void dna4_encode_avx512(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    int simd_len;
    
    memset(output, 0, byte_len);
    
    /* Process 64 characters at a time with AVX512 */
    simd_len = len & ~63;  /* Round down to multiple of 64 */
    
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
                /* bit_offset > 4: handle byte boundary crossing */
                int remaining_bits = 8 - bit_offset;
                output[byte_pos] |= (encoded >> (4 - remaining_bits));
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
                }
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
            /* bit_offset > 4: handle byte boundary crossing */
            int remaining_bits = 8 - bit_offset;
            output[byte_pos] |= (encoded >> (4 - remaining_bits));
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
            }
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
        /* 64 characters need 256 bits = 32 bytes */
        int bytes_needed = 32;
        uint8_t temp[32];
        
        /* Load required bytes from correct position */
        for (int b = 0; b < bytes_needed && (i * 4 / 8 + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[i * 4 / 8 + b];
        }
        
        /* Decode each 4-bit nibble using relative positioning */
        for (int j = 0; j < 64; j++) {
            int rel_bit_pos = j * 4;  /* Relative bit position in temp array */
            int byte_pos = rel_bit_pos / 8;
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded = 0;
            
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                /* bit_offset > 4: handle byte boundary crossing */
                int remaining_bits = 8 - bit_offset;
                encoded = (temp[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
                }
            }
            
            if (encoded >= 16) {
                encoded = 0; /* Default to '?' */
            }
            
            output[i + j] = kmersearch_dna4_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded = 0;
        
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            /* bit_offset > 4: handle byte boundary crossing */
            int remaining_bits = 8 - bit_offset;
            encoded = (input[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
            if (byte_pos + 1 < (len * 4 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
            }
        }
        
        if (encoded >= 16) {
            encoded = 0; /* Default to '?' */
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}
#endif

#ifdef __aarch64__
/* NEON implementations */
__attribute__((target("+simd")))
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
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 6) {
                output[byte_pos] |= (temp[j] << (6 - bit_offset));
            } else {
                /* bit_offset == 7: split across two bytes */
                output[byte_pos] |= (temp[j] >> 1);
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (temp[j] & 0x1) << 7;
                }
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 6) {
            output[byte_pos] |= (encoded << (6 - bit_offset));
        } else {
            /* bit_offset == 7: split across two bytes */
            output[byte_pos] |= (encoded >> 1);
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded & 0x1) << 7;
            }
        }
    }
}

__attribute__((target("+simd")))
static void dna2_decode_neon(const uint8_t* input, char* output, int len)
{
    /* Process 16 characters at a time with NEON */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    for (int i = 0; i < simd_len; i += 16) {
        /* Fixed: Load correct number of bytes for 16 characters (32 bits = 4 bytes) */
        int bytes_needed = 4;  /* 16 * 2 bits = 32 bits = 4 bytes */
        uint8_t temp[4];
        
        /* Fixed: Load from correct position in input array */
        for (int b = 0; b < bytes_needed && (i * 2 / 8 + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[i * 2 / 8 + b];
        }
        
        /* Decode each 2-bit pair using relative positions in temp array */
        for (int j = 0; j < 16; j++) {
            int rel_bit_pos = j * 2;  /* Relative position within temp array */
            int byte_pos = rel_bit_pos / 8;  /* Byte position within temp array */
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded;
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 6) {
                encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            } else {
                /* bit_offset == 7: bits span across two bytes */
                encoded = (temp[byte_pos] & 0x1) << 1;
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> 7) & 0x1;
                }
            }
            
            /* Range check */
            if (encoded >= 4) {
                encoded = 0; /* Default to 'A' */
            }
            
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded;
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 6) {
            encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        } else {
            /* bit_offset == 7: bits span across two bytes */
            encoded = (input[byte_pos] & 0x1) << 1;
            if (byte_pos + 1 < (len * 2 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> 7) & 0x1;
            }
        }
        
        /* Range check */
        if (encoded >= 4) {
            encoded = 0; /* Default to 'A' */
        }
        
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

__attribute__((target("+simd")))
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
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 4) {
                output[byte_pos] |= (encoded << (4 - bit_offset));
            } else {
                /* bit_offset > 4: split across two bytes */
                int remaining_bits = 8 - bit_offset;
                output[byte_pos] |= (encoded >> (4 - remaining_bits));
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
                }
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 4) {
            output[byte_pos] |= (encoded << (4 - bit_offset));
        } else {
            /* bit_offset > 4: split across two bytes */
            int remaining_bits = 8 - bit_offset;
            output[byte_pos] |= (encoded >> (4 - remaining_bits));
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
            }
        }
    }
}

__attribute__((target("+simd")))
static void dna4_decode_neon(const uint8_t* input, char* output, int len)
{
    /* Process 16 characters at a time with NEON */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    for (int i = 0; i < simd_len; i += 16) {
        /* Fixed: Load correct number of bytes for 16 characters (64 bits = 8 bytes) */
        int bytes_needed = 8;  /* 16 * 4 bits = 64 bits = 8 bytes */
        uint8_t temp[8];
        
        /* Fixed: Load from correct position in input array */
        for (int b = 0; b < bytes_needed && (i * 4 / 8 + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[i * 4 / 8 + b];
        }
        
        /* Decode each 4-bit nibble using relative positions in temp array */
        for (int j = 0; j < 16; j++) {
            int rel_bit_pos = j * 4;  /* Relative position within temp array */
            int byte_pos = rel_bit_pos / 8;  /* Byte position within temp array */
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded;
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                /* bit_offset > 4: bits span across two bytes */
                int remaining_bits = 8 - bit_offset;
                encoded = (temp[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
                }
                encoded &= 0xF;
            }
            
            /* Range check */
            if (encoded >= 16) {
                encoded = 0; /* Default to '?' */
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
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            /* bit_offset > 4: bits span across two bytes */
            int remaining_bits = 8 - bit_offset;
            encoded = (input[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
            if (byte_pos + 1 < (len * 4 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
            }
            encoded &= 0xF;
        }
        
        /* Range check */
        if (encoded >= 16) {
            encoded = 0; /* Default to '?' */
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}

/* SVE implementations */
__attribute__((target("+sve")))
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
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 6) {
                output[byte_pos] |= (temp[j] << (6 - bit_offset));
            } else {
                /* bit_offset == 7: split across two bytes */
                output[byte_pos] |= (temp[j] >> 1);
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (temp[j] & 0x1) << 7;
                }
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna2_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 6) {
            output[byte_pos] |= (encoded << (6 - bit_offset));
        } else {
            /* bit_offset == 7: split across two bytes */
            output[byte_pos] |= (encoded >> 1);
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded & 0x1) << 7;
            }
        }
    }
}

__attribute__((target("+sve")))
static void dna2_decode_sve(const uint8_t* input, char* output, int len)
{
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    
    for (int i = 0; i < simd_len; i += sve_len) {
        /* Fixed: Load correct number of bytes for sve_len characters */
        int bytes_needed = (sve_len * 2 + 7) / 8;
        uint8_t temp[bytes_needed];
        
        /* Fixed: Load from correct position in input array */
        for (int b = 0; b < bytes_needed && (i * 2 / 8 + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[i * 2 / 8 + b];
        }
        
        /* Decode each 2-bit pair using relative positions in temp array */
        for (int j = 0; j < sve_len && (i + j) < len; j++) {
            int rel_bit_pos = j * 2;  /* Relative position within temp array */
            int byte_pos = rel_bit_pos / 8;  /* Byte position within temp array */
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded;
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 6) {
                encoded = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
            } else {
                /* bit_offset == 7: bits span across two bytes */
                encoded = (temp[byte_pos] & 0x1) << 1;
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> 7) & 0x1;
                }
            }
            
            /* Range check */
            if (encoded >= 4) {
                encoded = 0; /* Default to 'A' */
            }
            
            output[i + j] = kmersearch_dna2_decode_table[encoded];
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8_t encoded;
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 6) {
            encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        } else {
            /* bit_offset == 7: bits span across two bytes */
            encoded = (input[byte_pos] & 0x1) << 1;
            if (byte_pos + 1 < (len * 2 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> 7) & 0x1;
            }
        }
        
        /* Range check */
        if (encoded >= 4) {
            encoded = 0; /* Default to 'A' */
        }
        
        output[i] = kmersearch_dna2_decode_table[encoded];
    }
    output[len] = '\0';
}

__attribute__((target("+sve")))
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
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 4) {
                output[byte_pos] |= (encoded << (4 - bit_offset));
            } else {
                /* bit_offset > 4: split across two bytes */
                int remaining_bits = 8 - bit_offset;
                output[byte_pos] |= (encoded >> (4 - remaining_bits));
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
                }
            }
        }
    }
    
    /* Handle remaining characters with scalar */
    for (int i = simd_len; i < len; i++) {
        uint8_t encoded = kmersearch_dna4_encode_table[(unsigned char)input[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 4) {
            output[byte_pos] |= (encoded << (4 - bit_offset));
        } else {
            /* bit_offset > 4: split across two bytes */
            int remaining_bits = 8 - bit_offset;
            output[byte_pos] |= (encoded >> (4 - remaining_bits));
            if (byte_pos + 1 < byte_len) {
                output[byte_pos + 1] |= (encoded << (4 + remaining_bits));
            }
        }
    }
}

__attribute__((target("+sve")))
static void dna4_decode_sve(const uint8_t* input, char* output, int len)
{
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    
    for (int i = 0; i < simd_len; i += sve_len) {
        /* Fixed: Load correct number of bytes for sve_len characters */
        int bytes_needed = (sve_len * 4 + 7) / 8;
        uint8_t temp[bytes_needed];
        
        /* Fixed: Load from correct position in input array */
        for (int b = 0; b < bytes_needed && (i * 4 / 8 + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[i * 4 / 8 + b];
        }
        
        /* Decode each 4-bit nibble using relative positions in temp array */
        for (int j = 0; j < sve_len && (i + j) < len; j++) {
            int rel_bit_pos = j * 4;  /* Relative position within temp array */
            int byte_pos = rel_bit_pos / 8;  /* Byte position within temp array */
            int bit_offset = rel_bit_pos % 8;
            uint8_t encoded;
            
            /* Fixed: Handle byte boundary crossing properly */
            if (bit_offset <= 4) {
                encoded = (temp[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                /* bit_offset > 4: bits span across two bytes */
                int remaining_bits = 8 - bit_offset;
                encoded = (temp[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
                if (byte_pos + 1 < bytes_needed) {
                    encoded |= (temp[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
                }
                encoded &= 0xF;
            }
            
            /* Range check */
            if (encoded >= 16) {
                encoded = 0; /* Default to '?' */
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
        
        /* Fixed: Handle byte boundary crossing properly */
        if (bit_offset <= 4) {
            encoded = (input[byte_pos] >> (4 - bit_offset)) & 0xF;
        } else {
            /* bit_offset > 4: bits span across two bytes */
            int remaining_bits = 8 - bit_offset;
            encoded = (input[byte_pos] & ((1 << remaining_bits) - 1)) << (4 - remaining_bits);
            if (byte_pos + 1 < (len * 4 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> (4 + remaining_bits)) & 0xF;
            }
            encoded &= 0xF;
        }
        
        /* Range check */
        if (encoded >= 16) {
            encoded = 0; /* Default to '?' */
        }
        
        output[i] = kmersearch_dna4_decode_table[encoded];
    }
    output[len] = '\0';
}
#endif

/*
 * AVX2 K-mer Processing Functions
 */

#ifdef __x86_64__
/* AVX2 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("avx2")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_avx2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with AVX2-optimized bit extraction where possible */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
    /* AVX2-optimized batch processing for aligned k-mer extraction */
    for (i = 0; i < simd_batch; i += 8)
    {
        /* Extract 8 k-mers in parallel using AVX2 */
        __m256i positions = _mm256_setr_epi32(i, i+1, i+2, i+3, i+4, i+5, i+6, i+7);
        
        /* Process each k-mer in the batch */
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit *kmer_key;
            
            /* Check bounds before extraction to avoid invalid access */
            int last_bit_pos = (pos + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;  /* Out of bounds, skip */
            }
            
            /* Create k-mer key (without occurrence count) */
            kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, pos, k);
            if (kmer_key == NULL)
                continue;  /* Skip if key creation failed */
                
            keys[key_count++] = PointerGetDatum(kmer_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit *kmer_key;
        
        /* Check bounds before extraction to avoid invalid access */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;  /* Out of bounds, skip */
        }
        
        /* Create k-mer key (without occurrence count) */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/* AVX2 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("avx2")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with AVX2-optimized expansion where possible */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
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
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}

/* AVX2 optimized version of kmersearch_count_matching_kmer_fast */
__attribute__((target("avx2")))
static int
kmersearch_count_matching_kmer_fast_avx2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* AVX2 optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    /* Safety check: ensure we have valid query keys */
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);  /* Use data size, not total size */
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Creating hash with keysize=%zu, query_nkeys=%d", 
         (size_t)VARBITBYTES(query_keys[0]), query_nkeys);
    
    query_hash = hash_create("QueryKmerHashAVX2", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* Insert all query k-mers into hash table using content as key */
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: NULL query key at index %d", i);
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Check each sequence k-mer against hash table */
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: NULL seq key at index %d", i);
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Size mismatch seq[%d]=%zu vs query[0]=%zu", 
                 i, (size_t)VARBITBYTES(seq_keys[i]), (size_t)VARBITBYTES(query_keys[0]));
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Returning match_count=%d", match_count);
    return match_count;
}

/* AVX512 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("avx512f,avx512bw")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_avx512(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with AVX512-optimized bit extraction (16 k-mers at a time) */
    simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
    
    /* AVX512-optimized batch processing */
    for (i = 0; i < simd_batch; i += 16)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 16 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit *kmer_key;
            
            /* Check bounds before extraction to avoid invalid access */
            int last_bit_pos = (pos + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;  /* Out of bounds, skip */
            }
            
            /* Create k-mer key (without occurrence count) */
            kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, pos, k);
            if (kmer_key == NULL)
                continue;  /* Skip if key creation failed */
                
            keys[key_count++] = PointerGetDatum(kmer_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit *kmer_key;
        
        /* Check bounds before extraction to avoid invalid access */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;  /* Out of bounds, skip */
        }
        
        /* Create k-mer key (without occurrence count) */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/* AVX512 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("avx512f,avx512bw")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with AVX512-optimized expansion (16 k-mers at a time) */
    simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
    
    /* AVX512-optimized batch processing for k-mer expansion */
    for (i = 0; i < simd_batch; i += 16)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 16 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers */
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}

/* AVX512 optimized version of kmersearch_count_matching_kmer_fast */
__attribute__((target("avx512f,avx512bw")))
static int
kmersearch_count_matching_kmer_fast_avx512(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0) {
        return 0;
    }
    
    /* AVX512 optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        return 0;
    }
    
    /* Validate the VarBit structure */
    if (VARSIZE(query_keys[0]) < VARHDRSZ + 4) {
        elog(ERROR, "Invalid VarBit size: %zu", (size_t)VARSIZE(query_keys[0]));
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
/* NEON optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+simd")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_neon(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with NEON-optimized bit extraction (4 k-mers at a time) */
    simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
    
    /* NEON-optimized batch processing */
    for (i = 0; i < simd_batch; i += 4)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 4 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit *kmer_key;
            
            /* Check bounds before extraction to avoid invalid access */
            int last_bit_pos = (pos + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;  /* Out of bounds, skip */
            }
            
            /* Create k-mer key (without occurrence count) */
            kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, pos, k);
            if (kmer_key == NULL)
                continue;  /* Skip if key creation failed */
                
            keys[key_count++] = PointerGetDatum(kmer_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit *kmer_key;
        
        /* Check bounds before extraction to avoid invalid access */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;  /* Out of bounds, skip */
        }
        
        /* Create k-mer key (without occurrence count) */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/* NEON optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+simd")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_neon(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with NEON-optimized expansion (4 k-mers at a time) */
    simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
    
    /* NEON-optimized batch processing for k-mer expansion */
    for (i = 0; i < simd_batch; i += 4)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 4 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers */
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}

/* NEON optimized version of kmersearch_count_matching_kmer_fast */
__attribute__((target("+simd")))
static int
kmersearch_count_matching_kmer_fast_neon(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* NEON optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmer_fast_neon: NULL query key detected");
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
/* SVE optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+sve")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_sve(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with SVE-optimized bit extraction (variable vector width) */
    /* SVE vector length is runtime determined, so we use a conservative batch size */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time (conservative) */
    
    /* SVE-optimized batch processing */
    for (i = 0; i < simd_batch; i += 8)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit *kmer_key;
            
            /* Check bounds before extraction to avoid invalid access */
            int last_bit_pos = (pos + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;  /* Out of bounds, skip */
            }
            
            /* Create k-mer key (without occurrence count) */
            kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, pos, k);
            if (kmer_key == NULL)
                continue;  /* Skip if key creation failed */
                
            keys[key_count++] = PointerGetDatum(kmer_key);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit *kmer_key;
        
        /* Check bounds before extraction to avoid invalid access */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;  /* Out of bounds, skip */
        }
        
        /* Create k-mer key (without occurrence count) */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/* SVE optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+sve")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_sve(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with SVE-optimized expansion (8 k-mers at a time, conservative) */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time (conservative) */
    
    /* SVE-optimized batch processing for k-mer expansion */
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
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}

/* SVE optimized version of kmersearch_count_matching_kmer_fast */
__attribute__((target("+sve")))
static int
kmersearch_count_matching_kmer_fast_sve(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* SVE optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmer_fast_sve: NULL query key detected");
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

/*
 * Get appropriate uint size for k-mer length
 */
size_t
kmersearch_get_kmer_uint_size(int k)
{
    if (k <= 8) return sizeof(uint16);
    if (k <= 16) return sizeof(uint32);
    if (k <= 32) return sizeof(uint64);
    ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                   errmsg("k-mer length must be between 4 and 32")));
}


/*
 * Extract DNA4 base and get all possible expansions directly as uint values
 */
static void
kmersearch_expand_dna4_base_direct(uint8 dna4_base, uint64 *expansions, int *count)
{
    /* DNA4 base expansion mapping - direct to 2-bit values */
    static const uint8 dna4_expansions[][4] = {
        {0},           /* A=0001 -> A(00) */
        {1},           /* C=0010 -> C(01) */
        {0, 1},        /* M=0011 -> A(00), C(01) */
        {2},           /* G=0100 -> G(10) */
        {0, 2},        /* R=0101 -> A(00), G(10) */
        {1, 2},        /* S=0110 -> C(01), G(10) */
        {0, 1, 2},     /* V=0111 -> A(00), C(01), G(10) */
        {3},           /* T=1000 -> T(11) */
        {0, 3},        /* W=1001 -> A(00), T(11) */
        {1, 3},        /* Y=1010 -> C(01), T(11) */
        {0, 1, 3},     /* H=1011 -> A(00), C(01), T(11) */
        {2, 3},        /* K=1100 -> G(10), T(11) */
        {0, 2, 3},     /* D=1101 -> A(00), G(10), T(11) */
        {1, 2, 3},     /* B=1110 -> C(01), G(10), T(11) */
        {0, 1, 2, 3}   /* N=1111 -> A(00), C(01), G(10), T(11) */
    };
    
    static const int dna4_expansion_counts[] = {1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
    
    if (dna4_base == 0 || dna4_base > 15) {
        *count = 0;
        return;
    }
    
    *count = dna4_expansion_counts[dna4_base - 1];
    for (int i = 0; i < *count; i++) {
        expansions[i] = dna4_expansions[dna4_base - 1][i];
    }
}

/*
 * Expand DNA4 k-mer to DNA2 k-mers directly as uint values
 * Based on kmersearch_expand_dna4_kmer2_to_dna2_direct() but returns uint arrays
 */
void *
kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[32][4];  /* Max k=32, max 4 expansions per base */
    int base_counts[32];
    int total_combinations = 1;
    void *results;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    int i, combo;
    
    *expansion_count = 0;
    
    /* Check if expansion will exceed limit */
    if (kmersearch_will_exceed_degenerate_limit_dna4_bits(dna4_seq, start_pos, k))
        return NULL;
    
    /* Extract expansion info for each base using the same logic as the original function */
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
    
    /* Allocate result array for the appropriate uint type */
    results = palloc(total_combinations * element_size);
    
    /* Generate all combinations directly as uint values */
    for (combo = 0; combo < total_combinations; combo++)
    {
        uint64 kmer_value = 0;
        int temp_combo = combo;
        
        /* Generate this combination */
        for (i = 0; i < k; i++)
        {
            int base_idx = temp_combo % base_counts[i];
            uint8 dna2_base = base_expansions[i][base_idx];
            kmer_value = (kmer_value << 2) | dna2_base;
            temp_combo /= base_counts[i];
        }
        
        /* Store in appropriate uint type */
        if (k <= 8) {
            uint16 *uint_results = (uint16 *) results;
            uint_results[combo] = (uint16) kmer_value;
        } else if (k <= 16) {
            uint32 *uint_results = (uint32 *) results;
            uint_results[combo] = (uint32) kmer_value;
        } else {
            uint64 *uint_results = (uint64 *) results;
            uint_results[combo] = kmer_value;
        }
    }
    
    *expansion_count = total_combinations;
    return results;
}

/*
 * DNA2 AVX2 implementation - DIRECT BIT MANIPULATION
 */
#ifdef __x86_64__
__attribute__((target("avx2")))
static void
kmersearch_extract_dna2_kmer2_as_uint_direct_avx2(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    bits8 *src_data = VARBITS(seq);
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    result = palloc(max_kmers * element_size);
    *output = result;
    
    /* AVX2 optimized k-mer extraction */
    if (k <= 8) {
        uint16 *output_array = (uint16 *) result;
        int simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
        
        /* AVX2 batch processing */
        for (i = 0; i < simd_batch; i += 16) {
            /* Process 16 k-mers in parallel */
            for (int batch_idx = 0; batch_idx < 16; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                /* Extract k-mer directly from VarBit data */
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint16) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint16) kmer_value;
            (*nkeys)++;
        }
    } else if (k <= 16) {
        uint32 *output_array = (uint32 *) result;
        int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
        
        /* AVX2 batch processing */
        for (i = 0; i < simd_batch; i += 8) {
            for (int batch_idx = 0; batch_idx < 8; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint32) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint32) kmer_value;
            (*nkeys)++;
        }
    } else {
        uint64 *output_array = (uint64 *) result;
        int simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
        
        /* AVX2 batch processing */
        for (i = 0; i < simd_batch; i += 4) {
            for (int batch_idx = 0; batch_idx < 4; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = kmer_value;
            (*nkeys)++;
        }
    }
}

/*
 * DNA2 AVX512 implementation - DIRECT BIT MANIPULATION
 */
__attribute__((target("avx512bw")))
static void
kmersearch_extract_dna2_kmer2_as_uint_direct_avx512(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    bits8 *src_data = VARBITS(seq);
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    result = palloc(max_kmers * element_size);
    *output = result;
    
    /* AVX512 optimized k-mer extraction */
    if (k <= 8) {
        uint16 *output_array = (uint16 *) result;
        int simd_batch = max_kmers & ~31;  /* Process 32 k-mers at a time */
        
        /* AVX512 batch processing */
        for (i = 0; i < simd_batch; i += 32) {
            for (int batch_idx = 0; batch_idx < 32; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint16) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint16) kmer_value;
            (*nkeys)++;
        }
    } else if (k <= 16) {
        uint32 *output_array = (uint32 *) result;
        int simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
        
        /* AVX512 batch processing */
        for (i = 0; i < simd_batch; i += 16) {
            for (int batch_idx = 0; batch_idx < 16; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint32) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint32) kmer_value;
            (*nkeys)++;
        }
    } else {
        uint64 *output_array = (uint64 *) result;
        int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
        
        /* AVX512 batch processing */
        for (i = 0; i < simd_batch; i += 8) {
            for (int batch_idx = 0; batch_idx < 8; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = kmer_value;
            (*nkeys)++;
        }
    }
}
#endif

/*
 * DNA2 NEON implementation - DIRECT BIT MANIPULATION
 */
#ifdef __aarch64__
__attribute__((target("+simd")))
static void
kmersearch_extract_dna2_kmer2_as_uint_direct_neon(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    bits8 *src_data = VARBITS(seq);
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    result = palloc(max_kmers * element_size);
    *output = result;
    
    /* NEON optimized k-mer extraction */
    if (k <= 8) {
        uint16 *output_array = (uint16 *) result;
        int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
        
        /* NEON batch processing */
        for (i = 0; i < simd_batch; i += 8) {
            for (int batch_idx = 0; batch_idx < 8; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint16) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint16) kmer_value;
            (*nkeys)++;
        }
    } else if (k <= 16) {
        uint32 *output_array = (uint32 *) result;
        int simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
        
        /* NEON batch processing */
        for (i = 0; i < simd_batch; i += 4) {
            for (int batch_idx = 0; batch_idx < 4; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint32) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint32) kmer_value;
            (*nkeys)++;
        }
    } else {
        uint64 *output_array = (uint64 *) result;
        int simd_batch = max_kmers & ~1;  /* Process 2 k-mers at a time */
        
        /* NEON batch processing */
        for (i = 0; i < simd_batch; i += 2) {
            for (int batch_idx = 0; batch_idx < 2; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = kmer_value;
            (*nkeys)++;
        }
    }
}

/*
 * DNA2 SVE implementation - DIRECT BIT MANIPULATION
 */
__attribute__((target("+sve")))
static void
kmersearch_extract_dna2_kmer2_as_uint_direct_sve(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    bits8 *src_data = VARBITS(seq);
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    result = palloc(max_kmers * element_size);
    *output = result;
    
    /* SVE optimized k-mer extraction */
    if (k <= 8) {
        uint16 *output_array = (uint16 *) result;
        int simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
        
        /* SVE batch processing */
        for (i = 0; i < simd_batch; i += 16) {
            for (int batch_idx = 0; batch_idx < 16; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint16) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint16) kmer_value;
            (*nkeys)++;
        }
    } else if (k <= 16) {
        uint32 *output_array = (uint32 *) result;
        int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
        
        /* SVE batch processing */
        for (i = 0; i < simd_batch; i += 8) {
            for (int batch_idx = 0; batch_idx < 8; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = (uint32) kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint32) kmer_value;
            (*nkeys)++;
        }
    } else {
        uint64 *output_array = (uint64 *) result;
        int simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
        
        /* SVE batch processing */
        for (i = 0; i < simd_batch; i += 4) {
            for (int batch_idx = 0; batch_idx < 4; batch_idx++) {
                uint64 kmer_value = 0;
                int pos = i + batch_idx;
                
                for (int j = 0; j < k; j++) {
                    int bit_pos = (pos + j) * 2;
                    int byte_pos = bit_pos / 8;
                    int bit_offset = bit_pos % 8;
                    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                    kmer_value = (kmer_value << 2) | base_bits;
                }
                
                output_array[*nkeys] = kmer_value;
                (*nkeys)++;
            }
        }
        
        /* Handle remaining k-mers */
        for (i = simd_batch; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            
            for (int j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = kmer_value;
            (*nkeys)++;
        }
    }
}
#endif

/*
 * DNA2 scalar implementation - DIRECT BIT MANIPULATION
 */
static void
kmersearch_extract_dna2_kmer2_as_uint_direct_scalar(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    bits8 *src_data = VARBITS(seq);
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    result = palloc(max_kmers * element_size);
    *output = result;
    
    /* CRITICAL: Direct bit manipulation - NO VarBit intermediate objects */
    if (k <= 8) {
        uint16 *output_array = (uint16 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            int j;
            
            /* Extract k-mer directly from VarBit data */
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint16) kmer_value;
            (*nkeys)++;
        }
    } else if (k <= 16) {
        uint32 *output_array = (uint32 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            int j;
            
            /* Extract k-mer directly from VarBit data */
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint32) kmer_value;
            (*nkeys)++;
        }
    } else {
        uint64 *output_array = (uint64 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            int j;
            
            /* Extract k-mer directly from VarBit data */
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = kmer_value;
            (*nkeys)++;
        }
    }
}

/*
 * DNA4 AVX2 implementation - DIRECT DEGENERATE EXPANSION
 */
#ifdef __x86_64__
__attribute__((target("avx2")))
static void
kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx2(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    /* Allocate result array with room for expansions (max 10 expansions per k-mer position) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    /* AVX2 optimized k-mer extraction with expansion */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
    /* AVX2 batch processing */
    for (i = 0; i < simd_batch; i += 8) {
        for (int batch_idx = 0; batch_idx < 8; batch_idx++) {
            int pos = i + batch_idx;
            void *expanded_uints;
            int expansion_count;
            int j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers directly as uint values */
            expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_uints || expansion_count == 0)
                continue;
            
            /* Copy expanded uint values to result array */
            if (k <= 8) {
                uint16 *src_array = (uint16 *) expanded_uints;
                uint16 *dst_array = (uint16 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else if (k <= 16) {
                uint32 *src_array = (uint32 *) expanded_uints;
                uint32 *dst_array = (uint32 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else {
                uint64 *src_array = (uint64 *) expanded_uints;
                uint64 *dst_array = (uint64 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            }
            
            pfree(expanded_uints);
        }
    }
    
    /* Handle remaining k-mers */
    for (i = simd_batch; i <= seq_bases - k; i++) {
        void *expanded_uints;
        int expansion_count;
        int j;
        
        expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_uints || expansion_count == 0)
            continue;
        
        if (k <= 8) {
            uint16 *src_array = (uint16 *) expanded_uints;
            uint16 *dst_array = (uint16 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else if (k <= 16) {
            uint32 *src_array = (uint32 *) expanded_uints;
            uint32 *dst_array = (uint32 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else {
            uint64 *src_array = (uint64 *) expanded_uints;
            uint64 *dst_array = (uint64 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        }
        
        pfree(expanded_uints);
    }
    
    *nkeys = key_count;
}

/*
 * DNA4 AVX512 implementation - DIRECT DEGENERATE EXPANSION
 */
__attribute__((target("avx512bw")))
static void
kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx512(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    /* Allocate result array with room for expansions (max 10 expansions per k-mer position) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    /* AVX512 optimized k-mer extraction with expansion */
    simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
    
    /* AVX512 batch processing */
    for (i = 0; i < simd_batch; i += 16) {
        for (int batch_idx = 0; batch_idx < 16; batch_idx++) {
            int pos = i + batch_idx;
            void *expanded_uints;
            int expansion_count;
            int j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers directly as uint values */
            expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_uints || expansion_count == 0)
                continue;
            
            /* Copy expanded uint values to result array */
            if (k <= 8) {
                uint16 *src_array = (uint16 *) expanded_uints;
                uint16 *dst_array = (uint16 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else if (k <= 16) {
                uint32 *src_array = (uint32 *) expanded_uints;
                uint32 *dst_array = (uint32 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else {
                uint64 *src_array = (uint64 *) expanded_uints;
                uint64 *dst_array = (uint64 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            }
            
            pfree(expanded_uints);
        }
    }
    
    /* Handle remaining k-mers */
    for (i = simd_batch; i <= seq_bases - k; i++) {
        void *expanded_uints;
        int expansion_count;
        int j;
        
        expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_uints || expansion_count == 0)
            continue;
        
        if (k <= 8) {
            uint16 *src_array = (uint16 *) expanded_uints;
            uint16 *dst_array = (uint16 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else if (k <= 16) {
            uint32 *src_array = (uint32 *) expanded_uints;
            uint32 *dst_array = (uint32 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else {
            uint64 *src_array = (uint64 *) expanded_uints;
            uint64 *dst_array = (uint64 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        }
        
        pfree(expanded_uints);
    }
    
    *nkeys = key_count;
}
#endif

/*
 * DNA4 NEON implementation - DIRECT DEGENERATE EXPANSION
 */
#ifdef __aarch64__
__attribute__((target("+simd")))
static void
kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_neon(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    /* Allocate result array with room for expansions (max 10 expansions per k-mer position) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    /* NEON optimized k-mer extraction with expansion */
    int simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
    
    /* NEON batch processing */
    for (i = 0; i < simd_batch; i += 4) {
        for (int batch_idx = 0; batch_idx < 4; batch_idx++) {
            int pos = i + batch_idx;
            void *expanded_uints;
            int expansion_count;
            int j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers directly as uint values */
            expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_uints || expansion_count == 0)
                continue;
            
            /* Copy expanded uint values to result array */
            if (k <= 8) {
                uint16 *src_array = (uint16 *) expanded_uints;
                uint16 *dst_array = (uint16 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else if (k <= 16) {
                uint32 *src_array = (uint32 *) expanded_uints;
                uint32 *dst_array = (uint32 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else {
                uint64 *src_array = (uint64 *) expanded_uints;
                uint64 *dst_array = (uint64 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            }
            
            pfree(expanded_uints);
        }
    }
    
    /* Handle remaining k-mers */
    for (i = simd_batch; i <= seq_bases - k; i++) {
        void *expanded_uints;
        int expansion_count;
        int j;
        
        expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_uints || expansion_count == 0)
            continue;
        
        if (k <= 8) {
            uint16 *src_array = (uint16 *) expanded_uints;
            uint16 *dst_array = (uint16 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else if (k <= 16) {
            uint32 *src_array = (uint32 *) expanded_uints;
            uint32 *dst_array = (uint32 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else {
            uint64 *src_array = (uint64 *) expanded_uints;
            uint64 *dst_array = (uint64 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        }
        
        pfree(expanded_uints);
    }
    
    *nkeys = key_count;
}

/*
 * DNA4 SVE implementation - DIRECT DEGENERATE EXPANSION
 */
__attribute__((target("+sve")))
static void
kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_sve(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    /* Allocate result array with room for expansions (max 10 expansions per k-mer position) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    /* SVE optimized k-mer extraction with expansion */
    int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
    /* SVE batch processing */
    for (i = 0; i < simd_batch; i += 8) {
        for (int batch_idx = 0; batch_idx < 8; batch_idx++) {
            int pos = i + batch_idx;
            void *expanded_uints;
            int expansion_count;
            int j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers directly as uint values */
            expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_uints || expansion_count == 0)
                continue;
            
            /* Copy expanded uint values to result array */
            if (k <= 8) {
                uint16 *src_array = (uint16 *) expanded_uints;
                uint16 *dst_array = (uint16 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else if (k <= 16) {
                uint32 *src_array = (uint32 *) expanded_uints;
                uint32 *dst_array = (uint32 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            } else {
                uint64 *src_array = (uint64 *) expanded_uints;
                uint64 *dst_array = (uint64 *) result;
                for (j = 0; j < expansion_count; j++) {
                    dst_array[key_count++] = src_array[j];
                }
            }
            
            pfree(expanded_uints);
        }
    }
    
    /* Handle remaining k-mers */
    for (i = simd_batch; i <= seq_bases - k; i++) {
        void *expanded_uints;
        int expansion_count;
        int j;
        
        expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_uints || expansion_count == 0)
            continue;
        
        if (k <= 8) {
            uint16 *src_array = (uint16 *) expanded_uints;
            uint16 *dst_array = (uint16 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else if (k <= 16) {
            uint32 *src_array = (uint32 *) expanded_uints;
            uint32 *dst_array = (uint32 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else {
            uint64 *src_array = (uint64 *) expanded_uints;
            uint64 *dst_array = (uint64 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        }
        
        pfree(expanded_uints);
    }
    
    *nkeys = key_count;
}
#endif

/*
 * DNA4 scalar implementation - DIRECT DEGENERATE EXPANSION
 * Now using kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct() for direct uint output
 */
static void
kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_scalar(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    /* Allocate result array with room for expansions (max 10 expansions per k-mer position) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    /* Extract k-mers using direct uint expansion */
    for (i = 0; i <= seq_bases - k; i++)
    {
        void *expanded_uints;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers directly as uint values */
        expanded_uints = kmersearch_expand_dna4_kmer2_as_uint_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_uints || expansion_count == 0)
            continue;
        
        /* Copy expanded uint values to result array */
        if (k <= 8) {
            uint16 *src_array = (uint16 *) expanded_uints;
            uint16 *dst_array = (uint16 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else if (k <= 16) {
            uint32 *src_array = (uint32 *) expanded_uints;
            uint32 *dst_array = (uint32 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        } else {
            uint64 *src_array = (uint64 *) expanded_uints;
            uint64 *dst_array = (uint64 *) result;
            for (j = 0; j < expansion_count; j++) {
                dst_array[key_count++] = src_array[j];
            }
        }
        
        /* Free the expanded uint array */
        pfree(expanded_uints);
    }
    
    *nkeys = key_count;
}

/*
 * DNA2 function dispatch
 */
void
kmersearch_extract_dna2_kmer2_as_uint_direct(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        /* AVX512 implementation */
        kmersearch_extract_dna2_kmer2_as_uint_direct_avx512(seq, k, output, nkeys);
        return;
    }
    if (simd_capability >= SIMD_AVX2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        /* AVX2 implementation */
        kmersearch_extract_dna2_kmer2_as_uint_direct_avx2(seq, k, output, nkeys);
        return;
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        /* SVE implementation */
        kmersearch_extract_dna2_kmer2_as_uint_direct_sve(seq, k, output, nkeys);
        return;
    }
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        /* NEON implementation */
        kmersearch_extract_dna2_kmer2_as_uint_direct_neon(seq, k, output, nkeys);
        return;
    }
#endif
    kmersearch_extract_dna2_kmer2_as_uint_direct_scalar(seq, k, output, nkeys);
}

/*
 * DNA4 function dispatch
 */
void
kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        /* AVX512 implementation */
        kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx512(seq, k, output, nkeys);
        return;
    }
    if (simd_capability >= SIMD_AVX2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        /* AVX2 implementation */
        kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx2(seq, k, output, nkeys);
        return;
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        /* SVE implementation */
        kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_sve(seq, k, output, nkeys);
        return;
    }
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        /* NEON implementation */
        kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_neon(seq, k, output, nkeys);
        return;
    }
#endif
    kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_scalar(seq, k, output, nkeys);
}
