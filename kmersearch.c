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

/* Global SIMD capability */
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
int kmersearch_highfreq_analysis_batch_size = 10000;  /* Default batch size for high-frequency k-mer analysis */

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
/* Functions moved to kmersearch_util.c */
/* Function moved to kmersearch_freq.c */
/* Moved to kmersearch_freq.c */
/* Moved to kmersearch_kmer.c */

/* Function moved to kmersearch_cache.c */

/* Analysis dshash functions (defined in kmersearch_freq.c) */
/* B-2: Other functions */
Datum *kmersearch_extract_dna2_kmer2_direct(VarBit *seq, int k, int *nkeys);
Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct(VarBit *seq, int k, int *nkeys);
Datum *kmersearch_extract_dna2_ngram_key2_direct(VarBit *seq, int k, int *nkeys);
Datum *kmersearch_extract_dna4_ngram_key2_with_expansion_direct(VarBit *seq, int k, int *nkeys);
/* Moved to kmersearch_kmer.c */
/* Function moved to kmersearch_gin.c */

/* Actual min score cache functions */
/* Function moved to kmersearch_gin.c */

/* New parallel analysis functions */
/* Moved to kmersearch_freq.c */
/* Moved to kmersearch_freq.c */
/* Function moved to kmersearch_freq.c */

/* New memory-efficient k-mer functions */
/* Function moved to kmersearch_gin.c */
/* Buffer functions moved to kmersearch_freq.c */
/* Function moved to kmersearch_freq.c */
/* Function moved to kmersearch_freq.c */
/* Function moved to kmersearch_freq.c */

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

/* SIMD implementation functions */

/* Scalar versions */
static Datum *kmersearch_extract_dna2_kmer2_direct_scalar(VarBit *seq, int k, int *nkeys);

/* SIMD function declarations are now in kmersearch.h */

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
        free_query_pattern_cache_manager(&query_pattern_cache_manager);
    
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
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
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
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
        free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
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

    DefineCustomIntVariable("kmersearch.highfreq_analysis_batch_size",
                           "Batch size for high-frequency k-mer analysis",
                           "Controls the number of rows processed in each batch during analysis",
                           &kmersearch_highfreq_analysis_batch_size,
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
#ifdef __x86_64__
static simd_capability_t detect_cpu_capabilities(void)
{
    unsigned int eax, ebx, ecx, edx;
    bool has_avx2 = false;
    bool has_bmi2 = false;
    bool has_avx512f = false;
    bool has_avx512bw = false;
    bool has_avx512vbmi = false;
    bool has_avx512vbmi2 = false;
    
    /* Check CPUID support level */
    unsigned int max_leaf = __get_cpuid_max(0, NULL);
    
    /* Check for AVX support (required for AVX2/AVX512) */
    if (max_leaf >= 1) {
        bool has_xsave;
        bool has_osxsave;
        
        __cpuid(1, eax, ebx, ecx, edx);
        has_xsave = (ecx & (1 << 27)) != 0;
        has_osxsave = (ecx & (1 << 28)) != 0;
        
        if (has_xsave && has_osxsave) {
            /* Check OS support for YMM/ZMM registers */
            unsigned long long xcr0 = 0;
            bool ymm_enabled;
            bool zmm_enabled;
            
            __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx");
            ymm_enabled = (xcr0 & 0x6) == 0x6;
            zmm_enabled = (xcr0 & 0xe6) == 0xe6;
            
            /* Check extended features */
            if (max_leaf >= 7) {
                __cpuid_count(7, 0, eax, ebx, ecx, edx);
                
                /* AVX2 and BMI2 */
                if (ymm_enabled) {
                    has_avx2 = (ebx & (1 << 5)) != 0;
                    has_bmi2 = (ebx & (1 << 8)) != 0;
                }
                
                /* AVX512 features */
                if (zmm_enabled) {
                    has_avx512f = (ebx & (1 << 16)) != 0;
                    has_avx512bw = (ebx & (1 << 30)) != 0;
                    has_avx512vbmi = (ecx & (1 << 1)) != 0;
                    has_avx512vbmi2 = (ecx & (1 << 6)) != 0;
                }
            }
        }
    }
    
    /* Return highest supported capability */
    if (has_avx512f && has_avx512bw && has_avx512vbmi && has_avx512vbmi2)
        return SIMD_AVX512VBMI2;
    if (has_avx512f && has_avx512bw && has_avx512vbmi)
        return SIMD_AVX512VBMI;
    if (has_avx512f && has_avx512bw)
        return SIMD_AVX512BW;
    if (has_avx512f)
        return SIMD_AVX512F;
    if (has_avx2 && has_bmi2)
        return SIMD_BMI2;
    if (has_avx2)
        return SIMD_AVX2;
    
    return SIMD_NONE;
}
#elif defined(__aarch64__)
static sigjmp_buf jmpbuf;

static void sigill_handler(int sig) {
    siglongjmp(jmpbuf, 1);
}

__attribute__((target("+sve2,+sve,+simd")))
static simd_capability_t detect_cpu_capabilities(void)
{
    struct sigaction sa, old_sa;
    sa.sa_handler = sigill_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGILL, &sa, &old_sa);
    
    /* Check for SVE2 support */
#ifdef __ARM_FEATURE_SVE2
    if (sigsetjmp(jmpbuf, 1) == 0) {
        /* Try SVE2-specific instruction */
        size_t vl = svcntd();
        (void)vl;
        sigaction(SIGILL, &old_sa, NULL);
        return SIMD_SVE2;
    }
#endif
    
    /* Check for SVE support */
#ifdef __ARM_FEATURE_SVE
    if (sigsetjmp(jmpbuf, 1) == 0) {
        size_t vl = svcntb();
        (void)vl;
        sigaction(SIGILL, &old_sa, NULL);
        return SIMD_SVE;
    }
#endif
    
    /* Check for NEON support */
#ifdef __ARM_NEON
    if (sigsetjmp(jmpbuf, 1) == 0) {
        volatile uint8x8_t a = vdup_n_u8(1);
        (void)a;
        sigaction(SIGILL, &old_sa, NULL);
        return SIMD_NEON;
    }
#endif
    
    return SIMD_NONE;
}
#endif

/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using bit operations (scalar implementation)
 */
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
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

#ifdef __x86_64__
/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using AVX2+BMI2
 */
__attribute__((target("avx2,bmi2")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_avx2(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
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
    
    /* Phase 1: Extract expansion info for each base using BMI2 */
    for (i = 0; i < k; i += 8)
    {
        int bases_to_process = (k - i) < 8 ? (k - i) : 8;
        int j;
        
        /* Extract 8 DNA4 bases (32 bits) at once if possible */
        if (bases_to_process == 8)
        {
            int bit_pos = (start_pos + i) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint64_t data64 = 0;
            uint64_t mask = 0x0F0F0F0F0F0F0F0FULL;
            uint64_t extracted;
            
            /* Load 64 bits to handle unaligned access */
            memcpy(&data64, &data[byte_pos], 8);
            data64 = data64 >> bit_offset;
            
            /* Use BMI2 PEXT to extract 4-bit values efficiently */
            extracted = _pext_u64(data64, mask);
            
            /* Process each extracted base */
            for (j = 0; j < 8; j++)
            {
                uint8 encoded = (extracted >> (j * 4)) & 0xF;
                int exp_count = kmersearch_dna4_to_dna2_table[encoded][0];
                int m;
                
                base_counts[i + j] = exp_count;
                for (m = 0; m < exp_count; m++)
                {
                    base_expansions[i + j][m] = kmersearch_dna4_to_dna2_table[encoded][m + 1];
                }
                total_combinations *= exp_count;
            }
        }
        else
        {
            /* Process remaining bases using scalar method */
            for (j = 0; j < bases_to_process; j++)
            {
                int bit_pos = (start_pos + i + j) * 4;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 encoded;
                int exp_count;
                int m;
                
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
                base_counts[i + j] = exp_count;
                
                for (m = 0; m < exp_count; m++)
                {
                    base_expansions[i + j][m] = kmersearch_dna4_to_dna2_table[encoded][m + 1];
                }
                
                total_combinations *= exp_count;
            }
        }
    }
    
    /* Allocate result array */
    results = (VarBit **) palloc(total_combinations * sizeof(VarBit *));
    
    /* Phase 2: Generate all combinations (using AVX2 for bit packing when beneficial) */
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
        /* TODO: Optimize this part with AVX2 when processing multiple bases */
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
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using AVX512+VBMI2
 */
__attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_avx512(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    /* TODO: Implement AVX512+VBMI2 optimized version */
    /* For now, fall back to scalar implementation */
    return kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(dna4_seq, start_pos, k, expansion_count);
}
#endif /* __x86_64__ */

#ifdef __aarch64__
/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using NEON
 */
__attribute__((target("+simd")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_neon(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    /* TODO: Implement NEON optimized version */
    /* For now, fall back to scalar implementation */
    return kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(dna4_seq, start_pos, k, expansion_count);
}

/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using SVE+NEON
 */
__attribute__((target("+sve,+simd")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_sve(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    /* TODO: Implement SVE+NEON optimized version */
    /* For now, fall back to scalar implementation */
    return kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(dna4_seq, start_pos, k, expansion_count);
}

/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using SVE2
 */
__attribute__((target("+sve2")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_sve2(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    /* TODO: Implement SVE2 optimized version */
    /* For now, fall back to scalar implementation */
    return kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(dna4_seq, start_pos, k, expansion_count);
}
#endif /* __aarch64__ */

/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using bit operations (dispatch function)
 */
VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    int seq_bits = VARBITLEN(dna4_seq);
    
    /* Quick check for small k values where SIMD overhead isn't worth it */
    if (k < 8) {
        return kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(dna4_seq, start_pos, k, expansion_count);
    }

#ifdef __x86_64__
    /* AVX512 with VBMI2 - best x86 performance */
    if (simd_capability >= SIMD_AVX512VBMI2 && seq_bits >= SIMD_DNA4KMER2_AVX512_THRESHOLD) {
        return kmersearch_expand_dna4_kmer2_to_dna2_direct_avx512(dna4_seq, start_pos, k, expansion_count);
    }
    /* AVX2 with BMI2 - good x86 performance */
    if (simd_capability >= SIMD_BMI2 && seq_bits >= SIMD_DNA4KMER2_AVX2_THRESHOLD) {
        return kmersearch_expand_dna4_kmer2_to_dna2_direct_avx2(dna4_seq, start_pos, k, expansion_count);
    }
#elif defined(__aarch64__)
    /* SVE2 - best ARM performance */
    if (simd_capability >= SIMD_SVE2 && seq_bits >= SIMD_DNA4KMER2_SVE_THRESHOLD) {
        return kmersearch_expand_dna4_kmer2_to_dna2_direct_sve2(dna4_seq, start_pos, k, expansion_count);
    }
    /* SVE with NEON assist */
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_DNA4KMER2_SVE_THRESHOLD) {
        return kmersearch_expand_dna4_kmer2_to_dna2_direct_sve(dna4_seq, start_pos, k, expansion_count);
    }
    /* Pure NEON */
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_DNA4KMER2_NEON_THRESHOLD) {
        return kmersearch_expand_dna4_kmer2_to_dna2_direct_neon(dna4_seq, start_pos, k, expansion_count);
    }
#endif
    
    /* Fallback to scalar implementation */
    return kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(dna4_seq, start_pos, k, expansion_count);
}



/*
 * Extract k-mers only (without occurrence counts) from DNA2 bit sequence
 * This function is used for frequency analysis phase
 */

/*
 * Encode k-mer-only VarBit into compact KmerData (ignoring occurrence count bits)
 */




/*
 * Cache management functions
 */
/* Cache management functions moved to kmersearch_cache.c */




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
    int shared_count = 0;
    int i, j;
    
    /* Extract k-mers from DNA2 sequence (no degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna2_ngram_key2_direct(sequence, kmersearch_kmer_size, &seq_nkeys);
    if (seq_datum_keys != NULL && seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    /* Extract k-mers from query as ngram_key2 format */
    query_keys = kmersearch_extract_query_ngram_key2(query_string, kmersearch_kmer_size, &query_nkeys);
    
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
    int shared_count = 0;
    int i, j;
    
    /* Extract k-mers from DNA4 sequence (with degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna4_ngram_key2_with_expansion_direct(sequence, kmersearch_kmer_size, &seq_nkeys);
    if (seq_datum_keys != NULL && seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    /* Extract k-mers from query as ngram_key2 format */
    query_keys = kmersearch_extract_query_ngram_key2(query_string, kmersearch_kmer_size, &query_nkeys);
    
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

/* Function moved to kmersearch_gin.c */

/* Function moved to kmersearch_gin.c */

/*
 * K-mer based matching for DNA4 sequences
 */

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
    if (query_len < kmersearch_kmer_size) {
        ereport(ERROR, (errmsg("Query sequence must be at least %d bases long", kmersearch_kmer_size)));
    }
    
    /* Extract k-mers from DNA2 sequence (no degenerate expansion) */
    elog(LOG, "DNA2 Cache: Starting k-mer extraction from sequence");
    seq_datum_keys = kmersearch_extract_dna2_ngram_key2_direct(sequence, kmersearch_kmer_size, &result.seq_nkeys);
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
        query_keys = get_cached_query_kmer(query_string, kmersearch_kmer_size, &result.query_nkeys);
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
    if (query_len < kmersearch_kmer_size) {
        ereport(ERROR, (errmsg("Query sequence must be at least %d bases long", kmersearch_kmer_size)));
    }
    
    /* Extract k-mers from DNA4 sequence (with degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna4_ngram_key2_with_expansion_direct(sequence, kmersearch_kmer_size, &result.seq_nkeys);
    if (seq_datum_keys != NULL && result.seq_nkeys > 0) {
        seq_keys = (VarBit **) palloc(result.seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < result.seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    
    if (seq_keys != NULL && result.seq_nkeys > 0) {
        /* Extract k-mers from query (with degenerate expansion) */
        query_keys = get_cached_query_kmer(query_string, kmersearch_kmer_size, &result.query_nkeys);
        
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


/* Function moved to kmersearch_gin.c */



/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/*
 * Determine the optimal number of parallel workers for k-mer analysis
 */
static int
kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation)
{
    int max_workers;
    int table_size_factor;
    int auto_workers;
    BlockNumber total_blocks;
    
    /* Use the smaller of the two limits */
    max_workers = Min(max_parallel_workers, max_parallel_maintenance_workers);
    
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

/* Function moved to kmersearch_freq.c */
/*
 * Helper functions implementation
 */

/* Functions moved to kmersearch_util.c */



/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_kmer.c */


/* Function moved to kmersearch_freq.c */





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
    kmersearch_highfreq_kmer_cache_free_internal();
}

/*
 * Parallel high-frequency k-mer cache internal functions
 */






/* Function moved to kmersearch_cache.c */

/*
 * Lookup entry in parallel high-frequency k-mer cache
 */





/*
 * SIMD Implementation Functions
 */

/* Scalar implementations (fallback) */




#ifdef __x86_64__
/* AVX2 implementations */



/* AVX512 implementations */



#endif

#ifdef __aarch64__
/* NEON implementations */




/* SVE implementations */



#endif

/*
 * AVX2 K-mer Processing Functions
 */

#ifdef __x86_64__


/* AVX2 optimized version of kmersearch_count_matching_kmer_fast */



/* AVX512 optimized version of kmersearch_count_matching_kmer_fast */
#endif

#ifdef __aarch64__


/* NEON optimized version of kmersearch_count_matching_kmer_fast */
#endif

#ifdef __aarch64__


/* SVE optimized version of kmersearch_count_matching_kmer_fast */

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

/*
 * DNA2 AVX512 implementation - DIRECT BIT MANIPULATION
 */
#endif

/*
 * DNA2 NEON implementation - DIRECT BIT MANIPULATION
 */
#ifdef __x86_64__
/*
 * DNA2 AVX2 implementation - DIRECT BIT MANIPULATION
 */
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
    
    /* Debug logging */
    elog(DEBUG3, "kmersearch_extract_dna2_kmer2_as_uint_direct_scalar: seq_bits=%d, seq_bases=%d, k=%d, max_kmers=%d",
         seq_bits, seq_bases, k, max_kmers);
    
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
        /* Debug first few bytes of data */
        if (seq_bases > 0) {
            elog(DEBUG3, "kmersearch_extract_dna2: First 4 bytes of src_data: %02x %02x %02x %02x",
                 src_data[0], src_data[1], src_data[2], src_data[3]);
        }
        
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
                
                /* Debug first k-mer extraction */
                if (i == 0 && j < k) {
                    elog(DEBUG3, "kmersearch_extract_dna2: pos %d, byte_pos=%d, bit_offset=%d, src_byte=0x%02x, base_bits=%u, kmer_value=%lu",
                         j, byte_pos, bit_offset, src_data[byte_pos], base_bits, (unsigned long)kmer_value);
                }
            }
            
            /* Debug first few k-mers */
            if (i < 3) {
                elog(DEBUG3, "kmersearch_extract_dna2: k-mer[%d] = %u", i, (uint16)kmer_value);
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

/*
 * DNA4 AVX512 implementation - DIRECT DEGENERATE EXPANSION
 */
#endif

/*
 * DNA4 NEON implementation - DIRECT DEGENERATE EXPANSION
 */
#ifdef __aarch64__

/*
 * DNA4 SVE implementation - DIRECT DEGENERATE EXPANSION
 */
#endif

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
    int simd_batch, i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    /* Allocate result array with room for expansions (max 10 expansions per k-mer position) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    /* NEON optimized k-mer extraction with expansion */
    simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
    
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
    int simd_batch, i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return;
    }
    
    /* Allocate result array with room for expansions (max 10 expansions per k-mer position) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    /* SVE optimized k-mer extraction with expansion */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
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
