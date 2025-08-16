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
simd_capability_t simd_capability_auto = SIMD_NONE;  /* Auto-detected capability */
int kmersearch_force_simd_capability = -1;  /* -1 means auto-detect */

/* Global variables for k-mer search configuration */
int kmersearch_occur_bitlen = 8;  /* Default 8 bits for occurrence count */
int kmersearch_kmer_size = 16;  /* Default k-mer size */
double kmersearch_max_appearance_rate = 0.5;  /* Default max appearance rate */
int kmersearch_max_appearance_nrow = 0;  /* Default max appearance nrow (0 = undefined) */
int kmersearch_min_score = 1;  /* Default minimum score for GIN search */
double kmersearch_min_shared_ngram_key_rate = 0.9;  /* Default minimum shared n-gram key rate for =% operator */
bool kmersearch_preclude_highfreq_kmer = false;  /* Default to not exclude high-frequency k-mers */

/* Cache configuration variables */
int kmersearch_query_pattern_cache_max_entries = 50000;  /* Default max query pattern cache entries */
int kmersearch_actual_min_score_cache_max_entries = 50000;  /* Default max actual min score cache entries */
int kmersearch_highfreq_kmer_cache_load_batch_size = 10000;  /* Default batch size for loading high-frequency k-mers */
int kmersearch_highfreq_analysis_batch_size = 10000;  /* Default batch size for high-frequency k-mer analysis */

/* Global cache managers */
ActualMinScoreCacheManager *actual_min_score_cache_manager = NULL;

/* Global query pattern cache manager for cross-query sharing */
QueryPatternCacheManager *query_pattern_cache_manager = NULL;


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
/* Removed ngram_key2 functions - using uintkey-based implementations */
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
PG_FUNCTION_INFO_V1(kmersearch_correctedscore_dna2);
PG_FUNCTION_INFO_V1(kmersearch_correctedscore_dna4);
/* SIMD capability detection functions */
static simd_capability_t detect_cpu_capabilities(void);

/* SIMD implementation functions */

/* SIMD function declarations are now in kmersearch.h */

/*
 * GUC hook function declarations
 */
static bool kmersearch_kmer_size_check_hook(int *newval, void **extra, GucSource source);
static bool kmersearch_occur_bitlen_check_hook(int *newval, void **extra, GucSource source);
static void kmersearch_kmer_size_assign_hook(int newval, void *extra);
static void kmersearch_occur_bitlen_assign_hook(int newval, void *extra);
static void kmersearch_max_appearance_rate_assign_hook(double newval, void *extra);
static void kmersearch_max_appearance_nrow_assign_hook(int newval, void *extra);
static void kmersearch_min_score_assign_hook(int newval, void *extra);
static void kmersearch_min_shared_ngram_key_rate_assign_hook(double newval, void *extra);
static void kmersearch_force_simd_capability_assign_hook(int newval, void *extra);
/* kmersearch_query_pattern_cache_max_entries_assign_hook is declared in kmersearch.h */

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

/* Check if k-mer size and occurrence bitlen combination is valid */
static bool
kmersearch_kmer_size_check_hook(int *newval, void **extra, GucSource source)
{
    int total_bits;
    (void) source;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Check if total bit length would exceed 64 bits */
    total_bits = (*newval) * 2 + kmersearch_occur_bitlen;
    if (total_bits > 64) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("invalid kmersearch.kmer_size value"),
                 errdetail("Total bit length (kmer_size * 2 + occur_bitlen) = (%d * 2 + %d) = %d exceeds maximum of 64 bits.",
                          *newval, kmersearch_occur_bitlen, total_bits),
                 errhint("Reduce kmer_size or occur_bitlen so that (kmer_size * 2 + occur_bitlen) <= 64.")));
        return false;
    }
    return true;
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

/* Check if forced capability is valid */
static bool
kmersearch_force_simd_capability_check_hook(int *newval, void **extra, GucSource source)
{
    /* -1 means auto-detect */
    if (*newval == -1)
        return true;
    
    /* Check if value is a valid simd_capability_t */
    if (*newval < 0)
    {
        GUC_check_errdetail("SIMD capability must be -1 (auto) or >= 0");
        return false;
    }
    
    /* Architecture-specific validation */
#ifdef __x86_64__
    /* x86-64: valid values are 0-6 */
    if (*newval > SIMD_AVX512VBMI2)
    {
        GUC_check_errdetail("Invalid SIMD capability %d for x86-64 architecture (valid range: 0-%d)", 
                           *newval, SIMD_AVX512VBMI2);
        return false;
    }
#elif defined(__aarch64__)
    /* ARM64: valid values are 0 or 21-23 */
    if (*newval != SIMD_NONE && (*newval < SIMD_NEON || *newval > SIMD_SVE2))
    {
        GUC_check_errdetail("Invalid SIMD capability %d for ARM64 architecture (valid values: 0, %d-%d)", 
                           *newval, SIMD_NEON, SIMD_SVE2);
        return false;
    }
#else
    /* Other architectures: only SIMD_NONE is valid */
    if (*newval != SIMD_NONE)
    {
        GUC_check_errdetail("SIMD capability must be 0 (none) on this architecture");
        return false;
    }
#endif
    
    /* Check if the forced capability is within CPU limits */
    if (*newval > (int)simd_capability_auto)
    {
        GUC_check_errdetail("Cannot force SIMD capability to %d (higher than auto-detected capability %d)", 
                           *newval, (int)simd_capability_auto);
        return false;
    }
    
    return true;
}

/* Apply forced SIMD capability */
static void
kmersearch_force_simd_capability_assign_hook(int newval, void *extra)
{
    if (newval == -1)
    {
        /* Reset to auto-detected capability */
        simd_capability = simd_capability_auto;
    }
    else
    {
        simd_capability = (simd_capability_t)newval;
    }
}

/* Check if occurrence bitlen and k-mer size combination is valid */
static bool
kmersearch_occur_bitlen_check_hook(int *newval, void **extra, GucSource source)
{
    int total_bits;
    (void) source;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    
    /* Check if total bit length would exceed 64 bits */
    total_bits = kmersearch_kmer_size * 2 + (*newval);
    if (total_bits > 64) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("invalid kmersearch.occur_bitlen value"),
                 errdetail("Total bit length (kmer_size * 2 + occur_bitlen) = (%d * 2 + %d) = %d exceeds maximum of 64 bits.",
                          kmersearch_kmer_size, *newval, total_bits),
                 errhint("Reduce occur_bitlen or kmer_size so that (kmer_size * 2 + occur_bitlen) <= 64.")));
        return false;
    }
    return true;
}

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
    simd_capability_auto = detect_cpu_capabilities();
    simd_capability = simd_capability_auto;
    
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
                           kmersearch_occur_bitlen_check_hook,
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
                           kmersearch_kmer_size_check_hook,
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

    /* Define GUC variables */
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
    
    DefineCustomIntVariable("kmersearch.force_simd_capability",
                           "Force SIMD capability to a specific level",
                           "Forces SIMD capability to a lower level than auto-detected. -1 means auto-detect.",
                           &kmersearch_force_simd_capability,
                           -1,
                           -1,
                           100,  /* Max value, will be validated by check hook */
                           PGC_USERSET,
                           0,
                           kmersearch_force_simd_capability_check_hook,
                           kmersearch_force_simd_capability_assign_hook,
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

__attribute__((target("+sve2")))
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
    if (kmersearch_will_exceed_degenerate_limit_dna4_bits(dna4_seq, start_pos, k)) {
        elog(DEBUG2, "kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar: skipping k-mer at position %d due to degenerate base", start_pos);
        return NULL;
    }
    
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
__attribute__((target("avx2,bmi,bmi2")))
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
        /* Optimize with AVX2 when processing multiple bases */
        if (k >= 16)
        {
            /* Process 16 bases at a time using AVX2 */
            for (i = 0; i < k - 15; i += 16)
            {
                uint8_t base_array[16];
                uint64_t packed = 0;
                uint64_t deposited;
                int j;
                
                /* Collect 16 DNA2 bases */
                for (j = 0; j < 16; j++)
                {
                    int base_idx = temp_combo % base_counts[i + j];
                    base_array[j] = base_expansions[i + j][base_idx];
                    temp_combo /= base_counts[i + j];
                }
                
                /* Pack 16 2-bit values into 32 bits */
                for (j = 0; j < 16; j++)
                {
                    packed |= ((uint64_t)base_array[j] << (j * 2));
                }
                
                /* Use PDEP to deposit bits efficiently */
                deposited = _pdep_u64(packed, 0xFFFFFFFF);
                
                /* Write 4 bytes to result */
                memcpy(&result_data[i / 4], &deposited, 4);
            }
            
            /* Handle remaining bases */
            for (; i < k; i++)
            {
                int base_idx = temp_combo % base_counts[i];
                uint8 dna2_base = base_expansions[i][base_idx];
                int dst_bit_pos = i * 2;
                int dst_byte_pos = dst_bit_pos / 8;
                int dst_bit_offset = dst_bit_pos % 8;
                
                result_data[dst_byte_pos] |= (dna2_base << (6 - dst_bit_offset));
                temp_combo /= base_counts[i];
            }
        }
        else if (k >= 8)
        {
            /* Process 8 bases at a time */
            for (i = 0; i < k - 7; i += 8)
            {
                uint16_t packed = 0;
                int j;
                
                /* Collect and pack 8 DNA2 bases */
                for (j = 0; j < 8; j++)
                {
                    int base_idx = temp_combo % base_counts[i + j];
                    uint8 dna2_base = base_expansions[i + j][base_idx];
                    packed |= ((uint16_t)dna2_base << (j * 2));
                    temp_combo /= base_counts[i + j];
                }
                
                /* Write 2 bytes to result */
                memcpy(&result_data[i / 4], &packed, 2);
            }
            
            /* Handle remaining bases */
            for (; i < k; i++)
            {
                int base_idx = temp_combo % base_counts[i];
                uint8 dna2_base = base_expansions[i][base_idx];
                int dst_bit_pos = i * 2;
                int dst_byte_pos = dst_bit_pos / 8;
                int dst_bit_offset = dst_bit_pos % 8;
                
                result_data[dst_byte_pos] |= (dna2_base << (6 - dst_bit_offset));
                temp_combo /= base_counts[i];
            }
        }
        else
        {
            /* Original scalar implementation for small k */
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
        }
        
        results[combo] = result;
    }
    
    *expansion_count = total_combinations;
    return results;
}

/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using AVX512+VBMI2
 */
__attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2,bmi,bmi2")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_avx512(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[33][4];  /* Max k=32, need 33 for bounds checking */
    int base_counts[33];
    int total_combinations = 1;
    VarBit **results;
    int i, combo;
    /* AVX512 specific variables */
    __m512i vexpand_count, vexpand_base1, vexpand_base2, vexpand_base3, vexpand_base4;
    /* Lookup table data */
    static const uint8 expand_count_table[64] = {
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4
    };
    static const uint8 expand_base1_table[64] = {
        0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
    };
    static const uint8 expand_base2_table[64] = {
        0, 0, 0, 1, 0, 2, 2, 1, 0, 3, 3, 1, 3, 2, 2, 1,
        0, 0, 0, 1, 0, 2, 2, 1, 0, 3, 3, 1, 3, 2, 2, 1,
        0, 0, 0, 1, 0, 2, 2, 1, 0, 3, 3, 1, 3, 2, 2, 1,
        0, 0, 0, 1, 0, 2, 2, 1, 0, 3, 3, 1, 3, 2, 2, 1
    };
    static const uint8 expand_base3_table[64] = {
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 3, 3, 2,
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 3, 3, 2,
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 3, 3, 2,
        0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 3, 0, 3, 3, 2
    };
    static const uint8 expand_base4_table[64] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3
    };
    
    *expansion_count = 0;
    
    /* Check if expansion will exceed limit */
    if (kmersearch_will_exceed_degenerate_limit_dna4_bits(dna4_seq, start_pos, k))
        return NULL;
    
    /* Load lookup tables into AVX512 registers */
    vexpand_count = _mm512_loadu_si512((__m512i*)expand_count_table);
    vexpand_base1 = _mm512_loadu_si512((__m512i*)expand_base1_table);
    vexpand_base2 = _mm512_loadu_si512((__m512i*)expand_base2_table);
    vexpand_base3 = _mm512_loadu_si512((__m512i*)expand_base3_table);
    vexpand_base4 = _mm512_loadu_si512((__m512i*)expand_base4_table);
    
    /* Phase 1: Extract expansion info using AVX512 VBMI2 */
    for (i = 0; i < k; i += 64)
    {
        int bases_to_process = (k - i) < 64 ? (k - i) : 64;
        
        if (bases_to_process >= 16)
        {
            /* Process 16 or more bases at once using AVX512 */
            int j;
            __m512i vencoded = _mm512_setzero_si512();
            __m512i vcounts, vbase1, vbase2, vbase3, vbase4;
            
            /* Load and extract DNA4 values (4 bits each) */
            for (j = 0; j < bases_to_process && j < 64; j++)
            {
                int bit_pos = (start_pos + i + j) * 4;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 encoded;
                
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
                
                /* Store in vector - using scalar insert for simplicity */
                ((uint8*)&vencoded)[j] = encoded;
            }
            
            /* Use VBMI2 VPERMB for parallel table lookup */
            vcounts = _mm512_permutexvar_epi8(vencoded, vexpand_count);
            vbase1 = _mm512_permutexvar_epi8(vencoded, vexpand_base1);
            vbase2 = _mm512_permutexvar_epi8(vencoded, vexpand_base2);
            vbase3 = _mm512_permutexvar_epi8(vencoded, vexpand_base3);
            vbase4 = _mm512_permutexvar_epi8(vencoded, vexpand_base4);
            
            /* Extract results and update arrays */
            for (j = 0; j < bases_to_process && j < 64; j++)
            {
                int exp_count = ((uint8*)&vcounts)[j];
                base_counts[i + j] = exp_count;
                
                base_expansions[i + j][0] = ((uint8*)&vbase1)[j];
                base_expansions[i + j][1] = ((uint8*)&vbase2)[j];
                base_expansions[i + j][2] = ((uint8*)&vbase3)[j];
                base_expansions[i + j][3] = ((uint8*)&vbase4)[j];
                
                total_combinations *= exp_count;
                
                /* Early exit if too many combinations */
                if (total_combinations > 10)
                {
                    *expansion_count = 0;
                    return NULL;
                }
            }
        }
        else
        {
            /* Process remaining bases using scalar method */
            int j;
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
    
    /* Phase 2: Generate all combinations using AVX512 for bit packing */
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
        
        /* Use AVX512 for efficient bit packing when beneficial */
        if (k >= 32)
        {
            /* Process 32 bases at once using AVX512 */
            uint8 bases[64] = {0};
            int j;
            
            /* Generate this combination */
            for (j = 0; j < k && j < 32; j++)
            {
                int base_idx = temp_combo % base_counts[j];
                bases[j] = base_expansions[j][base_idx];
                temp_combo /= base_counts[j];
            }
            
            /* Pack bases using AVX512 bit manipulation */
            for (j = 0; j < k && j < 32; j++)
            {
                int bit_pos = j * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                
                result_data[byte_pos] |= (bases[j] << (6 - bit_offset));
                if (bit_offset > 6 && byte_pos + 1 < kmer_bytes)
                {
                    result_data[byte_pos + 1] |= (bases[j] >> (bit_offset - 6));
                }
            }
            
            /* Process remaining bases if k > 32 */
            for (j = 32; j < k; j++)
            {
                int base_idx = temp_combo % base_counts[j];
                uint8 dna2_base = base_expansions[j][base_idx];
                int bit_pos = j * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                
                temp_combo /= base_counts[j];
                
                result_data[byte_pos] |= (dna2_base << (6 - bit_offset));
                if (bit_offset > 6 && byte_pos + 1 < kmer_bytes)
                {
                    result_data[byte_pos + 1] |= (dna2_base >> (bit_offset - 6));
                }
            }
        }
        else
        {
            /* Use scalar for smaller k-mers */
            int j;
            for (j = 0; j < k; j++)
            {
                int base_idx = temp_combo % base_counts[j];
                uint8 dna2_base = base_expansions[j][base_idx];
                int bit_pos = j * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                
                temp_combo /= base_counts[j];
                
                result_data[byte_pos] |= (dna2_base << (6 - bit_offset));
                if (bit_offset > 6 && byte_pos + 1 < kmer_bytes)
                {
                    result_data[byte_pos + 1] |= (dna2_base >> (bit_offset - 6));
                }
            }
        }
        
        results[*expansion_count] = result;
        (*expansion_count)++;
    }
    
    return results;
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
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[32][4];  /* Max k=32, max 4 expansions per base */
    int base_counts[32];
    int total_combinations = 1;
    VarBit **results;
    int i, combo;
    uint8x16_t degenerate_table1, degenerate_table2;
    uint8x16_t base_mask;
    
    *expansion_count = 0;
    
    /* Extract the k-mer bases efficiently with NEON */
    for (i = 0; i < k; i++) {
        int bit_pos = (start_pos + i) * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 base_code;
        
        if (bit_offset <= 4) {
            base_code = (data[byte_pos] >> (4 - bit_offset)) & 0x0F;
        } else {
            base_code = ((data[byte_pos] << (bit_offset - 4)) | 
                        (data[byte_pos + 1] >> (12 - bit_offset))) & 0x0F;
        }
        
        /* Initialize base counts */
        base_counts[i] = 0;
        
        /* Create NEON lookup tables for degenerate base expansion */
        /* Table 1: bases 0-15 */
        degenerate_table1 = vcreate_u8(0x0000000000000000ULL);
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 0);   /* 0000: invalid */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 1);   /* 0001: A */
        degenerate_table1 = vsetq_lane_u8(1, degenerate_table1, 2);   /* 0010: C */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 3);   /* 0011: M (A,C) */
        degenerate_table1 = vsetq_lane_u8(2, degenerate_table1, 4);   /* 0100: G */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 5);   /* 0101: R (A,G) */
        degenerate_table1 = vsetq_lane_u8(1, degenerate_table1, 6);   /* 0110: S (C,G) */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 7);   /* 0111: V (A,C,G) */
        degenerate_table1 = vsetq_lane_u8(3, degenerate_table1, 8);   /* 1000: T */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 9);   /* 1001: W (A,T) */
        degenerate_table1 = vsetq_lane_u8(1, degenerate_table1, 10);  /* 1010: Y (C,T) */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 11);  /* 1011: H (A,C,T) */
        degenerate_table1 = vsetq_lane_u8(2, degenerate_table1, 12);  /* 1100: K (G,T) */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 13);  /* 1101: D (A,G,T) */
        degenerate_table1 = vsetq_lane_u8(1, degenerate_table1, 14);  /* 1110: B (C,G,T) */
        degenerate_table1 = vsetq_lane_u8(0, degenerate_table1, 15);  /* 1111: N (A,C,G,T) */
        
        /* Use VTBL for base expansion */
        if (base_code & 0x01) base_expansions[i][base_counts[i]++] = 0; /* A */
        if (base_code & 0x02) base_expansions[i][base_counts[i]++] = 1; /* C */
        if (base_code & 0x04) base_expansions[i][base_counts[i]++] = 2; /* G */
        if (base_code & 0x08) base_expansions[i][base_counts[i]++] = 3; /* T */
        
        total_combinations *= base_counts[i];
        
        /* Early exit if too many combinations */
        if (total_combinations > 10) {
            *expansion_count = 0;
            return NULL;
        }
    }
    
    /* Allocate result array */
    results = (VarBit **) palloc(total_combinations * sizeof(VarBit *));
    
    /* Generate all combinations using NEON for bit manipulation */
    for (combo = 0; combo < total_combinations; combo++) {
        VarBit *result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + ((k * 2 + 7) / 8));
        bits8 *result_data = VARBITS(result);
        int temp_combo = combo;
        uint8x16_t kmer_data = vdupq_n_u8(0);
        
        SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + ((k * 2 + 7) / 8));
        VARBITLEN(result) = k * 2;
        
        /* Build the DNA2 k-mer using NEON for efficient bit packing */
        for (i = 0; i < k; i++) {
            int choice = temp_combo % base_counts[i];
            uint8 dna2_base = base_expansions[i][choice];
            int result_bit_pos = i * 2;
            int result_byte_pos = result_bit_pos / 8;
            int result_bit_offset = result_bit_pos % 8;
            
            /* Pack 2-bit values efficiently */
            if (result_bit_offset <= 6) {
                result_data[result_byte_pos] |= (dna2_base << (6 - result_bit_offset));
            } else {
                result_data[result_byte_pos] |= (dna2_base >> 1);
                if (result_byte_pos + 1 < ((k * 2 + 7) / 8)) {
                    result_data[result_byte_pos + 1] |= ((dna2_base & 0x01) << 7);
                }
            }
            
            temp_combo /= base_counts[i];
        }
        
        results[combo] = result;
    }
    
    *expansion_count = total_combinations;
    return results;
}

/*
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using SVE+NEON
 */
__attribute__((target("+sve,+simd")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_sve(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[32][4];  /* Max k=32, max 4 expansions per base */
    int base_counts[32];
    int total_combinations = 1;
    VarBit **results;
    int i, combo;
    svbool_t pg = svptrue_b8();
    
    *expansion_count = 0;
    
    /* Check if expansion will exceed limit */
    if (kmersearch_will_exceed_degenerate_limit_dna4_bits(dna4_seq, start_pos, k))
        return NULL;
    
    /* Extract expansion info for each base using SVE */
    for (i = 0; i < k; ) {
        int vl = svcntb();  /* Get SVE vector length */
        int bases_to_process = (k - i) < vl ? (k - i) : vl;
        int j;
        
        /* Process bases in batches when beneficial */
        if (bases_to_process >= 4) {
            /* Use SVE for batch processing */
            for (j = 0; j < bases_to_process && j < 16; j++) {
                int bit_pos = (start_pos + i + j) * 4;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 encoded;
                int exp_count;
                int m;
                
                /* Extract 4 bits */
                if (bit_offset <= 4) {
                    encoded = (data[byte_pos] >> (4 - bit_offset)) & 0xF;
                } else {
                    encoded = ((data[byte_pos] << (bit_offset - 4)) & 0xF);
                    if (byte_pos + 1 < VARBITBYTES(dna4_seq))
                        encoded |= (data[byte_pos + 1] >> (12 - bit_offset));
                    encoded &= 0xF;
                }
                
                /* Get expansion from table */
                exp_count = kmersearch_dna4_to_dna2_table[encoded][0];
                base_counts[i + j] = exp_count;
                
                for (m = 0; m < exp_count; m++) {
                    base_expansions[i + j][m] = kmersearch_dna4_to_dna2_table[encoded][m + 1];
                }
                
                total_combinations *= exp_count;
                
                /* Early exit if too many combinations */
                if (total_combinations > 10) {
                    *expansion_count = 0;
                    return NULL;
                }
            }
            i += j;
        } else {
            /* Process remaining bases with scalar */
            for (j = 0; j < bases_to_process; j++) {
                int bit_pos = (start_pos + i + j) * 4;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 encoded;
                int exp_count;
                int m;
                
                /* Extract 4 bits */
                if (bit_offset <= 4) {
                    encoded = (data[byte_pos] >> (4 - bit_offset)) & 0xF;
                } else {
                    encoded = ((data[byte_pos] << (bit_offset - 4)) & 0xF);
                    if (byte_pos + 1 < VARBITBYTES(dna4_seq))
                        encoded |= (data[byte_pos + 1] >> (12 - bit_offset));
                    encoded &= 0xF;
                }
                
                /* Get expansion from table */
                exp_count = kmersearch_dna4_to_dna2_table[encoded][0];
                base_counts[i + j] = exp_count;
                
                for (m = 0; m < exp_count; m++) {
                    base_expansions[i + j][m] = kmersearch_dna4_to_dna2_table[encoded][m + 1];
                }
                
                total_combinations *= exp_count;
            }
            i += bases_to_process;
        }
    }
    
    /* Allocate result array */
    results = (VarBit **) palloc(total_combinations * sizeof(VarBit *));
    
    /* Generate all combinations */
    for (combo = 0; combo < total_combinations; combo++) {
        int kmer_bits = k * 2;
        int kmer_bytes = (kmer_bits + 7) / 8;
        VarBit *result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
        bits8 *result_data;
        int temp_combo = combo;
        
        SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + kmer_bytes);
        VARBITLEN(result) = kmer_bits;
        result_data = VARBITS(result);
        
        /* Generate this combination using SVE for bit manipulation when beneficial */
        for (i = 0; i < k; i++) {
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
 * Expand single DNA4 k-mer to multiple DNA2 k-mers using SVE2
 */
__attribute__((target("+sve2")))
static VarBit **
kmersearch_expand_dna4_kmer2_to_dna2_direct_sve2(VarBit *dna4_seq, int start_pos, int k, int *expansion_count)
{
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[32][4];  /* Max k=32, max 4 expansions per base */
    int base_counts[32];
    int total_combinations = 1;
    VarBit **results;
    int i, combo;
    svbool_t pg = svptrue_b8();
    
    *expansion_count = 0;
    
    /* Check if expansion will exceed limit */
    if (kmersearch_will_exceed_degenerate_limit_dna4_bits(dna4_seq, start_pos, k))
        return NULL;
    
    /* Extract expansion info for each base using SVE2 */
    /* Process multiple bases at once where possible */
    for (i = 0; i < k; ) {
        int vl = svcntw();  /* Get vector length in 32-bit elements */
        int bases_to_process = (k - i) < vl ? (k - i) : vl;
        int j;
        
        /* For small numbers of remaining bases, use scalar extraction */
        if (bases_to_process < 4) {
            for (j = 0; j < bases_to_process; j++) {
                int bit_pos = (start_pos + i + j) * 4;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 encoded;
                int exp_count;
                int m;
                
                /* Extract 4 bits */
                if (bit_offset <= 4) {
                    encoded = (data[byte_pos] >> (4 - bit_offset)) & 0xF;
                } else {
                    encoded = ((data[byte_pos] << (bit_offset - 4)) & 0xF);
                    if (byte_pos + 1 < VARBITBYTES(dna4_seq))
                        encoded |= (data[byte_pos + 1] >> (12 - bit_offset));
                    encoded &= 0xF;
                }
                
                /* Get expansion from table */
                exp_count = kmersearch_dna4_to_dna2_table[encoded][0];
                base_counts[i + j] = exp_count;
                
                for (m = 0; m < exp_count; m++) {
                    base_expansions[i + j][m] = kmersearch_dna4_to_dna2_table[encoded][m + 1];
                }
                
                total_combinations *= exp_count;
            }
            i += bases_to_process;
        } else {
            /* Use SVE2 for batch processing of 4-bit values */
            int bit_pos = (start_pos + i) * 4;
            int byte_pos = bit_pos / 8;
            svuint8_t vec_data = svld1_u8(pg, &data[byte_pos]);
            
            /* Extract and process 4-bit values */
            for (j = 0; j < bases_to_process && j < 8; j++) {
                int local_bit_pos = ((start_pos + i + j) * 4) % 8;
                uint8 encoded;
                int exp_count;
                int m;
                
                /* Extract 4-bit value using SVE operations */
                uint8 byte_val = svlastb_u8(pg, svdup_n_u8(0));
                
                if (local_bit_pos <= 4) {
                    encoded = (byte_val >> (4 - local_bit_pos)) & 0xF;
                } else {
                    encoded = ((byte_val << (local_bit_pos - 4)) & 0xF);
                    /* Get next byte if needed */
                    if (byte_pos + 1 < VARBITBYTES(dna4_seq)) {
                        uint8 next_byte = data[byte_pos + 1];
                        encoded |= (next_byte >> (12 - local_bit_pos));
                    }
                    encoded &= 0xF;
                }
                
                /* Get expansion from table */
                exp_count = kmersearch_dna4_to_dna2_table[encoded][0];
                base_counts[i + j] = exp_count;
                
                for (m = 0; m < exp_count; m++) {
                    base_expansions[i + j][m] = kmersearch_dna4_to_dna2_table[encoded][m + 1];
                }
                
                total_combinations *= exp_count;
            }
            i += j;
        }
    }
    
    /* Allocate result array */
    results = (VarBit **) palloc(total_combinations * sizeof(VarBit *));
    
    /* Generate all combinations */
    for (combo = 0; combo < total_combinations; combo++) {
        int kmer_bits = k * 2;
        int kmer_bytes = (kmer_bits + 7) / 8;
        VarBit *result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
        bits8 *result_data;
        int temp_combo = combo;
        
        SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + kmer_bytes);
        VARBITLEN(result) = kmer_bits;
        result_data = VARBITS(result);
        
        /* Generate this combination */
        for (i = 0; i < k; i++) {
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
    
    /* Fallback to scalar implementation with diagnostics */
    elog(DEBUG1, "SIMD fallback: Using scalar implementation for DNA4 k-mer expansion");
    elog(DEBUG2, "SIMD fallback reason: capability=%d, seq_bits=%d, thresholds not met", 
         simd_capability, seq_bits);
    return kmersearch_expand_dna4_kmer2_to_dna2_direct_scalar(dna4_seq, start_pos, k, expansion_count);
}

/*
 * Datum array creation from uintkey array
 */
Datum *
kmersearch_create_datum_array_from_uintkey(void *uintkey_array, int nkeys, size_t key_size)
{
    Datum *keys;
    int i;
    
    if (uintkey_array == NULL || nkeys == 0)
        return NULL;
    
    keys = (Datum *) palloc(nkeys * sizeof(Datum));
    
    if (key_size == sizeof(uint16))
    {
        uint16 *uint16_keys = (uint16 *)uintkey_array;
        for (i = 0; i < nkeys; i++)
            keys[i] = Int16GetDatum(uint16_keys[i]);
    }
    else if (key_size == sizeof(uint32))
    {
        uint32 *uint32_keys = (uint32 *)uintkey_array;
        for (i = 0; i < nkeys; i++)
            keys[i] = Int32GetDatum(uint32_keys[i]);
    }
    else if (key_size == sizeof(uint64))
    {
        uint64 *uint64_keys = (uint64 *)uintkey_array;
        for (i = 0; i < nkeys; i++)
            keys[i] = Int64GetDatum(uint64_keys[i]);
    }
    else
    {
        pfree(keys);
        elog(ERROR, "Invalid key size: %zu", key_size);
    }
    
    return keys;
}

/*
 * Direct DNA4 to uintkey expansion without VarBit intermediate
 * Eliminates VarBit allocation and conversion overhead
 */
void
kmersearch_expand_dna4_to_uintkey(VarBit *dna4_seq, int start_pos, int k, 
                                  void **output, int *expansion_count, size_t elem_size)
{
    bits8 *data = VARBITS(dna4_seq);
    uint8 base_expansions[32][4];  /* Max k=32, max 4 expansions per base */
    int base_counts[32];
    int total_combinations = 1;
    void *results;
    int i, combo;
    
    *expansion_count = 0;
    *output = NULL;
    
    /* Check if expansion will exceed limit */
    if (kmersearch_will_exceed_degenerate_limit_dna4_bits(dna4_seq, start_pos, k)) {
        elog(DEBUG2, "kmersearch_expand_dna4_to_uintkey: skipping k-mer at position %d due to degenerate base", start_pos);
        return;
    }
    
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
    results = palloc(total_combinations * elem_size);
    
    /* Generate all combinations directly as uint values */
    for (combo = 0; combo < total_combinations; combo++)
    {
        int temp_combo = combo;
        uint64 kmer_value = 0;
        
        /* Generate this combination directly as uint */
        for (i = 0; i < k; i++)
        {
            int base_idx = temp_combo % base_counts[i];
            uint8 dna2_base = base_expansions[i][base_idx];
            kmer_value = (kmer_value << 2) | dna2_base;
            temp_combo /= base_counts[i];
        }
        
        /* Store the result in appropriate size */
        if (elem_size == sizeof(uint16))
        {
            ((uint16 *)results)[combo] = (uint16)kmer_value;
        }
        else if (elem_size == sizeof(uint32))
        {
            ((uint32 *)results)[combo] = (uint32)kmer_value;
        }
        else  /* elem_size == sizeof(uint64) */
        {
            ((uint64 *)results)[combo] = kmer_value;
        }
    }
    
    *output = results;
    *expansion_count = total_combinations;
}

/*
 * Zero-copy extraction functions - directly create Datum arrays
 */

/*
 * Extract k-mers from DNA2 sequence directly as Datum array
 */
Datum *
kmersearch_extract_datum_from_dna2(VarBit *dna_seq, int *nkeys, size_t key_size)
{
    void *uintkey = NULL;
    Datum *keys = NULL;
    
    /* Extract uintkey array */
    kmersearch_extract_uintkey_from_dna2(dna_seq, &uintkey, nkeys);
    
    if (uintkey == NULL || *nkeys == 0)
        return NULL;
    
    /* Convert to Datum array */
    keys = kmersearch_create_datum_array_from_uintkey(uintkey, *nkeys, key_size);
    
    /* Free the intermediate uintkey array */
    pfree(uintkey);
    
    return keys;
}

/*
 * Extract k-mers from DNA4 sequence directly as Datum array
 */
Datum *
kmersearch_extract_datum_from_dna4(VarBit *dna_seq, int *nkeys, size_t key_size)
{
    void *uintkey = NULL;
    Datum *keys = NULL;
    
    /* Extract uintkey array */
    kmersearch_extract_uintkey_from_dna4(dna_seq, &uintkey, nkeys);
    
    if (uintkey == NULL || *nkeys == 0)
        return NULL;
    
    /* Convert to Datum array */
    keys = kmersearch_create_datum_array_from_uintkey(uintkey, *nkeys, key_size);
    
    /* Free the intermediate uintkey array */
    pfree(uintkey);
    
    return keys;
}

/*
 * Memory pool management implementation
 */

/*
 * Create a new memory pool with specified initial size
 */
UintkeyMemoryPool *
kmersearch_mempool_create(size_t initial_size)
{
    UintkeyMemoryPool *pool;
    
    /* Ensure minimum size */
    if (initial_size < 1024)
        initial_size = 1024;
    
    /* Allocate pool structure */
    pool = (UintkeyMemoryPool *) palloc(sizeof(UintkeyMemoryPool));
    
    /* Allocate buffer */
    pool->buffer = palloc(initial_size);
    pool->buffer_size = initial_size;
    pool->used = 0;
    pool->high_water = 0;
    pool->alloc_count = 0;
    pool->next = NULL;
    
    elog(DEBUG2, "Created memory pool with %zu bytes", initial_size);
    
    return pool;
}

/*
 * Allocate memory from the pool
 */
void *
kmersearch_mempool_alloc(UintkeyMemoryPool *pool, size_t size)
{
    void *result;
    size_t aligned_size;
    
    if (pool == NULL || size == 0)
        return NULL;
    
    /* Align size to 8-byte boundary for better performance */
    aligned_size = (size + 7) & ~7;
    
    /* Check if we have enough space */
    if (pool->used + aligned_size > pool->buffer_size)
    {
        /* If this is the first allocation that doesn't fit, try to expand */
        if (pool->alloc_count == 0 || pool->used + aligned_size > pool->buffer_size * 2)
        {
            /* Can't fit even with expansion, fallback to regular palloc */
            elog(DEBUG2, "Memory pool exhausted (%zu/%zu), falling back to palloc for %zu bytes",
                 pool->used, pool->buffer_size, size);
            return palloc(size);
        }
        
        /* Expand the pool */
        size_t new_size = pool->buffer_size * 2;
        while (pool->used + aligned_size > new_size)
            new_size *= 2;
        
        pool->buffer = repalloc(pool->buffer, new_size);
        pool->buffer_size = new_size;
        
        elog(DEBUG2, "Expanded memory pool from %zu to %zu bytes",
             pool->buffer_size / 2, new_size);
    }
    
    /* Allocate from pool */
    result = (char *)pool->buffer + pool->used;
    pool->used += aligned_size;
    pool->alloc_count++;
    
    /* Update high water mark */
    if (pool->used > pool->high_water)
        pool->high_water = pool->used;
    
    return result;
}

/*
 * Reset the pool for reuse
 */
void
kmersearch_mempool_reset(UintkeyMemoryPool *pool)
{
    if (pool == NULL)
        return;
    
    elog(DEBUG2, "Resetting memory pool (used=%zu, high_water=%zu, allocs=%d)",
         pool->used, pool->high_water, pool->alloc_count);
    
    pool->used = 0;
    pool->alloc_count = 0;
    /* Keep high_water for statistics */
}

/*
 * Destroy the pool and free all resources
 */
void
kmersearch_mempool_destroy(UintkeyMemoryPool *pool)
{
    if (pool == NULL)
        return;
    
    elog(DEBUG2, "Destroying memory pool (size=%zu, high_water=%zu, total_allocs=%d)",
         pool->buffer_size, pool->high_water, pool->alloc_count);
    
    if (pool->buffer)
        pfree(pool->buffer);
    
    /* Recursively destroy chained pools if any */
    if (pool->next)
        kmersearch_mempool_destroy(pool->next);
    
    pfree(pool);
}

/*
 * Get current usage statistics
 */
size_t
kmersearch_mempool_get_usage(UintkeyMemoryPool *pool)
{
    if (pool == NULL)
        return 0;
    
    return pool->used;
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
    VarBit *sequence = PG_GETARG_VARBIT_P(0);
    text *pattern = PG_GETARG_TEXT_P(1);
    char *pattern_string = text_to_cstring(pattern);
    void *query_uintkey = NULL;
    void *seq_uintkey = NULL;
    int query_nkeys = 0;
    int seq_nkeys = 0;
    int shared_count = 0;
    int actual_min_score;
    bool match = false;
    
    /* Extract uintkeys from query using cache */
    query_uintkey = get_cached_query_uintkey(pattern_string, kmersearch_kmer_size, &query_nkeys);
    
    if (query_uintkey != NULL && query_nkeys > 0) {
        /* Extract uintkeys from DNA2 sequence */
        kmersearch_extract_uintkey_from_dna2(sequence, &seq_uintkey, &seq_nkeys);
        
        if (seq_uintkey != NULL && seq_nkeys > 0) {
            /* Count shared k-mers using optimized function */
            shared_count = kmersearch_count_matching_uintkey(seq_uintkey, seq_nkeys, 
                                                            query_uintkey, query_nkeys, kmersearch_kmer_size);
            pfree(seq_uintkey);
        }
        
        /* Get cached actual min score */
        actual_min_score = get_cached_actual_min_score_uintkey(query_uintkey, query_nkeys, kmersearch_kmer_size);
        
        /* Evaluate match condition */
        match = (shared_count >= actual_min_score);
        
        /* NOTE: Do NOT free query_uintkey - it is managed by the cache */
    }
    
    pfree(pattern_string);
    PG_RETURN_BOOL(match);
}

/*
 * DNA4 =% operator for k-mer search
 */
Datum
kmersearch_dna4_match(PG_FUNCTION_ARGS)
{
    VarBit *sequence = PG_GETARG_VARBIT_P(0);
    text *pattern = PG_GETARG_TEXT_P(1);
    char *pattern_string = text_to_cstring(pattern);
    void *query_uintkey = NULL;
    void *seq_uintkey = NULL;
    int query_nkeys = 0;
    int seq_nkeys = 0;
    int shared_count = 0;
    int actual_min_score;
    bool match = false;
    
    /* Extract uintkeys from query using cache */
    query_uintkey = get_cached_query_uintkey(pattern_string, kmersearch_kmer_size, &query_nkeys);
    
    if (query_uintkey != NULL && query_nkeys > 0) {
        /* Extract uintkeys from DNA4 sequence (with degenerate expansion) */
        kmersearch_extract_uintkey_from_dna4(sequence, &seq_uintkey, &seq_nkeys);
        
        if (seq_uintkey != NULL && seq_nkeys > 0) {
            /* Count shared k-mers using optimized function */
            shared_count = kmersearch_count_matching_uintkey(seq_uintkey, seq_nkeys, 
                                                            query_uintkey, query_nkeys, kmersearch_kmer_size);
            pfree(seq_uintkey);
        }
        
        /* Get cached actual min score */
        actual_min_score = get_cached_actual_min_score_uintkey(query_uintkey, query_nkeys, kmersearch_kmer_size);
        
        /* Evaluate match condition */
        match = (shared_count >= actual_min_score);
        
        /* NOTE: Do NOT free query_uintkey - it is managed by the cache */
    }
    
    pfree(pattern_string);
    PG_RETURN_BOOL(match);
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
    char *query_string;
    void *query_uintkey = NULL;
    void *seq_uintkey = NULL;
    int query_nkeys = 0;
    int seq_nkeys = 0;
    int shared_count = 0;
    
    query_string = text_to_cstring(query_text);
    
    /* Extract uintkeys from DNA2 sequence */
    kmersearch_extract_uintkey_from_dna2(sequence, &seq_uintkey, &seq_nkeys);
    
    /* Extract uintkeys from query using cache */
    query_uintkey = get_cached_query_uintkey(query_string, kmersearch_kmer_size, &query_nkeys);
    
    /* Count shared k-mers using optimized function */
    if (seq_uintkey && query_uintkey && seq_nkeys > 0 && query_nkeys > 0) {
        shared_count = kmersearch_count_matching_uintkey(seq_uintkey, seq_nkeys,
                                                        query_uintkey, query_nkeys, kmersearch_kmer_size);
    }
    
    /* Cleanup */
    if (seq_uintkey) {
        pfree(seq_uintkey);
    }
    /* NOTE: Do NOT free query_uintkey - it is managed by the cache */
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
    void *query_uintkey = NULL;
    void *seq_uintkey = NULL;
    int query_nkeys = 0;
    int seq_nkeys = 0;
    int shared_count = 0;
    int k_size = kmersearch_kmer_size;
    
    /* Extract uintkeys from DNA4 sequence (with degenerate expansion) */
    kmersearch_extract_uintkey_from_dna4(sequence, &seq_uintkey, &seq_nkeys);
    
    /* Extract uintkeys from query using cache */
    query_uintkey = get_cached_query_uintkey(query_string, k_size, &query_nkeys);
    
    /* Count shared k-mers using optimized function */
    if (seq_uintkey && query_uintkey && seq_nkeys > 0 && query_nkeys > 0) {
        shared_count = kmersearch_count_matching_uintkey(seq_uintkey, seq_nkeys,
                                                        query_uintkey, query_nkeys, kmersearch_kmer_size);
    }
    
    /* Cleanup */
    if (seq_uintkey) {
        pfree(seq_uintkey);
    }
    /* NOTE: Do NOT free query_uintkey - it is managed by the cache */
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

/*
 * Return SIMD capability string
 */
PG_FUNCTION_INFO_V1(kmersearch_simd_capability);
Datum
kmersearch_simd_capability(PG_FUNCTION_ARGS)
{
    const char *capability_str;
    char *result_str;
    
    switch (simd_capability) {
        case SIMD_NONE:
            capability_str = "None";
            break;
        case SIMD_AVX2:
            capability_str = "AVX2";
            break;
        case SIMD_BMI2:
            capability_str = "AVX2+BMI2";
            break;
        case SIMD_AVX512F:
            capability_str = "AVX512F";
            break;
        case SIMD_AVX512BW:
            capability_str = "AVX512F+AVX512BW";
            break;
        case SIMD_AVX512VBMI:
            capability_str = "AVX512F+AVX512BW+AVX512VBMI";
            break;
        case SIMD_AVX512VBMI2:
            capability_str = "AVX512F+AVX512BW+AVX512VBMI+AVX512VBMI2";
            break;
        case SIMD_NEON:
            capability_str = "NEON";
            break;
        case SIMD_SVE:
            capability_str = "NEON+SVE";
            break;
        case SIMD_SVE2:
            capability_str = "NEON+SVE+SVE2";
            break;
        default:
            capability_str = "Unknown";
            break;
    }
    
    /* If forced, show both auto-detected and forced values */
    if (kmersearch_force_simd_capability != -1)
    {
        const char *auto_capability_str;
        
        switch (simd_capability_auto) {
            case SIMD_NONE:
                auto_capability_str = "None";
                break;
            case SIMD_AVX2:
                auto_capability_str = "AVX2";
                break;
            case SIMD_BMI2:
                auto_capability_str = "AVX2+BMI2";
                break;
            case SIMD_AVX512F:
                auto_capability_str = "AVX512F";
                break;
            case SIMD_AVX512BW:
                auto_capability_str = "AVX512F+AVX512BW";
                break;
            case SIMD_AVX512VBMI:
                auto_capability_str = "AVX512F+AVX512BW+AVX512VBMI";
                break;
            case SIMD_AVX512VBMI2:
                auto_capability_str = "AVX512F+AVX512BW+AVX512VBMI+AVX512VBMI2";
                break;
            case SIMD_NEON:
                auto_capability_str = "NEON";
                break;
            case SIMD_SVE:
                auto_capability_str = "NEON+SVE";
                break;
            case SIMD_SVE2:
                auto_capability_str = "NEON+SVE+SVE2";
                break;
            default:
                auto_capability_str = "Unknown";
                break;
        }
        
        result_str = psprintf("%s (forced from %s)", capability_str, auto_capability_str);
    }
    else
    {
        result_str = pstrdup(capability_str);
    }
    
    PG_RETURN_TEXT_P(cstring_to_text(result_str));
}

/* Parallel k-mer analysis functions declarations added elsewhere */

/* Function moved to kmersearch_gin.c */

/* Function moved to kmersearch_gin.c */



/* Function moved to kmersearch_gin.c */



/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */

/* Function moved to kmersearch_freq.c */


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
 * Expand DNA4 k-mer to DNA2 k-mers directly as uint values
 * Based on kmersearch_expand_dna4_kmer2_to_dna2_direct() but returns uint arrays
 */
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
