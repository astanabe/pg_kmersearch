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
double kmersearch_min_shared_kmer_rate = 0.5;  /* Default minimum shared k-mer rate for =% operator */
bool kmersearch_preclude_highfreq_kmer = false;  /* Default to not exclude high-frequency k-mers */

/* Cache configuration variables */
int kmersearch_query_kmer_cache_max_entries = 50000;  /* Default max query-kmer cache entries */
int kmersearch_actual_min_score_cache_max_entries = 50000;  /* Default max actual min score cache entries */
int kmersearch_highfreq_kmer_cache_load_batch_size = 10000;  /* Default batch size for loading high-frequency k-mers */
int kmersearch_highfreq_analysis_batch_size = 10000;  /* Default batch size for high-frequency k-mer analysis */
int kmersearch_highfreq_analysis_hashtable_size = 1000000;  /* Default hash table size for high-frequency k-mer analysis */

/* Global cache managers */
ActualMinScoreCacheManager *actual_min_score_cache_manager = NULL;

/* Global query-kmer cache manager for cross-query sharing */
QueryKmerCacheManager *query_kmer_cache_manager = NULL;


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

void _PG_init(void);

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
PG_FUNCTION_INFO_V1(kmersearch_dna2_match);
PG_FUNCTION_INFO_V1(kmersearch_dna4_match);

PG_FUNCTION_INFO_V1(kmersearch_matchscore_dna2);
PG_FUNCTION_INFO_V1(kmersearch_matchscore_dna4);
static simd_capability_t detect_cpu_capabilities(void);


static bool kmersearch_kmer_size_check_hook(int *newval, void **extra, GucSource source);
static bool kmersearch_occur_bitlen_check_hook(int *newval, void **extra, GucSource source);
static void kmersearch_kmer_size_assign_hook(int newval, void *extra);
static void kmersearch_occur_bitlen_assign_hook(int newval, void *extra);
static void kmersearch_max_appearance_rate_assign_hook(double newval, void *extra);
static void kmersearch_max_appearance_nrow_assign_hook(int newval, void *extra);
static void kmersearch_min_score_assign_hook(int newval, void *extra);
static void kmersearch_min_shared_kmer_rate_assign_hook(double newval, void *extra);
static void kmersearch_force_simd_capability_assign_hook(int newval, void *extra);

static void 
clear_highfreq_cache_with_warning(void)
{
    bool had_valid_cache = global_highfreq_cache.is_valid;
    
    kmersearch_highfreq_kmer_cache_free_internal();
    
    if (had_valid_cache)
    {
        elog(WARNING, "High-frequency k-mer cache has been cleared. "
                      "You may need to manually execute kmersearch_highfreq_kmer_cache_load() "
                      "to reload the cache if needed.");
    }
}

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
    
    
    /* Clear query-kmer cache */
    if (query_kmer_cache_manager)
        kmersearch_free_query_kmer_cache_manager(&query_kmer_cache_manager);
    
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        kmersearch_free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
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
        kmersearch_free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
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
        kmersearch_free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
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
        kmersearch_free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
}

/* Min shared k-mer rate change affects actual min score cache */
static void
kmersearch_min_shared_kmer_rate_assign_hook(double newval, void *extra)
{
    (void) newval;  /* Suppress unused parameter warning */
    (void) extra;   /* Suppress unused parameter warning */
    /* Clear actual min score cache */
    if (actual_min_score_cache_manager)
        kmersearch_free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
}

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

/* Occurrence bit length change affects high-freq caches */
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
static bool guc_variables_initialized = false;

/*
 * Check if GUC variables are properly initialized
 * If not, report an error with instructions for shared_preload_libraries
 */
void
kmersearch_check_guc_initialization(void)
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
    
    DefineCustomRealVariable("kmersearch.min_shared_kmer_rate",
                            "Minimum shared k-mer rate for =% operator matching",
                            "Minimum ratio of shared k-mers between query and target sequence (0.0-1.0)",
                            &kmersearch_min_shared_kmer_rate,
                            0.5,
                            0.0,
                            1.0,
                            PGC_USERSET,
                            0,
                            NULL,
                            kmersearch_min_shared_kmer_rate_assign_hook,
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

    DefineCustomIntVariable("kmersearch.query_kmer_cache_max_entries",
                           "Maximum number of entries in query-kmer cache",
                           "Controls the maximum number of cached query k-mer extraction results",
                           &kmersearch_query_kmer_cache_max_entries,
                           50000,
                           1000,
                           10000000,
                           PGC_USERSET,
                           0,
                           NULL,
                           kmersearch_query_kmer_cache_max_entries_assign_hook,
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

    DefineCustomIntVariable("kmersearch.highfreq_analysis_hashtable_size",
                           "Hash table size for high-frequency k-mer analysis",
                           "Initial size of the hash table used during analysis",
                           &kmersearch_highfreq_analysis_hashtable_size,
                           1000000,
                           10000,
                           100000000,
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
        else
        {
            /* Expand the pool */
            size_t new_size = pool->buffer_size * 2;
            while (pool->used + aligned_size > new_size)
                new_size *= 2;
        
            pool->buffer = repalloc(pool->buffer, new_size);
            pool->buffer_size = new_size;
            
            elog(DEBUG2, "Expanded memory pool from %zu to %zu bytes",
                 pool->buffer_size / 2, new_size);
        }
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
    query_uintkey = kmersearch_get_cached_query_uintkey(pattern_string, kmersearch_kmer_size, &query_nkeys);
    
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
        actual_min_score = kmersearch_get_cached_actual_min_score_uintkey(query_uintkey, query_nkeys, kmersearch_kmer_size);
        
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
    query_uintkey = kmersearch_get_cached_query_uintkey(pattern_string, kmersearch_kmer_size, &query_nkeys);
    
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
        actual_min_score = kmersearch_get_cached_actual_min_score_uintkey(query_uintkey, query_nkeys, kmersearch_kmer_size);
        
        /* Evaluate match condition */
        match = (shared_count >= actual_min_score);
        
        /* NOTE: Do NOT free query_uintkey - it is managed by the cache */
    }
    
    pfree(pattern_string);
    PG_RETURN_BOOL(match);
}

/*
 * Match score functions - calculate similarity scores
 * 
 * These functions calculate similarity scores by counting all shared k-mers between
 * sequence and query without any filtering or high-frequency k-mer exclusion.
 */
Datum
kmersearch_matchscore_dna2(PG_FUNCTION_ARGS)
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
    query_uintkey = kmersearch_get_cached_query_uintkey(query_string, kmersearch_kmer_size, &query_nkeys);
    
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
kmersearch_matchscore_dna4(PG_FUNCTION_ARGS)
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
    query_uintkey = kmersearch_get_cached_query_uintkey(query_string, k_size, &query_nkeys);
    
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

/*
 * Module cleanup function
 */
void
_PG_fini(void)
{
    /* Free query-kmer cache manager on module unload (uses TopMemoryContext - needs manual cleanup) */
    /* DNA2/DNA4 cache managers are now local and automatically freed with QueryContext */
    kmersearch_free_query_kmer_cache_manager(&query_kmer_cache_manager);
    
    /* Free actual min score cache manager on module unload */
    kmersearch_free_actual_min_score_cache_manager(&actual_min_score_cache_manager);
    
    /* Free high-frequency k-mer cache on module unload */
    kmersearch_highfreq_kmer_cache_free_internal();
}
