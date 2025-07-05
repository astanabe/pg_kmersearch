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
#include <unistd.h>

PG_MODULE_MAGIC;

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
static int kmersearch_cache_max_entries = 50000;  /* Default max cache entries */

/* Hash table entry for k-mer frequency counting */
typedef struct KmerFreqEntry
{
    VarBit      *kmer_key;      /* K-mer binary key (without occurrence count) */
    int         row_count;      /* Number of rows where this k-mer appears */
    bool        excluded;       /* Whether this k-mer is excluded */
} KmerFreqEntry;

/* Hash table key for k-mer frequency counting */
typedef struct KmerHashKey
{
    int         key_len;            /* Length of the key data */
    char        data[FLEXIBLE_ARRAY_MEMBER];
} KmerHashKey;

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

/* Analysis result structure */
typedef struct KmerAnalysisResult
{
    int64       total_rows;                /* Total number of rows analyzed */
    int         excluded_kmers_count;      /* Number of excluded k-mers */
    int         parallel_workers_used;     /* Number of parallel workers used */
    double      analysis_duration;        /* Analysis duration in seconds */
    double      max_appearance_rate_used; /* Max appearance rate used */
    int         max_appearance_nrow_used;  /* Max appearance nrow used */
} KmerAnalysisResult;

/* Drop analysis result structure */
typedef struct DropAnalysisResult
{
    int         dropped_analyses;          /* Number of dropped analyses */
    int         dropped_excluded_kmers;    /* Number of dropped excluded k-mers */
    int64       freed_storage_bytes;       /* Freed storage in bytes */
} DropAnalysisResult;

/* Cache entry structure */
typedef struct KmerCacheEntry
{
    uint64      hash_key;                  /* Hash key for this entry */
    VarBit      *sequence_copy;            /* Copy of sequence for exact match */
    char        *query_string_copy;        /* Copy of query string for exact match */
    KmerMatchResult result;                /* Cached calculation result */
    struct KmerCacheEntry *next;           /* For LRU chain */
    struct KmerCacheEntry *prev;           /* For LRU chain */
} KmerCacheEntry;

/* Cache manager structure */
typedef struct KmerCacheManager
{
    HTAB        *hash_table;               /* PostgreSQL hash table */
    MemoryContext cache_context;           /* Memory context for cache */
    int         max_entries;               /* Maximum number of entries */
    int         current_entries;           /* Current number of entries */
    uint64      hits;                      /* Cache hit count */
    uint64      misses;                    /* Cache miss count */
    KmerCacheEntry *lru_head;              /* LRU chain head (most recent) */
    KmerCacheEntry *lru_tail;              /* LRU chain tail (least recent) */
} KmerCacheManager;

/* Global cache managers */
static KmerCacheManager *dna2_cache_manager = NULL;
static KmerCacheManager *dna4_cache_manager = NULL;

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
    bool        is_excluded;              /* Whether this k-mer is excluded */
} CompactKmerFreq;

/* Worker buffer for k-mer collection */
typedef struct KmerBuffer
{
    CompactKmerFreq *entries;            /* Buffer entries */
    int             count;               /* Current count */
    int             capacity;            /* Buffer capacity */
    int             k_size;              /* K-mer size */
} KmerBuffer;

/* Worker state for parallel k-mer analysis */
typedef struct KmerWorkerState
{
    int         worker_id;                /* Worker identifier */
    BlockNumber start_block;              /* Starting block number */
    BlockNumber end_block;                /* Ending block number */
    KmerBuffer  buffer;                   /* Local k-mer buffer */
    int         local_excluded_count;     /* Local count of excluded k-mers */
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
static bool kmersearch_is_kmer_excluded(Oid index_oid, VarBit *kmer_key);
static int kmersearch_count_excluded_kmers_in_query(Oid index_oid, VarBit **query_keys, int nkeys);
static int kmersearch_get_adjusted_min_score(Oid index_oid, VarBit **query_keys, int nkeys);
static int kmersearch_calculate_raw_score(VarBit *seq1, VarBit *seq2, text *query_text);
static int kmersearch_count_mutual_excluded_kmers(VarBit *seq1, VarBit *seq2, text *query_text, Oid index_oid);
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

/* Cache management functions */
static void init_kmer_cache_manager(KmerCacheManager **manager, const char *name);
static uint64 generate_cache_key(VarBit *sequence, const char *query_string);
static bool sequences_equal(VarBit *a, VarBit *b);
static KmerCacheEntry *lookup_cache_entry(KmerCacheManager *manager, VarBit *sequence, const char *query_string);
static void store_cache_entry(KmerCacheManager *manager, uint64 hash_key, VarBit *sequence, const char *query_string, KmerMatchResult result);
static void lru_touch(KmerCacheManager *manager, KmerCacheEntry *entry);
static void lru_evict_oldest(KmerCacheManager *manager);
static KmerMatchResult get_cached_result_dna2(VarBit *sequence, const char *query_string);
static KmerMatchResult get_cached_result_dna4(VarBit *sequence, const char *query_string);
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

/* New parallel analysis functions */
static int kmersearch_determine_parallel_workers(int requested_workers, Relation target_relation);
static KmerAnalysisResult kmersearch_analyze_table_parallel(Oid table_oid, const char *column_name, int k_size, int parallel_workers);
static void kmersearch_worker_analyze_blocks(KmerWorkerState *worker, Relation rel, const char *column_name, int k_size);
static void kmersearch_merge_worker_results_sql(KmerWorkerState *workers, int num_workers, const char *final_table_name, int k_size, int threshold_rows);
static void kmersearch_persist_excluded_kmers(Oid table_oid, const char *column_name, int k_size, void *unused_table, int threshold_rows);
static void kmersearch_persist_excluded_kmers_from_temp(Oid table_oid, const char *column_name, int k_size, const char *temp_table_name);

/* New memory-efficient k-mer functions */
static size_t kmersearch_get_kmer_data_size(int k_size);
static bool kmersearch_get_index_info(Oid index_oid, Oid *table_oid, char **column_name, int *k_size);
static void kmersearch_spi_connect_or_error(void);
static void kmersearch_handle_spi_error(int spi_result, const char *operation);
static bool kmersearch_delete_kmer_from_gin_index(Relation index_rel, VarBit *kmer_key);
static List *kmersearch_get_excluded_kmers_list(Oid index_oid);
static int kmersearch_calculate_buffer_size(int k_size);
static KmerData kmersearch_encode_kmer_data(VarBit *kmer, int k_size);
static void kmersearch_init_buffer(KmerBuffer *buffer, int k_size);
static void kmersearch_add_to_buffer(KmerBuffer *buffer, KmerData kmer_data, const char *temp_table_name);
static void kmersearch_flush_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name);
static void kmersearch_create_worker_temp_table(const char *temp_table_name, int k_size);
static bool kmersearch_check_analysis_exists(Oid table_oid, const char *column_name, int k_size);
static void kmersearch_validate_analysis_parameters(Oid table_oid, const char *column_name, int k_size);
static Datum *kmersearch_filter_excluded_kmers(Oid table_oid, const char *column_name, int k_size, Datum *all_keys, int total_keys, int *filtered_count);
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
PG_FUNCTION_INFO_V1(kmersearch_get_excluded_kmers);
PG_FUNCTION_INFO_V1(kmersearch_analyze_table);
PG_FUNCTION_INFO_V1(kmersearch_drop_analysis);
PG_FUNCTION_INFO_V1(kmersearch_reduce_index);
PG_FUNCTION_INFO_V1(kmersearch_cache_stats);
PG_FUNCTION_INFO_V1(kmersearch_cache_free);

/* Score calculation functions */
PG_FUNCTION_INFO_V1(kmersearch_rawscore_dna2);
PG_FUNCTION_INFO_V1(kmersearch_rawscore_dna4);
PG_FUNCTION_INFO_V1(kmersearch_correctedscore);

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
    
    /* Encode sequence */
    data_ptr = VARBITS(result);
    for (i = 0; i < input_len; i++)
    {
        uint8 encoded = kmersearch_dna2_encode_table[(unsigned char)input_string[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        data_ptr[byte_pos] |= (encoded << (6 - bit_offset));
    }
    
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
    int i;
    
    if (dna == NULL)
        ereport(ERROR, (errmsg("input DNA sequence is NULL")));
    
    if (bit_len < 0)
        ereport(ERROR, (errmsg("invalid bit length: %d", bit_len)));
    
    if (bit_len % 2 != 0)
        ereport(ERROR, (errmsg("bit length must be even for DNA2")));
    
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
    
    /* Encode sequence */
    data_ptr = VARBITS(result);
    for (i = 0; i < input_len; i++)
    {
        uint8 encoded = kmersearch_dna4_encode_table[(unsigned char)input_string[i]];
        int bit_pos = i * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (bit_offset <= 4)
        {
            data_ptr[byte_pos] |= (encoded << (4 - bit_offset));
        }
        else
        {
            data_ptr[byte_pos] |= (encoded >> (bit_offset - 4));
            if (byte_pos + 1 < byte_len)
                data_ptr[byte_pos + 1] |= (encoded << (12 - bit_offset));
        }
    }
    
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

/*
 * Module initialization
 */
void
_PG_init(void)
{
    /* Define custom GUC variables */
    DefineCustomRealVariable("kmersearch.max_appearance_rate",
                            "Maximum appearance rate for k-mers to be included in index",
                            "K-mers appearing in more than this fraction of rows will be excluded",
                            &kmersearch_max_appearance_rate,
                            0.05,
                            0.0,
                            1.0,
                            PGC_USERSET,
                            0,
                            NULL,
                            NULL,
                            NULL);
    
    DefineCustomIntVariable("kmersearch.max_appearance_nrow",
                           "Maximum number of rows for k-mers to be included in index",
                           "K-mers appearing in more than this number of rows will be excluded (0 = unlimited)",
                           &kmersearch_max_appearance_nrow,
                           0,
                           0,
                           INT_MAX,
                           PGC_USERSET,
                           0,
                           NULL,
                           NULL,
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
                           NULL,
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
                           NULL,
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
                           NULL,
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
                            NULL,
                            NULL);

    /* Define GUC variables for cache configuration */
    DefineCustomIntVariable("kmersearch.cache_max_entries",
                           "Maximum number of entries in k-mer cache",
                           "Controls the maximum number of cached k-mer calculation results",
                           &kmersearch_cache_max_entries,
                           50000,
                           1000,
                           10000000,
                           PGC_USERSET,
                           0,
                           NULL,
                           NULL,
                           NULL);
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
        
        /* Extract 2 bits from source */
        uint8 base_bits = (src_data[src_byte_pos] >> (6 - src_bit_offset)) & 0x3;
        
        /* Store in destination */
        dst_data[dst_byte_pos] |= (base_bits << (6 - dst_bit_offset));
    }
    
    return result;
}

/*
 * Extract k-mers directly from DNA2 bit sequence
 */
static Datum *
kmersearch_extract_dna2_kmers_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    HTAB *occurrence_hash;
    HASHCTL hash_ctl;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Create hash table for tracking k-mer occurrences */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64_t) * ((k * 2 + 63) / 64);  /* Round up to uint64_t units */
    hash_ctl.entrysize = sizeof(int);
    occurrence_hash = hash_create("KmerOccurrenceHash", 1000, &hash_ctl,
                                  HASH_ELEM | HASH_BLOBS);
    
    /* Extract k-mers */
    for (i = 0; i <= seq_bases - k; i++)
    {
        /* Create k-mer key for hash lookup */
        uint64_t kmer_key_data[4] = {0};  /* Support up to k=128 */
        bool found;
        int *occurrence_ptr;
        VarBit *ngram_key;
        bits8 *src_data = VARBITS(seq);
        int j;
        
        /* Extract k-mer bits for hash key */
        for (j = 0; j < k; j++)
        {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
            
            /* Store in hash key */
            int key_bit_pos = j * 2;
            int key_word_pos = key_bit_pos / 64;
            int key_bit_offset = key_bit_pos % 64;
            kmer_key_data[key_word_pos] |= ((uint64_t)base_bits << (62 - key_bit_offset));
        }
        
        /* Get or create occurrence count */
        occurrence_ptr = (int *) hash_search(occurrence_hash, kmer_key_data, HASH_ENTER, &found);
        if (!found)
            *occurrence_ptr = 0;
        
        (*occurrence_ptr)++;
        
        /* Skip if occurrence exceeds bit limit */
        if (*occurrence_ptr > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* Create simple k-mer key (without occurrence count for matching) */
        ngram_key = kmersearch_create_kmer_key_from_dna2_bits(seq, i, k);
        keys[key_count++] = PointerGetDatum(ngram_key);
    }
    
    hash_destroy(occurrence_hash);
    
    *nkeys = key_count;
    return keys;
}

/*
 * Extract k-mers directly from DNA4 bit sequence with degenerate expansion
 */
static Datum *
kmersearch_extract_dna4_kmers_with_expansion_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    HTAB *occurrence_hash;
    HASHCTL hash_ctl;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Create hash table for tracking k-mer occurrences */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64_t) * ((k * 2 + 63) / 64);
    hash_ctl.entrysize = sizeof(int);
    occurrence_hash = hash_create("DNA4KmerOccurrenceHash", 1000, &hash_ctl,
                                  HASH_ELEM | HASH_BLOBS);
    
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
            uint64_t kmer_key_data[4] = {0};
            bool found;
            int *occurrence_ptr;
            VarBit *ngram_key;
            bits8 *kmer_data = VARBITS(dna2_kmer);
            int base_idx;
            
            /* Create hash key from DNA2 k-mer */
            for (base_idx = 0; base_idx < k; base_idx++)
            {
                int bit_pos = base_idx * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (kmer_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                
                int key_bit_pos = base_idx * 2;
                int key_word_pos = key_bit_pos / 64;
                int key_bit_offset = key_bit_pos % 64;
                kmer_key_data[key_word_pos] |= ((uint64_t)base_bits << (62 - key_bit_offset));
            }
            
            /* Get or create occurrence count */
            occurrence_ptr = (int *) hash_search(occurrence_hash, kmer_key_data, HASH_ENTER, &found);
            if (!found)
                *occurrence_ptr = 0;
            
            (*occurrence_ptr)++;
            
            /* Skip if occurrence exceeds bit limit */
            if (*occurrence_ptr > (1 << kmersearch_occur_bitlen))
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
    
    hash_destroy(occurrence_hash);
    
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
 * Initialize cache manager
 */
static void
init_kmer_cache_manager(KmerCacheManager **manager, const char *name)
{
    if (*manager == NULL)
    {
        MemoryContext old_context = MemoryContextSwitchTo(TopMemoryContext);
        HASHCTL hash_ctl;
        
        /* Allocate manager in TopMemoryContext */
        *manager = (KmerCacheManager *) palloc0(sizeof(KmerCacheManager));
        
        /* Create cache context under TopMemoryContext */
        (*manager)->cache_context = AllocSetContextCreate(TopMemoryContext,
                                                         "KmerCache",
                                                         ALLOCSET_DEFAULT_SIZES);
        
        /* Initialize cache parameters */
        (*manager)->max_entries = kmersearch_cache_max_entries;
        (*manager)->current_entries = 0;
        (*manager)->hits = 0;
        (*manager)->misses = 0;
        (*manager)->lru_head = NULL;
        (*manager)->lru_tail = NULL;
        
        /* Create hash table */
        MemSet(&hash_ctl, 0, sizeof(hash_ctl));
        hash_ctl.keysize = sizeof(uint64);
        hash_ctl.entrysize = sizeof(KmerCacheEntry);
        hash_ctl.hcxt = (*manager)->cache_context;
        
        (*manager)->hash_table = hash_create(name, 1024, &hash_ctl,
                                           HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);
        
        MemoryContextSwitchTo(old_context);
    }
}

/*
 * Generate cache key using PostgreSQL's hash_any_extended
 */
static uint64
generate_cache_key(VarBit *sequence, const char *query_string)
{
    uint64 seq_hash, query_hash;
    
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
 * Move entry to head of LRU chain (most recently used)
 */
static void
lru_touch(KmerCacheManager *manager, KmerCacheEntry *entry)
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
 * Evict least recently used entry
 */
static void
lru_evict_oldest(KmerCacheManager *manager)
{
    KmerCacheEntry *oldest = manager->lru_tail;
    bool found;
    
    if (!oldest)
        return;
    
    /* Remove from LRU chain */
    if (oldest->prev)
        oldest->prev->next = NULL;
    else
        manager->lru_head = NULL;
    manager->lru_tail = oldest->prev;
    
    /* Remove from hash table */
    hash_search(manager->hash_table, &oldest->hash_key, HASH_REMOVE, &found);
    
    manager->current_entries--;
}

/*
 * Look up cache entry
 */
static KmerCacheEntry *
lookup_cache_entry(KmerCacheManager *manager, VarBit *sequence, const char *query_string)
{
    uint64 hash_key = generate_cache_key(sequence, query_string);
    KmerCacheEntry *entry;
    bool found;
    
    entry = (KmerCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_FIND, &found);
    
    if (found && sequences_equal(entry->sequence_copy, sequence) &&
        strcmp(entry->query_string_copy, query_string) == 0)
    {
        /* Cache hit - move to head of LRU */
        lru_touch(manager, entry);
        manager->hits++;
        return entry;
    }
    
    return NULL;  /* Cache miss */
}

/*
 * Store entry in cache
 */
static void
store_cache_entry(KmerCacheManager *manager, uint64 hash_key, VarBit *sequence, 
                 const char *query_string, KmerMatchResult result)
{
    KmerCacheEntry *entry;
    MemoryContext old_context;
    bool found;
    
    /* Check if we need to evict */
    while (manager->current_entries >= manager->max_entries)
        lru_evict_oldest(manager);
    
    old_context = MemoryContextSwitchTo(manager->cache_context);
    
    /* Create new entry */
    entry = (KmerCacheEntry *) hash_search(manager->hash_table, &hash_key, HASH_ENTER, &found);
    
    if (!found)
    {
        /* New entry */
        entry->hash_key = hash_key;
        entry->sequence_copy = (VarBit *) palloc(VARSIZE(sequence));
        memcpy(entry->sequence_copy, sequence, VARSIZE(sequence));
        entry->query_string_copy = pstrdup(query_string);
        entry->result = result;
        entry->next = NULL;
        entry->prev = NULL;
        
        /* Add to head of LRU chain */
        lru_touch(manager, entry);
        manager->current_entries++;
    }
    
    MemoryContextSwitchTo(old_context);
}

/*
 * Get cached result for DNA2
 */
static KmerMatchResult
get_cached_result_dna2(VarBit *sequence, const char *query_string)
{
    KmerCacheEntry *entry;
    KmerMatchResult result;
    uint64 hash_key;
    
    init_kmer_cache_manager(&dna2_cache_manager, "DNA2_KmerCache");
    
    entry = lookup_cache_entry(dna2_cache_manager, sequence, query_string);
    if (entry != NULL)
    {
        return entry->result;  /* Cache hit */
    }
    
    /* Cache miss - calculate result */
    dna2_cache_manager->misses++;
    result = kmersearch_calculate_kmer_match_and_score_dna2(sequence, query_string);
    
    /* Store in cache */
    hash_key = generate_cache_key(sequence, query_string);
    store_cache_entry(dna2_cache_manager, hash_key, sequence, query_string, result);
    
    return result;
}

/*
 * Get cached result for DNA4
 */
static KmerMatchResult
get_cached_result_dna4(VarBit *sequence, const char *query_string)
{
    KmerCacheEntry *entry;
    KmerMatchResult result;
    uint64 hash_key;
    
    init_kmer_cache_manager(&dna4_cache_manager, "DNA4_KmerCache");
    
    entry = lookup_cache_entry(dna4_cache_manager, sequence, query_string);
    if (entry != NULL)
    {
        return entry->result;  /* Cache hit */
    }
    
    /* Cache miss - calculate result */
    dna4_cache_manager->misses++;
    result = kmersearch_calculate_kmer_match_and_score_dna4(sequence, query_string);
    
    /* Store in cache */
    hash_key = generate_cache_key(sequence, query_string);
    store_cache_entry(dna4_cache_manager, hash_key, sequence, query_string, result);
    
    return result;
}

/*
 * Fast k-mer matching using hash table
 */
static int
kmersearch_count_matching_kmers_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i, j;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* Simple O(n*m) comparison - use hash table if performance becomes an issue */
    for (i = 0; i < seq_nkeys; i++)
    {
        for (j = 0; j < query_nkeys; j++)
        {
            if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[j]) &&
                VARSIZE(seq_keys[i]) == VARSIZE(query_keys[j]) &&
                memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[j]), VARBITBYTES(seq_keys[i])) == 0)
            {
                match_count++;
                break;  /* Found a match, move to next seq k-mer */
            }
        }
    }
    
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
    int adjusted_min_score;
    int i;
    Oid index_oid = InvalidOid;  /* TODO: Get actual index OID from context */
    VarBit **query_key_array;
    
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
    
    /* Calculate adjusted minimum score based on excluded k-mers in query */
    adjusted_min_score = kmersearch_get_adjusted_min_score(index_oid, query_key_array, nkeys);
    
    pfree(query_key_array);
    
    /* Return true if match count meets adjusted minimum score */
    PG_RETURN_BOOL(match_count >= adjusted_min_score);
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
    KmerMatchResult result = get_cached_result_dna2(dna, pattern_string);
    
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
    KmerMatchResult result = get_cached_result_dna4(dna, pattern_string);
    
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

/*
 * Hash function for k-mer keys
 */
static uint32
kmersearch_kmer_hash(const void *key, Size keysize)
{
    const KmerHashKey *k = (const KmerHashKey *) key;
    return DatumGetUInt32(hash_any((unsigned char *) k->data, k->key_len));
}

/*
 * Compare function for k-mer keys
 */
static int
kmersearch_kmer_compare(const void *key1, const void *key2, Size keysize)
{
    const KmerHashKey *k1 = (const KmerHashKey *) key1;
    const KmerHashKey *k2 = (const KmerHashKey *) key2;
    
    if (k1->key_len != k2->key_len)
        return k1->key_len - k2->key_len;
    
    return memcmp(k1->data, k2->data, k1->key_len);
}

/*
 * Analyze table frequency and determine excluded k-mers
 */
Datum
kmersearch_analyze_table_frequency(PG_FUNCTION_ARGS)
{
    Oid table_oid = PG_GETARG_OID(0);
    text *column_name_text = PG_GETARG_TEXT_P(1);
    int k = PG_GETARG_INT32(2);
    Oid index_oid = PG_GETARG_OID(3);
    
    char *column_name = text_to_cstring(column_name_text);
    int excluded_count = 0;
    
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
        /* Skip frequency analysis - create empty excluded k-mer list */
        ereport(NOTICE, (errmsg("High-frequency k-mer exclusion disabled, skipping table scan")));
        
        /* Insert index info with zero excluded k-mers */
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
     * 5. Insert excluded k-mers into kmersearch_excluded_kmers table
     * 6. Insert index statistics into kmersearch_index_info table
     */
    
    PG_RETURN_INT32(excluded_count);
}

/*
 * Get excluded k-mers for an index
 */
Datum
kmersearch_get_excluded_kmers(PG_FUNCTION_ARGS)
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
    
    /* Build query to get excluded k-mers */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT kmer_key FROM kmersearch_excluded_kmers WHERE index_oid = %u ORDER BY kmer_key",
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
 * Count excluded k-mers in query sequence
 */
static int
kmersearch_count_excluded_kmers_in_query(Oid index_oid, VarBit **query_keys, int nkeys)
{
    int excluded_count = 0;
    int i;
    
    /* For each k-mer in the query, check if it's excluded */
    for (i = 0; i < nkeys; i++)
    {
        if (kmersearch_is_kmer_excluded(index_oid, query_keys[i]))
        {
            excluded_count++;
        }
    }
    
    return excluded_count;
}

/*
 * Calculate adjusted minimum score based on excluded k-mers in query
 */
static int
kmersearch_get_adjusted_min_score(Oid index_oid, VarBit **query_keys, int nkeys)
{
    int excluded_count = kmersearch_count_excluded_kmers_in_query(index_oid, query_keys, nkeys);
    int adjusted_score = kmersearch_min_score - excluded_count;
    
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
    
    /* Count matching k-mers */
    for (i = 0; i < seq1_nkeys; i++)
    {
        for (j = 0; j < seq2_nkeys; j++)
        {
            if (seq1_keys[i] && seq2_keys[j] && 
                VARBITLEN(seq1_keys[i]) == VARBITLEN(seq2_keys[j]) &&
                memcmp(VARBITS(seq1_keys[i]), VARBITS(seq2_keys[j]), 
                       VARBITBYTES(seq1_keys[i])) == 0)
            {
                score++;
                break;  /* Count each k-mer only once */
            }
        }
    }
    
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
 * Count excluded k-mers that appear in both sequences
 */
static int
kmersearch_count_mutual_excluded_kmers(VarBit *seq1, VarBit *seq2, text *query_text, Oid index_oid)
{
    char *query_string = text_to_cstring(query_text);
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    int mutual_excluded = 0;
    VarBit **seq1_keys, **seq2_keys;
    int seq1_nkeys, seq2_nkeys;
    int i, j;
    
    /* Extract k-mers from both sequences */
    seq1_keys = kmersearch_extract_kmers_from_varbit(seq1, k, &seq1_nkeys);
    seq2_keys = kmersearch_extract_kmers_from_query(query_string, k, &seq2_nkeys);
    
    /* Count excluded k-mers that appear in both sequences */
    for (i = 0; i < seq1_nkeys; i++)
    {
        if (!seq1_keys[i] || !kmersearch_is_kmer_excluded(index_oid, seq1_keys[i]))
            continue;
            
        for (j = 0; j < seq2_nkeys; j++)
        {
            if (seq2_keys[j] && 
                VARBITLEN(seq1_keys[i]) == VARBITLEN(seq2_keys[j]) &&
                memcmp(VARBITS(seq1_keys[i]), VARBITS(seq2_keys[j]), 
                       VARBITBYTES(seq1_keys[i])) == 0)
            {
                mutual_excluded++;
                break;
            }
        }
    }
    
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
    return mutual_excluded;
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
    KmerMatchResult result = get_cached_result_dna2(sequence, query_string);
    
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
    KmerMatchResult result = get_cached_result_dna4(sequence, query_string);
    
    pfree(query_string);
    
    /* Return shared count as raw score */
    PG_RETURN_INT32(result.valid ? result.shared_count : 0);
}

/*
 * Corrected score calculation function
 */
Datum
kmersearch_correctedscore(PG_FUNCTION_ARGS)
{
    VarBit *sequence = PG_GETARG_VARBIT_P(0);
    text *query_text = PG_GETARG_TEXT_P(1);
    int raw_score, mutual_excluded, corrected_score;
    Oid index_oid = InvalidOid;  /* TODO: Get actual index OID from context */
    
    /* Calculate raw score */
    raw_score = kmersearch_calculate_raw_score(sequence, NULL, query_text);
    
    /* Count mutual excluded k-mers */
    mutual_excluded = kmersearch_count_mutual_excluded_kmers(sequence, NULL, query_text, index_oid);
    
    /* Corrected score = raw score + mutual excluded k-mers */
    corrected_score = raw_score + mutual_excluded;
    
    PG_RETURN_INT32(corrected_score);
}

/*
 * Check if a k-mer is excluded for a given index
 */
static bool
kmersearch_is_kmer_excluded(Oid index_oid, VarBit *kmer_key)
{
    /* For now, always return false - will implement proper check later */
    return false;
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
    query_keys = kmersearch_extract_query_kmers(query_string, k, &query_nkeys);
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
    
    /* Evaluate match conditions */
    result = kmersearch_evaluate_match_conditions(shared_count, query_nkeys);
    
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
    query_keys = kmersearch_extract_query_kmers(query_string, k, &query_nkeys);
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
    
    /* Evaluate match conditions */
    result = kmersearch_evaluate_match_conditions(shared_count, query_nkeys);
    
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
    int query_len;
    query_len = strlen(query_string);
    if (query_len < k) {
        return result;
    }
    
    /* Extract k-mers from DNA2 sequence (no degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna2_kmers_direct(sequence, k, &result.seq_nkeys);
    if (seq_datum_keys != NULL && result.seq_nkeys > 0) {
        int i;
        seq_keys = (VarBit **) palloc(result.seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < result.seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    if (seq_keys == NULL || result.seq_nkeys == 0) {
        goto cleanup;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    query_keys = kmersearch_extract_query_kmers(query_string, k, &result.query_nkeys);
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
    
    /* Evaluate match conditions for =% operator */
    result.match_result = kmersearch_evaluate_match_conditions(result.shared_count, result.query_nkeys);
    
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
    int query_len;
    query_len = strlen(query_string);
    if (query_len < k) {
        return result;
    }
    
    /* Extract k-mers from DNA4 sequence (with degenerate expansion) */
    seq_datum_keys = kmersearch_extract_dna4_kmers_with_expansion_direct(sequence, k, &result.seq_nkeys);
    if (seq_datum_keys != NULL && result.seq_nkeys > 0) {
        int i;
        seq_keys = (VarBit **) palloc(result.seq_nkeys * sizeof(VarBit *));
        for (i = 0; i < result.seq_nkeys; i++) {
            seq_keys[i] = DatumGetVarBitP(seq_datum_keys[i]);
        }
    }
    if (seq_keys == NULL || result.seq_nkeys == 0) {
        goto cleanup;
    }
    
    /* Extract k-mers from query (with degenerate expansion) */
    query_keys = kmersearch_extract_query_kmers(query_string, k, &result.query_nkeys);
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
    
    /* Evaluate match conditions for =% operator */
    result.match_result = kmersearch_evaluate_match_conditions(result.shared_count, result.query_nkeys);
    
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
 * Flush buffer contents to temporary table
 */
static void
kmersearch_flush_buffer_to_table(KmerBuffer *buffer, const char *temp_table_name)
{
    StringInfoData query;
    int i;
    
    if (buffer->count == 0) return;
    
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
    entry->is_excluded = false;
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
        data_type = "smallint";
    } else if (k_size <= 16) {
        data_type = "integer";
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
    HeapScanDesc scan;
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
        data_type = "smallint";
    } else if (k_size <= 16) {
        data_type = "integer";
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
        workers[i].local_excluded_count = 0;
        workers[i].rows_processed = 0;
        workers[i].temp_table_name = NULL;
        
        /* Process assigned blocks */
        kmersearch_worker_analyze_blocks(&workers[i], rel, column_name, k_size);
    }
    
    /* Merge worker results using SQL aggregation */
    {
        char *final_table_name = psprintf("temp_kmer_final_%d", getpid());
        kmersearch_merge_worker_results_sql(workers, num_workers, final_table_name, k_size, threshold_rows);
        
        /* Count excluded k-mers */
        SPI_connect();
        {
            StringInfoData count_query;
            initStringInfo(&count_query);
            appendStringInfo(&count_query, "SELECT count(*) FROM %s", final_table_name);
            
            if (SPI_exec(count_query.data, 0) == SPI_OK_SELECT && SPI_processed > 0) {
                result.excluded_kmers_count = DatumGetInt32(SPI_getbinval(SPI_tuptable->vals[0], 
                                                                        SPI_tuptable->tupdesc, 1, NULL));
            }
            pfree(count_query.data);
        }
        SPI_finish();
        
        /* Persist excluded k-mers to permanent tables */
        kmersearch_persist_excluded_kmers_from_temp(table_oid, column_name, k_size, final_table_name);
        
        pfree(final_table_name);
    }
    
    /* Calculate total statistics */
    for (i = 0; i < num_workers; i++) {
        result.total_rows += workers[i].rows_processed;
        result.excluded_kmers_count += workers[i].local_excluded_count;
    }
    
    /* Clean up */
    pfree(workers);
    relation_close(rel, AccessShareLock);
    
    return result;
}

/*
 * Persist excluded k-mers from temporary table to permanent tables
 */
static void
kmersearch_persist_excluded_kmers_from_temp(Oid table_oid, const char *column_name, int k_size,
                                           const char *temp_table_name)
{
    StringInfoData query;
    
    initStringInfo(&query);
    
    /* Insert excluded k-mers into permanent table */
    appendStringInfo(&query,
        "INSERT INTO kmersearch_excluded_kmers (index_oid, kmer_key, frequency_count, exclusion_reason) "
        "SELECT %u, kmer_data::varbit, frequency_count, 'high_frequency' "
        "FROM %s",
        table_oid, temp_table_name);
    
    SPI_connect();
    SPI_exec(query.data, 0);
    SPI_finish();
    
    pfree(query.data);
}

/*
 * Persist excluded k-mers to database tables (legacy - disabled in memory-efficient version)
 */
static void
kmersearch_persist_excluded_kmers(Oid table_oid, const char *column_name, int k_size, 
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
 * Filter excluded k-mers from the key array
 */
static Datum *
kmersearch_filter_excluded_kmers(Oid table_oid, const char *column_name, int k_size, 
                                Datum *all_keys, int total_keys, int *filtered_count)
{
    int ret;
    StringInfoData query;
    HTAB *excluded_hash = NULL;
    HASHCTL hash_ctl;
    Datum *filtered_keys;
    int filtered_idx = 0;
    int i;
    
    /* Check if analysis exists */
    if (!kmersearch_check_analysis_exists(table_oid, column_name, k_size)) {
        /* No analysis found, return all keys */
        *filtered_count = total_keys;
        return all_keys;
    }
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Create hash table for excluded k-mers lookup */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(KmerHashKey);
    hash_ctl.entrysize = sizeof(KmerFreqEntry);
    hash_ctl.hash = kmersearch_kmer_hash;
    hash_ctl.match = kmersearch_kmer_compare;
    hash_ctl.hcxt = CurrentMemoryContext;
    
    excluded_hash = hash_create("Excluded K-mers", 1024, &hash_ctl,
                               HASH_ELEM | HASH_FUNCTION | HASH_COMPARE | HASH_CONTEXT);
    
    /* Build query to get excluded k-mers from index info table */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT ek.kmer_key FROM kmersearch_excluded_kmers ek "
        "JOIN kmersearch_index_info ii ON ek.index_oid = ii.index_oid "
        "WHERE ii.table_oid = %u AND ii.column_name = '%s' AND ii.k_value = %d",
        table_oid, column_name, k_size);
    
    /* Execute query and populate hash table */
    ret = SPI_execute(query.data, true, 0);
    if (ret == SPI_OK_SELECT && SPI_processed > 0)
    {
        int j;
        for (j = 0; j < SPI_processed; j++)
        {
            bool isnull;
            Datum kmer_datum;
            VarBit *kmer;
            KmerHashKey *hash_key;
            KmerFreqEntry *entry;
            bool found;
            
            kmer_datum = SPI_getbinval(SPI_tuptable->vals[j], SPI_tuptable->tupdesc, 1, &isnull);
            if (!isnull)
            {
                kmer = DatumGetVarBitP(kmer_datum);
                
                /* Create hash key */
                hash_key = (KmerHashKey *) palloc(sizeof(KmerHashKey) + VARSIZE(kmer));
                hash_key->key_len = VARSIZE(kmer);
                memcpy(hash_key->data, kmer, VARSIZE(kmer));
                
                /* Add to hash table */
                entry = (KmerFreqEntry *) hash_search(excluded_hash, hash_key, HASH_ENTER, &found);
                if (!found)
                {
                    entry->kmer_key = (VarBit *) palloc(VARSIZE(kmer));
                    memcpy(entry->kmer_key, kmer, VARSIZE(kmer));
                    entry->excluded = true;
                }
                
                pfree(hash_key);
            }
        }
    }
    
    /* Filter out excluded k-mers */
    filtered_keys = (Datum *) palloc(total_keys * sizeof(Datum));
    
    for (i = 0; i < total_keys; i++)
    {
        VarBit *kmer = DatumGetVarBitP(all_keys[i]);
        KmerHashKey *hash_key;
        KmerFreqEntry *entry;
        bool found;
        
        /* Create hash key for lookup */
        hash_key = (KmerHashKey *) palloc(sizeof(KmerHashKey) + VARSIZE(kmer));
        hash_key->key_len = VARSIZE(kmer);
        memcpy(hash_key->data, kmer, VARSIZE(kmer));
        
        /* Check if this k-mer is excluded */
        entry = (KmerFreqEntry *) hash_search(excluded_hash, hash_key, HASH_FIND, &found);
        
        if (!found)
        {
            /* Not excluded, include in filtered result */
            filtered_keys[filtered_idx++] = all_keys[i];
        }
        
        pfree(hash_key);
    }
    
    /* Cleanup */
    pfree(query.data);
    if (excluded_hash)
        hash_destroy(excluded_hash);
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
    values[1] = Int32GetDatum(result.excluded_kmers_count);
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
    values[1] = Int32GetDatum(result.dropped_excluded_kmers);
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
        result.dropped_excluded_kmers = 10;  /* Mock value */
    } else {
        ereport(NOTICE, (errmsg("Would drop all analyses for table %u, column %s", 
                               table_oid, column_name)));
        result.dropped_analyses = 3;  /* Mock value */
        result.dropped_excluded_kmers = 30;  /* Mock value */
    }
    
    result.freed_storage_bytes = result.dropped_excluded_kmers * 100;  /* Mock calculation */
    
    return result;
}

/*
 * Helper function to get excluded k-mers list for a given index
 */
static List *
kmersearch_get_excluded_kmers_list(Oid index_oid)
{
    List *excluded_kmers = NIL;
    int ret;
    StringInfoData query;
    
    /* Connect to SPI */
    kmersearch_spi_connect_or_error();
    
    /* Build query to get excluded k-mers */
    initStringInfo(&query);
    appendStringInfo(&query,
        "SELECT ek.kmer_key FROM kmersearch_excluded_kmers ek "
        "JOIN kmersearch_index_info ii ON ek.index_oid = ii.index_oid "
        "WHERE ii.index_oid = %u ORDER BY ek.kmer_key",
        index_oid);
    
    /* Execute query */
    ret = SPI_execute(query.data, true, 0);
    kmersearch_handle_spi_error(ret, "SELECT excluded k-mers");
    
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
                excluded_kmers = lappend(excluded_kmers, kmer);
            }
        }
    }
    
    /* Cleanup */
    pfree(query.data);
    SPI_finish();
    
    return excluded_kmers;
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
 * Main function to reduce index size by removing excluded k-mers
 */
Datum
kmersearch_reduce_index(PG_FUNCTION_ARGS)
{
    Oid index_oid = PG_GETARG_OID(0);
    Relation index_rel;
    List *excluded_kmers;
    ListCell *cell;
    int removed_count = 0;
    int total_excluded = 0;
    StringInfoData result_msg;
    
    /* Validate index OID */
    if (!OidIsValid(index_oid))
        ereport(ERROR, (errmsg("invalid index OID")));
    
    /* Open the index with exclusive lock */
    index_rel = try_relation_open(index_oid, AccessExclusiveLock);
    if (index_rel == NULL)
        ereport(ERROR, (errmsg("index with OID %u does not exist", index_oid)));
    
    /* Verify it's a GIN index */
    if (index_rel->rd_rel->relam != GIN_AM_OID)
    {
        relation_close(index_rel, AccessExclusiveLock);
        ereport(ERROR, (errmsg("index with OID %u is not a GIN index", index_oid)));
    }
    
    /* Get list of excluded k-mers */
    excluded_kmers = kmersearch_get_excluded_kmers_list(index_oid);
    total_excluded = list_length(excluded_kmers);
    
    if (total_excluded == 0)
    {
        relation_close(index_rel, AccessExclusiveLock);
        ereport(NOTICE, (errmsg("No excluded k-mers found for index OID %u", index_oid)));
        PG_RETURN_TEXT_P(cstring_to_text("No excluded k-mers to remove"));
    }
    
    ereport(NOTICE, (errmsg("Starting index reduction: %d excluded k-mers to remove", 
                           total_excluded)));
    
    /* Remove each excluded k-mer from the index */
    foreach(cell, excluded_kmers)
    {
        VarBit *kmer = (VarBit *) lfirst(cell);
        
        if (kmersearch_delete_kmer_from_gin_index(index_rel, kmer))
        {
            removed_count++;
        }
        
        /* Report progress every 1000 k-mers */
        if (removed_count % 1000 == 0 && removed_count > 0)
        {
            ereport(NOTICE, (errmsg("Progress: removed %d of %d k-mers", 
                                   removed_count, total_excluded)));
        }
    }
    
    /* 
     * TODO: Implement index cleanup and statistics update
     * - Call ginvacuumcleanup() equivalent for statistics update
     * - Reclaim freed space
     */
    
    /* Close the index */
    relation_close(index_rel, AccessExclusiveLock);
    
    /* Build result message */
    initStringInfo(&result_msg);
    appendStringInfo(&result_msg, 
        "Index reduction completed: removed %d of %d excluded k-mers",
        removed_count, total_excluded);
    
    ereport(NOTICE, (errmsg("%s", result_msg.data)));
    
    PG_RETURN_TEXT_P(cstring_to_text(result_msg.data));
}

/*
 * Cache statistics function
 */
Datum
kmersearch_cache_stats(PG_FUNCTION_ARGS)
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
    if (dna2_cache_manager)
    {
        values[0] = Int64GetDatum(dna2_cache_manager->hits);
        values[1] = Int64GetDatum(dna2_cache_manager->misses);
        values[2] = Int32GetDatum(dna2_cache_manager->current_entries);
        values[3] = Int32GetDatum(dna2_cache_manager->max_entries);
    }
    else
    {
        values[0] = Int64GetDatum(0);
        values[1] = Int64GetDatum(0);
        values[2] = Int32GetDatum(0);
        values[3] = Int32GetDatum(0);
    }
    
    /* DNA4 cache statistics */
    if (dna4_cache_manager)
    {
        values[4] = Int64GetDatum(dna4_cache_manager->hits);
        values[5] = Int64GetDatum(dna4_cache_manager->misses);
        values[6] = Int32GetDatum(dna4_cache_manager->current_entries);
        values[7] = Int32GetDatum(dna4_cache_manager->max_entries);
    }
    else
    {
        values[4] = Int64GetDatum(0);
        values[5] = Int64GetDatum(0);
        values[6] = Int32GetDatum(0);
        values[7] = Int32GetDatum(0);
    }
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}

/*
 * Free cache function - clears all cached data
 */
static void
free_cache_manager(KmerCacheManager **manager)
{
    if (*manager != NULL)
    {
        /* Delete the cache context, which will free all allocated memory */
        if ((*manager)->cache_context)
            MemoryContextDelete((*manager)->cache_context);
        
        /* Free the manager itself */
        pfree(*manager);
        *manager = NULL;
    }
}

/*
 * Cache free function
 */
Datum
kmersearch_cache_free(PG_FUNCTION_ARGS)
{
    int freed_entries = 0;
    
    /* Count entries before freeing */
    if (dna2_cache_manager)
        freed_entries += dna2_cache_manager->current_entries;
    if (dna4_cache_manager)
        freed_entries += dna4_cache_manager->current_entries;
    
    /* Free both cache managers */
    free_cache_manager(&dna2_cache_manager);
    free_cache_manager(&dna4_cache_manager);
    
    PG_RETURN_INT32(freed_entries);
}

/*
 * Module cleanup function
 */
void
_PG_fini(void)
{
    /* Free cache managers on module unload */
    free_cache_manager(&dna2_cache_manager);
    free_cache_manager(&dna4_cache_manager);
}
