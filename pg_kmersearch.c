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
#include <ctype.h>
#include <string.h>
#include <math.h>

PG_MODULE_MAGIC;

/* Global variables for k-mer search configuration */
static int kmersearch_occur_bitlen = 8;  /* Default 8 bits for occurrence count */
static int kmersearch_kmer_size = 8;  /* Default k-mer size */
static double kmersearch_max_appearance_rate = 0.05;  /* Default max appearance rate */
static int kmersearch_max_appearance_nrow = 0;  /* Default max appearance nrow (0 = undefined) */
static int kmersearch_min_score = 1;  /* Default minimum score for GIN search */
static double kmersearch_min_shared_ngram_key_rate = 0.9;  /* Default minimum shared n-gram key rate for =% operator */

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
static Datum *kmersearch_extract_dna2_kmers_direct(VarBit *seq, int k, int *nkeys);
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
        
        /* Create n-gram key */
        ngram_key = kmersearch_create_ngram_key_from_dna2_bits(seq, i, k, *occurrence_ptr);
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
            {
                pfree(dna2_kmer);
                continue;
            }
            
            /* Create n-gram key with occurrence count */
            ngram_key = kmersearch_create_ngram_key_with_occurrence_from_dna2(dna2_kmer, k, *occurrence_ptr);
            keys[key_count++] = PointerGetDatum(ngram_key);
            
            pfree(dna2_kmer);
        }
        
        pfree(expanded_kmers);
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
    int kmer_bytes = (kmer_bits + 7) / 8;
    memcpy(dst_data, src_data, kmer_bytes);
    
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
 * Fast k-mer matching using hash table
 */
static int
kmersearch_count_matching_kmers_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    HTAB *query_hash;
    HASHCTL hash_ctl;
    int match_count = 0;
    int i;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* Create hash table for query k-mers */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(VarBit *);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    query_hash = hash_create("QueryKmerHash", query_nkeys, &hash_ctl,
                             HASH_ELEM | HASH_FUNCTION);
    
    /* Insert all query k-mers into hash table */
    for (i = 0; i < query_nkeys; i++)
    {
        bool found;
        hash_search(query_hash, &query_keys[i], HASH_ENTER, &found);
    }
    
    /* Count matches in sequence k-mers */
    for (i = 0; i < seq_nkeys; i++)
    {
        bool found;
        hash_search(query_hash, &seq_keys[i], HASH_FIND, &found);
        if (found)
            match_count++;
    }
    
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
    
    /* Use common calculation function */
    KmerMatchResult result = kmersearch_calculate_kmer_match_and_score_dna2(dna, pattern_string);
    
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
    
    /* Use common calculation function */
    KmerMatchResult result = kmersearch_calculate_kmer_match_and_score_dna4(dna, pattern_string);
    
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
    
    /* Implementation to retrieve excluded k-mers from system table */
    /* This would return an array of excluded k-mer keys */
    
    PG_RETURN_NULL(); /* Placeholder */
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
    
    /* Use common calculation function */
    KmerMatchResult result = kmersearch_calculate_kmer_match_and_score_dna2(sequence, query_string);
    
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
    
    /* Use common calculation function */
    KmerMatchResult result = kmersearch_calculate_kmer_match_and_score_dna4(sequence, query_string);
    
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
    seq_keys = (VarBit **)kmersearch_extract_dna2_kmers_direct(sequence, k, &result.seq_nkeys);
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
    seq_keys = (VarBit **)kmersearch_extract_dna4_kmers_with_expansion_direct(sequence, k, &result.seq_nkeys);
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
    
    return result;
}
