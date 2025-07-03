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
PG_FUNCTION_INFO_V1(kmersearch_extract_value);
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
PG_FUNCTION_INFO_V1(kmersearch_rawscore);
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
kmersearch_extract_value(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *dna = (kmersearch_dna2 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    char *sequence;
    int seq_len;
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    sequence = kmersearch_dna2_to_string((VarBit *)dna);
    seq_len = strlen(sequence);
    
    keys = kmersearch_extract_kmers(sequence, seq_len, k, nkeys);
    
    pfree(sequence);
    
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
    char *sequence;
    int seq_len;
    int k = kmersearch_kmer_size;  /* k-mer length from GUC variable */
    
    if (k < 4 || k > 64)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 64")));
    
    sequence = kmersearch_dna4_to_string((VarBit *)dna);
    seq_len = strlen(sequence);
    
    keys = kmersearch_extract_kmers_with_degenerate(sequence, seq_len, k, nkeys);
    
    pfree(sequence);
    
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
    
    char *dna_string = kmersearch_dna2_to_string(dna);
    char *pattern_string = text_to_cstring(pattern);
    
    /* Simple pattern matching - in reality this would use k-mer similarity scoring */
    bool result = (strstr(dna_string, pattern_string) != NULL);
    
    pfree(dna_string);
    pfree(pattern_string);
    
    PG_RETURN_BOOL(result);
}

/*
 * DNA4 =% operator for k-mer search
 */
Datum
kmersearch_dna4_match(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    text *pattern = PG_GETARG_TEXT_P(1);
    
    char *dna_string = kmersearch_dna4_to_string(dna);
    char *pattern_string = text_to_cstring(pattern);
    
    /* Simple pattern matching - in reality this would use k-mer similarity scoring */
    bool result = (strstr(dna_string, pattern_string) != NULL);
    
    pfree(dna_string);
    pfree(pattern_string);
    
    PG_RETURN_BOOL(result);
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
 * Raw score calculation function
 */
Datum
kmersearch_rawscore(PG_FUNCTION_ARGS)
{
    VarBit *sequence = PG_GETARG_VARBIT_P(0);
    text *query_text = PG_GETARG_TEXT_P(1);
    int score;
    
    /* Calculate raw score between sequence and query */
    score = kmersearch_calculate_raw_score(sequence, NULL, query_text);
    
    PG_RETURN_INT32(score);
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

