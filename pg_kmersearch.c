#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/builtins.h"
#include "utils/varbit.h"
#include "utils/memutils.h"
#include "utils/varlena.h"
#include "varatt.h"
#include <ctype.h>
#include <string.h>

PG_MODULE_MAGIC;

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
    kmersearch_dna2 *result;
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
    result = (kmersearch_dna2 *) palloc0(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    
    /* Set bit length */
    *((int32 *) result->data) = bit_len;
    
    /* Encode sequence */
    data_ptr = result->data + sizeof(int32);
    for (i = 0; i < input_len; i++)
    {
        uint8 encoded = kmersearch_dna2_encode_table[(unsigned char)input_string[i]];
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        data_ptr[byte_pos] |= (encoded << (6 - bit_offset));
    }
    
    PG_RETURN_POINTER(result);
}

/*
 * DNA2 output function
 */
Datum
kmersearch_dna2_out(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *dna = (kmersearch_dna2 *) PG_GETARG_POINTER(0);
    int32 bit_len = *((int32 *) dna->data);
    int char_len = bit_len / 2;
    char *result;
    bits8 *data_ptr = dna->data + sizeof(int32);
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
    PG_RETURN_CSTRING(result);
}

/*
 * DNA2 receive function (binary input)
 */
Datum
kmersearch_dna2_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    kmersearch_dna2 *result;
    int32 bit_len;
    int byte_len;
    
    bit_len = pq_getmsgint(buf, sizeof(int32));
    byte_len = (bit_len + 7) / 8;
    
    result = (kmersearch_dna2 *) palloc(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    
    *((int32 *) result->data) = bit_len;
    pq_copymsgbytes(buf, (char *) (result->data + sizeof(int32)), byte_len);
    
    PG_RETURN_POINTER(result);
}

/*
 * DNA2 send function (binary output)
 */
Datum
kmersearch_dna2_send(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *dna = (kmersearch_dna2 *) PG_GETARG_POINTER(0);
    StringInfoData buf;
    int32 bit_len = *((int32 *) dna->data);
    int byte_len = (bit_len + 7) / 8;
    
    pq_begintypsend(&buf);
    pq_sendint32(&buf, bit_len);
    pq_sendbytes(&buf, (char *) (dna->data + sizeof(int32)), byte_len);
    
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
    kmersearch_dna4 *result;
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
    result = (kmersearch_dna4 *) palloc0(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    
    /* Set bit length */
    *((int32 *) result->data) = bit_len;
    
    /* Encode sequence */
    data_ptr = result->data + sizeof(int32);
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
    
    PG_RETURN_POINTER(result);
}

/*
 * DNA4 output function
 */
Datum
kmersearch_dna4_out(PG_FUNCTION_ARGS)
{
    kmersearch_dna4 *dna = (kmersearch_dna4 *) PG_GETARG_POINTER(0);
    int32 bit_len = *((int32 *) dna->data);
    int char_len = bit_len / 4;
    char *result;
    bits8 *data_ptr = dna->data + sizeof(int32);
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
    PG_RETURN_CSTRING(result);
}

/*
 * DNA4 receive function (binary input)
 */
Datum
kmersearch_dna4_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    kmersearch_dna4 *result;
    int32 bit_len;
    int byte_len;
    
    bit_len = pq_getmsgint(buf, sizeof(int32));
    byte_len = (bit_len + 7) / 8;
    
    result = (kmersearch_dna4 *) palloc(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + byte_len);
    
    *((int32 *) result->data) = bit_len;
    pq_copymsgbytes(buf, (char *) (result->data + sizeof(int32)), byte_len);
    
    PG_RETURN_POINTER(result);
}

/*
 * DNA4 send function (binary output)
 */
Datum
kmersearch_dna4_send(PG_FUNCTION_ARGS)
{
    kmersearch_dna4 *dna = (kmersearch_dna4 *) PG_GETARG_POINTER(0);
    StringInfoData buf;
    int32 bit_len = *((int32 *) dna->data);
    int byte_len = (bit_len + 7) / 8;
    
    pq_begintypsend(&buf);
    pq_sendint32(&buf, bit_len);
    pq_sendbytes(&buf, (char *) (dna->data + sizeof(int32)), byte_len);
    
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}