/*
 * kmersearch_datatype.c - DNA2/DNA4 datatype input/output and type conversion functions
 *
 * This module contains all DNA2 and DNA4 datatype-related functionality including:
 * - Input/output functions (in, out, recv, send)
 * - Type conversion and utility functions
 * - Length calculation functions  
 * - String conversion functions
 * - Validation helper functions
 * - Encoding/decoding tables and lookup functions
 */

#include "kmersearch.h"

/* PostgreSQL function info declarations for datatype functions */
PG_FUNCTION_INFO_V1(kmersearch_dna2_in);
PG_FUNCTION_INFO_V1(kmersearch_dna2_out);
PG_FUNCTION_INFO_V1(kmersearch_dna2_recv);
PG_FUNCTION_INFO_V1(kmersearch_dna2_send);
PG_FUNCTION_INFO_V1(kmersearch_dna4_in);
PG_FUNCTION_INFO_V1(kmersearch_dna4_out);
PG_FUNCTION_INFO_V1(kmersearch_dna4_recv);
PG_FUNCTION_INFO_V1(kmersearch_dna4_send);
PG_FUNCTION_INFO_V1(kmersearch_dna2_eq);
PG_FUNCTION_INFO_V1(kmersearch_dna4_eq);
PG_FUNCTION_INFO_V1(kmersearch_dna2_char_length);
PG_FUNCTION_INFO_V1(kmersearch_dna4_char_length);

/*
 * DNA4 to DNA2 expansion table
 * Each entry contains: [count, base1, base2, base3, base4]
 * For DNA4 degenerate codes that expand to multiple DNA2 bases
 */
static const uint8 kmersearch_dna4_to_dna2_table[16][5] = {
    {1, 0, 0, 0, 0},     /* 0000 - A */
    {1, 1, 0, 0, 0},     /* 0001 - C */
    {2, 0, 1, 0, 0},     /* 0010 - M (A,C) */
    {1, 2, 0, 0, 0},     /* 0011 - G */
    {2, 0, 2, 0, 0},     /* 0100 - R (A,G) */
    {2, 1, 2, 0, 0},     /* 0101 - S (C,G) */
    {3, 0, 1, 2, 0},     /* 0110 - V (A,C,G) */
    {1, 3, 0, 0, 0},     /* 0111 - T */
    {2, 0, 3, 0, 0},     /* 1000 - W (A,T) */
    {2, 1, 3, 0, 0},     /* 1001 - Y (C,T) */
    {3, 0, 1, 3, 0},     /* 1010 - H (A,C,T) */
    {2, 2, 3, 0, 0},     /* 1011 - K (G,T) */
    {3, 0, 2, 3, 0},     /* 1100 - D (A,G,T) */
    {3, 1, 2, 3, 0},     /* 1101 - B (C,G,T) */
    {4, 0, 1, 2, 3}      /* 1110 - N (A,C,G,T) */
};

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
    
    /* Encode sequence using SIMD dispatch */
    data_ptr = VARBITS(result);
    simd_dispatch.dna2_encode(input_string, (uint8_t*)data_ptr, input_len);
    
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
    
    if (dna == NULL)
        ereport(ERROR, (errmsg("input DNA sequence is NULL")));
    
    if (bit_len < 0)
        ereport(ERROR, (errmsg("invalid bit length: %d", bit_len)));
    
    if (bit_len % 2 != 0)
        ereport(ERROR, (errmsg("bit length must be even for DNA2")));
    
    result = (char *) palloc(char_len + 1);
    
    /* Decode sequence using SIMD dispatch */
    simd_dispatch.dna2_decode((uint8_t*)data_ptr, result, char_len);
    
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
    
    /* Encode sequence using SIMD dispatch */
    data_ptr = VARBITS(result);
    simd_dispatch.dna4_encode(input_string, (uint8_t*)data_ptr, input_len);
    
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
    
    if (dna == NULL)
        ereport(ERROR, (errmsg("input DNA sequence is NULL")));
    
    if (bit_len < 0)
        ereport(ERROR, (errmsg("invalid bit length: %d", bit_len)));
    
    if (bit_len % 4 != 0)
        ereport(ERROR, (errmsg("bit length must be multiple of 4 for DNA4")));
    
    result = (char *) palloc(char_len + 1);
    
    /* Decode sequence using SIMD dispatch */
    simd_dispatch.dna4_decode((uint8_t*)data_ptr, result, char_len);
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
 * DNA2 character length function
 */
Datum
kmersearch_dna2_char_length(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int32 char_len = bit_len / 2;  /* 2 bits per character */
    
    PG_RETURN_INT32(char_len);
}

/*
 * DNA4 character length function
 */
Datum
kmersearch_dna4_char_length(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int32 char_len = bit_len / 4;  /* 4 bits per character */
    
    PG_RETURN_INT32(char_len);
}

/*
 * Convert DNA2 VarBit to string (used internally)
 */
char *
kmersearch_dna2_to_string(VarBit *dna)
{
    int32 bit_len = VARBITLEN(dna);
    int char_len = bit_len / 2;
    char *result;
    bits8 *data_ptr = VARBITS(dna);
    
    if (dna == NULL)
        ereport(ERROR, (errmsg("input DNA sequence is NULL")));
    
    if (bit_len < 0)
        ereport(ERROR, (errmsg("invalid bit length: %d", bit_len)));
    
    if (bit_len % 2 != 0)
        ereport(ERROR, (errmsg("bit length must be even for DNA2")));
    
    result = (char *) palloc(char_len + 1);
    
    /* Decode sequence using SIMD dispatch */
    simd_dispatch.dna2_decode((uint8_t*)data_ptr, result, char_len);
    
    return result;
}

/*
 * Convert DNA4 VarBit to string (used internally)
 */
char *
kmersearch_dna4_to_string(VarBit *dna)
{
    int32 bit_len = VARBITLEN(dna);
    int char_len = bit_len / 4;
    char *result;
    bits8 *data_ptr = VARBITS(dna);
    
    if (dna == NULL)
        ereport(ERROR, (errmsg("input DNA sequence is NULL")));
    
    if (bit_len < 0)
        ereport(ERROR, (errmsg("invalid bit length: %d", bit_len)));
    
    if (bit_len % 4 != 0)
        ereport(ERROR, (errmsg("bit length must be multiple of 4 for DNA4")));
    
    result = (char *) palloc(char_len + 1);
    
    /* Decode sequence using SIMD dispatch */
    simd_dispatch.dna4_decode((uint8_t*)data_ptr, result, char_len);
    
    return result;
}