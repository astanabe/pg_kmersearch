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

/* Forward declaration for comparison function */
static int dna_compare_simd(const uint8_t* a, const uint8_t* b, int bit_len);

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
PG_FUNCTION_INFO_V1(kmersearch_dna2_to_bytea);
PG_FUNCTION_INFO_V1(kmersearch_dna4_to_bytea);

/* BTree comparison functions */
PG_FUNCTION_INFO_V1(kmersearch_dna2_cmp);
PG_FUNCTION_INFO_V1(kmersearch_dna4_cmp);
PG_FUNCTION_INFO_V1(kmersearch_dna2_lt);
PG_FUNCTION_INFO_V1(kmersearch_dna2_le);
PG_FUNCTION_INFO_V1(kmersearch_dna2_gt);
PG_FUNCTION_INFO_V1(kmersearch_dna2_ge);
PG_FUNCTION_INFO_V1(kmersearch_dna2_ne);
PG_FUNCTION_INFO_V1(kmersearch_dna4_lt);
PG_FUNCTION_INFO_V1(kmersearch_dna4_le);
PG_FUNCTION_INFO_V1(kmersearch_dna4_gt);
PG_FUNCTION_INFO_V1(kmersearch_dna4_ge);
PG_FUNCTION_INFO_V1(kmersearch_dna4_ne);
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
    VarBit *a = PG_GETARG_VARBIT_P(0);
    VarBit *b = PG_GETARG_VARBIT_P(1);
    int32 bit_len_a = VARBITLEN(a);
    int32 bit_len_b = VARBITLEN(b);
    int32 byte_len;
    bool result;
    
    /* Compare bit lengths first */
    if (bit_len_a != bit_len_b)
        PG_RETURN_BOOL(false);
    
    /* Compare actual bit data, not the entire varbit structure */
    byte_len = (bit_len_a + 7) / 8;
    result = (byte_len == 0) || (memcmp(VARBITS(a), VARBITS(b), byte_len) == 0);
    PG_RETURN_BOOL(result);
}

/*
 * DNA4 equality operator
 */
Datum
kmersearch_dna4_eq(PG_FUNCTION_ARGS)
{
    VarBit *a = PG_GETARG_VARBIT_P(0);
    VarBit *b = PG_GETARG_VARBIT_P(1);
    int32 bit_len_a = VARBITLEN(a);
    int32 bit_len_b = VARBITLEN(b);
    int32 byte_len;
    bool result;
    
    /* Compare bit lengths first */
    if (bit_len_a != bit_len_b)
        PG_RETURN_BOOL(false);
    
    /* Compare actual bit data, not the entire varbit structure */
    byte_len = (bit_len_a + 7) / 8;
    result = (byte_len == 0) || (memcmp(VARBITS(a), VARBITS(b), byte_len) == 0);
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

/*
 * SIMD-optimized comparison functions
 */

/* Scalar comparison function (fallback) */
int
dna_compare_scalar(const uint8_t* a, const uint8_t* b, int bit_len)
{
    int byte_len = (bit_len + 7) / 8;
    return memcmp(a, b, byte_len);
}

/*
 * SIMD comparison functions - enabled when proper compiler flags are set
 */

#ifdef __x86_64__
/* AVX2 comparison function */
__attribute__((target("avx2")))
int
dna_compare_avx2(const uint8_t* a, const uint8_t* b, int bit_len)
{
    int byte_len = (bit_len + 7) / 8;
    int simd_len = byte_len & ~31;  /* Process 32 bytes at a time */
    
    for (int i = 0; i < simd_len; i += 32) {
        __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
        __m256i cmp = _mm256_cmpeq_epi8(va, vb);
        
        if (_mm256_movemask_epi8(cmp) != 0xFFFFFFFF) {
            /* Found difference, need detailed comparison */
            for (int j = 0; j < 32; j++) {
                if (a[i + j] != b[i + j]) {
                    return (a[i + j] < b[i + j]) ? -1 : 1;
                }
            }
        }
    }
    
    /* Handle remaining bytes with scalar comparison */
    for (int i = simd_len; i < byte_len; i++) {
        if (a[i] != b[i]) {
            return (a[i] < b[i]) ? -1 : 1;
        }
    }
    
    return 0;
}

/* AVX512 comparison function */
__attribute__((target("avx512f,avx512bw")))
int
dna_compare_avx512(const uint8_t* a, const uint8_t* b, int bit_len)
{
    int byte_len = (bit_len + 7) / 8;
    int simd_len = byte_len & ~63;  /* Process 64 bytes at a time */
    int remaining;
    
    for (int i = 0; i < simd_len; i += 64) {
        __m512i va = _mm512_loadu_si512((const __m512i*)(a + i));
        __m512i vb = _mm512_loadu_si512((const __m512i*)(b + i));
        __mmask64 cmp = _mm512_cmpeq_epi8_mask(va, vb);
        
        if (cmp != 0xFFFFFFFFFFFFFFFF) {
            /* Found difference, need detailed comparison */
            for (int j = 0; j < 64; j++) {
                if (a[i + j] != b[i + j]) {
                    return (a[i + j] < b[i + j]) ? -1 : 1;
                }
            }
        }
    }
    
    /* Handle remaining bytes with smaller SIMD or scalar */
    remaining = byte_len - simd_len;
    if (remaining >= 32) {
        return dna_compare_avx2(a + simd_len, b + simd_len, remaining * 8);
    } else {
        return dna_compare_scalar(a + simd_len, b + simd_len, remaining * 8);
    }
}
#endif

#ifdef __aarch64__
/* NEON comparison function */
__attribute__((target("+simd")))
int
dna_compare_neon(const uint8_t* a, const uint8_t* b, int bit_len)
{
    int byte_len = (bit_len + 7) / 8;
    int simd_len = byte_len & ~15;  /* Process 16 bytes at a time */
    
    for (int i = 0; i < simd_len; i += 16) {
        uint8x16_t va = vld1q_u8(a + i);
        uint8x16_t vb = vld1q_u8(b + i);
        uint8x16_t cmp = vceqq_u8(va, vb);
        
        /* Check if all elements are equal */
        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        if (vgetq_lane_u64(cmp64, 0) != 0xFFFFFFFFFFFFFFFF ||
            vgetq_lane_u64(cmp64, 1) != 0xFFFFFFFFFFFFFFFF) {
            /* Found difference, need detailed comparison */
            for (int j = 0; j < 16; j++) {
                if (a[i + j] != b[i + j]) {
                    return (a[i + j] < b[i + j]) ? -1 : 1;
                }
            }
        }
    }
    
    /* Handle remaining bytes with scalar comparison */
    for (int i = simd_len; i < byte_len; i++) {
        if (a[i] != b[i]) {
            return (a[i] < b[i]) ? -1 : 1;
        }
    }
    
    return 0;
}

/* SVE comparison function */
__attribute__((target("+sve")))
int
dna_compare_sve(const uint8_t* a, const uint8_t* b, int bit_len)
{
    int byte_len = (bit_len + 7) / 8;
    int vector_len = svcntb();
    
    for (int i = 0; i < byte_len; i += vector_len) {
        svbool_t pg = svwhilelt_b8(i, byte_len);
        svuint8_t va = svld1_u8(pg, a + i);
        svuint8_t vb = svld1_u8(pg, b + i);
        svbool_t cmp = svcmpeq_u8(pg, va, vb);
        
        if (!svptest_any(svptrue_b8(), cmp)) {
            /* Found difference, need detailed comparison */
            int remaining = (byte_len - i < vector_len) ? byte_len - i : vector_len;
            for (int j = 0; j < remaining; j++) {
                if (a[i + j] != b[i + j]) {
                    return (a[i + j] < b[i + j]) ? -1 : 1;
                }
            }
        }
    }
    
    return 0;
}
#endif

/* Main SIMD dispatch function */
static int
dna_compare_simd(const uint8_t* a, const uint8_t* b, int bit_len)
{
    /* Use the dispatch table for SIMD function selection */
    return simd_dispatch.dna_compare(a, b, bit_len);
}

/*
 * BTree comparison functions
 */

/*
 * DNA2 comparison function
 */
Datum
kmersearch_dna2_cmp(PG_FUNCTION_ARGS)
{
    VarBit *a = PG_GETARG_VARBIT_P(0);
    VarBit *b = PG_GETARG_VARBIT_P(1);
    int32 bit_len_a = VARBITLEN(a);
    int32 bit_len_b = VARBITLEN(b);
    int result;
    
    /* Compare bit lengths first */
    if (bit_len_a < bit_len_b)
        PG_RETURN_INT32(-1);
    if (bit_len_a > bit_len_b)
        PG_RETURN_INT32(1);
    
    /* Same bit length, compare bit data using SIMD */
    result = dna_compare_simd(VARBITS(a), VARBITS(b), bit_len_a);
    PG_RETURN_INT32(result);
}

/*
 * DNA4 comparison function
 */
Datum
kmersearch_dna4_cmp(PG_FUNCTION_ARGS)
{
    VarBit *a = PG_GETARG_VARBIT_P(0);
    VarBit *b = PG_GETARG_VARBIT_P(1);
    int32 bit_len_a = VARBITLEN(a);
    int32 bit_len_b = VARBITLEN(b);
    int result;
    
    /* Compare bit lengths first */
    if (bit_len_a < bit_len_b)
        PG_RETURN_INT32(-1);
    if (bit_len_a > bit_len_b)
        PG_RETURN_INT32(1);
    
    /* Same bit length, compare bit data using SIMD */
    result = dna_compare_simd(VARBITS(a), VARBITS(b), bit_len_a);
    PG_RETURN_INT32(result);
}

/*
 * DNA2 to BYTEA conversion (bit length + bit data for hash uniqueness)
 */
Datum
kmersearch_dna2_to_bytea(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int byte_len = (bit_len + 7) / 8;
    bytea *result;
    int32 net_bit_len;
    
    /* Allocate space for bit length (4 bytes) + bit data */
    result = (bytea *) palloc(VARHDRSZ + 4 + byte_len);
    SET_VARSIZE(result, VARHDRSZ + 4 + byte_len);
    
    /* Store bit length in network byte order for consistency */
    net_bit_len = htonl(bit_len);
    memcpy(VARDATA(result), &net_bit_len, 4);
    
    /* Store bit data */
    if (byte_len > 0) {
        memcpy(VARDATA(result) + 4, VARBITS(dna), byte_len);
    }
    
    PG_RETURN_BYTEA_P(result);
}

/*
 * DNA4 to BYTEA conversion (bit length + bit data for consistency)
 */
Datum
kmersearch_dna4_to_bytea(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int byte_len = (bit_len + 7) / 8;
    bytea *result;
    int32 net_bit_len;
    
    /* Allocate space for bit length (4 bytes) + bit data */
    result = (bytea *) palloc(VARHDRSZ + 4 + byte_len);
    SET_VARSIZE(result, VARHDRSZ + 4 + byte_len);
    
    /* Store bit length in network byte order for consistency */
    net_bit_len = htonl(bit_len);
    memcpy(VARDATA(result), &net_bit_len, 4);
    
    /* Store bit data */
    if (byte_len > 0) {
        memcpy(VARDATA(result) + 4, VARBITS(dna), byte_len);
    }
    
    PG_RETURN_BYTEA_P(result);
}

/*
 * DNA2 comparison operators
 */
Datum
kmersearch_dna2_lt(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna2_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp < 0);
}

Datum
kmersearch_dna2_le(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna2_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp <= 0);
}

Datum
kmersearch_dna2_gt(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna2_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp > 0);
}

Datum
kmersearch_dna2_ge(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna2_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp >= 0);
}

Datum
kmersearch_dna2_ne(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna2_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp != 0);
}

/*
 * DNA4 comparison operators
 */
Datum
kmersearch_dna4_lt(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna4_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp < 0);
}

Datum
kmersearch_dna4_le(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna4_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp <= 0);
}

Datum
kmersearch_dna4_gt(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna4_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp > 0);
}

Datum
kmersearch_dna4_ge(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna4_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp >= 0);
}

Datum
kmersearch_dna4_ne(PG_FUNCTION_ARGS)
{
    int32 cmp = DatumGetInt32(DirectFunctionCall2(kmersearch_dna4_cmp,
                                                  PG_GETARG_DATUM(0),
                                                  PG_GETARG_DATUM(1)));
    PG_RETURN_BOOL(cmp != 0);
}

/*
 * Hash functions for DNA2 and DNA4 types
 * Required for hash-based operations like GROUP BY
 */

PG_FUNCTION_INFO_V1(kmersearch_dna2_hash);
PG_FUNCTION_INFO_V1(kmersearch_dna4_hash);
PG_FUNCTION_INFO_V1(kmersearch_dna2_hash_extended);
PG_FUNCTION_INFO_V1(kmersearch_dna4_hash_extended);

Datum
kmersearch_dna2_hash(PG_FUNCTION_ARGS)
{
    VarBit *data = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(data);
    int32 byte_len = (bit_len + 7) / 8;
    uint32 hash;
    
    /* Use PostgreSQL's hash_any function for consistent hashing */
    hash = DatumGetUInt32(hash_any(VARBITS(data), byte_len));
    
    PG_RETURN_INT32(hash);
}

Datum
kmersearch_dna4_hash(PG_FUNCTION_ARGS)
{
    VarBit *data = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(data);
    int32 byte_len = (bit_len + 7) / 8;
    uint32 hash;
    
    /* Use PostgreSQL's hash_any function for consistent hashing */
    hash = DatumGetUInt32(hash_any(VARBITS(data), byte_len));
    
    PG_RETURN_INT32(hash);
}

Datum
kmersearch_dna2_hash_extended(PG_FUNCTION_ARGS)
{
    VarBit *data = PG_GETARG_VARBIT_P(0);
    uint64 seed = PG_GETARG_INT64(1);
    int32 bit_len = VARBITLEN(data);
    int32 byte_len = (bit_len + 7) / 8;
    uint64 hash;
    
    /* Use PostgreSQL's hash_any_extended function for 64-bit hashing with seed */
    hash = DatumGetUInt64(hash_any_extended(VARBITS(data), byte_len, seed));
    
    PG_RETURN_INT64(hash);
}

Datum
kmersearch_dna4_hash_extended(PG_FUNCTION_ARGS)
{
    VarBit *data = PG_GETARG_VARBIT_P(0);
    uint64 seed = PG_GETARG_INT64(1);
    int32 bit_len = VARBITLEN(data);
    int32 byte_len = (bit_len + 7) / 8;
    uint64 hash;
    
    /* Use PostgreSQL's hash_any_extended function for 64-bit hashing with seed */
    hash = DatumGetUInt64(hash_any_extended(VARBITS(data), byte_len, seed));
    
    PG_RETURN_INT64(hash);
}

/*
 * SIMD encoding/decoding functions
 */

/* Scalar implementation of DNA2 encoding */
void dna2_encode_scalar(const char* input, uint8_t* output, int len)
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

#ifdef __x86_64__
__attribute__((target("avx2")))
void dna2_encode_avx2(const char* input, uint8_t* output, int len)
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

__attribute__((target("avx512f,avx512bw")))
void dna2_encode_avx512(const char* input, uint8_t* output, int len)
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
#endif /* __x86_64__ */

#ifdef __aarch64__
__attribute__((target("+simd")))
void dna2_encode_neon(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 2 + 7) / 8;
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    int i, j;
    uint8x16_t chars, mask_C, mask_G, mask_T, mask_U, encoded;
    uint8_t temp[16];
    
    memset(output, 0, byte_len);
    
    /* Process 16 characters at a time with NEON */
    for (i = 0; i < simd_len; i += 16) {
        chars = vld1q_u8((uint8_t*)(input + i));
        
        /* Create comparison masks for each DNA base */
        /* A=00, so no mask needed */
        mask_C = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('C')),
                          vceqq_u8(chars, vdupq_n_u8('c')));
        mask_G = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('G')),
                          vceqq_u8(chars, vdupq_n_u8('g')));
        mask_T = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('T')),
                          vceqq_u8(chars, vdupq_n_u8('t')));
        mask_U = vorrq_u8(vceqq_u8(chars, vdupq_n_u8('U')),
                          vceqq_u8(chars, vdupq_n_u8('u')));
        
        /* Combine T and U masks */
        mask_T = vorrq_u8(mask_T, mask_U);
        
        /* Generate 2-bit encoded values */
        encoded = vdupq_n_u8(0);
        encoded = vorrq_u8(encoded, vandq_u8(mask_C, vdupq_n_u8(1)));
        encoded = vorrq_u8(encoded, vandq_u8(mask_G, vdupq_n_u8(2)));
        encoded = vorrq_u8(encoded, vandq_u8(mask_T, vdupq_n_u8(3)));
        
        /* Store and pack bits */
        vst1q_u8(temp, encoded);
        
        for (j = 0; j < 16; j++) {
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
    for (i = simd_len; i < len; i++) {
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
void dna2_encode_sve(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 2 + 7) / 8;
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);
    int i, j;
    svbool_t pg, mask_C, mask_G, mask_T, mask_U;
    svuint8_t chars, encoded;
    uint8_t temp[sve_len];
    
    memset(output, 0, byte_len);
    
    for (i = 0; i < simd_len; i += sve_len) {
        pg = svwhilelt_b8_s32(i, len);
        chars = svld1_u8(pg, (uint8_t*)(input + i));
        
        /* Create comparison masks for each DNA base */
        /* A=00, so no mask needed */
        mask_C = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('C')),
                               svcmpeq_u8(pg, chars, svdup_n_u8('c')));
        mask_G = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('G')),
                               svcmpeq_u8(pg, chars, svdup_n_u8('g')));
        mask_T = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('T')),
                               svcmpeq_u8(pg, chars, svdup_n_u8('t')));
        mask_U = svorr_b_z(pg, svcmpeq_u8(pg, chars, svdup_n_u8('U')),
                               svcmpeq_u8(pg, chars, svdup_n_u8('u')));
        
        /* Combine T and U masks */
        mask_T = svorr_b_z(pg, mask_T, mask_U);
        
        /* Generate 2-bit encoded values */
        encoded = svdup_n_u8(0);
        encoded = svorr_u8_m(mask_C, encoded, svdup_n_u8(1));
        encoded = svorr_u8_m(mask_G, encoded, svdup_n_u8(2));
        encoded = svorr_u8_m(mask_T, encoded, svdup_n_u8(3));
        
        /* Store and pack bits */
        svst1_u8(pg, temp, encoded);
        
        for (j = 0; j < sve_len && (i + j) < len; j++) {
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
    for (i = simd_len; i < len; i++) {
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
#endif /* __aarch64__ */

/* Scalar implementation of DNA2 decoding */
void dna2_decode_scalar(const uint8_t* input, char* output, int len)
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

#ifdef __x86_64__
__attribute__((target("avx2")))
void dna2_decode_avx2(const uint8_t* input, char* output, int len)
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

__attribute__((target("avx512f,avx512bw")))
void dna2_decode_avx512(const uint8_t* input, char* output, int len)
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
#endif /* __x86_64__ */

#ifdef __aarch64__
__attribute__((target("+simd")))
void dna2_decode_neon(const uint8_t* input, char* output, int len)
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

__attribute__((target("+sve")))
void dna2_decode_sve(const uint8_t* input, char* output, int len)
{
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    int i, j, b;
    
    for (i = 0; i < simd_len; i += sve_len) {
        /* Fixed: Load correct number of bytes for sve_len characters */
        int bytes_needed = (sve_len * 2 + 7) / 8;
        uint8_t temp[bytes_needed];
        
        /* Fixed: Load from correct position in input array */
        for (b = 0; b < bytes_needed && (i * 2 / 8 + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[i * 2 / 8 + b];
        }
        
        /* Decode each 2-bit pair using relative positions in temp array */
        for (j = 0; j < sve_len && (i + j) < len; j++) {
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
    for (i = simd_len; i < len; i++) {
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
#endif /* __aarch64__ */

/*
 * DNA4 Encoding Functions
 */

/* Scalar implementation for DNA4 encoding */
void dna4_encode_scalar(const char* input, uint8_t* output, int len)
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

/* Scalar implementation for DNA4 decoding */
void dna4_decode_scalar(const uint8_t* input, char* output, int len)
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
/* AVX2 implementation for DNA4 encoding */
__attribute__((target("avx2")))
void dna4_encode_avx2(const char* input, uint8_t* output, int len)
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

/* AVX2 implementation for DNA4 decoding */
__attribute__((target("avx2")))
void dna4_decode_avx2(const uint8_t* input, char* output, int len)
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

/* AVX512 implementation for DNA4 encoding */
__attribute__((target("avx512bw")))
void dna4_encode_avx512(const char* input, uint8_t* output, int len)
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

/* AVX512 implementation for DNA4 decoding */
__attribute__((target("avx512bw")))
void dna4_decode_avx512(const uint8_t* input, char* output, int len)
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
#endif /* __x86_64__ */

#ifdef __aarch64__
/* NEON implementation for DNA4 encoding */
__attribute__((target("+simd")))
void dna4_encode_neon(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    int i, j;
    uint8x16_t chars;
    uint8_t temp[16];
    
    memset(output, 0, byte_len);
    
    /* Process 16 characters at a time with NEON */
    for (i = 0; i < simd_len; i += 16) {
        chars = vld1q_u8((uint8_t*)(input + i));
        
        /* Use lookup table approach for DNA4 encoding */
        vst1q_u8(temp, chars);
        
        for (j = 0; j < 16; j++) {
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
    for (i = simd_len; i < len; i++) {
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

/* NEON implementation for DNA4 decoding */
__attribute__((target("+simd")))
void dna4_decode_neon(const uint8_t* input, char* output, int len)
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

/* SVE implementation for DNA4 encoding */
__attribute__((target("+sve")))
void dna4_encode_sve(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    int i, j;
    uint8_t temp[sve_len];
    
    memset(output, 0, byte_len);
    
    for (i = 0; i < simd_len; i += sve_len) {
        svbool_t pg = svwhilelt_b8_s32(i, len);
        svuint8_t chars = svld1_u8(pg, (uint8_t*)(input + i));
        
        /* Use lookup table approach for DNA4 encoding */
        svst1_u8(pg, temp, chars);
        
        for (j = 0; j < sve_len && (i + j) < len; j++) {
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
    for (i = simd_len; i < len; i++) {
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

/* SVE implementation for DNA4 decoding */
__attribute__((target("+sve")))
void dna4_decode_sve(const uint8_t* input, char* output, int len)
{
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    int i, j, b;
    
    for (i = 0; i < simd_len; i += sve_len) {
        /* Fixed: Load correct number of bytes for sve_len characters */
        int bytes_needed = (sve_len * 4 + 7) / 8;
        uint8_t temp[bytes_needed];
        
        /* Fixed: Load from correct position in input array */
        for (b = 0; b < bytes_needed && (i * 4 / 8 + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[i * 4 / 8 + b];
        }
        
        /* Decode each 4-bit nibble using relative positions in temp array */
        for (j = 0; j < sve_len && (i + j) < len; j++) {
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
    for (i = simd_len; i < len; i++) {
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
#endif /* __aarch64__ */
