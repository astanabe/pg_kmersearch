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
    dna2_encode(input_string, (uint8_t*)data_ptr, input_len);
    
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
    dna2_decode((uint8_t*)data_ptr, result, char_len);
    
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
    dna4_encode(input_string, (uint8_t*)data_ptr, input_len);
    
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
    dna4_decode((uint8_t*)data_ptr, result, char_len);
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
    dna2_decode((uint8_t*)data_ptr, result, char_len);
    
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
    dna4_decode((uint8_t*)data_ptr, result, char_len);
    
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

/*
 * Main dispatch functions with threshold-based SIMD selection
 */
void
dna2_encode(const char* input, uint8_t* output, int len)
{
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && len >= SIMD_ENCODE_AVX512_THRESHOLD) {
        dna2_encode_avx512(input, output, len);
        return;
    }
    if (simd_capability >= SIMD_AVX2 && len >= SIMD_ENCODE_AVX2_THRESHOLD) {
        dna2_encode_avx2(input, output, len);
        return;
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && len >= SIMD_ENCODE_SVE_THRESHOLD) {
        dna2_encode_sve(input, output, len);
        return;
    }
    if (simd_capability >= SIMD_NEON && len >= SIMD_ENCODE_NEON_THRESHOLD) {
        dna2_encode_neon(input, output, len);
        return;
    }
#endif
    dna2_encode_scalar(input, output, len);
}

void
dna2_decode(const uint8_t* input, char* output, int bit_len)
{
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && bit_len >= SIMD_DECODE_AVX512_THRESHOLD) {
        dna2_decode_avx512(input, output, bit_len);
        return;
    }
    if (simd_capability >= SIMD_AVX2 && bit_len >= SIMD_DECODE_AVX2_THRESHOLD) {
        dna2_decode_avx2(input, output, bit_len);
        return;
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE2 && bit_len >= SIMD_DECODE_SVE_THRESHOLD) {
        dna2_decode_sve2(input, output, bit_len);
        return;
    }
    if (simd_capability >= SIMD_SVE && bit_len >= SIMD_DECODE_SVE_THRESHOLD) {
        dna2_decode_sve(input, output, bit_len);
        return;
    }
    if (simd_capability >= SIMD_NEON && bit_len >= SIMD_DECODE_NEON_THRESHOLD) {
        dna2_decode_neon(input, output, bit_len);
        return;
    }
#endif
    dna2_decode_scalar(input, output, bit_len);
}

void
dna4_encode(const char* input, uint8_t* output, int len)
{
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && len >= SIMD_ENCODE_AVX512_THRESHOLD) {
        dna4_encode_avx512(input, output, len);
        return;
    }
    if (simd_capability >= SIMD_AVX2 && len >= SIMD_ENCODE_AVX2_THRESHOLD) {
        dna4_encode_avx2(input, output, len);
        return;
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && len >= SIMD_ENCODE_SVE_THRESHOLD) {
        dna4_encode_sve(input, output, len);
        return;
    }
    if (simd_capability >= SIMD_NEON && len >= SIMD_ENCODE_NEON_THRESHOLD) {
        dna4_encode_neon(input, output, len);
        return;
    }
#endif
    dna4_encode_scalar(input, output, len);
}

void
dna4_decode(const uint8_t* input, char* output, int bit_len)
{
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && bit_len >= SIMD_DECODE_AVX512_THRESHOLD) {
        dna4_decode_avx512(input, output, bit_len);
        return;
    }
    if (simd_capability >= SIMD_AVX2 && bit_len >= SIMD_DECODE_AVX2_THRESHOLD) {
        dna4_decode_avx2(input, output, bit_len);
        return;
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && bit_len >= SIMD_DECODE_SVE_THRESHOLD) {
        dna4_decode_sve(input, output, bit_len);
        return;
    }
    if (simd_capability >= SIMD_NEON && bit_len >= SIMD_DECODE_NEON_THRESHOLD) {
        dna4_decode_neon(input, output, bit_len);
        return;
    }
#endif
    dna4_decode_scalar(input, output, bit_len);
}

int
dna_compare(const uint8_t* a, const uint8_t* b, int bit_len)
{
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && bit_len >= SIMD_COMPARE_AVX512_THRESHOLD) {
        return dna_compare_avx512(a, b, bit_len);
    }
    if (simd_capability >= SIMD_AVX2 && bit_len >= SIMD_COMPARE_AVX2_THRESHOLD) {
        return dna_compare_avx2(a, b, bit_len);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && bit_len >= SIMD_COMPARE_SVE_THRESHOLD) {
        return dna_compare_sve(a, b, bit_len);
    }
    if (simd_capability >= SIMD_NEON && bit_len >= SIMD_COMPARE_NEON_THRESHOLD) {
        return dna_compare_neon(a, b, bit_len);
    }
#endif
    return dna_compare_scalar(a, b, bit_len);
}

/* Main SIMD dispatch function - now renamed to dna_compare() */
static int
dna_compare_simd(const uint8_t* a, const uint8_t* b, int bit_len)
{
    /* This function is deprecated - use dna_compare() instead */
    return dna_compare(a, b, bit_len);
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
__attribute__((target("avx2,bmi2")))
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
        
        /* Extract 2-bit values and pack using PDEP */
        _mm256_storeu_si256((__m256i*)temp, encoded);
        
        /* Process 32 bases (64 bits) using PDEP for efficient bit packing */
        for (int j = 0; j < 32; j += 32) {
            uint64_t base_bits0 = 0, base_bits1 = 0;
            uint64_t deposited0, deposited1;
            int byte_offset = (i + j) / 4;
            
            /* Collect 2-bit values into 64-bit integers */
            for (int k = 0; k < 16 && (j + k) < 32; k++) {
                base_bits0 |= ((uint64_t)temp[j + k] << (k * 2));
            }
            for (int k = 0; k < 16 && (j + k + 16) < 32; k++) {
                base_bits1 |= ((uint64_t)temp[j + k + 16] << (k * 2));
            }
            
            /* Use PDEP to deposit bits efficiently */
            deposited0 = _pdep_u64(base_bits0, 0xFFFFFFFF);
            deposited1 = _pdep_u64(base_bits1, 0xFFFFFFFF);
            
            /* Write packed data to output */
            if (byte_offset < byte_len) {
                memcpy(&output[byte_offset], &deposited0, 
                       (byte_offset + 4 <= byte_len) ? 4 : (byte_len - byte_offset));
            }
            if (byte_offset + 4 < byte_len) {
                memcpy(&output[byte_offset + 4], &deposited1, 
                       (byte_offset + 8 <= byte_len) ? 4 : (byte_len - byte_offset - 4));
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
__attribute__((target("avx2,bmi2")))
void dna2_decode_avx2(const uint8_t* input, char* output, int len)
{
    /* Process 32 characters at a time with AVX2 SIMD optimizations */
    int simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    /* Create lookup table for VPSHUFB */
    const __m256i decode_lut = _mm256_setr_epi8(
        'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
        'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
        'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T',
        'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'
    );
    
    const __m256i mask_2bits = _mm256_set1_epi8(0x03);
    
    for (int i = 0; i < simd_len; i += 32) {
        /* Calculate how many bytes needed for 32 characters (32 * 2 bits = 64 bits = 8 bytes) */
        int byte_offset = (i * 2) / 8;
        int bytes_to_load = 9; /* Load extra byte for boundary cases */
        uint8_t temp[16] __attribute__((aligned(16))) = {0};
        
        /* Load required bytes */
        for (int b = 0; b < bytes_to_load && (byte_offset + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[byte_offset + b];
        }
        
        /* Use PEXT to extract 2-bit values efficiently (process 8 chars at a time) */
        for (int chunk = 0; chunk < 4; chunk++) {
            int char_offset = chunk * 8;
            int bit_offset = (i % 4) * 2;  /* Bit offset within byte */
            uint64_t src_bits;
            uint16_t extracted = 0;
            
            /* Load 64 bits starting from appropriate position */
            memcpy(&src_bits, &temp[char_offset * 2 / 8], 8);
            
            /* Extract 16 bits (8 * 2-bit values) using bit manipulation */
            for (int j = 0; j < 8; j++) {
                int total_bit_pos = bit_offset + j * 2;
                int byte_pos = total_bit_pos / 8;
                int bit_pos = total_bit_pos % 8;
                
                if (byte_pos < 8) {
                    uint8_t val = (temp[byte_pos] >> (6 - bit_pos)) & 0x3;
                    if (bit_pos == 7 && byte_pos + 1 < 8) {
                        val = ((temp[byte_pos] & 0x1) << 1) | ((temp[byte_pos + 1] >> 7) & 0x1);
                    }
                    extracted |= (val << (j * 2));
                }
            }
            
            /* Create vector from extracted values */
            
            /* Extract individual 2-bit values and decode using PSHUFB */
            {
                uint8_t decoded_chars[8];
                int j;
                for (j = 0; j < 8; j++) {
                    uint8_t two_bit_val = (extracted >> (j * 2)) & 0x3;
                    decoded_chars[j] = (two_bit_val < 4) ? "ACGT"[two_bit_val] : 'A';
            }
            
                /* Store decoded characters */
                memcpy(&output[i + char_offset], decoded_chars, 8);
            }
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

__attribute__((target("avx512f,avx512bw,avx512vbmi")))
void dna2_decode_avx512(const uint8_t* input, char* output, int len)
{
    /* Process 64 characters at a time with AVX512 VBMI optimizations */
    int simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    /* Create lookup table for VPERMB (64-byte table) */
    const __m512i decode_lut = _mm512_set_epi8(
        'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
        'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
        'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A',
        'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'A'
    );
    
    const __m512i mask_2bits = _mm512_set1_epi8(0x03);
    
    for (int i = 0; i < simd_len; i += 64) {
        /* 64 characters need 128 bits = 16 bytes */
        int byte_offset = (i * 2) / 8;
        uint8_t temp[32] __attribute__((aligned(32))) = {0};
        
        /* Load required bytes */
        int bytes_to_load = 17; /* Extra byte for boundary cases */
        for (int b = 0; b < bytes_to_load && (byte_offset + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[byte_offset + b];
        }
        
        /* Process in 16-character chunks for better efficiency */
        for (int chunk = 0; chunk < 4; chunk++) {
            __m128i extracted;
            uint8_t extracted_vals[16];
            __m512i expanded, masked, decoded;
            
            extracted = _mm_setzero_si128();
            
            /* Extract 2-bit values for 16 characters */
            for (int j = 0; j < 16; j++) {
                int char_idx = chunk * 16 + j;
                int total_bit_pos = char_idx * 2;
                int byte_pos = total_bit_pos / 8;
                int bit_pos = total_bit_pos % 8;
                uint8_t val = 0;
                
                if (byte_pos < 17) {
                    if (bit_pos <= 6) {
                        val = (temp[byte_pos] >> (6 - bit_pos)) & 0x3;
                    } else if (byte_pos + 1 < 17) {
                        val = ((temp[byte_pos] & 0x1) << 1) | ((temp[byte_pos + 1] >> 7) & 0x1);
                    }
                }
                
                extracted_vals[j] = val;
            }
            extracted = _mm_loadu_si128((__m128i*)extracted_vals);
            
            /* Expand to 512-bit for VPERMB */
            expanded = _mm512_zextsi128_si512(extracted);
            expanded = _mm512_maskz_permutexvar_epi8((__mmask64)0xFFFF, expanded, expanded);
            
            /* Use VPERMB for lookup */
            masked = _mm512_and_si512(expanded, mask_2bits);
            decoded = _mm512_permutexvar_epi8(masked, decode_lut);
            
            /* Store 16 decoded characters */
            _mm_storeu_si128((__m128i*)&output[i + chunk * 16], _mm512_castsi512_si128(decoded));
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
    /* Process 16 characters at a time with NEON VTBL optimizations */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    /* Create lookup table for VTBL */
    const uint8x16_t decode_lut = vld1q_u8((const uint8_t*)"ACGTACGTACGTACGT");
    const uint8x16_t mask_2bits = vdupq_n_u8(0x03);
    
    for (int i = 0; i < simd_len; i += 16) {
        int byte_offset = (i * 2) / 8;
        uint8_t temp[8] __attribute__((aligned(8))) = {0};
        
        /* Load required bytes (4 bytes for 16 chars + extra for boundaries) */
        int bytes_to_load = 5;
        for (int b = 0; b < bytes_to_load && (byte_offset + b) < (len * 2 + 7) / 8; b++) {
            temp[b] = input[byte_offset + b];
        }
        
        /* Extract 2-bit values for 16 characters */
        uint8x16_t extracted = vdupq_n_u8(0);
        uint8_t extracted_vals[16];
        
        for (int j = 0; j < 16; j++) {
            int bit_pos = j * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t val = 0;
            
            if (byte_pos < 5) {
                if (bit_offset <= 6) {
                    val = (temp[byte_pos] >> (6 - bit_offset)) & 0x3;
                } else if (byte_pos + 1 < 5) {
                    val = ((temp[byte_pos] & 0x1) << 1) | ((temp[byte_pos + 1] >> 7) & 0x1);
                }
            }
            
            extracted_vals[j] = val;
        }
        
        /* Load extracted values into NEON register */
        extracted = vld1q_u8(extracted_vals);
        
        /* Use VTBL for lookup */
        uint8x16_t masked = vandq_u8(extracted, mask_2bits);
        uint8x16_t decoded = vqtbl1q_u8(decode_lut, masked);
        
        /* Store result */
        vst1q_u8((uint8_t*)&output[i], decoded);
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

__attribute__((target("+sve2")))
void dna2_decode_sve2(const uint8_t* input, char* output, int len)
{
    /* Get SVE vector length */
    int sve_len = svcntb();
    int simd_len = len & ~(sve_len - 1);  /* Round down to SVE vector multiple */
    int i;
    svbool_t pg = svptrue_b8();
    
    /* Process multiple characters at once using SVE2 */
    for (i = 0; i < simd_len; i += sve_len) {
        int bytes_needed = (sve_len * 2 + 7) / 8;
        svuint8_t vec_data;
        svuint8_t extracted_bits;
        svuint8_t decoded_chars;
        int j;
        
        /* Load compressed data */
        vec_data = svld1_u8(pg, &input[i * 2 / 8]);
        
        /* Extract 2-bit values for each character position */
        /* SVE2 provides more efficient bit manipulation than SVE */
        for (j = 0; j < sve_len && (i + j) < len; j++) {
            int bit_pos = (i + j) * 2;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t encoded;
            
            /* Extract 2-bit value */
            if (bit_offset <= 6) {
                encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
            } else {
                /* bit_offset == 7: bits span across two bytes */
                encoded = (input[byte_pos] & 0x1) << 1;
                if (byte_pos + 1 < (len * 2 + 7) / 8) {
                    encoded |= (input[byte_pos + 1] >> 7) & 0x1;
                }
            }
            
            /* Range check and decode */
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
        
        /* Extract 2-bit value */
        if (bit_offset <= 6) {
            encoded = (input[byte_pos] >> (6 - bit_offset)) & 0x3;
        } else {
            /* bit_offset == 7: bits span across two bytes */
            encoded = (input[byte_pos] & 0x1) << 1;
            if (byte_pos + 1 < (len * 2 + 7) / 8) {
                encoded |= (input[byte_pos + 1] >> 7) & 0x1;
            }
        }
        
        /* Range check and decode */
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
__attribute__((target("avx2,bmi2")))
void dna4_encode_avx2(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    int simd_len;
    /* Create comparison vectors for efficient encoding */
    const __m256i vec_A = _mm256_set1_epi8('A');
    const __m256i vec_C = _mm256_set1_epi8('C');
    const __m256i vec_G = _mm256_set1_epi8('G');
    const __m256i vec_T = _mm256_set1_epi8('T');
    const __m256i vec_a = _mm256_set1_epi8('a');
    const __m256i vec_c = _mm256_set1_epi8('c');
    const __m256i vec_g = _mm256_set1_epi8('g');
    const __m256i vec_t = _mm256_set1_epi8('t');
    /* Degenerate base vectors */
    const __m256i vec_M = _mm256_set1_epi8('M');
    const __m256i vec_R = _mm256_set1_epi8('R');
    const __m256i vec_W = _mm256_set1_epi8('W');
    const __m256i vec_S = _mm256_set1_epi8('S');
    const __m256i vec_Y = _mm256_set1_epi8('Y');
    const __m256i vec_K = _mm256_set1_epi8('K');
    const __m256i vec_V = _mm256_set1_epi8('V');
    const __m256i vec_H = _mm256_set1_epi8('H');
    const __m256i vec_D = _mm256_set1_epi8('D');
    const __m256i vec_B = _mm256_set1_epi8('B');
    const __m256i vec_N = _mm256_set1_epi8('N');
    
    memset(output, 0, byte_len);
    
    /* Process 32 characters at a time with AVX2 */
    simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    for (int i = 0; i < simd_len; i += 32) {
        __m256i chars = _mm256_loadu_si256((__m256i*)(input + i));
        
        /* Create masks for each base type */
        __m256i mask_A = _mm256_or_si256(_mm256_cmpeq_epi8(chars, vec_A), _mm256_cmpeq_epi8(chars, vec_a));
        __m256i mask_C = _mm256_or_si256(_mm256_cmpeq_epi8(chars, vec_C), _mm256_cmpeq_epi8(chars, vec_c));
        __m256i mask_G = _mm256_or_si256(_mm256_cmpeq_epi8(chars, vec_G), _mm256_cmpeq_epi8(chars, vec_g));
        __m256i mask_T = _mm256_or_si256(_mm256_cmpeq_epi8(chars, vec_T), _mm256_cmpeq_epi8(chars, vec_t));
        
        /* Degenerate bases */
        __m256i mask_M = _mm256_cmpeq_epi8(chars, vec_M);
        __m256i mask_R = _mm256_cmpeq_epi8(chars, vec_R);
        __m256i mask_W = _mm256_cmpeq_epi8(chars, vec_W);
        __m256i mask_S = _mm256_cmpeq_epi8(chars, vec_S);
        __m256i mask_Y = _mm256_cmpeq_epi8(chars, vec_Y);
        __m256i mask_K = _mm256_cmpeq_epi8(chars, vec_K);
        __m256i mask_V = _mm256_cmpeq_epi8(chars, vec_V);
        __m256i mask_H = _mm256_cmpeq_epi8(chars, vec_H);
        __m256i mask_D = _mm256_cmpeq_epi8(chars, vec_D);
        __m256i mask_B = _mm256_cmpeq_epi8(chars, vec_B);
        __m256i mask_N = _mm256_cmpeq_epi8(chars, vec_N);
        
        /* Generate 4-bit encoded values using SIMD */
        __m256i encoded = _mm256_setzero_si256();
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_A, _mm256_set1_epi8(0x01)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_C, _mm256_set1_epi8(0x02)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_G, _mm256_set1_epi8(0x04)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_T, _mm256_set1_epi8(0x08)));
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_M, _mm256_set1_epi8(0x03))); /* A|C */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_R, _mm256_set1_epi8(0x05))); /* A|G */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_W, _mm256_set1_epi8(0x09))); /* A|T */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_S, _mm256_set1_epi8(0x06))); /* C|G */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_Y, _mm256_set1_epi8(0x0A))); /* C|T */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_K, _mm256_set1_epi8(0x0C))); /* G|T */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_V, _mm256_set1_epi8(0x07))); /* A|C|G */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_H, _mm256_set1_epi8(0x0B))); /* A|C|T */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_D, _mm256_set1_epi8(0x0D))); /* A|G|T */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_B, _mm256_set1_epi8(0x0E))); /* C|G|T */
        encoded = _mm256_or_si256(encoded, _mm256_and_si256(mask_N, _mm256_set1_epi8(0x0F))); /* A|C|G|T */
        
        /* Extract and pack encoded values */
        {
            uint8_t temp[32];
            int j;
        
        _mm256_storeu_si256((__m256i*)temp, encoded);
        
        /* Pack 4-bit values into output */
        for (j = 0; j < 32; j++) {
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 4) {
                output[byte_pos] |= (temp[j] << (4 - bit_offset));
            } else {
                /* bit_offset > 4: split across two bytes */
                int remaining_bits = 8 - bit_offset;
                output[byte_pos] |= (temp[j] >> (4 - remaining_bits));
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (temp[j] << (4 + remaining_bits));
                }
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
    /* Process 32 characters at a time with AVX2 SIMD optimizations */
    int simd_len = len & ~31;  /* Round down to multiple of 32 */
    
    /* Create lookup table for VPSHUFB - split into two halves for 4-bit lookups */
    const __m256i decode_lut_lo = _mm256_setr_epi8(
        '?', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N',
        '?', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
    );
    
    const __m256i mask_nibble = _mm256_set1_epi8(0x0F);
    
    for (int i = 0; i < simd_len; i += 32) {
        int byte_offset = (i * 4) / 8;
        uint8_t temp[17] __attribute__((aligned(32))) = {0}; /* Extra byte for boundary cases */
        
        /* Load required bytes */
        int bytes_to_load = 16 + 1; /* 32 * 4 bits = 128 bits = 16 bytes + 1 for boundary */
        for (int b = 0; b < bytes_to_load && (byte_offset + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[byte_offset + b];
        }
        
        /* Process in 16-character chunks for better efficiency */
        for (int chunk = 0; chunk < 2; chunk++) {
            __m128i nibbles = _mm_setzero_si128();
            uint8_t extracted[16];
            __m256i nibbles_256, masked, decoded;
            
            /* Extract 4-bit values for 16 characters */
            for (int j = 0; j < 16; j++) {
                int char_idx = chunk * 16 + j;
                int total_bit_pos = char_idx * 4;
                int byte_pos = total_bit_pos / 8;
                int bit_pos = total_bit_pos % 8;
                uint8_t val = 0;
                
                if (byte_pos < 17) {
                    if (bit_pos <= 4) {
                        val = (temp[byte_pos] >> (4 - bit_pos)) & 0x0F;
                    } else if (byte_pos + 1 < 17) {
                        int remaining = 8 - bit_pos;
                        val = ((temp[byte_pos] & ((1 << remaining) - 1)) << (4 - remaining)) |
                              ((temp[byte_pos + 1] >> (4 + remaining)) & 0x0F);
                    }
                }
                
                extracted[j] = val;
            }
            
            /* Load extracted nibbles */
            
            nibbles = _mm_loadu_si128((__m128i*)extracted);
            
            /* Use PSHUFB for lookup - process 16 at once */
            nibbles_256 = _mm256_zextsi128_si256(nibbles);
            nibbles_256 = _mm256_permute2x128_si256(nibbles_256, nibbles_256, 0);
            masked = _mm256_and_si256(nibbles_256, mask_nibble);
            decoded = _mm256_shuffle_epi8(decode_lut_lo, masked);
            
            /* Store 16 decoded characters */
            _mm_storeu_si128((__m128i*)&output[i + chunk * 16], _mm256_castsi256_si128(decoded));
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
__attribute__((target("avx512bw,avx512vbmi")))
void dna4_encode_avx512(const char* input, uint8_t* output, int len)
{
    int byte_len = (len * 4 + 7) / 8;
    int simd_len;
    /* Create comparison vectors for efficient encoding */
    const __m512i vec_A = _mm512_set1_epi8('A');
    const __m512i vec_C = _mm512_set1_epi8('C');
    const __m512i vec_G = _mm512_set1_epi8('G');
    const __m512i vec_T = _mm512_set1_epi8('T');
    const __m512i vec_a = _mm512_set1_epi8('a');
    const __m512i vec_c = _mm512_set1_epi8('c');
    const __m512i vec_g = _mm512_set1_epi8('g');
    const __m512i vec_t = _mm512_set1_epi8('t');
    
    memset(output, 0, byte_len);
    
    /* Process 64 characters at a time with AVX512 */
    simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    for (int i = 0; i < simd_len; i += 64) {
        __m512i chars = _mm512_loadu_si512((__m512i*)(input + i));
        
        /* Create masks for each base type using AVX512 mask registers */
        __mmask64 mask_A = _mm512_cmpeq_epi8_mask(chars, vec_A) | _mm512_cmpeq_epi8_mask(chars, vec_a);
        __mmask64 mask_C = _mm512_cmpeq_epi8_mask(chars, vec_C) | _mm512_cmpeq_epi8_mask(chars, vec_c);
        __mmask64 mask_G = _mm512_cmpeq_epi8_mask(chars, vec_G) | _mm512_cmpeq_epi8_mask(chars, vec_g);
        __mmask64 mask_T = _mm512_cmpeq_epi8_mask(chars, vec_T) | _mm512_cmpeq_epi8_mask(chars, vec_t);
        
        /* Use lookup table approach for DNA4 encoding */
        uint8_t temp[64];
        int j;
        
        _mm512_storeu_si512((__m512i*)temp, chars);
        
        for (j = 0; j < 64; j++) {
            temp[j] = kmersearch_dna4_encode_table[(unsigned char)temp[j]];
        }
        
        /* Pack 4-bit values into output */
        for (j = 0; j < 64; j++) {
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 4) {
                output[byte_pos] |= (temp[j] << (4 - bit_offset));
            } else {
                /* bit_offset > 4: handle byte boundary crossing */
                int remaining_bits = 8 - bit_offset;
                output[byte_pos] |= (temp[j] >> (4 - remaining_bits));
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (temp[j] << (4 + remaining_bits));
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
__attribute__((target("avx512bw,avx512vbmi")))
void dna4_decode_avx512(const uint8_t* input, char* output, int len)
{
    /* Process 64 characters at a time with AVX512 VBMI optimizations */
    int simd_len = len & ~63;  /* Round down to multiple of 64 */
    
    /* Create 64-byte lookup table for VPERMB */
    const __m512i decode_lut = _mm512_set_epi8(
        'N', 'B', 'D', 'K', 'H', 'Y', 'W', 'T',
        'V', 'S', 'R', 'G', 'M', 'C', 'A', '?',
        'N', 'B', 'D', 'K', 'H', 'Y', 'W', 'T',
        'V', 'S', 'R', 'G', 'M', 'C', 'A', '?',
        'N', 'B', 'D', 'K', 'H', 'Y', 'W', 'T',
        'V', 'S', 'R', 'G', 'M', 'C', 'A', '?',
        'N', 'B', 'D', 'K', 'H', 'Y', 'W', 'T',
        'V', 'S', 'R', 'G', 'M', 'C', 'A', '?'
    );
    
    const __m512i mask_nibble = _mm512_set1_epi8(0x0F);
    
    for (int i = 0; i < simd_len; i += 64) {
        int byte_offset = (i * 4) / 8;
        uint8_t temp[33] __attribute__((aligned(64))) = {0}; /* Extra byte for boundary cases */
        uint8_t nibbles[64] __attribute__((aligned(64)));
        __m512i nibbles_vec, masked, decoded;
        
        /* Load required bytes */
        int bytes_to_load = 32 + 1; /* 64 * 4 bits = 256 bits = 32 bytes + 1 for boundary */
        for (int b = 0; b < bytes_to_load && (byte_offset + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[byte_offset + b];
        }
        
        /* Process in 16-character chunks for extraction */
        
        for (int j = 0; j < 64; j++) {
            int total_bit_pos = j * 4;
            int byte_pos = total_bit_pos / 8;
            int bit_pos = total_bit_pos % 8;
            uint8_t val = 0;
            
            if (byte_pos < 33) {
                if (bit_pos <= 4) {
                    val = (temp[byte_pos] >> (4 - bit_pos)) & 0x0F;
                } else if (byte_pos + 1 < 33) {
                    int remaining = 8 - bit_pos;
                    val = ((temp[byte_pos] & ((1 << remaining) - 1)) << (4 - remaining)) |
                          ((temp[byte_pos + 1] >> (4 + remaining)) & 0x0F);
                }
            }
            
            nibbles[j] = val;
        }
        
        /* Load all nibbles into AVX512 register */
        nibbles_vec = _mm512_loadu_si512((__m512i*)nibbles);
        
        /* Use VPERMB for 64-byte lookup */
        masked = _mm512_and_si512(nibbles_vec, mask_nibble);
        decoded = _mm512_permutexvar_epi8(masked, decode_lut);
        
        /* Store 64 decoded characters */
        _mm512_storeu_si512((__m512i*)&output[i], decoded);
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
    
    memset(output, 0, byte_len);
    
    /* Create comparison vectors for efficient encoding */
    const uint8x16_t vec_A = vdupq_n_u8('A');
    const uint8x16_t vec_C = vdupq_n_u8('C');
    const uint8x16_t vec_G = vdupq_n_u8('G');
    const uint8x16_t vec_T = vdupq_n_u8('T');
    const uint8x16_t vec_a = vdupq_n_u8('a');
    const uint8x16_t vec_c = vdupq_n_u8('c');
    const uint8x16_t vec_g = vdupq_n_u8('g');
    const uint8x16_t vec_t = vdupq_n_u8('t');
    
    /* Process 16 characters at a time with NEON */
    for (i = 0; i < simd_len; i += 16) {
        uint8x16_t chars = vld1q_u8((uint8_t*)(input + i));
        
        /* Create masks for each base type */
        uint8x16_t mask_A = vorrq_u8(vceqq_u8(chars, vec_A), vceqq_u8(chars, vec_a));
        uint8x16_t mask_C = vorrq_u8(vceqq_u8(chars, vec_C), vceqq_u8(chars, vec_c));
        uint8x16_t mask_G = vorrq_u8(vceqq_u8(chars, vec_G), vceqq_u8(chars, vec_g));
        uint8x16_t mask_T = vorrq_u8(vceqq_u8(chars, vec_T), vceqq_u8(chars, vec_t));
        
        /* Generate 4-bit encoded values using SIMD */
        uint8x16_t encoded = vdupq_n_u8(0);
        encoded = vorrq_u8(encoded, vandq_u8(mask_A, vdupq_n_u8(0x01)));
        encoded = vorrq_u8(encoded, vandq_u8(mask_C, vdupq_n_u8(0x02)));
        encoded = vorrq_u8(encoded, vandq_u8(mask_G, vdupq_n_u8(0x04)));
        encoded = vorrq_u8(encoded, vandq_u8(mask_T, vdupq_n_u8(0x08)));
        
        /* Extract encoded values and handle degenerate bases */
        uint8_t temp[16];
        vst1q_u8(temp, encoded);
        
        /* Check for degenerate bases and update encoding */
        uint8_t chars_temp[16];
        vst1q_u8(chars_temp, chars);
        
        for (j = 0; j < 16; j++) {
            uint8_t ch = chars_temp[j];
            /* Only update if not already encoded (i.e., degenerate base) */
            if (temp[j] == 0) {
                temp[j] = kmersearch_dna4_encode_table[(unsigned char)ch];
            }
        }
        
        /* Pack 4-bit values into output */
        for (j = 0; j < 16; j++) {
            int bit_pos = (i + j) * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            if (bit_offset <= 4) {
                output[byte_pos] |= (temp[j] << (4 - bit_offset));
            } else {
                /* bit_offset > 4: split across two bytes */
                int remaining_bits = 8 - bit_offset;
                output[byte_pos] |= (temp[j] >> (4 - remaining_bits));
                if (byte_pos + 1 < byte_len) {
                    output[byte_pos + 1] |= (temp[j] << (4 + remaining_bits));
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
    /* Process 16 characters at a time with NEON VTBL optimizations */
    int simd_len = len & ~15;  /* Round down to multiple of 16 */
    
    /* Create lookup table for VTBL */
    const uint8x16_t decode_lut = vld1q_u8((const uint8_t*)"?ACMGRSVTWYHDKBN");
    const uint8x16_t mask_nibble = vdupq_n_u8(0x0F);
    
    for (int i = 0; i < simd_len; i += 16) {
        int byte_offset = (i * 4) / 8;
        uint8_t temp[9] __attribute__((aligned(16))) = {0}; /* Extra byte for boundaries */
        
        /* Load required bytes */
        int bytes_to_load = 8 + 1; /* 16 * 4 bits = 64 bits = 8 bytes + 1 for boundary */
        for (int b = 0; b < bytes_to_load && (byte_offset + b) < (len * 4 + 7) / 8; b++) {
            temp[b] = input[byte_offset + b];
        }
        
        /* Extract 4-bit values for 16 characters */
        uint8_t nibbles[16] __attribute__((aligned(16)));
        
        for (int j = 0; j < 16; j++) {
            int bit_pos = j * 4;
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            uint8_t val = 0;
            
            if (byte_pos < 9) {
                if (bit_offset <= 4) {
                    val = (temp[byte_pos] >> (4 - bit_offset)) & 0x0F;
                } else if (byte_pos + 1 < 9) {
                    int remaining = 8 - bit_offset;
                    val = ((temp[byte_pos] & ((1 << remaining) - 1)) << (4 - remaining)) |
                          ((temp[byte_pos + 1] >> (4 + remaining)) & 0x0F);
                }
            }
            
            nibbles[j] = val;
        }
        
        /* Load nibbles into NEON register */
        uint8x16_t nibbles_vec = vld1q_u8(nibbles);
        
        /* Use VTBL for lookup */
        uint8x16_t masked = vandq_u8(nibbles_vec, mask_nibble);
        uint8x16_t decoded = vqtbl1q_u8(decode_lut, masked);
        
        /* Store result */
        vst1q_u8((uint8_t*)&output[i], decoded);
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
