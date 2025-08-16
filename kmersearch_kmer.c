#include "kmersearch.h"

/*
 * kmersearch_kmer.c
 * 
 * Simple k-mer utility functions for the pg_kmersearch extension.
 * 
 * This module contains basic utility functions for:
 * - k-mer bit operations
 * - Simple data conversion
 * - Basic helper functions
 * 
 * Note: Complex memory management functions remain in kmersearch.c for stability
 */

/* Helper function to validate DNA4 characters */
static inline bool
kmersearch_is_valid_dna4_char(char c)
{
    return kmersearch_dna4_encode_table[(unsigned char)c] != 0 || c == 'A' || c == 'a';
}

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

/* Forward declarations for static functions */
static int kmersearch_count_matching_kmer_fast_scalar_simple(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
static int kmersearch_count_matching_kmer_fast_scalar_hashtable(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
static int kmersearch_count_matching_uintkey_scalar(void *seq_keys, int seq_nkeys, void *query_keys, int query_nkeys, int k_size);

/* SIMD-optimized k-mer matching functions */
#ifdef __x86_64__
static int kmersearch_count_matching_kmer_fast_avx2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
static int kmersearch_count_matching_kmer_fast_avx512(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
#endif
#ifdef __aarch64__
static int kmersearch_count_matching_kmer_fast_neon(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
static int kmersearch_count_matching_kmer_fast_sve(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys);
#endif

/*
 * Simple utility function: Set a specific bit in a bit array
 */
void
kmersearch_set_bit_at(bits8 *data, int bit_pos, int value)
{
    int byte_pos = bit_pos / 8;
    int bit_offset = bit_pos % 8;
    
    if (value) {
        data[byte_pos] |= (1 << (7 - bit_offset));
    } else {
        data[byte_pos] &= ~(1 << (7 - bit_offset));
    }
}

/*
 * Simple utility function: Count degenerate combinations in a k-mer string
 */
int
kmersearch_count_degenerate_combinations(const char *kmer, int k)
{
    int total_combinations = 1;
    int i;
    
    for (i = 0; i < k; i++) {
        char base = kmer[i];
        int count = 1;
        
        switch (base) {
            case 'M': case 'm': /* A,C */
            case 'R': case 'r': /* A,G */
            case 'W': case 'w': /* A,T */
            case 'S': case 's': /* C,G */
            case 'Y': case 'y': /* C,T */
            case 'K': case 'k': /* G,T */
                count = 2;
                break;
            case 'V': case 'v': /* A,C,G */
            case 'H': case 'h': /* A,C,T */
            case 'D': case 'd': /* A,G,T */
            case 'B': case 'b': /* C,G,T */
                count = 3;
                break;
            case 'N': case 'n': /* A,C,G,T */
                count = 4;
                break;
            default:
                count = 1; /* Standard bases A,C,G,T */
                break;
        }
        
        total_combinations *= count;
        
        /* Prevent overflow - limit to reasonable number */
        if (total_combinations > 1000) {
            return 1000;
        }
    }
    
    return total_combinations;
}

/*
 * Simple utility function: Check if DNA4 k-mer will exceed degenerate limit
 * 
 * Rules:
 * - If only 1 MRWSYKVHDB character: expand (return false)
 * - If 2 or more MRWSYKVHDB characters: don't expand (return true)
 * - If 1 or more N or ? characters: don't expand (return true)
 */
bool
kmersearch_will_exceed_degenerate_limit_dna4_bits(VarBit *seq, int start_pos, int k)
{
    int degenerate_count = 0;  /* Count of MRWSYKVHDB (expansion 2 or 3) */
    bits8 *data = VARBITS(seq);
    int i;
    
    for (i = 0; i < k; i++)
    {
        int bit_pos = (start_pos + i) * 4;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 encoded;
        int expansion_count;
        
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
        expansion_count = kmersearch_dna4_to_dna2_table[encoded][0];
        
        /* N or ? (expansion count 4 or 0/invalid) - always exceed limit */
        if (expansion_count == 4 || expansion_count == 0)
            return true;
        
        /* MRWSYKVHDB (expansion count 2 or 3) */
        if (expansion_count == 2 || expansion_count == 3)
        {
            degenerate_count++;
            /* 2 or more MRWSYKVHDB characters - exceed limit */
            if (degenerate_count >= 2)
                return true;
        }
    }
    
    /* Only 0 or 1 MRWSYKVHDB character - within limit */
    return false;
}


/*
 * Create k-mer key without occurrence count (for frequency analysis)
 */
/*
 * Create k-mer key from DNA2 bits without occurrence count
 */
/*
 * Create n-gram key from DNA2 bits with occurrence count
 */
/*
 * Create n-gram key from DNA4 bits with occurrence count (by converting to DNA2 first)
 */
/*
 * Create n-gram key from existing DNA2 k-mer with occurrence count
 */
/*
 * Find or add k-mer occurrence in sorted array for k <= 8
 */
int
kmersearch_find_or_add_kmer_occurrence16(KmerOccurrence16 *occurrences, int *count, uint16 kmer_value, int max_count)
{
    int left = 0, right = *count - 1;
    int insert_pos = *count;
    
    /* Binary search for existing k-mer */
    while (left <= right)
    {
        int mid = (left + right) / 2;
        if (occurrences[mid].kmer_value == kmer_value)
        {
            return ++occurrences[mid].count;  /* Found, increment and return */
        }
        else if (occurrences[mid].kmer_value < kmer_value)
        {
            left = mid + 1;
        }
        else
        {
            right = mid - 1;
            insert_pos = mid;
        }
    }
    
    /* Not found, insert new entry if space available */
    if (*count >= max_count)
        return -1;  /* Array full */
    
    /* Shift elements to make room */
    if (insert_pos < *count)
    {
        memmove(&occurrences[insert_pos + 1], &occurrences[insert_pos],
                (*count - insert_pos) * sizeof(KmerOccurrence16));
    }
    
    /* Insert new entry */
    occurrences[insert_pos].kmer_value = kmer_value;
    occurrences[insert_pos].count = 1;
    (*count)++;
    
    return 1;
}

/*
 * Find or add k-mer occurrence in sorted array for k <= 16
 */
int
kmersearch_find_or_add_kmer_occurrence32(KmerOccurrence32 *occurrences, int *count, uint32 kmer_value, int max_count)
{
    int left = 0, right = *count - 1;
    int insert_pos = *count;
    
    /* Binary search for existing k-mer */
    while (left <= right)
    {
        int mid = (left + right) / 2;
        if (occurrences[mid].kmer_value == kmer_value)
        {
            return ++occurrences[mid].count;  /* Found, increment and return */
        }
        else if (occurrences[mid].kmer_value < kmer_value)
        {
            left = mid + 1;
        }
        else
        {
            right = mid - 1;
            insert_pos = mid;
        }
    }
    
    /* Not found, insert new entry if space available */
    if (*count >= max_count)
        return -1;  /* Array full */
    
    /* Shift elements to make room */
    if (insert_pos < *count)
    {
        memmove(&occurrences[insert_pos + 1], &occurrences[insert_pos],
                (*count - insert_pos) * sizeof(KmerOccurrence32));
    }
    
    /* Insert new entry */
    occurrences[insert_pos].kmer_value = kmer_value;
    occurrences[insert_pos].count = 1;
    (*count)++;
    
    return 1;
}

/*
 * Find or add k-mer occurrence in sorted array for k <= 32
 */
int
kmersearch_find_or_add_kmer_occurrence64(KmerOccurrence64 *occurrences, int *count, uint64 kmer_value, int max_count)
{
    int left = 0, right = *count - 1;
    int insert_pos = *count;
    
    /* Binary search for existing k-mer */
    while (left <= right)
    {
        int mid = (left + right) / 2;
        if (occurrences[mid].kmer_value == kmer_value)
        {
            return ++occurrences[mid].count;  /* Found, increment and return */
        }
        else if (occurrences[mid].kmer_value < kmer_value)
        {
            left = mid + 1;
        }
        else
        {
            right = mid - 1;
            insert_pos = mid;
        }
    }
    
    /* Not found, insert new entry if space available */
    if (*count >= max_count)
        return -1;  /* Array full */
    
    /* Shift elements to make room */
    if (insert_pos < *count)
    {
        memmove(&occurrences[insert_pos + 1], &occurrences[insert_pos],
                (*count - insert_pos) * sizeof(KmerOccurrence64));
    }
    
    /* Insert new entry */
    occurrences[insert_pos].kmer_value = kmer_value;
    occurrences[insert_pos].count = 1;
    (*count)++;
    
    return 1;
}


/*
 * Extract k-mers from query string as ngram_key2 format (with occurrence counts)
 * This function generates ngram_key2 compatible with sequence extraction using SIMD optimizations
 */
/*
 * Helper function to get bit at position from bit array
 */
uint8
kmersearch_get_bit_at(bits8 *data, int bit_pos)
{
    int byte_pos = bit_pos / 8;
    int bit_offset = bit_pos % 8;
    return (data[byte_pos] >> (7 - bit_offset)) & 1;
}

/*
 * Fast check if degenerate combinations exceed limit (11+)
 * Returns true if combinations will exceed 10
 */
bool
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
 * Removed: kmersearch_remove_occurrence_bits() function 
 * This function was eliminated because the correct architecture uses ngram_key2
 * directly without removing occurrence bits for high-frequency k-mer comparisons.
 */




/* Note: kmersearch_extract_dna2_ngram_key2_direct kept in kmersearch.c due to SIMD dependencies */




/*
 * Convert uint-based k-mer to VarBit k-mer format
 */
/*
 * Create complete ngram_key2 from uint k-mer with occurrence count
 */
/*
 * Create an ngram_key2 from a kmer2 uint64 value and occurrence count
 * This function is used during frequency analysis to create ngram keys
 * from already extracted k-mers
 */
/* Forward declarations for SIMD variants */
/*
 * Scalar version: Extract k-mers directly from DNA2 bit sequence
 */
#ifdef __x86_64__
/* AVX2 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("avx2,bmi,bmi2")))
/* AVX512 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2,bmi,bmi2")))
#endif /* __x86_64__ */

#ifdef __aarch64__
/* NEON optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+simd")))
/* SVE optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+sve,+simd")))
/* SVE2 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+sve2")))
#endif /* __aarch64__ */

/* Forward declaration for DNA4 expansion functions */
/*
 * Scalar version: Extract k-mers directly from DNA4 bit sequence with degenerate expansion
 */
/*
 * Extract k-mers directly from DNA2 bit sequence (with SIMD dispatch)
 */
#ifdef __x86_64__
/* AVX2 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("avx2,bmi,bmi2")))
/* AVX512 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2,bmi,bmi2")))
#endif /* __x86_64__ */

#ifdef __aarch64__
/* NEON optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+simd")))
/* SVE optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+sve,+simd")))
/* SVE2 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+sve2")))
#endif /* __aarch64__ */

/*
 * Extract k-mers directly from DNA4 bit sequence with degenerate expansion (with SIMD dispatch)
 */
/*
 * SIMD-optimized implementations for k-mer matching
 */

#ifdef __x86_64__
/*
 * AVX2 optimized version of k-mer matching using hash table
 */
__attribute__((target("avx2,bmi,bmi2")))
static int
kmersearch_count_matching_kmer_fast_avx2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i, j;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    int batch_size = 8;  /* Process 8 keys at a time for cache efficiency */
    int batch_end;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* AVX2 optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    /* Safety check: ensure we have valid query keys */
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);  /* Use data size, not total size */
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Creating hash with keysize=%zu, query_nkeys=%d", 
         (size_t)VARBITBYTES(query_keys[0]), query_nkeys);
    
    query_hash = hash_create("QueryKmerHashAVX2", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* Insert all query k-mers into hash table using content as key */
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: NULL query key at index %d", i);
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Process sequence keys in batches for better cache performance */
    for (i = 0; i < seq_nkeys; i += batch_size)
    {
        batch_end = (i + batch_size < seq_nkeys) ? i + batch_size : seq_nkeys;
        
        /* Prefetch batch data */
        for (j = i; j < batch_end && j < seq_nkeys; j++) {
            if (seq_keys[j] != NULL) {
                _mm_prefetch((const char*)VARBITS(seq_keys[j]), _MM_HINT_T0);
            }
        }
        
        /* Process batch */
        for (j = i; j < batch_end; j++)
        {
            if (seq_keys[j] == NULL) {
                elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: NULL seq key at index %d", j);
                continue;
            }
            
            if (VARBITBYTES(seq_keys[j]) != VARBITBYTES(query_keys[0])) {
                elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Size mismatch seq[%d]=%zu vs query[0]=%zu", 
                     j, (size_t)VARBITBYTES(seq_keys[j]), (size_t)VARBITBYTES(query_keys[0]));
                continue;
            }
            
            if (hash_search(query_hash, VARBITS(seq_keys[j]), HASH_FIND, NULL))
            {
                match_count++;
            }
        }
    }
    
    hash_destroy(query_hash);
    
    elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Returning match_count=%d", match_count);
    return match_count;
}

/*
 * AVX512 optimized version of k-mer matching using hash table
 */
__attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2,bmi,bmi2")))
static int
kmersearch_count_matching_kmer_fast_avx512(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i, j;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    int batch_size = 16;  /* Process 16 keys at a time for cache efficiency */
    int batch_end;
    
    if (seq_nkeys == 0 || query_nkeys == 0) {
        return 0;
    }
    
    /* AVX512 optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        return 0;
    }
    
    /* Validate the VarBit structure */
    if (VARSIZE(query_keys[0]) < VARHDRSZ + 4) {
        elog(ERROR, "Invalid VarBit size: %zu", (size_t)VARSIZE(query_keys[0]));
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    query_hash = hash_create("QueryKmerHashAVX512", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Process sequence keys in larger batches for AVX512 */
    for (i = 0; i < seq_nkeys; i += batch_size)
    {
        batch_end = (i + batch_size < seq_nkeys) ? i + batch_size : seq_nkeys;
        
        /* Prefetch batch data with AVX512 prefetch instructions */
        for (j = i; j < batch_end && j < seq_nkeys; j++) {
            if (seq_keys[j] != NULL) {
                _mm_prefetch((const char*)VARBITS(seq_keys[j]), _MM_HINT_T0);
                /* Prefetch next cache line if key is large */
                if (VARBITBYTES(seq_keys[j]) > 64) {
                    _mm_prefetch((const char*)VARBITS(seq_keys[j]) + 64, _MM_HINT_T0);
                }
            }
        }
        
        /* Process batch */
        for (j = i; j < batch_end; j++)
        {
            if (seq_keys[j] == NULL) {
                continue;
            }
            
            if (VARBITBYTES(seq_keys[j]) != VARBITBYTES(query_keys[0])) {
                continue;
            }
            
            if (hash_search(query_hash, VARBITS(seq_keys[j]), HASH_FIND, NULL))
            {
                match_count++;
            }
        }
    }
    
    hash_destroy(query_hash);
    
    return match_count;
}
#endif /* __x86_64__ */

#ifdef __aarch64__
/*
 * NEON optimized version of k-mer matching using hash table
 */
__attribute__((target("+simd")))
static int
kmersearch_count_matching_kmer_fast_neon(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i, j;
    int batch_size = 4;  /* Process 4 keys at a time for cache efficiency */
    int batch_end;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* NEON optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmer_fast_neon: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHashNEON", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Process sequence keys with prefetching for better cache performance */
    for (i = 0; i < seq_nkeys; i += batch_size)
    {
        batch_end = (i + batch_size < seq_nkeys) ? i + batch_size : seq_nkeys;
        
        /* Prefetch batch data */
        for (j = i; j < batch_end && j < seq_nkeys; j++) {
            if (seq_keys[j] != NULL) {
                __builtin_prefetch(VARBITS(seq_keys[j]), 0, 1);
            }
        }
        
        /* Process batch */
        for (j = i; j < batch_end; j++)
        {
            if (seq_keys[j] == NULL) {
                continue;
            }
            
            if (VARBITBYTES(seq_keys[j]) != VARBITBYTES(query_keys[0])) {
                continue;
            }
            
            if (hash_search(query_hash, VARBITS(seq_keys[j]), HASH_FIND, NULL))
            {
                match_count++;
            }
        }
    }
    
    hash_destroy(query_hash);
    
    return match_count;
}

/*
 * SVE optimized version of k-mer matching using hash table
 */
__attribute__((target("+sve,+simd")))
static int
kmersearch_count_matching_kmer_fast_sve(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* SVE optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmer_fast_sve: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHashSVE", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    return match_count;
}

/*
 * SVE2 optimized version of k-mer matching using hash table
 */
__attribute__((target("+sve2")))
static int
kmersearch_count_matching_kmer_fast_sve2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* SVE2 optimized hash table implementation */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (query_keys[0] == NULL) {
        elog(LOG, "kmersearch_count_matching_kmer_fast_sve2: NULL query key detected");
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHashSVE2", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* SVE2 optimized insertion with vectorized memory operations */
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* SVE2 optimized search with advanced comparison instructions */
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    return match_count;
}
#endif /* __aarch64__ */

/*
 * Count matching k-mers using simple O(n*m) algorithm
 * Suitable for small data sets
 */
static int
kmersearch_count_matching_kmer_fast_scalar_simple(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i, j;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;

    /* Simple O(n*m) algorithm - good for small data sets */
    for (i = 0; i < seq_nkeys; i++)
    {
        for (j = 0; j < query_nkeys; j++)
        {
            if (VARBITLEN(seq_keys[i]) == VARBITLEN(query_keys[j]) &&
                VARSIZE(seq_keys[i]) == VARSIZE(query_keys[j]) &&
                memcmp(VARBITS(seq_keys[i]), VARBITS(query_keys[j]), VARBITBYTES(seq_keys[i])) == 0)
            {
                match_count++;
                break; /* Each seq key matches at most once */
            }
        }
    }

    return match_count;
}

/*
 * Count matching k-mers using hash table for O(n+m) performance
 * Suitable for large data sets
 */
static int
kmersearch_count_matching_kmer_fast_scalar_hashtable(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* Create hash table using VarBit content as key */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    /* Safety check: ensure we have valid query keys */
    if (query_keys[0] == NULL) {
        return 0;
    }
    
    hash_ctl.keysize = VARBITBYTES(query_keys[0]);  /* Use data size, not total size */
    hash_ctl.entrysize = sizeof(bool);
    hash_ctl.hash = tag_hash;
    
    query_hash = hash_create("QueryKmerHash", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* Insert all query k-mers into hash table using content as key */
    for (i = 0; i < query_nkeys; i++)
    {
        if (query_keys[i] == NULL) {
            continue;
        }
        hash_search(query_hash, VARBITS(query_keys[i]), HASH_ENTER, &found);
    }
    
    /* Check each sequence k-mer against hash table */
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    /* Cleanup */
    hash_destroy(query_hash);
    
    return match_count;
}

/*
 * Remove occurrence count from ngram_key2 to get k-mer only
 */
/*
 * Convert VarBit k-mer to uint16 (for k <= 8)
 */
/*
 * Convert VarBit k-mer to uint32 (for k <= 16)
 */
/*
 * Convert VarBit k-mer to uint64 (for k <= 32)
 */
/*
 * Extract uint keys with occurrence counting from DNA2 sequence (scalar implementation)
 */
static void
kmersearch_extract_uintkey_from_dna2_scalar(VarBit *seq, void **output, int *nkeys)
{
    int k = kmersearch_kmer_size;
    int occur_bitlen = kmersearch_occur_bitlen;
    int seq_bits = VARBITLEN(seq);
    int seq_len = seq_bits / 2;
    int max_kmers = seq_len - k + 1;
    int kmer_bits = k * 2;
    int total_bits = kmer_bits + occur_bitlen;
    size_t elem_size;
    void *result;
    int result_count = 0;
    int i;
    
    /* Validate parameters */
    if (k < 4 || k > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid k-mer size: %d", k),
                 errdetail("k-mer size must be between 4 and 32")));
    }
    
    if (occur_bitlen < 1 || occur_bitlen > 16) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid occurrence bit length: %d", occur_bitlen),
                 errdetail("Occurrence bit length must be between 1 and 16")));
    }
    
    if (total_bits > 64) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Total bit length exceeds 64 bits"),
                 errdetail("k-mer bits: %d, occurrence bits: %d, total: %d", 
                          kmer_bits, occur_bitlen, total_bits)));
    }
    
    if (max_kmers <= 0) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Query sequence must be at least %d bases long", k)));
        return;
    }
    
    /* Determine element size based on total bits */
    if (total_bits <= 16) {
        elem_size = sizeof(uint16);
    } else if (total_bits <= 32) {
        elem_size = sizeof(uint32);
    } else {
        elem_size = sizeof(uint64);
    }
    
    /* Allocate result array */
    result = palloc(max_kmers * elem_size);
    
    /* Allocate occurrence tracking array based on element size */
    void *occurrences = NULL;
    int occurrence_count = 0;
    
    if (elem_size == sizeof(uint16)) {
        occurrences = palloc0(max_kmers * sizeof(KmerOccurrence16));
    } else if (elem_size == sizeof(uint32)) {
        occurrences = palloc0(max_kmers * sizeof(KmerOccurrence32));
    } else {
        occurrences = palloc0(max_kmers * sizeof(KmerOccurrence64));
    }
    
    /* Extract k-mers and track occurrences */
    for (i = 0; i <= seq_len - k; i++) {
        int j;
        int current_count;
        bits8 *data = VARBITS(seq);
        
        /* Extract k-mer bits and store based on element size */
        if (elem_size == sizeof(uint16)) {
            uint16 kmer_value = 0;
            uint16 final_value;
            
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 nucleotide = (data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | nucleotide;
            }
            
            current_count = kmersearch_find_or_add_kmer_occurrence16(
                (KmerOccurrence16 *)occurrences, &occurrence_count,
                kmer_value, max_kmers);
            if (current_count < 0) continue;
            if (current_count > (1 << occur_bitlen)) continue;
            
            final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
            ((uint16 *)result)[result_count++] = final_value;
            
        } else if (elem_size == sizeof(uint32)) {
            uint32 kmer_value = 0;
            uint32 final_value;
            
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 nucleotide = (data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | nucleotide;
            }
            
            current_count = kmersearch_find_or_add_kmer_occurrence32(
                (KmerOccurrence32 *)occurrences, &occurrence_count,
                kmer_value, max_kmers);
            if (current_count < 0) continue;
            if (current_count > (1 << occur_bitlen)) continue;
            
            final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
            ((uint32 *)result)[result_count++] = final_value;
            
        } else {
            uint64 kmer_value = 0;
            uint64 final_value;
            
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 nucleotide = (data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | nucleotide;
            }
            
            current_count = kmersearch_find_or_add_kmer_occurrence64(
                (KmerOccurrence64 *)occurrences, &occurrence_count,
                kmer_value, max_kmers);
            if (current_count < 0) continue;
            if (current_count > (1 << occur_bitlen)) continue;
            
            final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
            ((uint64 *)result)[result_count++] = final_value;
        }
    }
    
    /* Free occurrence tracking array */
    pfree(occurrences);
    
    /* Reallocate to actual size if needed */
    if (result_count < max_kmers) {
        void *new_result = palloc(result_count * elem_size);
        memcpy(new_result, result, result_count * elem_size);
        pfree(result);
        result = new_result;
    }
    
    *output = result;
    *nkeys = result_count;
}

/*
 * Extract uint keys with occurrence counting from DNA2 sequence (dispatch function)
 */
void
kmersearch_extract_uintkey_from_dna2(VarBit *seq, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
    /* For now, always use scalar implementation */
    /* SIMD implementations can be added later based on capability and threshold */
#ifdef __x86_64__
    /* Future: Add AVX2/AVX512 dispatch based on simd_capability */
#elif defined(__aarch64__)
    /* Future: Add NEON/SVE dispatch based on simd_capability */
#endif
    
    kmersearch_extract_uintkey_from_dna2_scalar(seq, output, nkeys);
}

/*
 * Extract uint keys with occurrence counting from DNA4 sequence (scalar implementation)
 */
static void
kmersearch_extract_uintkey_from_dna4_scalar(VarBit *seq, void **output, int *nkeys)
{
    int k = kmersearch_kmer_size;
    int occur_bitlen = kmersearch_occur_bitlen;
    int seq_bits = VARBITLEN(seq);
    int seq_len = seq_bits / 4;  /* DNA4 uses 4 bits per character */
    int max_kmers = seq_len - k + 1;
    int kmer_bits = k * 2;  /* Output is DNA2 format */
    int total_bits = kmer_bits + occur_bitlen;
    size_t elem_size;
    void *result;
    int result_count = 0;
    int result_capacity;
    int i;
    
    /* Validate parameters */
    if (k < 4 || k > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid k-mer size: %d", k),
                 errdetail("k-mer size must be between 4 and 32")));
    }
    
    if (occur_bitlen < 1 || occur_bitlen > 16) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid occurrence bit length: %d", occur_bitlen),
                 errdetail("Occurrence bit length must be between 1 and 16")));
    }
    
    if (total_bits > 64) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Total bit length exceeds 64 bits"),
                 errdetail("k-mer bits: %d, occurrence bits: %d, total: %d", 
                          kmer_bits, occur_bitlen, total_bits)));
    }
    
    if (max_kmers <= 0) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Query sequence must be at least %d bases long", k)));
        return;
    }
    
    /* Determine element size based on total bits */
    if (total_bits <= 16) {
        elem_size = sizeof(uint16);
    } else if (total_bits <= 32) {
        elem_size = sizeof(uint32);
    } else {
        elem_size = sizeof(uint64);
    }
    
    /* Allocate result array with extra space for degenerate expansions */
    result_capacity = max_kmers * 16;  /* Conservative estimate for degenerate bases */
    result = palloc(result_capacity * elem_size);
    
    /* Allocate occurrence tracking array based on element size */
    void *occurrences = NULL;
    int occurrence_count = 0;
    
    if (elem_size == sizeof(uint16)) {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence16));
    } else if (elem_size == sizeof(uint32)) {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence32));
    } else {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence64));
    }
    
    /* Process each k-mer position */
    for (i = 0; i <= seq_len - k; i++) {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0) {
            continue;
        }
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++) {
            VarBit *dna2_kmer = expanded_kmers[j];
            int current_count;
            
            if (!dna2_kmer) {
                continue;
            }
            
            /* Convert VarBit to appropriate uint type and store */
            if (elem_size == sizeof(uint16)) {
                uint16 kmer_value = 0;
                uint16 final_value;
                bits8 *data = VARBITS(dna2_kmer);
                int nbits = VARBITLEN(dna2_kmer);
                int nbytes = (nbits + 7) / 8;
                int bi;
                
                for (bi = 0; bi < nbytes && bi < 2; bi++) {
                    kmer_value = (kmer_value << 8) | data[bi];
                }
                
                if (nbits % 8 != 0) {
                    int shift = 8 - (nbits % 8);
                    kmer_value >>= shift;
                }
                
                current_count = kmersearch_find_or_add_kmer_occurrence16(
                    (KmerOccurrence16 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    pfree(dna2_kmer);
                    continue;
                }
                
                if (current_count > (1 << occur_bitlen)) {
                    pfree(dna2_kmer);
                    continue;
                }
                
                final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                ((uint16 *)result)[result_count++] = final_value;
                
            } else if (elem_size == sizeof(uint32)) {
                uint32 kmer_value = 0;
                uint32 final_value;
                bits8 *data = VARBITS(dna2_kmer);
                int nbits = VARBITLEN(dna2_kmer);
                int nbytes = (nbits + 7) / 8;
                int bi;
                
                for (bi = 0; bi < nbytes && bi < 4; bi++) {
                    kmer_value = (kmer_value << 8) | data[bi];
                }
                
                if (nbits % 8 != 0) {
                    int shift = 8 - (nbits % 8);
                    kmer_value >>= shift;
                }
                
                current_count = kmersearch_find_or_add_kmer_occurrence32(
                    (KmerOccurrence32 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    pfree(dna2_kmer);
                    continue;
                }
                
                if (current_count > (1 << occur_bitlen)) {
                    pfree(dna2_kmer);
                    continue;
                }
                
                final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                ((uint32 *)result)[result_count++] = final_value;
                
            } else {
                uint64 kmer_value = 0;
                uint64 final_value;
                bits8 *data = VARBITS(dna2_kmer);
                int nbits = VARBITLEN(dna2_kmer);
                int nbytes = (nbits + 7) / 8;
                int bi;
                
                for (bi = 0; bi < nbytes && bi < 8; bi++) {
                    kmer_value = (kmer_value << 8) | data[bi];
                }
                
                if (nbits % 8 != 0) {
                    int shift = 8 - (nbits % 8);
                    kmer_value >>= shift;
                }
                
                current_count = kmersearch_find_or_add_kmer_occurrence64(
                    (KmerOccurrence64 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    pfree(dna2_kmer);
                    continue;
                }
                
                if (current_count > (1 << occur_bitlen)) {
                    pfree(dna2_kmer);
                    continue;
                }
                
                final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                ((uint64 *)result)[result_count++] = final_value;
            }
            
            /* Free the expanded k-mer */
            pfree(dna2_kmer);
        }
        
        /* Free the expansion array */
        if (expanded_kmers) {
            pfree(expanded_kmers);
        }
    }
    
    /* Free occurrence tracking array */
    pfree(occurrences);
    
    /* Reallocate to actual size if needed */
    if (result_count < result_capacity) {
        void *new_result = palloc(result_count * elem_size);
        memcpy(new_result, result, result_count * elem_size);
        pfree(result);
        result = new_result;
    }
    
    *output = result;
    *nkeys = result_count;
}

/*
 * Extract uint keys with occurrence counting from DNA4 sequence (dispatch function)
 */
void
kmersearch_extract_uintkey_from_dna4(VarBit *seq, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
    /* For now, always use scalar implementation */
    /* SIMD implementations can be added later based on capability and threshold */
#ifdef __x86_64__
    /* Future: Add AVX2/AVX512 dispatch based on simd_capability */
#elif defined(__aarch64__)
    /* Future: Add NEON/SVE dispatch based on simd_capability */
#endif
    
    kmersearch_extract_uintkey_from_dna4_scalar(seq, output, nkeys);
}

/*
 * Count matching uintkeys - scalar implementation with hash table
 */
static int
kmersearch_count_matching_uintkey_scalar(void *seq_keys, int seq_nkeys, void *query_keys, int query_nkeys, int k_size)
{
    int shared_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    /* Use hash table for O(n+m) performance */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (k_size <= 8) {
        hash_ctl.keysize = sizeof(uint16);
        hash_ctl.entrysize = sizeof(uint16);
    }
    else if (k_size <= 16) {
        hash_ctl.keysize = sizeof(uint32);
        hash_ctl.entrysize = sizeof(uint32);
    }
    else {
        hash_ctl.keysize = sizeof(uint64);
        hash_ctl.entrysize = sizeof(uint64);
    }
    
    hash_ctl.hash = tag_hash;
    query_hash = hash_create("QueryUintkeyHash", query_nkeys * 2, &hash_ctl,
                            HASH_ELEM | HASH_FUNCTION | HASH_BLOBS);
    
    /* Insert all query uintkeys into hash table */
    if (k_size <= 8) {
        uint16 *query = (uint16 *)query_keys;
        for (i = 0; i < query_nkeys; i++) {
            hash_search(query_hash, &query[i], HASH_ENTER, &found);
        }
        
        /* Check each sequence uintkey against hash table */
        uint16 *seq = (uint16 *)seq_keys;
        for (i = 0; i < seq_nkeys; i++) {
            if (hash_search(query_hash, &seq[i], HASH_FIND, NULL)) {
                shared_count++;
            }
        }
    }
    else if (k_size <= 16) {
        uint32 *query = (uint32 *)query_keys;
        for (i = 0; i < query_nkeys; i++) {
            hash_search(query_hash, &query[i], HASH_ENTER, &found);
        }
        
        uint32 *seq = (uint32 *)seq_keys;
        for (i = 0; i < seq_nkeys; i++) {
            if (hash_search(query_hash, &seq[i], HASH_FIND, NULL)) {
                shared_count++;
            }
        }
    }
    else {
        uint64 *query = (uint64 *)query_keys;
        for (i = 0; i < query_nkeys; i++) {
            hash_search(query_hash, &query[i], HASH_ENTER, &found);
        }
        
        uint64 *seq = (uint64 *)seq_keys;
        for (i = 0; i < seq_nkeys; i++) {
            if (hash_search(query_hash, &seq[i], HASH_FIND, NULL)) {
                shared_count++;
            }
        }
    }
    
    /* Clean up hash table */
    hash_destroy(query_hash);
    
    return shared_count;
}

/*
 * Count matching uintkeys - dispatch function
 * This will dispatch to SIMD implementations in the future
 */
int
kmersearch_count_matching_uintkey(void *seq_keys, int seq_nkeys, void *query_keys, int query_nkeys, int k_size)
{
    /* Future: Add SIMD dispatch based on simd_capability */
    return kmersearch_count_matching_uintkey_scalar(seq_keys, seq_nkeys, query_keys, query_nkeys, k_size);
}

/*
 * Extract uint keys with occurrence counting from text sequence
 */
void
kmersearch_extract_uintkey_from_text(const char *text, void **output, int *nkeys)
{
    int text_len = strlen(text);
    int bit_len = text_len * 4;  /* 4 bits per character for DNA4 */
    int byte_len = (bit_len + 7) / 8;
    VarBit *dna4_seq;
    bits8 *data_ptr;
    int i;
    
    /* Validate input characters */
    for (i = 0; i < text_len; i++) {
        if (!kmersearch_is_valid_dna4_char(text[i])) {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("Invalid character '%c' for DNA4 sequence", text[i]),
                     errhint("DNA4 accepts A,C,G,T,U,M,R,W,S,Y,K,V,H,D,B,N characters")));
        }
    }
    
    /* Allocate DNA4 VarBit */
    dna4_seq = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(dna4_seq, VARHDRSZ + sizeof(int32) + byte_len);
    VARBITLEN(dna4_seq) = bit_len;
    
    /* Encode text to DNA4 using SIMD dispatch */
    data_ptr = VARBITS(dna4_seq);
    dna4_encode(text, (uint8_t*)data_ptr, text_len);
    
    /* Extract uint keys from DNA4 sequence */
    kmersearch_extract_uintkey_from_dna4(dna4_seq, output, nkeys);
    
    /* Free temporary DNA4 sequence */
    pfree(dna4_seq);
}