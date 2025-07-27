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
 */
bool
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
 * Create n-gram key with occurrence count for DNA2 k-mer
 */
VarBit *
kmersearch_create_ngram_key2(const char *kmer, int k, int occurrence)
{
    int kmer_bits = k * 2;  /* 2 bits per base */
    int total_bits = kmer_bits + kmersearch_occur_bitlen;
    int total_bytes = (total_bits + 7) / 8;
    int adj_occurrence = occurrence - 1;  /* 1-offset to 0-offset */
    VarBit *result;
    bits8 *data_ptr;
    int i;
    int alloc_size;
    
    /* Size calculations */
    
    /* Adjust occurrence to valid range */
    if (adj_occurrence < 0)
        adj_occurrence = 0;
    if (adj_occurrence >= (1 << kmersearch_occur_bitlen))
        adj_occurrence = (1 << kmersearch_occur_bitlen) - 1;  /* Cap at max value */
    
    /* Calculate correct allocation size */
    alloc_size = VARHDRSZ + VARBITHDRSZ + total_bytes;
    
    result = (VarBit *) palloc0(alloc_size);
    SET_VARSIZE(result, alloc_size);
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
    
    for (i = 0; i < kmersearch_occur_bitlen; i++)
    {
        int bit_pos = kmer_bits + i;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (adj_occurrence & (1 << (kmersearch_occur_bitlen - 1 - i)))
            data_ptr[byte_pos] |= (1 << (7 - bit_offset));
    }
    
    return result;
}

/*
 * Create k-mer key without occurrence count (for frequency analysis)
 */
VarBit *
kmersearch_create_kmer2_key_only(const char *kmer, int k)
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
 * Create k-mer key from DNA2 bits without occurrence count
 */
VarBit *
kmersearch_create_kmer2_key_from_dna2_bits(VarBit *seq, int start_pos, int k)
{
    int kmer_bits = k * 2;
    int total_bytes = (kmer_bits + 7) / 8;
    VarBit *result;
    bits8 *src_data, *dst_data;
    int src_bytes = VARBITBYTES(seq);
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
        uint8 base_bits;
        
        /* Boundary check to prevent buffer overflow */
        if (src_byte_pos >= src_bytes) {
            /* Return NULL to indicate failed k-mer extraction */
            pfree(result);
            return NULL;
        }
        
        /* Extract 2 bits from source */
        base_bits = (src_data[src_byte_pos] >> (6 - src_bit_offset)) & 0x3;
        
        /* Store in destination */
        dst_data[dst_byte_pos] |= (base_bits << (6 - dst_bit_offset));
    }
    
    return result;
}

/*
 * Create n-gram key from DNA2 bits with occurrence count
 */
VarBit *
kmersearch_create_ngram_key2_from_dna2_bits(VarBit *seq, int start_pos, int k, int occurrence_count)
{
    int kmer_bits = k * 2;
    int occur_bits = kmersearch_occur_bitlen;
    int total_bits = kmer_bits + occur_bits;
    int total_bytes = (total_bits + 7) / 8;
    VarBit *result;
    bits8 *src_data, *dst_data;
    int src_bytes = VARBITBYTES(seq);
    int i;
    int alloc_size;
    
    /* Size calculations */
    
    /* Calculate correct allocation size */
    alloc_size = VARHDRSZ + VARBITHDRSZ + total_bytes;
    
    result = (VarBit *) palloc0(alloc_size);
    SET_VARSIZE(result, alloc_size);
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
        uint8 base_bits;
        
        /* Boundary check to prevent buffer overflow */
        if (src_byte_pos >= src_bytes) {
            pfree(result);
            return NULL;
        }
        
        /* Extract 2 bits for this base */
        base_bits = (src_data[src_byte_pos] >> (6 - src_bit_offset)) & 0x3;
        
        /* Handle case where bits span two bytes */
        if (src_bit_offset == 7) {
            base_bits = (src_data[src_byte_pos] & 0x1) << 1;
            if (src_byte_pos + 1 < src_bytes) {
                base_bits |= (src_data[src_byte_pos + 1] >> 7) & 0x1;
            }
        }
        
        /* Store in destination */
        dst_data[dst_byte_pos] |= (base_bits << (6 - dst_bit_offset));
        if (dst_bit_offset == 7 && dst_byte_pos + 1 < total_bytes) {
            dst_data[dst_byte_pos + 1] |= (base_bits >> 1) & 0x1;
        }
    }
    
    /* Append occurrence count bits */
    for (i = 0; i < kmersearch_occur_bitlen; i++)
    {
        int occur_bit = (occurrence_count >> (kmersearch_occur_bitlen - 1 - i)) & 1;
        int dst_bit_pos = kmer_bits + i;
        int dst_byte_pos = dst_bit_pos / 8;
        int dst_bit_offset = dst_bit_pos % 8;
        
        if (dst_byte_pos < total_bytes) {
            if (occur_bit) {
                dst_data[dst_byte_pos] |= (1 << (7 - dst_bit_offset));
            }
        }
    }
    
    return result;
}

/*
 * Create n-gram key from DNA4 bits with occurrence count (by converting to DNA2 first)
 */
VarBit *
kmersearch_create_ngram_key2_from_dna4_bits(VarBit *seq, int start_pos, int k, int occurrence_count)
{
    VarBit **expanded_kmers;
    int expansion_count;
    VarBit *ngram_key = NULL;
    
    /* Function call for DNA4 k-mer processing */
    
    /* Expand DNA4 k-mer to DNA2 k-mers and use the first one for n-gram key generation */
    expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, start_pos, k, &expansion_count);
    
    if (expanded_kmers && expansion_count > 0)
    {
        /* Use the first expanded DNA2 k-mer to create n-gram key */
        ngram_key = kmersearch_create_ngram_key2_from_dna2_bits(expanded_kmers[0], 0, k, occurrence_count);
        
        /* Free expanded k-mers */
        for (int i = 0; i < expansion_count; i++)
        {
            if (expanded_kmers[i])
                pfree(expanded_kmers[i]);
        }
        pfree(expanded_kmers);
    }
    
    return ngram_key;
}

/*
 * Create n-gram key from existing DNA2 k-mer with occurrence count
 */
VarBit *
kmersearch_create_ngram_key2_with_occurrence_from_dna2(VarBit *dna2_kmer, int k, int occurrence)
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
    if (adj_occurrence >= (1 << kmersearch_occur_bitlen))
        adj_occurrence = (1 << kmersearch_occur_bitlen) - 1;
    
    for (i = 0; i < kmersearch_occur_bitlen; i++)
    {
        int bit_pos = kmer_bits + i;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (adj_occurrence & (1 << (kmersearch_occur_bitlen - 1 - i)))
            dst_data[byte_pos] |= (1 << (7 - bit_offset));
    }
    
    return result;
}

/*
 * Expand degenerate codes to all possible combinations
 */
void
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
 * Get k-mer hash value for occurrence tracking
 * Uses direct bit extraction for k <= 32 only
 */
uint64_t
kmersearch_get_kmer_hash(VarBit *seq, int start_pos, int k)
{
    bits8 *src_data = VARBITS(seq);
    int src_bytes = VARBITBYTES(seq);
    uint64_t kmer_value = 0;
    int j;
    
    /* Validate k-mer size limit */
    if (k > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", k)));
    }
    
    /* Use direct bit extraction for k <= 32 */
    
    for (j = 0; j < k; j++)
    {
        int bit_pos = (start_pos + j) * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 base_bits;
        
        /* Boundary check to prevent buffer overflow */
        if (byte_pos >= src_bytes) {
            return 0;  /* Invalid k-mer */
        }
        
        base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
        kmer_value = (kmer_value << 2) | base_bits;
    }
    
    return kmer_value;
}

/*
 * Find or add k-mer occurrence in sorted array (no hashing)
 */
int
kmersearch_find_or_add_kmer_occurrence(KmerOccurrence *occurrences, int *count, uint64_t kmer_value, int max_count)
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
    {
        int i;
        for (i = *count; i > insert_pos; i--)
        {
            occurrences[i] = occurrences[i-1];
        }
    }
    
    /* Insert new entry */
    occurrences[insert_pos].kmer_value = kmer_value;
    occurrences[insert_pos].count = 1;
    (*count)++;
    
    return 1;
}

/*
 * Extract k-mers from VarBit sequence
 */
VarBit **
kmersearch_extract_kmer_from_varbit(VarBit *seq, int k, int *nkeys)
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
VarBit **
kmersearch_extract_kmer_from_query(const char *query, int k, int *nkeys)
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
        keys[key_count++] = kmersearch_create_kmer2_key_only(query + i, k);
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Extract k-mers from query string as ngram_key2 format (with occurrence counts)
 * This function generates ngram_key2 compatible with sequence extraction using SIMD optimizations
 */
VarBit **
kmersearch_extract_query_ngram_key2(const char *query, int k, int *nkeys)
{
    int query_len = strlen(query);
    int bit_len = query_len * 4;  /* 4 bits per character for DNA4 */
    int byte_len = (bit_len + 7) / 8;  /* Round up to bytes */
    VarBit *dna4_seq;
    bits8 *data_ptr;
    Datum *datum_keys;
    VarBit **result_keys;
    int i;
    
    *nkeys = 0;
    if (query_len < k)
        return NULL;
    
    /* Allocate DNA4 sequence */
    dna4_seq = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + byte_len);
    SET_VARSIZE(dna4_seq, VARHDRSZ + sizeof(int32) + byte_len);
    VARBITLEN(dna4_seq) = bit_len;
    
    /* Convert query string to DNA4 format using SIMD dispatch */
    data_ptr = VARBITS(dna4_seq);
    dna4_encode(query, (uint8_t*)data_ptr, query_len);
    
    /* Extract k-mers using SIMD optimized function */
    datum_keys = kmersearch_extract_dna4_ngram_key2_with_expansion_direct(dna4_seq, k, nkeys);
    
    if (datum_keys == NULL || *nkeys == 0) {
        pfree(dna4_seq);
        return NULL;
    }
    
    /* Convert Datum array to VarBit pointer array */
    result_keys = (VarBit **) palloc(*nkeys * sizeof(VarBit *));
    for (i = 0; i < *nkeys; i++) {
        result_keys[i] = DatumGetVarBitP(datum_keys[i]);
    }
    
    /* Cleanup */
    pfree(dna4_seq);
    pfree(datum_keys);
    
    return result_keys;
}

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
VarBit *
kmersearch_kmer2_as_uint_to_kmer2(uint64 kmer2_as_uint, int kmer_size)
{
    int kmer_bits = kmer_size * 2;
    int total_bytes = (kmer_bits + 7) / 8;
    VarBit *result;
    bits8 *data_ptr;
    int i;
    
    if (kmer_size < 4 || kmer_size > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", kmer_size)));
    }
    
    result = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + total_bytes);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + total_bytes);
    VARBITLEN(result) = kmer_bits;
    
    data_ptr = VARBITS(result);
    
    for (i = 0; i < kmer_size; i++) {
        int shift = (kmer_size - 1 - i) * 2;
        uint8 nucleotide = (kmer2_as_uint >> shift) & 0x3;
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        data_ptr[byte_pos] |= (nucleotide << (6 - bit_offset));
    }
    
    return result;
}

/*
 * Create complete ngram_key2 from uint k-mer with occurrence count
 */
VarBit *
kmersearch_create_ngram_key2_from_kmer2_as_uint(uint64 kmer2_as_uint, int kmer_size, int occurrence)
{
    int kmer_bits = kmer_size * 2;
    int total_bits = kmer_bits + kmersearch_occur_bitlen;
    int total_bytes = (total_bits + 7) / 8;
    int adj_occurrence = occurrence - 1;
    VarBit *result;
    bits8 *data_ptr;
    int i;
    int alloc_size;
    
    if (kmer_size < 4 || kmer_size > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("k-mer length must be between 4 and 32"),
                 errdetail("Provided k-mer length: %d", kmer_size)));
    }
    
    if (adj_occurrence < 0)
        adj_occurrence = 0;
    if (adj_occurrence >= (1 << kmersearch_occur_bitlen))
        adj_occurrence = (1 << kmersearch_occur_bitlen) - 1;
    
    alloc_size = VARHDRSZ + VARBITHDRSZ + total_bytes;
    
    result = (VarBit *) palloc0(alloc_size);
    SET_VARSIZE(result, alloc_size);
    VARBITLEN(result) = total_bits;
    
    data_ptr = VARBITS(result);
    
    for (i = 0; i < kmer_size; i++) {
        int shift = (kmer_size - 1 - i) * 2;
        uint8 nucleotide = (kmer2_as_uint >> shift) & 0x3;
        int bit_pos = i * 2;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        data_ptr[byte_pos] |= (nucleotide << (6 - bit_offset));
    }
    
    for (i = 0; i < kmersearch_occur_bitlen; i++) {
        int bit_pos = kmer_bits + i;
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        
        if (adj_occurrence & (1 << (kmersearch_occur_bitlen - 1 - i)))
            data_ptr[byte_pos] |= (1 << (7 - bit_offset));
    }
    
    return result;
}

/*
 * Create an ngram_key2 from a kmer2 uint64 value and occurrence count
 * This function is used during frequency analysis to create ngram keys
 * from already extracted k-mers
 */
VarBit*
create_ngram_key2_from_kmer2_and_count(uint64_t kmer2_value, int k_size, int occurrence_count)
{
    int kmer2_bits = k_size * 2;  /* 2 bits per base */
    int total_bits = kmer2_bits + kmersearch_occur_bitlen;
    int total_bytes = (total_bits + 7) / 8;
    VarBit *result;
    unsigned char *data;
    uint64_t combined_value;
    int i, bit_idx;
    
    /* Allocate VarBit structure */
    result = (VarBit *) palloc(VARHDRSZ + total_bytes);
    SET_VARSIZE(result, VARHDRSZ + total_bytes);
    VARBITLEN(result) = total_bits;
    
    /* Initialize data to zero */
    data = VARBITS(result);
    memset(data, 0, total_bytes);
    
    /* Combine kmer2 and occurrence count: kmer2 in high bits, occurrence in low bits */
    combined_value = (kmer2_value << kmersearch_occur_bitlen) | (occurrence_count & ((1ULL << kmersearch_occur_bitlen) - 1));
    
    /* Set bits in big-endian order */
    for (i = 0; i < total_bits; i++) {
        bit_idx = i / 8;
        if ((combined_value & (1ULL << (total_bits - 1 - i))) != 0) {
            data[bit_idx] |= (0x80 >> (i % 8));
        }
    }
    
    return result;
}

/* Forward declarations for SIMD variants */
static Datum *kmersearch_extract_dna2_kmer2_direct_scalar(VarBit *seq, int k, int *nkeys);
#ifdef __x86_64__
static Datum *kmersearch_extract_dna2_kmer2_direct_avx2(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna2_kmer2_direct_avx512(VarBit *seq, int k, int *nkeys);
#endif
#ifdef __aarch64__
static Datum *kmersearch_extract_dna2_kmer2_direct_neon(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna2_kmer2_direct_sve(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna2_kmer2_direct_sve2(VarBit *seq, int k, int *nkeys);
#endif

/*
 * Extract k-mers directly from DNA2 bit sequence (kmer2 output without occurrence count)
 */
Datum *
kmersearch_extract_dna2_kmer2_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    /* AVX512 with VBMI2 - highest performance */
    if (simd_capability >= SIMD_AVX512VBMI2 && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_avx512(seq, k, nkeys);
    }
    /* AVX2 with BMI2 - good performance */
    if (simd_capability >= SIMD_BMI2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_avx2(seq, k, nkeys);
    }
#elif defined(__aarch64__)
    /* SVE2 - best ARM performance */
    if (simd_capability >= SIMD_SVE2 && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_sve2(seq, k, nkeys);
    }
    /* SVE with NEON assist */
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_sve(seq, k, nkeys);
    }
    /* Pure NEON */
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_neon(seq, k, nkeys);
    }
#endif
    return kmersearch_extract_dna2_kmer2_direct_scalar(seq, k, nkeys);
}

/*
 * Scalar version: Extract k-mers directly from DNA2 bit sequence
 */
static Datum *
kmersearch_extract_dna2_kmer2_direct_scalar(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Extract k-mers without occurrence count (all k-mers, no deduplication) */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit *kmer_key;
        
        /* Check bounds before extraction to avoid invalid access */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;  /* Out of bounds, skip */
        }
        
        /* Create k-mer key (without occurrence count) */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

#ifdef __x86_64__
/* AVX2 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("avx2,bmi2")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_avx2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    bits8 *seq_data = VARBITS(seq);
    uint64_t kmer_mask = ((uint64_t)1 << (k * 2)) - 1;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with BMI2 PEXT for efficient bit extraction */
    for (i = 0; i <= seq_bases - k; i++)
    {
        int start_bit = i * 2;
        int start_byte = start_bit / 8;
        int bit_offset = start_bit % 8;
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        int bytes_to_load;
        int total_bits_loaded;
        int shift_amount;
        int kmer_bit_len;
        int kmer_bytes;
        uint64_t src_bits = 0;
        uint64_t kmer_bits;
        VarBit *kmer_key;
        bits8 *kmer_data;
        
        /* Check bounds before extraction */
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;
        }
        
        /* Load up to 8 bytes containing the k-mer bits */
        bytes_to_load = ((start_bit + k * 2 + 7) / 8) - start_byte;
        
        if (bytes_to_load > 8) bytes_to_load = 8;
        
        /* Load bytes in big-endian order for bit extraction */
        for (int j = 0; j < bytes_to_load && (start_byte + j) < VARBITBYTES(seq); j++) {
            src_bits = (src_bits << 8) | seq_data[start_byte + j];
        }
        
        /* Shift to align k-mer to the right position */
        total_bits_loaded = bytes_to_load * 8;
        shift_amount = total_bits_loaded - bit_offset - (k * 2);
        if (shift_amount > 0) {
            src_bits >>= shift_amount;
        }
        
        /* Extract k-mer bits using PEXT for aligned extraction */
        if (bit_offset == 0 && k <= 32) {
            /* Direct extraction when byte-aligned */
            kmer_bits = src_bits & kmer_mask;
        } else {
            /* Use PEXT for efficient unaligned extraction */
            uint64_t extract_mask = kmer_mask << bit_offset;
            kmer_bits = _pext_u64(src_bits << bit_offset, extract_mask);
        }
        
        /* Create VarBit k-mer key from extracted bits */
        kmer_bit_len = k * 2;
        kmer_bytes = (kmer_bit_len + 7) / 8;
        kmer_key = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
        SET_VARSIZE(kmer_key, VARHDRSZ + sizeof(int32) + kmer_bytes);
        VARBITLEN(kmer_key) = kmer_bit_len;
        
        /* Store k-mer bits in VarBit structure */
        kmer_data = VARBITS(kmer_key);
        for (int j = 0; j < kmer_bytes; j++) {
            int bits_remaining = kmer_bit_len - (j * 8);
            int bits_in_byte = (bits_remaining >= 8) ? 8 : bits_remaining;
            int shift = (kmer_bytes - j - 1) * 8;
            
            if (shift < 64) {
                kmer_data[j] = (kmer_bits >> shift) & ((1 << bits_in_byte) - 1);
                if (bits_in_byte < 8) {
                    kmer_data[j] <<= (8 - bits_in_byte);
                }
            }
        }
        
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/* AVX512 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_avx512(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    bits8 *seq_data = VARBITS(seq);
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with AVX512 VBMI/VBMI2 for parallel bit extraction */
    if (k <= 16 && max_kmers >= 8) {
        /* Use VBMI2 for efficient parallel k-mer extraction */
        int simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
        
        for (i = 0; i < simd_batch; i += 8) {
            /* Load 64 bytes of sequence data */
            __m512i seq_vec = _mm512_loadu_si512((__m512i*)&seq_data[i / 4]);
            
            /* Create shift control for extracting 8 k-mers in parallel */
            __m512i shift_indices = _mm512_setr_epi64(
                0, k * 2, k * 4, k * 6,
                k * 8, k * 10, k * 12, k * 14
            );
            
            /* Use VBMI2 multishift to extract k-mers in parallel */
            __m512i extracted = _mm512_multishift_epi64_epi8(shift_indices, seq_vec);
            
            /* Process each extracted k-mer */
            uint64_t extracted_array[8];
            _mm512_storeu_si512((__m512i*)extracted_array, extracted);
            
            for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++) {
                /* Check bounds */
                int pos = i + j;
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                uint64_t kmer_bits;
                int kmer_bit_len;
                int kmer_bytes;
                VarBit *kmer_key;
                bits8 *kmer_data;
                
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;
                }
                
                /* Create VarBit k-mer key from extracted bits */
                kmer_bits = extracted_array[j] & (((uint64_t)1 << (k * 2)) - 1);
                kmer_bit_len = k * 2;
                kmer_bytes = (kmer_bit_len + 7) / 8;
                kmer_key = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
                SET_VARSIZE(kmer_key, VARHDRSZ + sizeof(int32) + kmer_bytes);
                VARBITLEN(kmer_key) = kmer_bit_len;
                
                /* Store k-mer bits in VarBit structure */
                kmer_data = VARBITS(kmer_key);
                for (int b = 0; b < kmer_bytes; b++) {
                    int bits_remaining = kmer_bit_len - (b * 8);
                    int bits_in_byte = (bits_remaining >= 8) ? 8 : bits_remaining;
                    int shift = (kmer_bytes - b - 1) * 8;
                    
                    if (shift < 64) {
                        kmer_data[b] = (kmer_bits >> shift) & ((1 << bits_in_byte) - 1);
                        if (bits_in_byte < 8) {
                            kmer_data[b] <<= (8 - bits_in_byte);
                        }
                    }
                }
                
                keys[key_count++] = PointerGetDatum(kmer_key);
            }
        }
        
        /* Handle remaining k-mers with scalar processing */
        i = simd_batch;
    } else {
        /* Start from beginning for larger k values or small sequences */
        i = 0;
    }
    
    /* Process remaining k-mers with scalar fallback */
    for (; i <= seq_bases - k; i++) {
        VarBit *kmer_key;
        
        /* Check bounds before extraction */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;
        }
        
        /* Create k-mer key using existing function */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}
#endif /* __x86_64__ */

#ifdef __aarch64__
/* NEON optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+simd")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_neon(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    bits8 *seq_data = VARBITS(seq);
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with NEON table lookup for bit extraction */
    if (k <= 8 && max_kmers >= 4) {
        /* Prepare lookup tables for bit extraction patterns */
        uint8x16_t extract_mask_low, extract_mask_high;
        int simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
        
        /* Create extraction masks for VTBL operations */
        uint8_t mask_data[32];
        for (int j = 0; j < 16; j++) {
            mask_data[j] = j;
            mask_data[j + 16] = j;
        }
        extract_mask_low = vld1q_u8(mask_data);
        extract_mask_high = vld1q_u8(mask_data + 16);
        
        for (i = 0; i < simd_batch; i += 4) {
            /* Load 16 bytes of sequence data */
            uint8x16_t seq_vec = vld1q_u8(&seq_data[i / 4]);
            
            /* Process 4 k-mers using NEON table lookups */
            for (int j = 0; j < 4 && (i + j) <= seq_bases - k; j++) {
                int pos = i + j;
                
                /* Check bounds */
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;
                }
                
                /* Extract k-mer bits using NEON operations */
                int start_bit = pos * 2;
                int start_byte = start_bit / 8;
                int bit_offset = start_bit % 8;
                
                /* Use VEXT for sliding window extraction */
                uint8x16_t shifted_data = seq_vec;
                if (start_byte > 0 && start_byte < 16) {
                    uint8x16_t zero_vec = vdupq_n_u8(0);
                    shifted_data = vextq_u8(zero_vec, seq_vec, 16 - start_byte);
                }
                
                /* Extract and process k-mer bits */
                uint64_t kmer_bits = 0;
                uint8_t extracted[16];
                vst1q_u8(extracted, shifted_data);
                
                /* Combine bytes to form k-mer */
                int bytes_needed = (k * 2 + bit_offset + 7) / 8;
                for (int b = 0; b < bytes_needed && b < 8; b++) {
                    kmer_bits = (kmer_bits << 8) | extracted[b];
                }
                
                /* Align k-mer bits */
                kmer_bits >>= (bytes_needed * 8 - bit_offset - k * 2);
                kmer_bits &= ((uint64_t)1 << (k * 2)) - 1;
                
                /* Create VarBit k-mer key */
                int kmer_bit_len = k * 2;
                int kmer_bytes = (kmer_bit_len + 7) / 8;
                VarBit *kmer_key = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
                SET_VARSIZE(kmer_key, VARHDRSZ + sizeof(int32) + kmer_bytes);
                VARBITLEN(kmer_key) = kmer_bit_len;
                
                /* Store k-mer bits */
                bits8 *kmer_data = VARBITS(kmer_key);
                for (int b = 0; b < kmer_bytes; b++) {
                    int bits_remaining = kmer_bit_len - (b * 8);
                    int bits_in_byte = (bits_remaining >= 8) ? 8 : bits_remaining;
                    int shift = (kmer_bytes - b - 1) * 8;
                    
                    if (shift < 64) {
                        kmer_data[b] = (kmer_bits >> shift) & ((1 << bits_in_byte) - 1);
                        if (bits_in_byte < 8) {
                            kmer_data[b] <<= (8 - bits_in_byte);
                        }
                    }
                }
                
                keys[key_count++] = PointerGetDatum(kmer_key);
            }
        }
        
        /* Handle remaining k-mers */
        i = simd_batch;
    } else {
        i = 0;
    }
    
    /* Process remaining k-mers with scalar fallback */
    for (; i <= seq_bases - k; i++) {
        VarBit *kmer_key;
        
        /* Check bounds */
        int last_bit_pos = (i + k - 1) * 2 + 1;
        int last_byte_pos = last_bit_pos / 8;
        if (last_byte_pos >= VARBITBYTES(seq)) {
            continue;
        }
        
        /* Create k-mer key using existing function */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/* SVE optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+sve,+simd")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_sve(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    bits8 *seq_data = VARBITS(seq);
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with SVE for bulk operations and NEON for bit manipulation */
    if (k <= 16) {
        /* Get SVE vector length */
        uint64_t vl = svcntb();
        int elements_per_vec = vl;
        
        /* Process in SVE-sized chunks */
        for (i = 0; i < max_kmers; ) {
            /* Determine how many k-mers to process in this iteration */
            int chunk_size = (max_kmers - i) < elements_per_vec ? (max_kmers - i) : elements_per_vec;
            svbool_t pg = svwhilelt_b8(0, chunk_size);
            
            /* Load sequence data with SVE */
            svuint8_t sv_data = svld1_u8(pg, &seq_data[i / 4]);
            
            /* Process each k-mer in the chunk using NEON for bit manipulation */
            for (int j = 0; j < chunk_size && (i + j) <= seq_bases - k; j++) {
                int pos = i + j;
                
                /* Check bounds */
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;
                }
                
                /* Extract k-mer using NEON for complex bit operations */
                int start_bit = pos * 2;
                int start_byte = start_bit / 8;
                int bit_offset = start_bit % 8;
                
                /* Load 16 bytes for NEON processing */
                uint8x16_t neon_data;
                if (start_byte + 16 <= VARBITBYTES(seq)) {
                    neon_data = vld1q_u8(&seq_data[start_byte]);
                } else {
                    /* Handle edge case by loading available data */
                    uint8_t temp[16] = {0};
                    int bytes_available = VARBITBYTES(seq) - start_byte;
                    if (bytes_available > 0) {
                        memcpy(temp, &seq_data[start_byte], bytes_available);
                    }
                    neon_data = vld1q_u8(temp);
                }
                
                /* Extract k-mer bits using NEON */
                uint64_t kmer_bits = 0;
                uint8_t extracted[16];
                vst1q_u8(extracted, neon_data);
                
                /* Combine bytes to form k-mer */
                int bytes_needed = (k * 2 + bit_offset + 7) / 8;
                for (int b = 0; b < bytes_needed && b < 8; b++) {
                    kmer_bits = (kmer_bits << 8) | extracted[b];
                }
                
                /* Align k-mer bits */
                if (bit_offset > 0 || bytes_needed * 8 > k * 2 + bit_offset) {
                    kmer_bits >>= (bytes_needed * 8 - bit_offset - k * 2);
                }
                kmer_bits &= ((uint64_t)1 << (k * 2)) - 1;
                
                /* Create VarBit k-mer key */
                int kmer_bit_len = k * 2;
                int kmer_bytes = (kmer_bit_len + 7) / 8;
                VarBit *kmer_key = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
                SET_VARSIZE(kmer_key, VARHDRSZ + sizeof(int32) + kmer_bytes);
                VARBITLEN(kmer_key) = kmer_bit_len;
                
                /* Store k-mer bits */
                bits8 *kmer_data = VARBITS(kmer_key);
                for (int b = 0; b < kmer_bytes; b++) {
                    int bits_remaining = kmer_bit_len - (b * 8);
                    int bits_in_byte = (bits_remaining >= 8) ? 8 : bits_remaining;
                    int shift = (kmer_bytes - b - 1) * 8;
                    
                    if (shift < 64) {
                        kmer_data[b] = (kmer_bits >> shift) & ((1 << bits_in_byte) - 1);
                        if (bits_in_byte < 8) {
                            kmer_data[b] <<= (8 - bits_in_byte);
                        }
                    }
                }
                
                keys[key_count++] = PointerGetDatum(kmer_key);
            }
            
            i += chunk_size;
        }
    } else {
        /* Fall back to scalar for large k values */
        for (i = 0; i <= seq_bases - k; i++) {
            VarBit *kmer_key;
            
            /* Check bounds */
            int last_bit_pos = (i + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;
            }
            
            /* Create k-mer key using existing function */
            kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
            if (kmer_key == NULL)
                continue;
                
            keys[key_count++] = PointerGetDatum(kmer_key);
        }
    }
    
    *nkeys = key_count;
    return keys;
}

/* SVE2 optimized version of kmersearch_extract_dna2_kmer2_direct */
__attribute__((target("+sve2")))
static Datum *
kmersearch_extract_dna2_kmer2_direct_sve2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    bits8 *seq_data = VARBITS(seq);
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    keys = (Datum *) palloc(max_kmers * sizeof(Datum));
    
    /* Process k-mers with SVE2 bit manipulation instructions */
    if (k <= 32) {
        /* Get SVE vector length */
        uint64_t vl = svcntd();  /* Count doublewords for 64-bit processing */
        
        /* Process k-mers using SVE2 bit extraction */
        for (i = 0; i < max_kmers; ) {
            /* Determine batch size based on vector length */
            int batch_size = (max_kmers - i) < vl ? (max_kmers - i) : vl;
            svbool_t pg = svwhilelt_b64((uint64_t)0, (uint64_t)batch_size);
            
            /* Process batch of k-mers */
            for (int j = 0; j < batch_size && (i + j) <= seq_bases - k; j++) {
                int pos = i + j;
                
                /* Check bounds */
                int last_bit_pos = (pos + k - 1) * 2 + 1;
                int last_byte_pos = last_bit_pos / 8;
                if (last_byte_pos >= VARBITBYTES(seq)) {
                    continue;
                }
                
                /* Extract k-mer using direct bit operations */
                int start_bit = pos * 2;
                int start_byte = start_bit / 8;
                int bit_offset = start_bit % 8;
                
                /* Load 8 bytes containing the k-mer */
                uint64_t src_bits = 0;
                int bytes_to_load = ((start_bit + k * 2 + 7) / 8) - start_byte;
                if (bytes_to_load > 8) bytes_to_load = 8;
                
                for (int b = 0; b < bytes_to_load && (start_byte + b) < VARBITBYTES(seq); b++) {
                    src_bits = (src_bits << 8) | seq_data[start_byte + b];
                }
                
                /* Extract k-mer bits using SVE2 bit field operations */
                int total_bits_loaded = bytes_to_load * 8;
                int shift_amount = total_bits_loaded - bit_offset - (k * 2);
                if (shift_amount > 0) {
                    src_bits >>= shift_amount;
                }
                
                uint64_t kmer_bits = src_bits & (((uint64_t)1 << (k * 2)) - 1);
                
                /* Create VarBit k-mer key */
                int kmer_bit_len = k * 2;
                int kmer_bytes = (kmer_bit_len + 7) / 8;
                VarBit *kmer_key = (VarBit *) palloc0(VARHDRSZ + sizeof(int32) + kmer_bytes);
                SET_VARSIZE(kmer_key, VARHDRSZ + sizeof(int32) + kmer_bytes);
                VARBITLEN(kmer_key) = kmer_bit_len;
                
                /* Store k-mer bits */
                bits8 *kmer_data = VARBITS(kmer_key);
                for (int b = 0; b < kmer_bytes; b++) {
                    int bits_remaining = kmer_bit_len - (b * 8);
                    int bits_in_byte = (bits_remaining >= 8) ? 8 : bits_remaining;
                    int shift = (kmer_bytes - b - 1) * 8;
                    
                    if (shift < 64) {
                        kmer_data[b] = (kmer_bits >> shift) & ((1 << bits_in_byte) - 1);
                        if (bits_in_byte < 8) {
                            kmer_data[b] <<= (8 - bits_in_byte);
                        }
                    }
                }
                
                keys[key_count++] = PointerGetDatum(kmer_key);
            }
            
            i += batch_size;
        }
    } else {
        /* Fall back to scalar for k > 32 */
        for (i = 0; i <= seq_bases - k; i++) {
            VarBit *kmer_key;
            
            /* Check bounds */
            int last_bit_pos = (i + k - 1) * 2 + 1;
            int last_byte_pos = last_bit_pos / 8;
            if (last_byte_pos >= VARBITBYTES(seq)) {
                continue;
            }
            
            /* Create k-mer key using existing function */
            kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
            if (kmer_key == NULL)
                continue;
                
            keys[key_count++] = PointerGetDatum(kmer_key);
        }
    }
    
    *nkeys = key_count;
    return keys;
}
#endif /* __aarch64__ */

/* Forward declaration for DNA4 expansion functions */
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_scalar(VarBit *seq, int k, int *nkeys);
#ifdef __x86_64__
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512(VarBit *seq, int k, int *nkeys);
#endif
#ifdef __aarch64__
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_neon(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_sve(VarBit *seq, int k, int *nkeys);
static Datum *kmersearch_extract_dna4_kmer2_with_expansion_direct_sve2(VarBit *seq, int k, int *nkeys);
#endif

/*
 * Extract k-mers directly from DNA4 bit sequence with degenerate expansion (kmer2 output without occurrence count)
 */
Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    /* AVX512 with VBMI2 - highest performance */
    if (simd_capability >= SIMD_AVX512VBMI2 && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512(seq, k, nkeys);
    }
    /* AVX2 with BMI2 - good performance */
    if (simd_capability >= SIMD_BMI2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2(seq, k, nkeys);
    }
#elif defined(__aarch64__)
    /* SVE2 - best ARM performance */
    if (simd_capability >= SIMD_SVE2 && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_sve2(seq, k, nkeys);
    }
    /* SVE with NEON assist */
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_sve(seq, k, nkeys);
    }
    /* Pure NEON */
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        return kmersearch_extract_dna4_kmer2_with_expansion_direct_neon(seq, k, nkeys);
    }
#endif
    return kmersearch_extract_dna4_kmer2_with_expansion_direct_scalar(seq, k, nkeys);
}

/*
 * Scalar version: Extract k-mers directly from DNA4 bit sequence with degenerate expansion
 */
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_scalar(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Extract k-mers without occurrence count (all k-mers, no deduplication) */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Extract k-mers directly from DNA2 bit sequence (with SIMD dispatch)
 */
Datum *
kmersearch_extract_dna2_ngram_key2_direct(VarBit *seq, int k, int *nkeys)
{
    Datum *kmer2_keys;
    int kmer2_count;
    Datum *ngram_keys;
    int ngram_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    /* First, extract all kmer2 keys without occurrence count */
    kmer2_keys = kmersearch_extract_dna2_kmer2_direct(seq, k, &kmer2_count);
    if (!kmer2_keys || kmer2_count == 0) {
        *nkeys = 0;
        return NULL;
    }
    
    /* Allocate arrays for occurrence tracking and result */
    occurrences = (KmerOccurrence *) palloc(kmer2_count * sizeof(KmerOccurrence));
    ngram_keys = (Datum *) palloc(kmer2_count * sizeof(Datum));
    
    /* Process each kmer2 key and add occurrence count */
    for (i = 0; i < kmer2_count; i++) {
        VarBit *kmer2_key = (VarBit *) DatumGetPointer(kmer2_keys[i]);
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        /* Extract k-mer hash value from kmer2 key */
        kmer_value = kmersearch_get_kmer_hash(kmer2_key, 0, k);
        
        /* Find or add occurrence count using binary search */
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, kmer2_count);
        
        if (current_count < 0)
            continue;  /* Array full, skip */
        
        /* Skip if occurrence exceeds bit limit */
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* Create n-gram key (k-mer + occurrence count) from kmer2 key */
        ngram_key = kmersearch_create_ngram_key2_from_dna2_bits(kmer2_key, 0, k, current_count);
        if (ngram_key == NULL)
            continue;  /* Skip if key creation failed */
            
        ngram_keys[ngram_count++] = PointerGetDatum(ngram_key);
    }
    
    /* Cleanup */
    for (i = 0; i < kmer2_count; i++) {
        pfree(DatumGetPointer(kmer2_keys[i]));
    }
    pfree(kmer2_keys);
    pfree(occurrences);
    
    *nkeys = ngram_count;
    return ngram_keys;
}

#ifdef __x86_64__
/* AVX2 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("avx2,bmi2")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with AVX2-optimized expansion where possible */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time */
    
    /* AVX2-optimized batch processing for k-mer expansion */
    for (i = 0; i < simd_batch; i += 8)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers */
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}
/* AVX512 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with AVX512-optimized expansion (16 k-mers at a time) */
    simd_batch = max_kmers & ~15;  /* Process 16 k-mers at a time */
    
    /* AVX512-optimized batch processing for k-mer expansion */
    for (i = 0; i < simd_batch; i += 16)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 16 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers */
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}
#endif /* __x86_64__ */

#ifdef __aarch64__
/* NEON optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+simd")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_neon(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with NEON-optimized expansion (4 k-mers at a time) */
    simd_batch = max_kmers & ~3;  /* Process 4 k-mers at a time */
    
    /* NEON-optimized batch processing for k-mer expansion */
    for (i = 0; i < simd_batch; i += 4)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 4 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers */
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}
/* SVE optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+sve")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_sve(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    int simd_batch;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with SVE-optimized expansion (8 k-mers at a time, conservative) */
    simd_batch = max_kmers & ~7;  /* Process 8 k-mers at a time (conservative) */
    
    /* SVE-optimized batch processing for k-mer expansion */
    for (i = 0; i < simd_batch; i += 8)
    {
        /* Process each k-mer in the batch */
        for (int j = 0; j < 8 && (i + j) <= seq_bases - k; j++)
        {
            int pos = i + j;
            VarBit **expanded_kmers;
            int expansion_count;
            int exp_j;
            
            /* Expand DNA4 k-mer to DNA2 k-mers */
            expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, pos, k, &expansion_count);
            
            if (!expanded_kmers || expansion_count == 0)
                continue;
            
            /* Process each expanded k-mer */
            for (exp_j = 0; exp_j < expansion_count; exp_j++)
            {
                VarBit *dna2_kmer = expanded_kmers[exp_j];
                
                /* Add kmer2 key directly (without occurrence count) */
                if (dna2_kmer)
                    keys[key_count++] = PointerGetDatum(dna2_kmer);
            }
            
            /* Free only the array, not the individual kmers since we're using them */
            if (expanded_kmers)
                pfree(expanded_kmers);
        }
    }
    
    /* Handle remaining k-mers with scalar processing */
    for (i = simd_batch; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}

/* SVE2 optimized version of kmersearch_extract_dna4_kmer2_with_expansion_direct */
__attribute__((target("+sve2")))
static Datum *
kmersearch_extract_dna4_kmer2_with_expansion_direct_sve2(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    Datum *keys;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0)
        return NULL;
    
    /* Allocate keys array with room for expansions */
    keys = (Datum *) palloc(max_kmers * 10 * sizeof(Datum));  /* Max 10 expansions */
    
    /* Process k-mers with SVE2 optimized expansion */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit **expanded_kmers;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer to DNA2 k-mers */
        expanded_kmers = kmersearch_expand_dna4_kmer2_to_dna2_direct(seq, i, k, &expansion_count);
        
        if (!expanded_kmers || expansion_count == 0)
            continue;
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++)
        {
            VarBit *dna2_kmer = expanded_kmers[j];
            
            /* Add kmer2 key directly (without occurrence count) */
            if (dna2_kmer)
                keys[key_count++] = PointerGetDatum(dna2_kmer);
        }
        
        /* Free only the array, not the individual kmers since we're using them */
        if (expanded_kmers)
            pfree(expanded_kmers);
    }
    
    *nkeys = key_count;
    return keys;
}
#endif /* __aarch64__ */

/*
 * Extract k-mers directly from DNA4 bit sequence with degenerate expansion (with SIMD dispatch)
 */
Datum *
kmersearch_extract_dna4_ngram_key2_with_expansion_direct(VarBit *seq, int k, int *nkeys)
{
    Datum *kmer2_keys;
    int kmer2_count;
    Datum *ngram_keys;
    int ngram_count = 0;
    int i;
    KmerOccurrence *occurrences;
    int occurrence_count = 0;
    
    /* First, extract all kmer2 keys without occurrence count */
    kmer2_keys = kmersearch_extract_dna4_kmer2_with_expansion_direct(seq, k, &kmer2_count);
    if (!kmer2_keys || kmer2_count == 0) {
        *nkeys = 0;
        return NULL;
    }
    
    /* Allocate arrays for occurrence tracking and result */
    occurrences = (KmerOccurrence *) palloc(kmer2_count * sizeof(KmerOccurrence));
    ngram_keys = (Datum *) palloc(kmer2_count * sizeof(Datum));
    
    /* Process each kmer2 key and add occurrence count */
    for (i = 0; i < kmer2_count; i++) {
        VarBit *kmer2_key = (VarBit *) DatumGetPointer(kmer2_keys[i]);
        uint64_t kmer_value;
        int current_count;
        VarBit *ngram_key;
        
        /* Extract k-mer hash value from kmer2 key */
        kmer_value = kmersearch_get_kmer_hash(kmer2_key, 0, k);
        
        /* Find or add occurrence count using binary search */
        current_count = kmersearch_find_or_add_kmer_occurrence(occurrences, &occurrence_count, 
                                                              kmer_value, kmer2_count);
        
        if (current_count < 0)
            continue;  /* Array full, skip */
        
        /* Skip if occurrence exceeds bit limit */
        if (current_count > (1 << kmersearch_occur_bitlen))
            continue;
        
        /* Create n-gram key (k-mer + occurrence count) from kmer2 key */
        ngram_key = kmersearch_create_ngram_key2_from_dna2_bits(kmer2_key, 0, k, current_count);
        if (ngram_key == NULL)
            continue;  /* Skip if key creation failed */
            
        ngram_keys[ngram_count++] = PointerGetDatum(ngram_key);
    }
    
    /* Cleanup */
    for (i = 0; i < kmer2_count; i++) {
        pfree(DatumGetPointer(kmer2_keys[i]));
    }
    pfree(kmer2_keys);
    pfree(occurrences);
    
    *nkeys = ngram_count;
    return ngram_keys;
}

/*
 * SIMD-optimized implementations for k-mer matching
 */

#ifdef __x86_64__
/*
 * AVX2 optimized version of k-mer matching using hash table
 */
__attribute__((target("avx2")))
static int
kmersearch_count_matching_kmer_fast_avx2(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
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
    
    /* Check each sequence k-mer against hash table */
    for (i = 0; i < seq_nkeys; i++)
    {
        if (seq_keys[i] == NULL) {
            elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: NULL seq key at index %d", i);
            continue;
        }
        
        if (VARBITBYTES(seq_keys[i]) != VARBITBYTES(query_keys[0])) {
            elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Size mismatch seq[%d]=%zu vs query[0]=%zu", 
                 i, (size_t)VARBITBYTES(seq_keys[i]), (size_t)VARBITBYTES(query_keys[0]));
            continue;
        }
        
        if (hash_search(query_hash, VARBITS(seq_keys[i]), HASH_FIND, NULL))
        {
            match_count++;
        }
    }
    
    hash_destroy(query_hash);
    
    elog(LOG, "kmersearch_count_matching_kmer_fast_avx2: Returning match_count=%d", match_count);
    return match_count;
}

/*
 * AVX512 optimized version of k-mer matching using hash table
 */
__attribute__((target("avx512f,avx512bw")))
static int
kmersearch_count_matching_kmer_fast_avx512(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
    HTAB *query_hash;
    HASHCTL hash_ctl;
    bool found;
    
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
#endif /* __x86_64__ */

#ifdef __aarch64__
/*
 * NEON optimized version of k-mer matching using hash table
 */
static int
kmersearch_count_matching_kmer_fast_neon(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int match_count = 0;
    int i;
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
 * SVE optimized version of k-mer matching using hash table
 */
__attribute__((target("+sve")))
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
#endif /* __aarch64__ */

/*
 * Count matching k-mers between two key arrays (with SIMD dispatch)
 */
int
kmersearch_count_matching_kmer_fast(VarBit **seq_keys, int seq_nkeys, VarBit **query_keys, int query_nkeys)
{
    int key_combinations;
    
    /* Validate input parameters */
    if (!seq_keys || !query_keys) {
        return 0;
    }
    
    if (seq_nkeys <= 0 || query_nkeys <= 0) {
        return 0;
    }
    
    key_combinations = seq_nkeys * query_nkeys;
    
    /* For small datasets, O(n*m) might be faster than hash table overhead */
    if (key_combinations < 100) {
        return kmersearch_count_matching_kmer_fast_scalar_simple(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && key_combinations >= SIMD_KEYCOMB_AVX512_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_avx512(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    if (simd_capability >= SIMD_AVX2 && key_combinations >= SIMD_KEYCOMB_AVX2_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_avx2(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && key_combinations >= SIMD_KEYCOMB_SVE_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_sve(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
    if (simd_capability >= SIMD_NEON && key_combinations >= SIMD_KEYCOMB_NEON_THRESHOLD) {
        return kmersearch_count_matching_kmer_fast_neon(seq_keys, seq_nkeys, query_keys, query_nkeys);
    }
#endif
    return kmersearch_count_matching_kmer_fast_scalar_hashtable(seq_keys, seq_nkeys, query_keys, query_nkeys);
}

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