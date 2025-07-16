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
    int occur_bits = kmersearch_occur_bitlen;
    int total_bits = kmer_bits + occur_bits;
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
    if (adj_occurrence >= (1 << occur_bits))
        adj_occurrence = (1 << occur_bits) - 1;  /* Cap at max value */
    
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
    for (i = 0; i < occur_bits; i++)
    {
        int occur_bit = (occurrence_count >> (occur_bits - 1 - i)) & 1;
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
 * Uses direct bit extraction for k <= 32, PostgreSQL hash for k > 32
 */
uint64_t
kmersearch_get_kmer_hash(VarBit *seq, int start_pos, int k)
{
    bits8 *src_data = VARBITS(seq);
    int src_bytes = VARBITBYTES(seq);
    
    /* For k <= 32, use direct bit extraction for performance and deterministic values */
    if (k <= 32) {
        uint64_t kmer_value = 0;
        int j;
        
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
    
    /* For k > 32, use PostgreSQL's hash function */
    else {
        int kmer_bits = k * 2;
        int start_bit = start_pos * 2;
        int start_byte = start_bit / 8;
        int start_bit_offset = start_bit % 8;
        int kmer_bytes = (kmer_bits + start_bit_offset + 7) / 8;
        
        /* Boundary check to prevent buffer overflow */
        if (start_byte + kmer_bytes > src_bytes) {
            return 0;  /* Invalid k-mer */
        }
        
        /* Hash the k-mer data starting from the appropriate byte position */
        return hash_any_extended(src_data + start_byte, kmer_bytes, start_bit_offset);
    }
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
    simd_dispatch.dna4_encode(query, (uint8_t*)data_ptr, query_len);
    
    /* Extract k-mers using SIMD optimized function */
    datum_keys = kmersearch_extract_dna4_kmer2_with_expansion_direct(dna4_seq, k, nkeys);
    
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

/*
 * Extract k-mers from query string with degenerate code expansion
 */
VarBit **
kmersearch_extract_query_kmer_with_degenerate(const char *query, int k, int *nkeys)
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
                VarBit *kmer_key = kmersearch_create_kmer2_key_only(expanded[j], k);
                keys[key_count++] = kmer_key;
                pfree(expanded[j]);
            }
        }
        else
        {
            /* Simple case - no degenerate codes */
            VarBit *kmer_key = kmersearch_create_kmer2_key_only(kmer, k);
            keys[key_count++] = kmer_key;
        }
    }
    
    *nkeys = key_count;
    return keys;
}


/*
 * Convert VarBit to hexadecimal string representation
 */
char *
kmersearch_varbit_to_hex_string(VarBit *varbit)
{
    int bitlen = VARBITLEN(varbit);
    int bytelen = VARBITBYTES(varbit);
    unsigned char *bits = VARBITS(varbit);
    char *hex_string;
    int i;
    
    hex_string = palloc(bytelen * 2 + 1);
    
    for (i = 0; i < bytelen; i++) {
        sprintf(hex_string + i * 2, "%02x", bits[i]);
    }
    
    hex_string[bytelen * 2] = '\0';
    return hex_string;
}

/* Note: kmersearch_extract_dna2_kmer2_direct kept in kmersearch.c due to SIMD dependencies */

/*
 * Extract k-mers from DNA2 sequence without occurrence count (internal function)
 */
Datum *
kmersearch_extract_dna2_kmer2_only(VarBit *seq, int k, int *nkeys)
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
    
    /* Extract k-mers */
    for (i = 0; i <= seq_bases - k; i++)
    {
        VarBit *kmer_key;
        
        /* Create k-mer key (without occurrence count) */
        kmer_key = kmersearch_create_kmer2_key_from_dna2_bits(seq, i, k);
        if (kmer_key == NULL)
            continue;  /* Skip if key creation failed */
            
        keys[key_count++] = PointerGetDatum(kmer_key);
    }
    
    *nkeys = key_count;
    return keys;
}

/*
 * Encode k-mer-only VarBit into compact KmerData (ignoring occurrence count bits)
 */
KmerData
kmersearch_encode_kmer2_only_data(VarBit *kmer, int k_size)
{
    KmerData result;
    unsigned char *bits;
    int i, bit_offset;
    
    memset(&result, 0, sizeof(KmerData));
    bits = VARBITS(kmer);
    
    if (k_size <= 8) {
        result.k8_data = 0;
        for (i = 0; i < k_size; i++) {
            int byte_pos, bit_in_byte, nucleotide;
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k8_data |= (nucleotide << (2 * (k_size - 1 - i)));
        }
    } else if (k_size <= 16) {
        result.k16_data = 0;
        for (i = 0; i < k_size; i++) {
            int byte_pos, bit_in_byte, nucleotide;
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k16_data |= (nucleotide << (2 * (k_size - 1 - i)));
        }
    } else if (k_size <= 32) {
        result.k32_data = 0;
        for (i = 0; i < k_size; i++) {
            int byte_pos, bit_in_byte, nucleotide;
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k32_data |= ((uint64)nucleotide << (2 * (k_size - 1 - i)));
        }
    } else {
        /* For k > 32, split across high and low 64-bit values */
        result.k64_data.high = 0;
        result.k64_data.low = 0;
        for (i = 0; i < k_size; i++) {
            int byte_pos, bit_in_byte, nucleotide;
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            
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
 * Encode VarBit k-mer into compact KmerData
 */
KmerData
kmersearch_encode_kmer_data(VarBit *kmer, int k_size)
{
    KmerData result;
    unsigned char *bits;
    int i, bit_offset;
    int byte_pos, bit_in_byte, nucleotide;
    
    memset(&result, 0, sizeof(KmerData));
    bits = VARBITS(kmer);
    
    if (k_size <= 8) {
        result.k8_data = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k8_data |= (nucleotide << (2 * (k_size - 1 - i)));
        }
    } else if (k_size <= 16) {
        result.k16_data = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k16_data |= (nucleotide << (2 * (k_size - 1 - i)));
        }
    } else if (k_size <= 32) {
        result.k32_data = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            result.k32_data |= ((uint64)nucleotide << (2 * (k_size - 1 - i)));
        }
    } else {
        /* For k > 32, split across high and low 64-bit values */
        result.k64_data.high = 0;
        result.k64_data.low = 0;
        for (i = 0; i < k_size; i++) {
            bit_offset = i * 2;
            byte_pos = bit_offset / 8;
            bit_in_byte = bit_offset % 8;
            nucleotide = (bits[byte_pos] >> (6 - bit_in_byte)) & 0x3;
            
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
 * Extract k-mers from DNA sequence and create n-gram keys
 */
Datum *
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
                
                ngram_key = kmersearch_create_ngram_key2(expanded[j], k, occurrence);
                keys[key_count++] = PointerGetDatum(ngram_key);
                
                pfree(expanded[j]);
            }
        }
        else
        {
            /* Simple case - no degenerate codes */
            int occurrence = 1;
            VarBit *ngram_key = kmersearch_create_ngram_key2(kmer, k, occurrence);
            keys[key_count++] = PointerGetDatum(ngram_key);
        }
    }
    
    *nkeys = key_count;
    return keys;
}

