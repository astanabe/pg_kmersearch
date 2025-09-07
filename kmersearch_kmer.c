#include "kmersearch.h"

/*
 * kmersearch_kmer.c
 * Simple k-mer utility functions for the pg_kmersearch extension.
 */
/* DNA4 to DNA2 expansion table: [expansion_count, base1, base2, base3, base4] */
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

/* Direct text character to DNA2 expansion table: [expansion_count, base1, base2, base3, base4]
 * DNA2 encoding: A=0, C=1, G=2, T=3 */
static const uint8 kmersearch_text_to_dna2_table[256][5] = {
    /* Initialize all entries to invalid (0 expansions) */
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, /* 64 */
    /* A=65 */ {1, 0, 0, 0, 0},
    /* B=66 */ {3, 1, 2, 3, 0},  /* C,G,T */
    /* C=67 */ {1, 1, 0, 0, 0},
    /* D=68 */ {3, 0, 2, 3, 0},  /* A,G,T */
    {0, 0, 0, 0, 0}, /* E=69 */
    {0, 0, 0, 0, 0}, /* F=70 */
    /* G=71 */ {1, 2, 0, 0, 0},
    /* H=72 */ {3, 0, 1, 3, 0},  /* A,C,T */
    {0, 0, 0, 0, 0}, /* I=73 */
    {0, 0, 0, 0, 0}, /* J=74 */
    /* K=75 */ {2, 2, 3, 0, 0},  /* G,T */
    {0, 0, 0, 0, 0}, /* L=76 */
    /* M=77 */ {2, 0, 1, 0, 0},  /* A,C */
    /* N=78 */ {4, 0, 1, 2, 3},  /* A,C,G,T */
    {0, 0, 0, 0, 0}, /* O=79 */
    {0, 0, 0, 0, 0}, /* P=80 */
    {0, 0, 0, 0, 0}, /* Q=81 */
    /* R=82 */ {2, 0, 2, 0, 0},  /* A,G */
    /* S=83 */ {2, 1, 2, 0, 0},  /* C,G */
    /* T=84 */ {1, 3, 0, 0, 0},
    /* U=85 */ {1, 3, 0, 0, 0},  /* U as T */
    /* V=86 */ {3, 0, 1, 2, 0},  /* A,C,G */
    /* W=87 */ {2, 0, 3, 0, 0},  /* A,T */
    {0, 0, 0, 0, 0}, /* X=88 */
    /* Y=89 */ {2, 1, 3, 0, 0},  /* C,T */
    {0, 0, 0, 0, 0}, /* Z=90 */
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    /* a=97 */ {1, 0, 0, 0, 0},
    /* b=98 */ {3, 1, 2, 3, 0},  /* C,G,T */
    /* c=99 */ {1, 1, 0, 0, 0},
    /* d=100 */ {3, 0, 2, 3, 0},  /* A,G,T */
    {0, 0, 0, 0, 0}, /* e=101 */
    {0, 0, 0, 0, 0}, /* f=102 */
    /* g=103 */ {1, 2, 0, 0, 0},
    /* h=104 */ {3, 0, 1, 3, 0},  /* A,C,T */
    {0, 0, 0, 0, 0}, /* i=105 */
    {0, 0, 0, 0, 0}, /* j=106 */
    /* k=107 */ {2, 2, 3, 0, 0},  /* G,T */
    {0, 0, 0, 0, 0}, /* l=108 */
    /* m=109 */ {2, 0, 1, 0, 0},  /* A,C */
    /* n=110 */ {4, 0, 1, 2, 3},  /* A,C,G,T */
    {0, 0, 0, 0, 0}, /* o=111 */
    {0, 0, 0, 0, 0}, /* p=112 */
    {0, 0, 0, 0, 0}, /* q=113 */
    /* r=114 */ {2, 0, 2, 0, 0},  /* A,G */
    /* s=115 */ {2, 1, 2, 0, 0},  /* C,G */
    /* t=116 */ {1, 3, 0, 0, 0},
    /* u=117 */ {1, 3, 0, 0, 0},  /* U as T */
    /* v=118 */ {3, 0, 1, 2, 0},  /* A,C,G */
    /* w=119 */ {2, 0, 3, 0, 0},  /* A,T */
    {0, 0, 0, 0, 0}, /* x=120 */
    /* y=121 */ {2, 1, 3, 0, 0},  /* C,T */
    {0, 0, 0, 0, 0}, /* z=122 */
    /* Rest of the table initialized to {0, 0, 0, 0, 0} */
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}
};

static int kmersearch_count_matching_uintkey_scalar(void *seq_keys, int seq_nkeys, void *query_keys, int query_nkeys, int k_size);

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
 * Extract k-mers from query string as uintkey format (with occurrence counts)
 * This function generates uintkey compatible with sequence extraction using SIMD optimizations
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
 * Extract k-mers directly from DNA2 bit sequence (with SIMD dispatch)
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
    void *occurrences = NULL;
    int occurrence_count = 0;
    
    /* Validate parameters */
    if (k < 4 || k > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid k-mer size: %d", k),
                 errdetail("k-mer size must be between 4 and 32")));
    }
    
    if (occur_bitlen < 0 || occur_bitlen > 16) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid occurrence bit length: %d", occur_bitlen),
                 errdetail("Occurrence bit length must be between 0 and 16")));
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
            if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) continue;
            
            if (occur_bitlen == 0) {
                /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                if (current_count == 1) {
                    final_value = kmer_value;
                    ((uint16 *)result)[result_count++] = final_value;
                }
            } else {
                /* Include occurrence count in the uintkey */
                final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                ((uint16 *)result)[result_count++] = final_value;
            }
            
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
            if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) continue;
            
            if (occur_bitlen == 0) {
                /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                if (current_count == 1) {
                    final_value = kmer_value;
                    ((uint32 *)result)[result_count++] = final_value;
                }
            } else {
                /* Include occurrence count in the uintkey */
                final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                ((uint32 *)result)[result_count++] = final_value;
            }
            
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
            if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) continue;
            
            if (occur_bitlen == 0) {
                /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                if (current_count == 1) {
                    final_value = kmer_value;
                    ((uint64 *)result)[result_count++] = final_value;
                }
            } else {
                /* Include occurrence count in the uintkey */
                final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                ((uint64 *)result)[result_count++] = final_value;
            }
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
    
    /* SIMD capability check */
    
    /* SIMD implementations can be added later based on capability and threshold */
#ifdef __x86_64__
    /* Future: Add AVX2/AVX512 dispatch based on simd_capability */
    if (simd_capability >= SIMD_AVX2) {
        /* AVX2 available but not yet implemented for DNA2 extraction, using scalar */
    }
#elif defined(__aarch64__)
    /* Future: Add NEON/SVE dispatch based on simd_capability */
    if (simd_capability >= SIMD_NEON) {
        /* NEON available but not yet implemented for DNA2 extraction, using scalar */
    }
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
    void *occurrences = NULL;
    int occurrence_count = 0;
    
    /* Validate parameters */
    if (k < 4 || k > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid k-mer size: %d", k),
                 errdetail("k-mer size must be between 4 and 32")));
    }
    
    if (occur_bitlen < 0 || occur_bitlen > 16) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid occurrence bit length: %d", occur_bitlen),
                 errdetail("Occurrence bit length must be between 0 and 16")));
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
    
    if (elem_size == sizeof(uint16)) {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence16));
    } else if (elem_size == sizeof(uint32)) {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence32));
    } else {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence64));
    }
    
    /* Process each k-mer position */
    for (i = 0; i <= seq_len - k; i++) {
        void *expanded_uintkeys;
        int expansion_count;
        int j;
        
        /* Expand DNA4 k-mer directly to uintkey format */
        kmersearch_expand_dna4_to_uintkey(seq, i, k, &expanded_uintkeys, &expansion_count, elem_size);
        
        if (!expanded_uintkeys || expansion_count == 0) {
            continue;
        }
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++) {
            int current_count;
            
            /* Process directly as appropriate uint type */
            if (elem_size == sizeof(uint16)) {
                uint16 kmer_value = ((uint16 *)expanded_uintkeys)[j];
                uint16 final_value;
                
                current_count = kmersearch_find_or_add_kmer_occurrence16(
                    (KmerOccurrence16 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    continue;
                }
                
                if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) {
                    continue;
                }
                
                if (occur_bitlen == 0) {
                    /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                    if (current_count == 1) {
                        final_value = kmer_value;
                        ((uint16 *)result)[result_count++] = final_value;
                    }
                } else {
                    /* Include occurrence count in the uintkey */
                    final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                    ((uint16 *)result)[result_count++] = final_value;
                }
                
            } else if (elem_size == sizeof(uint32)) {
                uint32 kmer_value = ((uint32 *)expanded_uintkeys)[j];
                uint32 final_value;
                
                current_count = kmersearch_find_or_add_kmer_occurrence32(
                    (KmerOccurrence32 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    continue;
                }
                
                if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) {
                    continue;
                }
                
                if (occur_bitlen == 0) {
                    /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                    if (current_count == 1) {
                        final_value = kmer_value;
                        ((uint32 *)result)[result_count++] = final_value;
                    }
                } else {
                    /* Include occurrence count in the uintkey */
                    final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                    ((uint32 *)result)[result_count++] = final_value;
                }
                
            } else {
                uint64 kmer_value = ((uint64 *)expanded_uintkeys)[j];
                uint64 final_value;
                
                current_count = kmersearch_find_or_add_kmer_occurrence64(
                    (KmerOccurrence64 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    continue;
                }
                
                if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) {
                    continue;
                }
                
                if (occur_bitlen == 0) {
                    /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                    if (current_count == 1) {
                        final_value = kmer_value;
                        ((uint64 *)result)[result_count++] = final_value;
                    }
                } else {
                    /* Include occurrence count in the uintkey */
                    final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                    ((uint64 *)result)[result_count++] = final_value;
                }
            }
        }
        
        /* Free the expansion array */
        if (expanded_uintkeys) {
            pfree(expanded_uintkeys);
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
    
    /* SIMD capability check */
    
    /* SIMD implementations can be added later based on capability and threshold */
#ifdef __x86_64__
    /* Future: Add AVX2/AVX512 dispatch based on simd_capability */
    if (simd_capability >= SIMD_AVX2) {
        /* AVX2 available but not yet implemented for DNA4 extraction, using scalar */
    }
#elif defined(__aarch64__)
    /* Future: Add NEON/SVE dispatch based on simd_capability */
    if (simd_capability >= SIMD_NEON) {
        /* NEON available but not yet implemented for DNA4 extraction, using scalar */
    }
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
    int total_bits;
    
    if (seq_nkeys == 0 || query_nkeys == 0)
        return 0;
    
    total_bits = k_size * 2 + kmersearch_occur_bitlen;
    
    /* Use hash table for O(n+m) performance */
    memset(&hash_ctl, 0, sizeof(hash_ctl));
    
    if (total_bits <= 16) {
        hash_ctl.keysize = sizeof(uint16);
        hash_ctl.entrysize = sizeof(uint16);
    }
    else if (total_bits <= 32) {
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
    if (total_bits <= 16) {
        uint16 *query = (uint16 *)query_keys;
        uint16 *seq;
        
        for (i = 0; i < query_nkeys; i++) {
            hash_search(query_hash, &query[i], HASH_ENTER, &found);
        }
        
        /* Check each sequence uintkey against hash table */
        seq = (uint16 *)seq_keys;
        for (i = 0; i < seq_nkeys; i++) {
            if (hash_search(query_hash, &seq[i], HASH_FIND, NULL)) {
                shared_count++;
            }
        }
    }
    else if (total_bits <= 32) {
        uint32 *query = (uint32 *)query_keys;
        uint32 *seq;
        
        for (i = 0; i < query_nkeys; i++) {
            hash_search(query_hash, &query[i], HASH_ENTER, &found);
        }
        
        seq = (uint32 *)seq_keys;
        for (i = 0; i < seq_nkeys; i++) {
            if (hash_search(query_hash, &seq[i], HASH_FIND, NULL)) {
                shared_count++;
            }
        }
    }
    else {
        uint64 *query = (uint64 *)query_keys;
        uint64 *seq;
        
        for (i = 0; i < query_nkeys; i++) {
            hash_search(query_hash, &query[i], HASH_ENTER, &found);
        }
        
        seq = (uint64 *)seq_keys;
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
 * Helper function to check if text will exceed degenerate limit
 * 
 * Rules (matching kmersearch_will_exceed_degenerate_limit_dna4_bits):
 * - If only 1 MRWSYKVHDB character: expand (return false)
 * - If 2 or more MRWSYKVHDB characters: don't expand (return true)
 * - If 1 or more N characters: don't expand (return true)
 * 
 * This ensures expansion limit of 3 (max expansion from single V/H/D/B)
 */
static bool
kmersearch_text_will_exceed_degenerate_limit(const char *text, int start_pos, int k)
{
    int degenerate_count = 0;  /* Count of MRWSYKVHDB characters */
    int i;
    
    for (i = 0; i < k && start_pos + i < strlen(text); i++)
    {
        char c = toupper(text[start_pos + i]);
        
        /* N - always exceed limit (would expand to 4) */
        if (c == 'N')
        {
            return true;
        }
        /* MRWSYKVHDB - allow only 1 */
        else if (c == 'M' || c == 'R' || c == 'W' || c == 'S' || c == 'Y' || c == 'K' ||
                 c == 'V' || c == 'H' || c == 'D' || c == 'B')
        {
            degenerate_count++;
            /* 2 or more degenerate characters - exceed limit */
            if (degenerate_count >= 2)
                return true;
        }
    }
    
    /* Only 0 or 1 degenerate character - within limit */
    return false;
}

/*
 * Expand text k-mer directly to uint keys
 */
static void
kmersearch_expand_text_to_uintkey(const char *text, int start_pos, int k, 
                                  void **output, int *expansion_count, size_t elem_size)
{
    uint8 base_expansions[32][4];  /* Max k=32, max 4 expansions per base */
    int base_counts[32];
    int total_combinations = 1;
    void *results;
    int i, combo;
    
    *expansion_count = 0;
    *output = NULL;
    
    /* Check if expansion will exceed limit */
    if (kmersearch_text_will_exceed_degenerate_limit(text, start_pos, k)) {
        return;
    }
    
    /* Extract expansion info for each base */
    for (i = 0; i < k; i++)
    {
        unsigned char c = text[start_pos + i];
        const uint8 *expansion = kmersearch_text_to_dna2_table[c];
        int exp_count = expansion[0];
        int j;
        
        if (exp_count == 0) {
            /* Invalid character */
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("Invalid character '%c' for DNA sequence", c),
                     errhint("DNA accepts A,C,G,T,U,M,R,W,S,Y,K,V,H,D,B,N characters")));
        }
        
        base_counts[i] = exp_count;
        for (j = 0; j < exp_count; j++)
        {
            base_expansions[i][j] = expansion[j + 1];
        }
        
        total_combinations *= exp_count;
    }
    
    /* Allocate result array */
    results = palloc(total_combinations * elem_size);
    
    /* Generate all combinations directly as uint values */
    for (combo = 0; combo < total_combinations; combo++)
    {
        int temp_combo = combo;
        uint64 kmer_value = 0;
        
        /* Generate this combination directly as uint */
        for (i = 0; i < k; i++)
        {
            int base_idx = temp_combo % base_counts[i];
            uint8 dna2_base = base_expansions[i][base_idx];
            kmer_value = (kmer_value << 2) | dna2_base;
            temp_combo /= base_counts[i];
        }
        
        /* Store the result in appropriate size */
        if (elem_size == sizeof(uint16))
        {
            ((uint16 *)results)[combo] = (uint16)kmer_value;
        }
        else if (elem_size == sizeof(uint32))
        {
            ((uint32 *)results)[combo] = (uint32)kmer_value;
        }
        else  /* elem_size == sizeof(uint64) */
        {
            ((uint64 *)results)[combo] = kmer_value;
        }
    }
    
    *output = results;
    *expansion_count = total_combinations;
}

/*
 * Extract uint keys from text string containing DNA sequence
 * Directly processes text without DNA4 intermediate representation
 */
void
kmersearch_extract_uintkey_from_text(const char *text, void **output, int *nkeys)
{
    int text_len = strlen(text);
    int k = kmersearch_kmer_size;
    int occur_bitlen = kmersearch_occur_bitlen;
    int max_kmers = text_len - k + 1;
    int kmer_bits = k * 2;  /* DNA2 format: 2 bits per base */
    int total_bits = kmer_bits + occur_bitlen;
    size_t elem_size;
    void *result;
    int result_count = 0;
    int result_capacity;
    int i;
    void *occurrences = NULL;
    int occurrence_count = 0;
    
    /* Validate parameters */
    if (k < 4 || k > 32) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid k-mer size: %d", k),
                 errdetail("k-mer size must be between 4 and 32")));
    }
    
    if (occur_bitlen < 0 || occur_bitlen > 16) {
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("Invalid occurrence bit length: %d", occur_bitlen),
                 errdetail("Occurrence bit length must be between 0 and 16")));
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
        *output = NULL;
        *nkeys = 0;
        return;
    }
    
    /* Validate all input characters first */
    for (i = 0; i < text_len; i++) {
        unsigned char c = text[i];
        if (kmersearch_text_to_dna2_table[c][0] == 0) {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("Invalid character '%c' for DNA sequence", c),
                     errhint("DNA accepts A,C,G,T,U,M,R,W,S,Y,K,V,H,D,B,N characters")));
        }
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
    
    if (elem_size == sizeof(uint16)) {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence16));
    } else if (elem_size == sizeof(uint32)) {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence32));
    } else {
        occurrences = palloc0(result_capacity * sizeof(KmerOccurrence64));
    }
    
    /* Process each k-mer position */
    for (i = 0; i <= text_len - k; i++) {
        void *expanded_uintkeys;
        int expansion_count;
        int j;
        
        /* Expand text k-mer directly to uintkey format */
        kmersearch_expand_text_to_uintkey(text, i, k, &expanded_uintkeys, &expansion_count, elem_size);
        
        if (!expanded_uintkeys || expansion_count == 0) {
            continue;
        }
        
        /* Process each expanded k-mer */
        for (j = 0; j < expansion_count; j++) {
            int current_count;
            
            /* Process directly as appropriate uint type */
            if (elem_size == sizeof(uint16)) {
                uint16 kmer_value = ((uint16 *)expanded_uintkeys)[j];
                uint16 final_value;
                
                current_count = kmersearch_find_or_add_kmer_occurrence16(
                    (KmerOccurrence16 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    continue;
                }
                
                if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) {
                    continue;
                }
                
                if (occur_bitlen == 0) {
                    /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                    if (current_count == 1) {
                        final_value = kmer_value;
                        ((uint16 *)result)[result_count++] = final_value;
                    }
                } else {
                    /* Include occurrence count in the uintkey */
                    final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                    ((uint16 *)result)[result_count++] = final_value;
                }
                
            } else if (elem_size == sizeof(uint32)) {
                uint32 kmer_value = ((uint32 *)expanded_uintkeys)[j];
                uint32 final_value;
                
                current_count = kmersearch_find_or_add_kmer_occurrence32(
                    (KmerOccurrence32 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    continue;
                }
                
                if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) {
                    continue;
                }
                
                if (occur_bitlen == 0) {
                    /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                    if (current_count == 1) {
                        final_value = kmer_value;
                        ((uint32 *)result)[result_count++] = final_value;
                    }
                } else {
                    /* Include occurrence count in the uintkey */
                    final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                    ((uint32 *)result)[result_count++] = final_value;
                }
                
            } else {
                uint64 kmer_value = ((uint64 *)expanded_uintkeys)[j];
                uint64 final_value;
                
                current_count = kmersearch_find_or_add_kmer_occurrence64(
                    (KmerOccurrence64 *)occurrences, &occurrence_count,
                    kmer_value, result_capacity);
                if (current_count < 0) {
                    continue;
                }
                
                if (occur_bitlen > 0 && current_count > (1 << occur_bitlen)) {
                    continue;
                }
                
                if (occur_bitlen == 0) {
                    /* For occur_bitlen=0, only output unique k-mers (first occurrence) */
                    if (current_count == 1) {
                        final_value = kmer_value;
                        ((uint64 *)result)[result_count++] = final_value;
                    }
                } else {
                    /* Include occurrence count in the uintkey */
                    final_value = (kmer_value << occur_bitlen) | ((current_count - 1) & ((1 << occur_bitlen) - 1));
                    ((uint64 *)result)[result_count++] = final_value;
                }
            }
        }
        
        /* Free the expansion array */
        if (expanded_uintkeys) {
            pfree(expanded_uintkeys);
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