# Implementation Plan: K-mer Extraction with Native uint Output (REVISED)

## Overview

This document outlines the implementation plan for new k-mer extraction functions that output native unsigned integer types instead of VarBit structures, providing significant performance improvements for high-frequency k-mer analysis.

**CRITICAL PERFORMANCE REQUIREMENT**: These functions must avoid creating intermediate VarBit k-mer structures entirely. All k-mer extraction must be performed through direct bit manipulation on the input VarBit sequence data to achieve the target performance gains.

## Functions to Implement

### Primary Functions

1. **`kmersearch_extract_dna2_kmer2_as_uint_direct(VarBit *seq, int k, void **output, int *nkeys)`**
   - Extracts k-mers from DNA2 sequences as native uint types
   - Output type depends on k-mer size: uint16 (k≤8), uint32 (k≤16), uint64 (k≤32)
   - **CRITICAL**: Must extract k-mers directly from VarBit data without creating intermediate VarBit k-mer objects

2. **`kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(VarBit *seq, int k, void **output, int *nkeys)`**
   - Extracts k-mers from DNA4 sequences with degenerate expansion as native uint types
   - Output type depends on k-mer size: uint16 (k≤8), uint32 (k≤16), uint64 (k≤32)
   - **CRITICAL**: Must perform degenerate expansion directly to uint values without creating intermediate VarBit k-mer objects

## Data Structures

### Enhanced KmerData Union

The existing `KmerData` union (kmersearch.h:270-275) will be leveraged:

```c
typedef union KmerData
{
    uint16      k8_data;                 /* k <= 8: 16 bits */
    uint32      k16_data;                /* k <= 16: 32 bits */  
    uint64      k32_data;                /* k <= 32: 64 bits */
} KmerData;
```

### New Result Structure

```c
typedef struct KmerUintResult
{
    void        *data;                   /* Pointer to uint16/uint32/uint64 array */
    int         count;                   /* Number of k-mers extracted */
    size_t      element_size;            /* Size of each element (2, 4, or 8 bytes) */
    int         k_size;                  /* K-mer size used */
} KmerUintResult;
```

## Implementation Architecture

### SIMD Dispatch System

Following the existing pattern in kmersearch.c:858-873, both functions will use runtime SIMD capability detection:

```c
/* DNA2 function dispatch */
void *kmersearch_extract_dna2_kmer2_as_uint_direct(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_as_uint_direct_avx512(seq, k, output, nkeys);
    }
    if (simd_capability >= SIMD_AVX2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_as_uint_direct_avx2(seq, k, output, nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_as_uint_direct_sve(seq, k, output, nkeys);
    }
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_as_uint_direct_neon(seq, k, output, nkeys);
    }
#endif
    return kmersearch_extract_dna2_kmer2_as_uint_direct_scalar(seq, k, output, nkeys);
}
```

### SIMD Implementation Variants

#### x86_64 Implementations
- **AVX512**: 16 k-mers parallel processing
- **AVX2**: 8 k-mers parallel processing  
- **Scalar**: Fallback implementation

#### aarch64 Implementations
- **SVE**: Variable vector width (8+ k-mers)
- **NEON**: 4 k-mers parallel processing
- **Scalar**: Fallback implementation

## Detailed Implementation Plan

### Phase 1: Core Infrastructure

#### 1.1 Header Declarations (kmersearch.h)

Add function declarations:

```c
/* Native uint k-mer extraction functions */
void *kmersearch_extract_dna2_kmer2_as_uint_direct(VarBit *seq, int k, void **output, int *nkeys);
void *kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(VarBit *seq, int k, void **output, int *nkeys);

/* Direct bit manipulation helper functions */
size_t kmersearch_get_kmer_uint_size(int k);
uint64 kmersearch_extract_dna2_kmer_direct_bits(VarBit *seq, int start_pos, int k);
void kmersearch_extract_dna4_kmer_expansions_direct_bits(VarBit *seq, int start_pos, int k, uint64 *output, int *count);
```

#### 1.2 Utility Functions

```c
/* Get appropriate uint size for k-mer length */
size_t kmersearch_get_kmer_uint_size(int k)
{
    if (k <= 8) return sizeof(uint16);
    if (k <= 16) return sizeof(uint32);
    if (k <= 32) return sizeof(uint64);
    ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                   errmsg("k-mer length must be between 4 and 32")));
}

/* Extract DNA2 k-mer directly as uint64 - NO VarBit intermediate */
uint64 kmersearch_extract_dna2_kmer_direct_bits(VarBit *seq, int start_pos, int k)
{
    bits8 *src_data = VARBITS(seq);
    uint64 kmer_value = 0;
    int j;
    
    /* Direct bit extraction from VarBit data */
    for (j = 0; j < k; j++)
    {
        int bit_pos = (start_pos + j) * 2;  /* DNA2: 2 bits per base */
        int byte_pos = bit_pos / 8;
        int bit_offset = bit_pos % 8;
        uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
        kmer_value = (kmer_value << 2) | base_bits;
    }
    
    return kmer_value;
}
```

### Phase 2: Scalar Implementations

#### 2.1 DNA2 Scalar Implementation - DIRECT BIT MANIPULATION

```c
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_scalar(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    bits8 *src_data = VARBITS(seq);
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return NULL;
    }
    
    result = palloc(max_kmers * element_size);
    *output = result;
    
    /* CRITICAL: Direct bit manipulation - NO VarBit intermediate objects */
    if (k <= 8) {
        uint16 *output_array = (uint16 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            int j;
            
            /* Extract k-mer directly from VarBit data */
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint16) kmer_value;
            (*nkeys)++;
        }
    } else if (k <= 16) {
        uint32 *output_array = (uint32 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            int j;
            
            /* Extract k-mer directly from VarBit data */
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = (uint32) kmer_value;
            (*nkeys)++;
        }
    } else {
        uint64 *output_array = (uint64 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            uint64 kmer_value = 0;
            int j;
            
            /* Extract k-mer directly from VarBit data */
            for (j = 0; j < k; j++) {
                int bit_pos = (i + j) * 2;
                int byte_pos = bit_pos / 8;
                int bit_offset = bit_pos % 8;
                uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
                kmer_value = (kmer_value << 2) | base_bits;
            }
            
            output_array[*nkeys] = kmer_value;
            (*nkeys)++;
        }
    }
    
    return result;
}
```

#### 2.2 DNA4 Scalar Implementation - DIRECT DEGENERATE EXPANSION

```c
/* Extract DNA4 base and get all possible expansions directly as uint values */
static void kmersearch_expand_dna4_base_direct(uint8 dna4_base, uint64 *expansions, int *count)
{
    /* DNA4 base expansion mapping - direct to 2-bit values */
    static const uint8 dna4_expansions[][4] = {
        {0},           /* A=0001 -> A(00) */
        {1},           /* C=0010 -> C(01) */
        {0, 1},        /* M=0011 -> A(00), C(01) */
        {2},           /* G=0100 -> G(10) */
        {0, 2},        /* R=0101 -> A(00), G(10) */
        {1, 2},        /* S=0110 -> C(01), G(10) */
        {0, 1, 2},     /* V=0111 -> A(00), C(01), G(10) */
        {3},           /* T=1000 -> T(11) */
        {0, 3},        /* W=1001 -> A(00), T(11) */
        {1, 3},        /* Y=1010 -> C(01), T(11) */
        {0, 1, 3},     /* H=1011 -> A(00), C(01), T(11) */
        {2, 3},        /* K=1100 -> G(10), T(11) */
        {0, 2, 3},     /* D=1101 -> A(00), G(10), T(11) */
        {1, 2, 3},     /* B=1110 -> C(01), G(10), T(11) */
        {0, 1, 2, 3}   /* N=1111 -> A(00), C(01), G(10), T(11) */
    };
    
    static const int dna4_expansion_counts[] = {1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
    
    if (dna4_base == 0 || dna4_base > 15) {
        *count = 0;
        return;
    }
    
    *count = dna4_expansion_counts[dna4_base - 1];
    for (int i = 0; i < *count; i++) {
        expansions[i] = dna4_expansions[dna4_base - 1][i];
    }
}

static void *kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_scalar(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    bits8 *src_data = VARBITS(seq);
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return NULL;
    }
    
    /* Allocate for maximum possible expansions (4^k, but limited to reasonable size) */
    result = palloc(max_kmers * 1024 * element_size);  /* Conservative allocation */
    *output = result;
    
    /* CRITICAL: Direct degenerate expansion - NO VarBit intermediate objects */
    for (i = 0; i <= seq_bases - k; i++) {
        /* Extract DNA4 k-mer and generate all expansions directly */
        uint64 current_expansions[1024];  /* Stack allocation for current k-mer expansions */
        int current_count = 1;
        current_expansions[0] = 0;
        
        /* Build all possible expansions for this k-mer position */
        for (int pos = 0; pos < k; pos++) {
            int bit_pos = (i + pos) * 4;  /* DNA4: 4 bits per base */
            int byte_pos = bit_pos / 8;
            int bit_offset = bit_pos % 8;
            
            uint8 dna4_base;
            if (bit_offset <= 4) {
                dna4_base = (src_data[byte_pos] >> (4 - bit_offset)) & 0xF;
            } else {
                /* Base spans two bytes */
                uint8 high_bits = (src_data[byte_pos] & ((1 << (8 - bit_offset)) - 1)) << (bit_offset - 4);
                uint8 low_bits = (src_data[byte_pos + 1] >> (12 - bit_offset)) & ((1 << (bit_offset - 4)) - 1);
                dna4_base = high_bits | low_bits;
            }
            
            /* Get expansions for this base */
            uint64 base_expansions[4];
            int base_expansion_count;
            kmersearch_expand_dna4_base_direct(dna4_base, base_expansions, &base_expansion_count);
            
            /* Skip if expansion would exceed reasonable limits */
            if (current_count * base_expansion_count > 10) {
                current_count = 0;
                break;
            }
            
            /* Multiply current expansions by base expansions */
            int new_count = current_count * base_expansion_count;
            uint64 temp_expansions[1024];
            
            for (int c = 0; c < current_count; c++) {
                for (int b = 0; b < base_expansion_count; b++) {
                    temp_expansions[c * base_expansion_count + b] = 
                        (current_expansions[c] << 2) | base_expansions[b];
                }
            }
            
            /* Copy back to current_expansions */
            for (int t = 0; t < new_count; t++) {
                current_expansions[t] = temp_expansions[t];
            }
            current_count = new_count;
        }
        
        /* Store results in appropriate uint type */
        if (k <= 8) {
            uint16 *output_array = (uint16 *) result;
            for (int j = 0; j < current_count; j++) {
                output_array[key_count++] = (uint16) current_expansions[j];
            }
        } else if (k <= 16) {
            uint32 *output_array = (uint32 *) result;
            for (int j = 0; j < current_count; j++) {
                output_array[key_count++] = (uint32) current_expansions[j];
            }
        } else {
            uint64 *output_array = (uint64 *) result;
            for (int j = 0; j < current_count; j++) {
                output_array[key_count++] = current_expansions[j];
            }
        }
    }
    
    *nkeys = key_count;
    return result;
}
```

### Phase 3: SIMD Optimizations

#### 3.1 AVX2 Implementation Pattern

```c
__attribute__((target("avx2")))
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_avx2(VarBit *seq, int k, void **output, int *nkeys)
{
    /* AVX2-optimized parallel k-mer extraction */
    /* Process 8 k-mers simultaneously using 256-bit vectors */
    /* Direct bit manipulation to avoid VarBit overhead */
    /* Use _mm256_load_si256() and bit shifting operations */
}
```

#### 3.2 AVX512 Implementation Pattern

```c
__attribute__((target("avx512f,avx512bw")))
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_avx512(VarBit *seq, int k, void **output, int *nkeys)
{
    /* AVX512-optimized parallel k-mer extraction */
    /* Process 16 k-mers simultaneously using 512-bit vectors */
    /* Use _mm512_load_si512() and advanced bit manipulation */
}
```

#### 3.3 ARM NEON Implementation Pattern

```c
__attribute__((target("neon")))
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_neon(VarBit *seq, int k, void **output, int *nkeys)
{
    /* NEON-optimized parallel k-mer extraction */
    /* Process 4 k-mers simultaneously using 128-bit vectors */
    /* Use vld1q_u8() and bit manipulation instructions */
}
```

#### 3.4 ARM SVE Implementation Pattern

```c
__attribute__((target("sve")))
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_sve(VarBit *seq, int k, void **output, int *nkeys)
{
    /* SVE-optimized parallel k-mer extraction */
    /* Variable vector width processing */
    /* Use svld1() and scalable vector operations */
}
```

### Phase 4: Integration with High-Frequency Analysis

#### 4.1 Modified Analysis Function

Update `kmersearch_perform_highfreq_analysis` to use new functions:

```c
/* In kmersearch_freq.c, replace VarBit extraction with uint extraction */
if (is_dna4_type) {
    void *kmer_uints = kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(sequence, k_size, &output, &nkeys);
    /* Cast to appropriate uint type based on k_size */
    if (k_size <= 8) {
        uint16 *kmer_array = (uint16 *) output;
        /* Process uint16 array directly */
    } else if (k_size <= 16) {
        uint32 *kmer_array = (uint32 *) output;
        /* Process uint32 array directly */
    } else {
        uint64 *kmer_array = (uint64 *) output;
        /* Process uint64 array directly */
    }
} else {
    void *kmer_uints = kmersearch_extract_dna2_kmer2_as_uint_direct(sequence, k_size, &output, &nkeys);
    /* Cast to appropriate uint type based on k_size */
}
```

## Performance Optimizations

### Memory Layout Optimization

1. **Contiguous Arrays**: Direct uint arrays for better cache locality
2. **Type-Specific Processing**: Separate code paths for uint16/uint32/uint64
3. **SIMD-Friendly Alignment**: 32-byte aligned memory allocation for AVX2/AVX512
4. **Zero VarBit Intermediate Objects**: Complete elimination of temporary VarBit k-mer creation

### Algorithmic Improvements

1. **Batch Processing**: Process multiple k-mers per SIMD operation
2. **Reduced Function Calls**: Eliminate per-kmer function call overhead
3. **Direct Bit Manipulation**: Avoid VarBit structure construction entirely
4. **In-place Degenerate Expansion**: Generate all expansions directly to uint arrays

## Critical Performance Requirements

### Elimination of VarBit Intermediates

**MANDATORY**: The following patterns must be completely avoided:

❌ **FORBIDDEN**: Creating VarBit k-mer objects and then converting to uint
```c
VarBit *kmer = kmersearch_create_kmer2_key_from_dna2_bits(seq, pos, k);
uint64 hash = kmersearch_get_kmer_hash(kmer, 0, k);
pfree(kmer);
```

✅ **REQUIRED**: Direct bit manipulation to uint
```c
uint64 kmer_value = 0;
for (int j = 0; j < k; j++) {
    int bit_pos = (pos + j) * 2;
    int byte_pos = bit_pos / 8;
    int bit_offset = bit_pos % 8;
    uint8 base_bits = (src_data[byte_pos] >> (6 - bit_offset)) & 0x3;
    kmer_value = (kmer_value << 2) | base_bits;
}
```

### Memory Allocation Efficiency

- Single allocation for entire result array
- No per-kmer memory allocations
- Immediate deallocation of temporary expansion arrays

### Function Call Overhead Elimination

- Inline bit manipulation operations
- Avoid function calls within inner loops
- Direct array access patterns

## Testing Strategy

### Unit Tests

1. **Correctness Verification**: Compare output with existing VarBit-based functions
2. **Edge Cases**: Test boundary conditions (k=4, k=32, short sequences)
3. **SIMD Validation**: Verify all SIMD variants produce identical results
4. **Performance Validation**: Confirm significant speedup over VarBit-based methods

### Performance Tests

1. **Benchmark Comparison**: Measure performance vs. existing functions (target: 20-40% improvement)
2. **Memory Usage**: Verify reduced memory footprint (target: 50-75% reduction)
3. **Large Dataset Testing**: Test with real genomic data
4. **Memory Allocation Profiling**: Confirm elimination of VarBit intermediate allocations

### Regression Tests

Add to existing test suite:
- `sql/13_uint_extraction.sql`
- `expected/13_uint_extraction.out`

## Implementation Timeline

### Phase 1: Infrastructure (Week 1)
- Add header declarations
- Implement direct bit manipulation utility functions
- Create test framework

### Phase 2: Scalar Implementation (Week 2)
- Implement scalar DNA2 function with direct bit manipulation
- Implement scalar DNA4 function with direct degenerate expansion
- Basic unit tests and performance validation

### Phase 3: SIMD Optimization (Weeks 3-4)
- Implement x86_64 SIMD variants
- Implement aarch64 SIMD variants
- Performance testing

### Phase 4: Integration (Week 5)
- Integrate with high-frequency analysis
- Complete testing suite
- Documentation updates

## Files to Modify

### Core Implementation
- `kmersearch.c`: Add new function implementations
- `kmersearch.h`: Add function declarations and structures

### Integration Points
- `kmersearch_freq.c`: Update high-frequency analysis to use new functions
- `kmersearch_kmer.c`: Add helper functions

### Testing
- `sql/13_uint_extraction.sql`: New test cases
- `expected/13_uint_extraction.out`: Expected test results

## Expected Performance Gains

### Memory Efficiency
- **Reduction**: 50-75% memory usage (eliminate VarBit headers and intermediate objects)
- **Cache Performance**: Improved due to contiguous uint arrays
- **Allocation Overhead**: Eliminated per-kmer allocation costs

### Processing Speed
- **Small k-mers (k≤8)**: 25-35% improvement using uint16 + direct bit manipulation
- **Medium k-mers (k≤16)**: 30-40% improvement using uint32 + direct bit manipulation
- **Large k-mers (k≤32)**: 35-45% improvement using uint64 + direct bit manipulation
- **SIMD Environments**: Additional 10-20% improvement

### High-Frequency Analysis
- **Overall speedup**: 30-50% for typical genomic datasets
- **Parallel processing**: Better scalability with reduced memory pressure
- **Memory bandwidth**: Improved utilization due to smaller data structures

## Risk Mitigation

### Compatibility
- Maintain existing VarBit-based functions for backward compatibility
- Gradual migration path for existing code

### Portability
- Comprehensive fallback to scalar implementations
- Runtime SIMD detection prevents crashes on older hardware

### Testing
- Extensive cross-validation with existing implementations
- Performance regression testing on multiple architectures

## Future Enhancements

### Additional Optimizations
- GPU acceleration using CUDA/OpenCL
- Memory-mapped file processing for large datasets
- Adaptive batch size based on available memory

### Extended Functionality
- Direct integration with GIN index construction
- Streaming k-mer extraction for very large sequences
- Compressed output formats for space efficiency