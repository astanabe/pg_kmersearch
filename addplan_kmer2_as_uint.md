# Implementation Plan: K-mer Extraction with Native uint Output

## Overview

This document outlines the implementation plan for new k-mer extraction functions that output native unsigned integer types instead of VarBit structures, providing significant performance improvements for high-frequency k-mer analysis.

## Functions to Implement

### Primary Functions

1. **`kmersearch_extract_dna2_kmer2_as_uint_direct(VarBit *seq, int k, void **output, int *nkeys)`**
   - Extracts k-mers from DNA2 sequences as native uint types
   - Output type depends on k-mer size: uint16 (k≤8), uint32 (k≤16), uint64 (k≤32)

2. **`kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(VarBit *seq, int k, void **output, int *nkeys)`**
   - Extracts k-mers from DNA4 sequences with degenerate expansion as native uint types
   - Output type depends on k-mer size: uint16 (k≤8), uint32 (k≤16), uint64 (k≤32)

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

/* Helper functions */
size_t kmersearch_get_kmer_uint_size(int k);
KmerData kmersearch_extract_single_kmer2_as_uint(VarBit *seq, int start_pos, int k);
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
```

### Phase 2: Scalar Implementations

#### 2.1 DNA2 Scalar Implementation

```c
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_scalar(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 2;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return NULL;
    }
    
    result = palloc(max_kmers * element_size);
    *output = result;
    
    if (k <= 8) {
        uint16 *output_array = (uint16 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            output_array[*nkeys] = (uint16) kmersearch_get_kmer_hash_fast(seq, i, k);
            (*nkeys)++;
        }
    } else if (k <= 16) {
        uint32 *output_array = (uint32 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            output_array[*nkeys] = (uint32) kmersearch_get_kmer_hash_fast(seq, i, k);
            (*nkeys)++;
        }
    } else {
        uint64 *output_array = (uint64 *) result;
        for (i = 0; i <= seq_bases - k; i++) {
            output_array[*nkeys] = kmersearch_get_kmer_hash_fast(seq, i, k);
            (*nkeys)++;
        }
    }
    
    return result;
}
```

#### 2.2 DNA4 Scalar Implementation

```c
static void *kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_scalar(VarBit *seq, int k, void **output, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    int seq_bases = seq_bits / 4;
    int max_kmers = (seq_bases >= k) ? (seq_bases - k + 1) : 0;
    size_t element_size = kmersearch_get_kmer_uint_size(k);
    void *result;
    int key_count = 0;
    int i;
    
    *nkeys = 0;
    if (max_kmers <= 0) {
        *output = NULL;
        return NULL;
    }
    
    /* Allocate for maximum possible expansions (10x) */
    result = palloc(max_kmers * 10 * element_size);
    *output = result;
    
    for (i = 0; i <= seq_bases - k; i++) {
        uint64 *expanded_hashes;
        int expansion_count;
        int j;
        
        /* Get expanded k-mer hashes directly */
        expanded_hashes = kmersearch_expand_dna4_kmer_to_uint64_hashes(seq, i, k, &expansion_count);
        
        if (!expanded_hashes || expansion_count == 0)
            continue;
        
        /* Store results in appropriate uint type */
        if (k <= 8) {
            uint16 *output_array = (uint16 *) result;
            for (j = 0; j < expansion_count; j++) {
                output_array[key_count++] = (uint16) expanded_hashes[j];
            }
        } else if (k <= 16) {
            uint32 *output_array = (uint32 *) result;
            for (j = 0; j < expansion_count; j++) {
                output_array[key_count++] = (uint32) expanded_hashes[j];
            }
        } else {
            uint64 *output_array = (uint64 *) result;
            for (j = 0; j < expansion_count; j++) {
                output_array[key_count++] = expanded_hashes[j];
            }
        }
        
        pfree(expanded_hashes);
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
}
```

#### 3.2 AVX512 Implementation Pattern

```c
__attribute__((target("avx512f,avx512bw")))
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_avx512(VarBit *seq, int k, void **output, int *nkeys)
{
    /* AVX512-optimized parallel k-mer extraction */
    /* Process 16 k-mers simultaneously using 512-bit vectors */
}
```

#### 3.3 ARM NEON Implementation Pattern

```c
__attribute__((target("neon")))
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_neon(VarBit *seq, int k, void **output, int *nkeys)
{
    /* NEON-optimized parallel k-mer extraction */
    /* Process 4 k-mers simultaneously using 128-bit vectors */
}
```

#### 3.4 ARM SVE Implementation Pattern

```c
__attribute__((target("sve")))
static void *kmersearch_extract_dna2_kmer2_as_uint_direct_sve(VarBit *seq, int k, void **output, int *nkeys)
{
    /* SVE-optimized parallel k-mer extraction */
    /* Variable vector width processing */
}
```

### Phase 4: Integration with High-Frequency Analysis

#### 4.1 Modified Analysis Function

Update `kmersearch_perform_highfreq_analysis` to use new functions:

```c
/* In kmersearch_freq.c, replace VarBit extraction with uint extraction */
if (is_dna4_type) {
    void *kmer_uints = kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(sequence, k_size, &output, &nkeys);
    uint64 *kmer_array = (uint64 *) output;  /* Cast based on k_size */
} else {
    void *kmer_uints = kmersearch_extract_dna2_kmer2_as_uint_direct(sequence, k_size, &output, &nkeys);
    uint64 *kmer_array = (uint64 *) output;  /* Cast based on k_size */
}
```

## Performance Optimizations

### Memory Layout Optimization

1. **Contiguous Arrays**: Direct uint arrays for better cache locality
2. **Type-Specific Processing**: Separate code paths for uint16/uint32/uint64
3. **SIMD-Friendly Alignment**: 32-byte aligned memory allocation for AVX2/AVX512

### Algorithmic Improvements

1. **Batch Processing**: Process multiple k-mers per SIMD operation
2. **Reduced Function Calls**: Eliminate per-kmer function call overhead
3. **Direct Bit Manipulation**: Avoid VarBit structure construction

## Testing Strategy

### Unit Tests

1. **Correctness Verification**: Compare output with existing VarBit-based functions
2. **Edge Cases**: Test boundary conditions (k=4, k=32, short sequences)
3. **SIMD Validation**: Verify all SIMD variants produce identical results

### Performance Tests

1. **Benchmark Comparison**: Measure performance vs. existing functions
2. **Memory Usage**: Verify reduced memory footprint
3. **Large Dataset Testing**: Test with real genomic data

### Regression Tests

Add to existing test suite:
- `sql/13_uint_extraction.sql`
- `expected/13_uint_extraction.out`

## Implementation Timeline

### Phase 1: Infrastructure (Week 1)
- Add header declarations
- Implement utility functions
- Create test framework

### Phase 2: Scalar Implementation (Week 2)
- Implement scalar DNA2 function
- Implement scalar DNA4 function
- Basic unit tests

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
- **Reduction**: 50-75% memory usage (eliminate VarBit headers)
- **Cache Performance**: Improved due to contiguous uint arrays

### Processing Speed
- **Small k-mers (k≤8)**: 15-25% improvement using uint16
- **Medium k-mers (k≤16)**: 20-30% improvement using uint32
- **Large k-mers (k≤32)**: 25-35% improvement using uint64
- **SIMD Environments**: Additional 5-15% improvement

### High-Frequency Analysis
- **Overall speedup**: 20-40% for typical genomic datasets
- **Parallel processing**: Better scalability with reduced memory pressure

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