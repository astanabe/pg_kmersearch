# SIMD Functionality Enablement Plan

## Overview
This document outlines the plan to enable all SIMD (Single Instruction, Multiple Data) optimizations in pg_kmersearch. Currently, only encoding/decoding functions utilize SIMD, while comparison and k-mer extraction functions have incomplete or disabled SIMD implementations.

## CRITICAL ISSUE DISCOVERED
During refactoring, the SIMD dispatch mechanism was severely damaged:
1. **Data size thresholds were removed** - SIMD selection previously considered both CPU capabilities AND data size
2. **Dispatch table approach is fundamentally flawed** - Different functions require different size thresholds
3. **Original implementation in master branch** shows proper per-function threshold checking

## Current Status (as of refactor branch commit b4fbcbe)

### ✅ Fully Enabled SIMD Functions (Using Dispatch Table)
- `dna2_encode` (AVX2, AVX512, NEON, SVE)
- `dna2_decode` (AVX2, AVX512, NEON, SVE)
- `dna4_encode` (AVX2, AVX512, NEON, SVE)
- `dna4_decode` (AVX2, AVX512, NEON, SVE)
- `dna_compare` (AVX2, AVX512, NEON, SVE) - Phase 1 completed on 2025-07-27

### ❌ Disabled/Incomplete SIMD Functions
1. **K-mer Extraction (as_uint)** - Functions not implemented
2. **K-mer Extraction (direct)** - SIMD implementations exist but missing data size threshold checks
3. **K-mer Matching** - SIMD implementations exist but missing data size threshold checks

### ⚠️ Critical Issues Identified
1. **Dispatch table approach** - All SIMD functions currently use dispatch table without data size checks
2. **Missing thresholds** - K-mer extraction functions in `kmersearch_kmer.c` use only `simd_capability` without checking data size
3. **Inconsistent implementation** - Some functions have proper threshold checks (in master) but were changed to dispatch table

## SIMD Implementations Available in Master Branch (commit 517bb58)

### Fully Implemented SIMD Functions in Master
After investigating master branch commit 517bb58, the following SIMD implementations are already complete and can be retrieved:

#### 1. K-mer Extraction (as_uint variants) - FULLY IMPLEMENTED
- **DNA2 extraction with proper threshold checking**:
  - `kmersearch_extract_dna2_kmer2_as_uint_direct()` - Dispatch function with data size thresholds
  - `kmersearch_extract_dna2_kmer2_as_uint_direct_scalar()` - Scalar implementation
  - `kmersearch_extract_dna2_kmer2_as_uint_direct_avx2()` - AVX2 implementation
  - `kmersearch_extract_dna2_kmer2_as_uint_direct_avx512()` - AVX512 implementation
  - `kmersearch_extract_dna2_kmer2_as_uint_direct_neon()` - NEON implementation
  - `kmersearch_extract_dna2_kmer2_as_uint_direct_sve()` - SVE implementation

- **DNA4 extraction with proper threshold checking**:
  - `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()` - Dispatch with thresholds
  - `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_scalar()` - Scalar
  - `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx2()` - AVX2
  - `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx512()` - AVX512
  - `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_neon()` - NEON
  - `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_sve()` - SVE

#### 2. K-mer Matching - FULLY IMPLEMENTED
- `kmersearch_count_matching_kmer_fast()` - Dispatch function with key combination thresholds
- `kmersearch_count_matching_kmer_fast_scalar_simple()` - Simple O(n*m) for small datasets
- `kmersearch_count_matching_kmer_fast_scalar_hashtable()` - Hash table for larger datasets
- `kmersearch_count_matching_kmer_fast_avx2()` - AVX2 optimized hash table
- `kmersearch_count_matching_kmer_fast_avx512()` - AVX512 optimized
- `kmersearch_count_matching_kmer_fast_neon()` - NEON optimized
- `kmersearch_count_matching_kmer_fast_sve()` - SVE optimized

#### 3. Proper Threshold-Based Dispatch Examples from Master
All functions in master properly check both CPU capability AND data size:
```c
// Example from master: kmersearch_extract_dna2_kmer2_as_uint_direct()
if (simd_capability >= SIMD_AVX512BW && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
    kmersearch_extract_dna2_kmer2_as_uint_direct_avx512(seq, k, output, nkeys);
    return;
}
```

## Implementation Plan

### Phase 0: Remove Dispatch Table System [COMPLETED - 2025-07-27]
**Priority: CRITICAL** - Must be completed before any other work
**Status**: Successfully completed without errors or warnings

#### Actions to Take from Master Branch
1. **Copy threshold-based dispatch patterns** from master for all functions
2. **Retrieve complete as_uint implementations** - These are fully implemented with proper SIMD
3. **Retrieve k-mer matching SIMD implementations** - Complete set with proper thresholds
4. **Study the dispatch patterns** - Master shows correct per-function threshold checking

#### 0.1 Remove dispatch table from all functions [COMPLETED]
- **Files**: `kmersearch.c`, `kmersearch_datatype.c`, `kmersearch_kmer.c`
- **Tasks**:
  - [x] Remove `simd_dispatch_table_t` structure from `kmersearch.h`
  - [x] Remove `simd_dispatch` global variable
  - [x] Remove `init_simd_dispatch_table()` function
  - [x] Update all SIMD-enabled functions to use direct selection

#### 0.2 Restore per-function SIMD selection with data size thresholds
- **Pattern to follow** (from master branch):
```c
Datum *
kmersearch_extract_dna2_kmer2_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
    /* Use SIMD based on runtime capability and data size thresholds */
#ifdef __x86_64__
    if (simd_capability >= SIMD_AVX512BW && seq_bits >= SIMD_EXTRACT_AVX512_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_avx512(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_AVX2 && seq_bits >= SIMD_EXTRACT_AVX2_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_avx2(seq, k, nkeys);
    }
#elif defined(__aarch64__)
    if (simd_capability >= SIMD_SVE && seq_bits >= SIMD_EXTRACT_SVE_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_sve(seq, k, nkeys);
    }
    if (simd_capability >= SIMD_NEON && seq_bits >= SIMD_EXTRACT_NEON_THRESHOLD) {
        return kmersearch_extract_dna2_kmer2_direct_neon(seq, k, nkeys);
    }
#endif
    return kmersearch_extract_dna2_kmer2_direct_scalar(seq, k, nkeys);
}
```

- **Pattern for encode/decode functions**:
```c
void
dna2_encode_simd(const char *input, uint8_t *output, int len)
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
dna2_decode_simd(const uint8_t *input, char *output, int bit_len)
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
```

#### 0.3 Create new dispatch functions and update call sites [COMPLETED]
1. **Create new encoding/decoding functions** (`kmersearch_datatype.c`):
   - [x] Create `dna2_encode()` - Main function with threshold-based SIMD dispatch
   - [x] Create `dna2_decode()` - Main function with threshold-based SIMD dispatch
   - [x] Create `dna4_encode()` - Main function with threshold-based SIMD dispatch
   - [x] Create `dna4_decode()` - Main function with threshold-based SIMD dispatch
   
   **New thresholds to add to `kmersearch.h`**:
   ```c
   /* SIMD encoding thresholds (input character length)
    * Initially set to same values as SIMD_EXTRACT thresholds
    * TODO: These values need performance testing for optimal settings
    */
   #define SIMD_ENCODE_AVX2_THRESHOLD    512     /* 512 chars: Use AVX2 for encoding */
   #define SIMD_ENCODE_AVX512_THRESHOLD  1024    /* 1024 chars: Use AVX512 for encoding */
   #define SIMD_ENCODE_NEON_THRESHOLD    256     /* 256 chars: Use NEON for encoding */
   #define SIMD_ENCODE_SVE_THRESHOLD     512     /* 512 chars: Use SVE for encoding */
   
   /* SIMD decoding thresholds (bit length)
    * Initially set to same values as SIMD_EXTRACT thresholds
    * TODO: These values need performance testing for optimal settings
    */
   #define SIMD_DECODE_AVX2_THRESHOLD    512     /* 512 bits: Use AVX2 for decoding */
   #define SIMD_DECODE_AVX512_THRESHOLD  1024    /* 1024 bits: Use AVX512 for decoding */
   #define SIMD_DECODE_NEON_THRESHOLD    256     /* 256 bits: Use NEON for decoding */
   #define SIMD_DECODE_SVE_THRESHOLD     512     /* 512 bits: Use SVE for decoding */
   ```

2. **Rename and update comparison function** (`kmersearch_datatype.c`):
   - [x] Create `dna_compare()` - Main comparison function with threshold-based SIMD dispatch
   - [x] Keep `dna_compare_simd()` as deprecated wrapper calling `dna_compare()`
   - [x] Update function to use threshold-based SIMD dispatch instead of dispatch table

3. **Update all call sites**:
   - [x] Replace all `simd_dispatch.dna2_encode()` calls with `dna2_encode()`
   - [x] Replace all `simd_dispatch.dna2_decode()` calls with `dna2_decode()`
   - [x] Replace all `simd_dispatch.dna4_encode()` calls with `dna4_encode()`
   - [x] Replace all `simd_dispatch.dna4_decode()` calls with `dna4_decode()`
   - [x] Keep `dna_compare_simd()` calls working through wrapper function

3. **K-mer extraction functions** (`kmersearch.c`, `kmersearch_kmer.c`):
   - [x] `kmersearch_extract_dna2_kmer2_direct()` - Updated with threshold-based dispatch
   - [x] `kmersearch_extract_dna4_kmer2_with_expansion_direct()` - Updated with threshold-based dispatch
   - [ ] `kmersearch_extract_dna2_kmer2_as_uint_direct()` - To be imported from master
   - [ ] `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()` - To be imported from master

4. **K-mer matching functions** (`kmersearch_kmer.c`):
   - [x] `kmersearch_count_matching_kmer_fast()` - Updated with threshold-based dispatch (TODOs for SIMD implementations)

#### 0.4 Verify thresholds are appropriate
- **Current threshold values** in `kmersearch.h`:
  - **Comparison thresholds** (bit length) - Already exist and will be used:
    - `SIMD_COMPARE_AVX2_THRESHOLD = 128` (128 bits)
    - `SIMD_COMPARE_AVX512_THRESHOLD = 256` (256 bits)
    - `SIMD_COMPARE_NEON_THRESHOLD = 64` (64 bits)
    - `SIMD_COMPARE_SVE_THRESHOLD = 128` (128 bits)
  - **K-mer extraction thresholds** (sequence bit length):
    - `SIMD_EXTRACT_AVX2_THRESHOLD = 512` (512 bits)
    - `SIMD_EXTRACT_AVX512_THRESHOLD = 1024` (1024 bits)
    - `SIMD_EXTRACT_NEON_THRESHOLD = 256` (256 bits)
    - `SIMD_EXTRACT_SVE_THRESHOLD = 512` (512 bits)
  - **K-mer matching thresholds** (key combination count):
    - `SIMD_KEYCOMB_AVX2_THRESHOLD = 128` (128 combinations)
    - `SIMD_KEYCOMB_AVX512_THRESHOLD = 256` (256 combinations)
    - `SIMD_KEYCOMB_NEON_THRESHOLD = 64` (64 combinations)
    - `SIMD_KEYCOMB_SVE_THRESHOLD = 128` (128 combinations)
    
- **New thresholds added** [COMPLETED]:
  - **Encoding thresholds** (character count) - Initially same as EXTRACT thresholds:
    - `SIMD_ENCODE_AVX2_THRESHOLD = 512` (512 characters)
    - `SIMD_ENCODE_AVX512_THRESHOLD = 1024` (1024 characters)
    - `SIMD_ENCODE_NEON_THRESHOLD = 256` (256 characters)
    - `SIMD_ENCODE_SVE_THRESHOLD = 512` (512 characters)
  - **Decoding thresholds** (bit length) - Initially same as EXTRACT thresholds:
    - `SIMD_DECODE_AVX2_THRESHOLD = 512` (512 bits)
    - `SIMD_DECODE_AVX512_THRESHOLD = 1024` (1024 bits)
    - `SIMD_DECODE_NEON_THRESHOLD = 256` (256 bits)
    - `SIMD_DECODE_SVE_THRESHOLD = 512` (512 bits)
    
- **TODO**: All threshold values need comprehensive performance testing to determine optimal settings. The initial values for ENCODE/DECODE thresholds are copied from EXTRACT thresholds and may not be optimal.

### Phase 1: Enable DNA Comparison SIMD Dispatch [COMPLETED - 2025-07-27]
**Priority: High** - Implementation already exists, just needs integration
**Status**: Successfully built without errors or warnings
**NOTE**: This phase used dispatch table but needs to be refactored in Phase 0
**IMPORTANT**: The function `dna_compare_simd()` created in this phase will be renamed to `dna_compare()` during Phase 0 refactoring

#### 1.1 Update simd_dispatch_table_t initialization
- **File**: `kmersearch.c`
- **Function**: `init_simd_dispatch_table()`
- **Tasks**:
  - [x] Add `simd_dispatch.dna_compare = dna_compare_scalar;` to default initialization
  - [x] Add AVX2: `simd_dispatch.dna_compare = dna_compare_avx2;`
  - [x] Add AVX512: `simd_dispatch.dna_compare = dna_compare_avx512;`
  - [x] Add NEON: `simd_dispatch.dna_compare = dna_compare_neon;`
  - [x] Add SVE: `simd_dispatch.dna_compare = dna_compare_sve;`

#### 1.2 Update dna_compare_simd to use dispatch table
- **File**: `kmersearch_datatype.c`
- **Function**: `dna_compare_simd()`
- **Tasks**:
  - [x] Replace direct SIMD selection with `simd_dispatch.dna_compare(a, b, bit_len)`
  - [x] Remove redundant switch statement

#### 1.3 Export comparison functions in header
- **File**: `kmersearch.h`
- **Tasks**:
  - [x] Add function declarations for all dna_compare variants

### Phase 2: Import K-mer Extraction SIMD (as_uint variants) from Master
**Priority: High** - Complete implementations available in master branch
**Status**: Can be copied directly from master commit 517bb58

#### 2.1 Import DNA2 k-mer extraction SIMD functions
- **Source**: master branch `kmersearch.c`
- **Functions to import**:
  - [x] `kmersearch_extract_dna2_kmer2_as_uint_direct()` - Already has proper threshold dispatch
  - [x] `kmersearch_extract_dna2_kmer2_as_uint_direct_scalar()` - Complete implementation
  - [x] `kmersearch_extract_dna2_kmer2_as_uint_direct_avx2()` - Complete with SIMD optimization
  - [x] `kmersearch_extract_dna2_kmer2_as_uint_direct_avx512()` - Complete with SIMD optimization
  - [x] `kmersearch_extract_dna2_kmer2_as_uint_direct_neon()` - Complete with SIMD optimization
  - [x] `kmersearch_extract_dna2_kmer2_as_uint_direct_sve()` - Complete with SIMD optimization

#### 2.2 Import DNA4 k-mer extraction SIMD functions
- **Source**: master branch `kmersearch.c`
- **Functions to import**:
  - [x] `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()` - Has proper threshold dispatch
  - [x] `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_scalar()` - Complete
  - [x] `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx2()` - Complete
  - [x] `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx512()` - Complete
  - [x] `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_neon()` - Complete
  - [x] `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_sve()` - Complete

#### 2.3 Import K-mer matching SIMD functions
- **Source**: master branch `kmersearch.c`
- **Functions to import**:
  - [x] `kmersearch_count_matching_kmer_fast()` - Has proper key combination threshold dispatch
  - [x] `kmersearch_count_matching_kmer_fast_scalar_simple()` - For small datasets
  - [x] `kmersearch_count_matching_kmer_fast_scalar_hashtable()` - For larger datasets
  - [x] `kmersearch_count_matching_kmer_fast_avx2()` - AVX2 optimized
  - [x] `kmersearch_count_matching_kmer_fast_avx512()` - AVX512 optimized
  - [x] `kmersearch_count_matching_kmer_fast_neon()` - NEON optimized
  - [x] `kmersearch_count_matching_kmer_fast_sve()` - SVE optimized

### Phase 3: Complete K-mer Extraction SIMD (direct variants)
**Priority: Low** - Implementations exist but need optimization
**Status**: Basic implementations completed, but not truly SIMD-optimized

#### 3.1 AVX2 implementations
- **File**: `kmersearch_kmer.c`
- **Current Status**: 
  - `kmersearch_extract_dna2_kmer2_direct_avx2()` - Implemented but mostly uses scalar operations
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2()` - Stub calling scalar version
- **Tasks**:
  - [ ] Optimize `kmersearch_extract_dna2_kmer2_direct_avx2()` with true SIMD operations
  - [ ] Implement `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2()` with SIMD

#### 3.2 AVX512 implementations
- **File**: `kmersearch_kmer.c`
- **Current Status**: 
  - `kmersearch_extract_dna2_kmer2_direct_avx512()` - Implemented but mostly uses scalar operations
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512()` - Stub calling scalar version
- **Tasks**:
  - [ ] Optimize `kmersearch_extract_dna2_kmer2_direct_avx512()` with true SIMD operations
  - [ ] Implement `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512()` with SIMD

#### 3.3 NEON implementations
- **File**: `kmersearch_kmer.c`
- **Current Status**: 
  - `kmersearch_extract_dna2_kmer2_direct_neon()` - Implemented but mostly uses scalar operations
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_neon()` - Stub calling scalar version
- **Tasks**:
  - [ ] Optimize `kmersearch_extract_dna2_kmer2_direct_neon()` with true SIMD operations
  - [ ] Implement `kmersearch_extract_dna4_kmer2_with_expansion_direct_neon()` with SIMD

#### 3.4 SVE implementations
- **File**: `kmersearch_kmer.c`
- **Current Status**: 
  - `kmersearch_extract_dna2_kmer2_direct_sve()` - Implemented but mostly uses scalar operations
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_sve()` - Stub calling scalar version
- **Tasks**:
  - [ ] Optimize `kmersearch_extract_dna2_kmer2_direct_sve()` with true SIMD operations
  - [ ] Implement `kmersearch_extract_dna4_kmer2_with_expansion_direct_sve()` with SIMD

## Testing Strategy

### Unit Tests
1. **Comparison Tests**
   - [ ] Create test cases comparing scalar vs SIMD results
   - [ ] Test edge cases (different bit lengths, alignment)

2. **K-mer Extraction Tests**
   - [ ] Verify identical output between scalar and SIMD versions
   - [ ] Test performance improvements
   - [ ] Test degenerate code handling in DNA4

### Performance Benchmarks
- [ ] Create benchmark suite for each SIMD function
- [ ] Measure speedup vs scalar implementation
- [ ] Document optimal thresholds for SIMD usage

### Regression Tests
- [ ] Run full `make installcheck` after each phase
- [ ] Verify no functional changes

## Implementation Guidelines

### Code Standards
1. **Alignment**: Ensure proper memory alignment for SIMD operations
2. **Fallback**: Always maintain scalar fallback for unsupported CPUs
3. **Documentation**: Document SIMD algorithms and assumptions

### Safety Considerations
1. **CPU Detection**: Verify CPU capabilities before using SIMD instructions
2. **Bounds Checking**: Ensure SIMD operations don't read beyond allocated memory
3. **Error Handling**: Graceful degradation to scalar on SIMD failures

## Expected Benefits

### Performance Improvements
- **DNA Comparison**: 2-8x speedup for long sequences
- **K-mer Extraction**: 4-16x speedup for batch processing
- **Overall**: 30-50% improvement in GIN index creation time

### Resource Efficiency
- Reduced CPU usage for large datasets
- Better cache utilization
- Lower power consumption on supported hardware

## Risk Mitigation

### Potential Issues
1. **Compatibility**: Some PostgreSQL builds may not support certain SIMD instructions
2. **Debugging**: SIMD code is harder to debug
3. **Maintenance**: Requires platform-specific expertise

### Mitigation Strategies
1. Extensive testing on multiple platforms
2. Clear documentation and code comments
3. Maintain scalar implementations as reference

## Timeline Estimate (Revised)

- **Phase 0**: 1-2 days (CRITICAL - dispatch table removal and threshold restoration)
- **Phase 1**: COMPLETED
- **Phase 2**: 1-2 days (HIGH - import from master branch, no new implementation needed)
- **Phase 3**: 2-3 weeks (LOW - optimization of existing implementations)

Total estimated time: 1-2 weeks (reduced from 3-4 weeks due to master branch availability)

## Success Metrics

1. All SIMD functions pass regression tests
2. Performance benchmarks show measurable improvements
3. No increase in bug reports or stability issues
4. Code coverage remains at 100% for SIMD paths

## Summary of Current Implementation Status

### What's Working
1. **SIMD dispatch infrastructure** - CPU capability detection and dispatch table are functional
2. **Encoding/Decoding** - Full SIMD implementations for DNA2/DNA4 encode/decode
3. **Comparison** - DNA comparison with SIMD implementations (Phase 1 completed)
4. **Basic k-mer extraction** - SIMD function stubs exist and dispatch correctly

### What's Missing
1. **Data size thresholds** - All functions dispatch based only on CPU capability, ignoring data size
2. **True SIMD optimization** - K-mer extraction "SIMD" functions mostly use scalar operations
3. **K-mer matching SIMD** - `kmersearch_count_matching_kmer_fast()` lacks SIMD variants

### Critical Path Forward
1. **Phase 0 is mandatory** - Must restore per-function threshold checking before any optimization
2. **Phase 2 can proceed in parallel** - as_uint variants don't exist yet, so can be implemented correctly from start
3. **Phase 3 is lowest priority** - Current implementations work correctly, just not optimized

## Functions That Must Be Implemented (Not Available in Master)

After thorough investigation, both master and refactor branches use dispatch table for encode/decode functions. Since we are completely removing the dispatch table system, we need to create new dispatch functions:

### 1. New Dispatch Functions to Create
These functions will replace both dispatch table usage and any `_simd` wrapper functions:

#### Encoding Functions
- **`dna2_encode()`** - Main encoding function with threshold-based SIMD dispatch
  - Replaces: `simd_dispatch.dna2_encode()` calls
  - Checks: CPU capability AND input length against `SIMD_ENCODE_*_THRESHOLD`
  
- **`dna4_encode()`** - Main encoding function with threshold-based SIMD dispatch
  - Replaces: `simd_dispatch.dna4_encode()` calls
  - Checks: CPU capability AND input length against `SIMD_ENCODE_*_THRESHOLD`

#### Decoding Functions  
- **`dna2_decode()`** - Main decoding function with threshold-based SIMD dispatch
  - Replaces: `simd_dispatch.dna2_decode()` calls
  - Checks: CPU capability AND bit length against `SIMD_DECODE_*_THRESHOLD`
  
- **`dna4_decode()`** - Main decoding function with threshold-based SIMD dispatch
  - Replaces: `simd_dispatch.dna4_decode()` calls
  - Checks: CPU capability AND bit length against `SIMD_DECODE_*_THRESHOLD`

#### Comparison Function
- **`dna_compare()`** - Main comparison function with threshold-based SIMD dispatch
  - Replaces: `dna_compare_simd()` calls in refactor branch (rename required)
  - Replaces: `simd_dispatch.dna_compare()` calls if any
  - Checks: CPU capability AND bit length against `SIMD_COMPARE_*_THRESHOLD`
  - **Note**: This is a rename of `dna_compare_simd()` to follow consistent naming convention

### 2. Implementation Pattern
Each function will follow this pattern:
```c
void dna2_encode(const char *input, uint8_t *output, int len)
{
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
```

### 3. Already Available from Master
All SIMD implementations can be imported from master:
- All individual SIMD encode/decode/compare implementations (AVX2, AVX512, NEON, SVE)
- All as_uint k-mer extraction functions (complete with proper dispatch)
- All k-mer matching functions (complete with proper dispatch)

## Recommended Implementation Order

1. **First: Phase 0.1-0.2** - Remove dispatch table infrastructure
   - Remove `simd_dispatch_table_t` and related code
   - Remove `init_simd_dispatch_table()` function

2. **Second: Phase 0.3** - Create new dispatch functions
   - Create 5 new main functions: `dna2_encode()`, `dna2_decode()`, `dna4_encode()`, `dna4_decode()`, `dna_compare()`
   - Add new threshold definitions to `kmersearch.h`
   - Update all call sites to use new functions instead of dispatch table

3. **Third: Import from Master** - Copy all complete SIMD implementations:
   - All as_uint k-mer extraction functions (Phase 2.1-2.2)
   - All k-mer matching functions (Phase 2.3)
   - Ensure all individual SIMD implementations are present

4. **Fourth: Phase 0.4** - Update remaining functions to use threshold-based dispatch
   - Update k-mer extraction functions in `kmersearch_kmer.c`
   - Verify all functions follow the same pattern

5. **Fifth: Phase 3** - Optimize the existing direct k-mer extraction SIMD implementations

### Key Changes from Original Plan
- **No `_simd` wrapper functions needed** - We create main dispatch functions directly
- **Complete removal of dispatch table** - Both master and refactor branches use it, but we're eliminating it entirely
- **Simpler architecture** - Each function directly checks thresholds and dispatches to appropriate SIMD implementation

This approach eliminates the dispatch table completely while maintaining the same performance optimization through threshold-based SIMD selection.