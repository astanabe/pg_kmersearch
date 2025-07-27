# SIMD Optimization Plan for Direct K-mer Extraction

## Overview
This document outlines the remaining SIMD optimization work for pg_kmersearch, focusing on Phase 3: optimizing direct k-mer extraction functions with advanced SIMD instructions.

## SIMD Capability Hierarchy

### Updated simd_capability_t enum
```c
typedef enum {
    SIMD_NONE = 0,              /* No SIMD support */
    
    /* x86/x64 SIMD capabilities */
    SIMD_AVX2 = 1,              /* AVX2 support */
    SIMD_BMI2 = 2,              /* AVX2 + BMI2 support */
    SIMD_AVX512F = 3,           /* AVX512F support */
    SIMD_AVX512BW = 4,          /* AVX512BW support */
    SIMD_AVX512VBMI = 5,        /* AVX512 + VBMI support */
    SIMD_AVX512VBMI2 = 6,       /* AVX512 + VBMI + VBMI2 support */
    
    /* ARM SIMD capabilities */
    SIMD_NEON = 21,             /* ARM NEON support */
    SIMD_SVE = 22,              /* ARM SVE support */
    SIMD_SVE2 = 23,             /* ARM SVE2 support */
} simd_capability_t;
```

## Phase 3: Direct K-mer Extraction SIMD Optimizations

### Overview
Optimize `kmersearch_extract_dna2_kmer2_direct_*()` and `kmersearch_extract_dna4_kmer2_with_expansion_direct_*()` functions with advanced SIMD instructions.

### 3.1 x86/x64 Implementations

#### AVX2 Functions (Requires SIMD_BMI2)
- **Functions** (names remain unchanged):
  - `kmersearch_extract_dna2_kmer2_direct_avx2()`
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2()`
- **Requirements**: SIMD_BMI2 (AVX2 + BMI2)
- **Key Instructions**:
  - `_pext_u64()` - Parallel bit extraction
  - `_pdep_u64()` - Parallel bit deposit
  - AVX2 vector operations for bulk processing
- **Implementation Strategy**:
  ```c
  __attribute__((target("avx2,bmi2")))
  static Datum *
  kmersearch_extract_dna2_kmer2_direct_avx2(VarBit *seq, int k, int *nkeys)
  {
      /* Use BMI2 PEXT for efficient k-mer extraction */
      uint64_t mask = ((1ULL << (k * 2)) - 1);
      for (int i = 0; i < max_kmers; i++) {
          uint64_t src_bits = load_64bits_at_offset(seq, i * 2);
          uint64_t kmer = _pext_u64(src_bits, mask << (i * 2 % 64));
          /* Process extracted k-mer */
      }
  }
  ```

#### AVX512 Functions (Requires SIMD_AVX512VBMI2)
- **Functions** (names remain unchanged):
  - `kmersearch_extract_dna2_kmer2_direct_avx512()`
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512()`
- **Requirements**: SIMD_AVX512VBMI2
- **Key Instructions**:
  - `_mm512_multishift_epi64_epi8()` - Variable bit shift (VBMI)
  - `_mm512_shrdi_epi32/64()` - Double shift (VBMI2)
  - `_mm512_permutexvar_epi8()` - Byte permutation (VBMI)
- **Implementation Strategy**:
  ```c
  __attribute__((target("avx512f,avx512bw,avx512vbmi,avx512vbmi2")))
  static Datum *
  kmersearch_extract_dna2_kmer2_direct_avx512(VarBit *seq, int k, int *nkeys)
  {
      /* Use VBMI2 for parallel k-mer extraction */
      __m512i shift_control = _mm512_setr_epi64(0, k*2, k*4, k*6, ...);
      __m512i src_data = _mm512_loadu_si512(seq_data);
      __m512i kmers = _mm512_multishift_epi64_epi8(shift_control, src_data);
      /* Process 8 k-mers in parallel */
  }
  ```

### 3.2 ARM Implementations

#### NEON Functions (Optimized with bit manipulation extensions)
- **Functions**:
  - `kmersearch_extract_dna2_kmer2_direct_neon()`
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_neon()`
- **Requirements**: SIMD_NEON
- **Optimization Strategy**:
  - Use `VTBL/VTBX` for table-based bit extraction
  - Apply `VEXT` for sliding window operations
  - Leverage NEON bit manipulation capabilities
- **Implementation Approach**:
  ```c
  __attribute__((target("+simd")))
  static Datum *
  kmersearch_extract_dna2_kmer2_direct_neon(VarBit *seq, int k, int *nkeys)
  {
      /* Prepare lookup tables for bit extraction patterns */
      uint8x16_t extract_tables[4];
      prepare_extraction_tables(extract_tables, k);
      
      /* Use VTBL for flexible bit extraction */
      for (int i = 0; i < max_kmers; i += 8) {
          uint8x16_t data = vld1q_u8(&seq_data[i/4]);
          uint8x16_t extracted = vtbl4q_u8(extract_tables, data);
          /* Process extracted k-mers */
      }
  }
  ```

#### SVE Functions (Optimized with NEON bit manipulation assistance)
- **Functions**:
  - `kmersearch_extract_dna2_kmer2_direct_sve()`
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_sve()`
- **Requirements**: SIMD_SVE
- **Optimization Strategy**:
  - Primary: Use SVE for scalable vector operations
  - Secondary: Utilize NEON for fixed-size bit manipulation
  - Hybrid approach for SVE1 systems without SVE2
- **Implementation Approach**:
  ```c
  __attribute__((target("+sve")))
  static Datum *
  kmersearch_extract_dna2_kmer2_direct_sve(VarBit *seq, int k, int *nkeys)
  {
      /* SVE1 with NEON assistance for bit manipulation */
      int vl = svcntb();
      for (int i = 0; i < seq_len; i += vl) {
          svuint8_t sv_data = svld1_u8(svptrue_b8(), &seq_data[i]);
          /* Use NEON for complex bit extraction within chunks */
          /* Use SVE for bulk operations and predication */
      }
  }
  ```

#### SVE2 Functions (New)
- **Functions** (NEW):
  - `kmersearch_extract_dna2_kmer2_direct_sve2()`
  - `kmersearch_extract_dna4_kmer2_with_expansion_direct_sve2()`
- **Requirements**: SIMD_SVE2
- **Key Instructions**:
  - `svbext` - Bit extraction (similar to x86 PEXT)
  - `svbdep` - Bit deposit (similar to x86 PDEP)
  - SVE2 complex bit manipulation instructions
- **Implementation Strategy**:
  ```c
  __attribute__((target("+sve2")))
  static Datum *
  kmersearch_extract_dna2_kmer2_direct_sve2(VarBit *seq, int k, int *nkeys)
  {
      /* Use SVE2 bit manipulation for efficient extraction */
      svbool_t pg = svptrue_b64();
      svuint64_t positions = svindex_u64(0, k*2);
      svuint64_t mask = svdup_n_u64((1ULL << (k*2)) - 1);
      
      for (int i = 0; i < seq_len; i += svcntd()) {
          svuint64_t data = svld1_u64(pg, &seq_data[i]);
          svuint64_t kmers = svbext_u64_z(pg, data, positions, mask);
          /* Process extracted k-mers */
      }
  }
  ```

### 3.3 Infrastructure Updates

#### Dispatch Logic Updates

##### x86/x64 Dispatch
```c
Datum *
kmersearch_extract_dna2_kmer2_direct(VarBit *seq, int k, int *nkeys)
{
    int seq_bits = VARBITLEN(seq);
    
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
```

#### CPU Capability Detection
```c
static simd_capability_t
detect_cpu_capabilities(void)
{
#ifdef __x86_64__
    unsigned int eax, ebx, ecx, edx;
    
    /* Get CPUID information */
    __cpuid(0, eax, ebx, ecx, edx);
    unsigned int max_leaf = eax;
    
    if (max_leaf >= 7) {
        __cpuid_count(7, 0, eax, ebx, ecx, edx);
        bool has_avx2 = (ebx & (1 << 5)) != 0;
        bool has_bmi2 = (ebx & (1 << 8)) != 0;
        bool has_avx512f = (ebx & (1 << 16)) != 0;
        bool has_avx512bw = (ebx & (1 << 30)) != 0;
        bool has_avx512vbmi = (ecx & (1 << 1)) != 0;
        bool has_avx512vbmi2 = (ecx & (1 << 6)) != 0;
        
        /* Check OS support for AVX/AVX512 */
        __cpuid(1, eax, ebx, ecx, edx);
        bool has_xsave = (ecx & (1 << 27)) != 0;
        bool has_osxsave = (ecx & (1 << 28)) != 0;
        
        if (has_xsave && has_osxsave) {
            uint64_t xcr0 = _xgetbv(0);
            bool ymm_enabled = (xcr0 & 0x6) == 0x6;
            bool zmm_enabled = (xcr0 & 0xe6) == 0xe6;
            
            /* Return highest supported capability */
            if (zmm_enabled && has_avx512vbmi2) return SIMD_AVX512VBMI2;
            if (zmm_enabled && has_avx512vbmi) return SIMD_AVX512VBMI;
            if (zmm_enabled && has_avx512bw) return SIMD_AVX512BW;
            if (zmm_enabled && has_avx512f) return SIMD_AVX512F;
            if (ymm_enabled && has_avx2 && has_bmi2) return SIMD_BMI2;
            if (ymm_enabled && has_avx2) return SIMD_AVX2;
        }
    }
#elif defined(__aarch64__)
    #ifdef __ARM_FEATURE_SVE2
        return SIMD_SVE2;
    #elif defined(__ARM_FEATURE_SVE)
        return SIMD_SVE;
    #elif defined(__ARM_NEON)
        return SIMD_NEON;
    #endif
#endif
    return SIMD_NONE;
}
```

## Implementation Tasks

### Phase 3.1: Update simd_capability_t enum
- [x] Update enum definition in kmersearch.h with new values
- [x] Ensure proper spacing between x86 and ARM values

### Phase 3.2: x86/x64 Optimizations
- [x] Update AVX2 functions to require SIMD_BMI2
  - [x] Modify dispatch logic to check for SIMD_BMI2 instead of SIMD_AVX2
  - [x] Add BMI2 instructions (PEXT/PDEP) to `kmersearch_extract_dna2_kmer2_direct_avx2()`
  - [x] Add BMI2 instructions to `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx2()`
- [x] Update AVX512 functions to require SIMD_AVX512VBMI2
  - [x] Modify dispatch logic to check for SIMD_AVX512VBMI2
  - [x] Add VBMI/VBMI2 instructions to `kmersearch_extract_dna2_kmer2_direct_avx512()`
  - [x] Add VBMI/VBMI2 instructions to `kmersearch_extract_dna4_kmer2_with_expansion_direct_avx512()`

### Phase 3.3: ARM Optimizations
- [x] Optimize NEON implementations with bit manipulation techniques
  - [x] Implement table-based extraction for `kmersearch_extract_dna2_kmer2_direct_neon()`
  - [x] Implement table-based extraction for `kmersearch_extract_dna4_kmer2_with_expansion_direct_neon()`
- [x] Update SVE functions to use NEON for bit manipulation assistance
  - [x] Add NEON hybrid approach to `kmersearch_extract_dna2_kmer2_direct_sve()`
  - [x] Add NEON hybrid approach to `kmersearch_extract_dna4_kmer2_with_expansion_direct_sve()`
- [x] Create new SVE2 implementations
  - [x] Implement `kmersearch_extract_dna2_kmer2_direct_sve2()`
  - [x] Implement `kmersearch_extract_dna4_kmer2_with_expansion_direct_sve2()`
  - [x] Add function declarations in kmersearch_kmer.c

### Phase 3.4: Infrastructure Updates
- [x] Update dispatch functions with new requirements
- [x] Update CPU capability detection function
- [x] Add proper fallback chains
- [x] Update build system for SVE2 support

## Expected Performance Improvements

### x86/x64
- **AVX2+BMI2**: 2-3x speedup for bit extraction operations
- **AVX512+VBMI2**: 4-8x speedup with parallel k-mer extraction

### ARM
- **NEON (optimized)**: 2x speedup with table lookups
- **SVE+NEON**: 2-3x speedup with hybrid approach
- **SVE2**: 3-5x speedup with native bit manipulation

## Testing Strategy

1. **Functional Testing**
   - Verify identical results between scalar and SIMD implementations
   - Test edge cases (various k values, sequence lengths)
   - Test with and without required CPU features

2. **Performance Testing**
   - Benchmark each implementation variant
   - Measure speedup vs scalar baseline
   - Validate threshold values for SIMD selection

3. **Compatibility Testing**
   - Test on CPUs with different capability levels
   - Verify proper fallback behavior
   - Test build on different compiler versions

## Success Metrics

1. All SIMD implementations pass regression tests
2. Measured performance improvements match expectations
3. Proper fallback for unsupported CPU features
4. No performance regression for small data sizes