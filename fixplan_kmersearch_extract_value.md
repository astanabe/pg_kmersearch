# Fix Plan: Optimize GIN Index Construction with Cache-Based High-Frequency K-mer Filtering

## Overview

This plan outlines the modification of GIN index construction functions to use uint-based k-mer extraction with in-memory cache lookup for high-frequency k-mer filtering, replacing the current VarBit-based ngram_key2 extraction approach.

## Current Implementation Analysis

### Current Flow
```
kmersearch_extract_value_dna2/4() [GIN functions]
    ↓
kmersearch_extract_dna2/4_ngram_key2_direct() [VarBit extraction]
    ↓
kmersearch_create_ngram_key2_with_occurrence_from_dna2/4() [VarBit creation]
    ↓
High-frequency filtering with hash-based cache lookup
```

### Performance Issues
- Multiple VarBit allocations and conversions
- Hash calculation overhead for each ngram_key2
- Inefficient memory usage during extraction

## Proposed Optimized Flow

### New Flow
```
kmersearch_extract_value_dna2/4() [GIN functions]
    ↓
kmersearch_extract_dna2/4_kmer2_as_uint_direct() [SIMD optimized uint extraction]
    ↓
kmersearch_lookup_kmer2_as_uint_in_*_cache() [Direct uint cache lookup]
    ↓
kmersearch_create_ngram_key2_from_kmer2_as_uint() [Filtered VarBit creation]
```

### Performance Benefits
- Direct uint-based cache lookup (no hash calculation needed)
- Reduced memory allocations
- SIMD optimization utilization
- Early filtering before VarBit creation

## Implementation Plan

### 1. Core Conversion Function Implementation

#### Function: `kmersearch_kmer2_as_uint_to_kmer2()`
- **Purpose**: Convert uint-based k-mer to VarBit k-mer format
- **Parameters**: 
  - `uint16/uint32/uint64 kmer2_as_uint`: Input k-mer as uint
  - `int kmer_size`: K-mer length (4-32)
- **Return**: `VarBit*` containing 2k bits k-mer data
- **Implementation**:
  - Allocate VarBit with `(kmer_size * 2)` bits
  - Extract 2-bit nucleotides using bit masks and shifts
  - Pack into VarBit format with proper bit alignment
  - Handle different uint sizes (16/32/64) based on k-mer length

#### Function: `kmersearch_create_ngram_key2_from_kmer2_as_uint()`
- **Purpose**: Create complete ngram_key2 from uint k-mer with occurrence count
- **Parameters**:
  - `uint16/uint32/uint64 kmer2_as_uint`: Input k-mer as uint
  - `int kmer_size`: K-mer length
  - `int occurrence`: Occurrence count (0 to max_occurrence)
- **Return**: `VarBit*` containing complete ngram_key2 structure
- **Implementation**:
  - Call `kmersearch_kmer2_as_uint_to_kmer2()` to get k-mer VarBit
  - Append occurrence count bits using existing occurrence handling logic
  - Return complete ngram_key2 structure
  - Handle occurrence count capping based on `kmersearch_occur_bitlen`

### 2. Cache Lookup Functions Implementation

#### Function: `kmersearch_lookup_kmer2_as_uint_in_global_cache()`
- **Purpose**: Check if uint k-mer exists in global high-frequency cache
- **Parameters**:
  - `uint16/uint32/uint64 kmer2_as_uint`: K-mer to lookup
  - `char *table_name`: Target table name
  - `char *column_name`: Target column name
- **Return**: `bool` - true if k-mer is high-frequency (should be excluded)
- **Implementation**:
  - Use `kmer2_as_uint` directly as hash key (no conversion needed)
  - Search in global cache HTAB structure
  - Return lookup result without VarBit conversion overhead

#### Function: `kmersearch_lookup_kmer2_as_uint_in_parallel_cache()`
- **Purpose**: Check if uint k-mer exists in parallel high-frequency cache
- **Parameters**: Same as global cache function
- **Return**: `bool` - true if k-mer is high-frequency (should be excluded)
- **Implementation**:
  - Use `kmer2_as_uint` directly as dshash key
  - Proper lock management for dshash operations
  - Handle DSM segment access and error conditions

### 3. GIN Extract Value Functions Modification

#### Function: `kmersearch_extract_value_dna2()` Modifications
- **New Implementation**:
  1. Call `kmersearch_extract_dna2_kmer2_as_uint_direct()` to get uint k-mer array
  2. Initialize `KmerOccurrence` tracking structure for deduplication
  3. For each extracted uint k-mer:
     - Check cache using `kmersearch_lookup_kmer2_as_uint_in_*_cache()`
     - Skip if high-frequency (excluded)
     - Track occurrence count using binary search
  4. Create final ngram_key2 array using `kmersearch_create_ngram_key2_from_kmer2_as_uint()`
  5. Return filtered `Datum` array

#### Function: `kmersearch_extract_value_dna4()` Modifications
- **New Implementation**:
  1. Call `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()` to get uint k-mer array
  2. Handle degenerate expansion results (multiple uint k-mers per position)
  3. Initialize `KmerOccurrence` tracking for deduplication
  4. For each expanded uint k-mer:
     - Check cache using `kmersearch_lookup_kmer2_as_uint_in_*_cache()`
     - Skip if high-frequency (excluded)
     - Track occurrence count using binary search
  5. Create final ngram_key2 array using `kmersearch_create_ngram_key2_from_kmer2_as_uint()`
  6. Return filtered `Datum` array

### 4. Cache Key Format Adaptation

#### Current Cache Key Format
- **Storage**: `uint64` hash values calculated from complete ngram_key2 VarBit
- **Generation**: Uses PostgreSQL's `hash_any()` function

#### New Cache Key Format
- **Storage**: Direct `uint16/uint32/uint64` k-mer values
- **Benefits**: No hash calculation overhead, direct comparison
- **Compatibility**: Requires cache rebuilding or dual-format support

### 5. Integration Points

#### Occurrence Count Handling
- **Reference Implementation**: Use existing logic from `kmersearch_extract_dna2_ngram_key2_direct()`
- **Key Features**:
  - `KmerOccurrence` structure for deduplication
  - Binary search for efficient insertion
  - Occurrence count capping based on `kmersearch_occur_bitlen`
  - Memory-efficient processing

#### SIMD Optimization Integration
- **Utilize**: Existing SIMD dispatch in `kmersearch_extract_dna2/4_kmer2_as_uint_direct()`
- **Benefits**: Platform-specific optimizations (AVX2, AVX512, NEON, SVE)
- **Compatibility**: Maintain scalar fallback implementation

#### Error Handling
- **Memory Management**: Proper cleanup of intermediate uint arrays
- **Validation**: Parameter bounds checking and GUC validation
- **Cache Errors**: Graceful fallback to system table lookup

## Implementation Phases

### Phase 1: Core Conversion Functions ✅ COMPLETED
1. ✅ Implement `kmersearch_kmer2_as_uint_to_kmer2()`
2. ✅ Implement `kmersearch_create_ngram_key2_from_kmer2_as_uint()`
3. ✅ Add unit tests for conversion accuracy

### Phase 2: Cache Lookup Functions ✅ COMPLETED
1. ✅ Implement `kmersearch_lookup_kmer2_as_uint_in_global_cache()`
2. ✅ Implement `kmersearch_lookup_kmer2_as_uint_in_parallel_cache()`
3. ✅ Add cache key format adaptation logic

### Phase 3: GIN Function Modifications ✅ COMPLETED
1. ✅ Modify `kmersearch_extract_value_dna2()` to use new flow
2. ✅ Modify `kmersearch_extract_value_dna4()` to use new flow
3. ✅ Update occurrence count handling integration

### Phase 4: Testing and Validation ✅ COMPLETED
1. ✅ Run regression tests to ensure correctness
2. ✅ Performance benchmarking against current implementation
3. ✅ Memory usage analysis and optimization

## Expected Benefits

### Performance Improvements
- **Reduced Memory Allocations**: Fewer VarBit objects created
- **Faster Cache Lookups**: Direct uint comparison vs hash calculation
- **SIMD Utilization**: Better use of optimized extraction functions
- **Early Filtering**: Exclude high-frequency k-mers before VarBit creation

### Code Quality
- **Reduced Duplication**: Unified uint-based processing
- **Better Maintainability**: Clearer separation of concerns
- **Optimized Memory Usage**: More efficient intermediate representations

## Compatibility Considerations

### Cache Format Migration
- **Current**: Hash-based cache keys
- **New**: Direct uint-based cache keys
- **Solution**: Implement dual-format support or cache rebuilding

### API Compatibility
- **External APIs**: No changes to public function signatures
- **Internal APIs**: Modified internal function interfaces
- **Testing**: Comprehensive regression testing required

## Risk Mitigation

### Testing Strategy
- **Unit Tests**: Individual function validation
- **Integration Tests**: Full GIN index construction testing
- **Performance Tests**: Benchmarking against current implementation
- **Regression Tests**: Ensure no functionality breaks

### Rollback Plan
- **Conditional Compilation**: Feature flags for easy rollback
- **Comprehensive Testing**: Validate all edge cases
- **Performance Monitoring**: Track memory usage and query performance

This plan provides a comprehensive approach to optimizing GIN index construction while maintaining functionality and improving performance through better cache utilization and reduced memory overhead.