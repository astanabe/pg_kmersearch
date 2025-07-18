# K-mer Analysis Optimization Plan: From ngram_key2 to kmer2_as_uint Storage

## Overview

This plan optimizes the high-frequency k-mer analysis system by switching from storing ngram_key2 (k-mer + occurrence count) to kmer2_as_uint (k-mer only as integer) in the `kmersearch_highfreq_kmer` table. This reduces memory usage and simplifies the analysis pipeline.

**IMPORTANT**: The core kmer2_as_uint extraction functions are **already implemented** and ready for use:
- ✅ `kmersearch_extract_dna2_kmer2_as_uint_direct()`
- ✅ `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()`
- ✅ SIMD optimizations (AVX512, AVX2, NEON, SVE)
- ✅ Function dispatch tables

## Key Changes Summary

1. **Phase 2 Elimination**: Remove the current Phase 2 that extracts ngram_key2 from data
2. **Direct kmer2 Storage**: Store kmer2 as uint16/uint32/uint64 based on k-mer size
3. **Simplified Pipeline**: Single-pass analysis with direct database storage
4. **Updated Cache System**: Load kmer2_as_uint values into memory caches
5. **Modified GIN Extraction**: Extract kmer2_as_uint, compare with memory cache, then convert to ngram_key2 for GIN index storage

## Important Note: GIN Index Storage Format

**CRITICAL**: While the high-frequency k-mer analysis and caching system will be optimized to use kmer2_as_uint format, the **GIN index itself continues to store ngram_key2** (k-mer + occurrence count) format. The optimization affects only:

- **kmersearch_highfreq_kmer table**: Stores kmer2_as_uint instead of ngram_key
- **Memory caches**: Load and store kmer2_as_uint values
- **High-frequency filtering**: Compare kmer2_as_uint values during GIN index construction

The GIN index construction process will:
1. Extract kmer2_as_uint from DNA sequences
2. Filter high-frequency k-mers using kmer2_as_uint cache lookup
3. **Convert remaining kmer2_as_uint to ngram_key2 format** (add occurrence count)
4. **Store ngram_key2 in the GIN index** (existing format preserved)

## Data Type Mapping Summary

**Database Storage Types (PostgreSQL):**
- k ≤ 8: `smallint` (16-bit signed integer)
- k ≤ 16: `integer` (32-bit signed integer)  
- k ≤ 32: `bigint` (64-bit signed integer)

**Memory Storage Types (C):**
- k ≤ 8: `uint16` (16-bit unsigned integer)
- k ≤ 16: `uint32` (32-bit unsigned integer)
- k ≤ 32: `uint64` (64-bit unsigned integer)

**Conversion Logic:**
- Database → Memory: Use appropriate DatumGetInt16/32/64 and cast to uint16/32/64
- Memory → Database: Use appropriate Int16/32/64GetDatum with proper type conversion
- All values are stored as positive integers (k-mer bit patterns), so signed/unsigned conversion is safe

## 1. Database Schema Changes

### 1.1 kmersearch_highfreq_kmer Table Modification

**Current Schema:**
```sql
CREATE TABLE kmersearch_highfreq_kmer (
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    ngram_key varbit NOT NULL,              -- TO BE REMOVED
    detection_reason text,
    created_at timestamp with time zone DEFAULT now(),
    PRIMARY KEY (table_oid, column_name, ngram_key)
);
```

**New Schema:**
```sql
CREATE TABLE kmersearch_highfreq_kmer (
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    kmer2_as_uint bigint NOT NULL,          -- NEW: kmer2 as uint16/uint32/uint64
    detection_reason text,
    created_at timestamp with time zone DEFAULT now(),
    PRIMARY KEY (table_oid, column_name, kmer2_as_uint)
);
```

**Column Type Selection Logic:**
- k ≤ 8: `smallint` (PostgreSQL database storage) ↔ `uint16` (C memory storage)
- k ≤ 16: `integer` (PostgreSQL database storage) ↔ `uint32` (C memory storage)
- k ≤ 32: `bigint` (PostgreSQL database storage) ↔ `uint64` (C memory storage)

### 1.2 Index Updates

Update indexes to use the new `kmer2_as_uint` column instead of `ngram_key`.

## 2. Core Analysis Function Modifications

### 2.1 kmersearch_perform_highfreq_analysis()

**Location**: `kmersearch_freq.c`

**Current Architecture:**
- Phase 1: Parallel frequency analysis → temp_kmer_final table
- Phase 2: Parallel ngram_key2 extraction → kmersearch_highfreq_kmer table

**New Architecture:**
- Phase 1: Parallel frequency analysis → temp_kmer_final table
- Direct Storage: Insert kmer2_as_uint values directly into kmersearch_highfreq_kmer table

**Key Modifications:**

1. **Remove Phase 2 Infrastructure:**
   - Remove `kmersearch_perform_highfreq_analysis_phase2_parallel()`
   - Remove `kmersearch_perform_highfreq_analysis_phase2_worker()`
   - Remove temporary table creation for Phase 2

2. **Update Direct Storage Logic:**
   ```c
   // After Phase 1 completion, directly insert kmer2_as_uint
   // Instead of: INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, ngram_key)
   // Use: INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, kmer2_as_uint)
   
   // kmer2_as_uint values are obtained directly from:
   // - kmersearch_extract_dna2_kmer2_as_uint_direct() for DNA2 data
   // - kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct() for DNA4 data
   ```

### 2.2 kmersearch_undo_highfreq_analysis()

**Location**: `kmersearch_freq.c`

**Modifications:**
- Update SQL queries to reference `kmer2_as_uint` instead of `ngram_key`
- Update statistics calculation to count `kmer2_as_uint` entries
- Maintain same cleanup logic for metadata tables

## 3. Cache System Modifications

### 3.1 Global Cache Structure Updates

**Current Structure:**
```c
typedef struct HighfreqKmerCache {
    HighfreqCacheKey current_cache_key;
    MemoryContext cache_context;
    HTAB *highfreq_hash;                 // VarBit* → hash mapping
    VarBit **highfreq_kmers;             // Array of ngram_key2
    int highfreq_count;
    bool is_valid;
} HighfreqKmerCache;
```

**New Structure:**
```c
typedef struct HighfreqKmerCache {
    HighfreqCacheKey current_cache_key;
    MemoryContext cache_context;
    HTAB *highfreq_hash;                 // uint64 → boolean mapping
    uint64 *highfreq_kmers;              // Array of kmer2_as_uint (uint16/uint32/uint64)
    int highfreq_count;
    bool is_valid;
} HighfreqKmerCache;
```

**Memory Storage Types:**
- k ≤ 8: Store as `uint16` in memory (16-bit integers)
- k ≤ 16: Store as `uint32` in memory (32-bit integers)
- k ≤ 32: Store as `uint64` in memory (64-bit integers)

### 3.2 Parallel Cache Structure Updates

**Current Structure:**
```c
typedef struct ParallelHighfreqKmerCacheEntry {
    uint64 kmer_hash;                    // Hash of ngram_key2
    int32 frequency_count;
    HighfreqCacheKey cache_key;
} ParallelHighfreqKmerCacheEntry;
```

**New Structure:**
```c
typedef struct ParallelHighfreqKmerCacheEntry {
    uint64 kmer2_as_uint;                // Direct kmer2 as uint16/uint32/uint64
    int32 frequency_count;               // Keep for consistency
    HighfreqCacheKey cache_key;
} ParallelHighfreqKmerCacheEntry;
```

**Memory Storage Types (Same as Global Cache):**
- k ≤ 8: Store as `uint16` in memory (16-bit integers)
- k ≤ 16: Store as `uint32` in memory (32-bit integers)
- k ≤ 32: Store as `uint64` in memory (64-bit integers)

### 3.3 Cache Loading Functions

#### kmersearch_highfreq_kmer_cache_load()

**Location**: `kmersearch_cache.c`

**Current Implementation:**
```sql
SELECT DISTINCT hkm.ngram_key FROM kmersearch_highfreq_kmer hkm
WHERE hkm.table_oid = %u AND hkm.column_name = '%s'
```

**New Implementation:**
```sql
SELECT DISTINCT hkm.kmer2_as_uint FROM kmersearch_highfreq_kmer hkm
WHERE hkm.table_oid = %u AND hkm.column_name = '%s'
```

**Cache Storage Logic:**
```c
// Instead of: hash = kmersearch_ngram_key_to_hash(ngram_key);
// Use: direct storage of kmer2_as_uint from database
// Database → Memory type mapping:
// smallint → uint16, integer → uint32, bigint → uint64

uint64 kmer2_as_uint;
if (kmer_size <= 8) {
    kmer2_as_uint = (uint64)DatumGetInt16(values[0]);  // smallint → uint16
} else if (kmer_size <= 16) {
    kmer2_as_uint = (uint64)DatumGetInt32(values[0]);  // integer → uint32
} else {
    kmer2_as_uint = DatumGetInt64(values[0]);          // bigint → uint64
}
global_highfreq_cache.highfreq_kmers[i] = kmer2_as_uint;
```

#### kmersearch_parallel_highfreq_kmer_cache_load()

**Location**: `kmersearch_cache.c`

**Similar Changes:**
- Update SQL query to select `kmer2_as_uint`
- Store kmer2_as_uint directly in dshash entries with same type mapping:
  - k ≤ 8: smallint (database) → uint16 (memory)
  - k ≤ 16: integer (database) → uint32 (memory)
  - k ≤ 32: bigint (database) → uint64 (memory)
- Remove hash calculation overhead

## 4. GIN Index Extraction Function Modifications

### 4.1 Existing K-mer Extraction Functions

**IMPORTANT**: These functions are **already implemented** and available for use:

#### kmersearch_extract_dna2_kmer2_as_uint_direct()
```c
uint64 *kmersearch_extract_dna2_kmer2_as_uint_direct(
    const unsigned char *dna2_data,
    int32 dna2_len,
    int kmer_size,
    int *num_kmers
);
```

**Return Type**: Always returns `uint64 *` array, but actual values are:
- k ≤ 8: uint16 values (stored in uint64 array)
- k ≤ 16: uint32 values (stored in uint64 array)  
- k ≤ 32: uint64 values (stored in uint64 array)

#### kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()
```c
uint64 *kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(
    const unsigned char *dna4_data,
    int32 dna4_len,
    int kmer_size,
    int *num_kmers
);
```

**Return Type**: Same as DNA2 function - always returns `uint64 *` array with values sized according to k-mer length.

**STATUS**: **✅ ALREADY IMPLEMENTED** - These functions are ready for use in the new architecture.

### 4.2 GIN Extract Value Function Updates

#### kmersearch_extract_value_dna2()

**Location**: `kmersearch_gin.c:33-76`

**Current Flow:**
1. Extract ngram_key2 using `kmersearch_extract_dna2_ngram_key2_direct()`
2. Filter high-frequency k-mers using ngram_key2 comparison
3. Return filtered ngram_key2 array

**New Flow:**
1. Extract kmer2_as_uint using `kmersearch_extract_dna2_kmer2_as_uint_direct()`
2. Filter high-frequency k-mers using kmer2_as_uint comparison with cache
3. **Convert remaining kmer2_as_uint to VarBit format**
4. **Add occurrence count to create final ngram_key2** (follow existing logic in `kmersearch_extract_dna2_ngram_key2_direct()`)
5. **Return filtered ngram_key2 array for GIN index storage**

**IMPORTANT**: The GIN index continues to store ngram_key2 format. Only the high-frequency filtering process is optimized to use kmer2_as_uint.

**Implementation Changes:**
```c
// Phase 1: Extract kmer2_as_uint using EXISTING function
uint64 *kmer2_as_uint_keys = kmersearch_extract_dna2_kmer2_as_uint_direct(
    dna2_data, dna2_len, kmer_size, &num_kmers);

// Phase 2: Filter against high-frequency cache using kmer2_as_uint
for (int i = 0; i < num_kmers; i++) {
    if (kmersearch_lookup_kmer2_as_uint_in_cache(kmer2_as_uint_keys[i])) {
        // Skip high-frequency k-mer (filtering optimization)
        continue;
    }
    
    // Phase 3: Convert to VarBit and add occurrence count → ngram_key2
    // This maintains compatibility with existing GIN index format
    // IMPORTANT: Use the same occurrence counting logic as in existing functions:
    // - kmersearch_extract_dna2_ngram_key2_direct() (lines 787-846)
    // - kmersearch_extract_dna4_ngram_key2_with_expansion_direct() (lines 931-989)
    VarBit *ngram_key2 = kmersearch_create_ngram_key2_from_kmer2_as_uint(
        kmer2_as_uint_keys[i], occurrence_count);
    
    // Add ngram_key2 to result array (GIN index stores ngram_key2 format)
    result_array[result_count++] = ngram_key2;
}
```

#### kmersearch_extract_value_dna4()

**Location**: `kmersearch_gin.c:84-127`

**Similar modifications for DNA4 with degenerate code expansion using the EXISTING `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()` function.**

**IMPORTANT**: Both DNA2 and DNA4 extraction functions will:
1. **Use existing kmer2_as_uint extraction functions** for efficient high-frequency filtering
2. Convert filtered results back to ngram_key2 format for GIN index storage
3. Maintain full compatibility with existing GIN index infrastructure

**STATUS**: **✅ EXTRACTION FUNCTIONS READY** - Only need to integrate with cache lookup and ngram_key2 conversion.

### 4.3 Cache Lookup Functions

#### New Lookup Functions

**kmersearch_lookup_kmer2_as_uint_in_global_cache():**
```c
bool kmersearch_lookup_kmer2_as_uint_in_global_cache(uint64 kmer2_as_uint) {
    if (!global_highfreq_cache.is_valid || global_highfreq_cache.highfreq_count == 0)
        return false;
    
    bool found;
    hash_search(global_highfreq_cache.highfreq_hash, &kmer2_as_uint, HASH_FIND, &found);
    return found;
}
```

**kmersearch_lookup_kmer2_as_uint_in_parallel_cache():**
```c
bool kmersearch_lookup_kmer2_as_uint_in_parallel_cache(uint64 kmer2_as_uint) {
    if (!parallel_cache_hash)
        return false;
    
    ParallelHighfreqKmerCacheEntry *entry = dshash_find(parallel_cache_hash, &kmer2_as_uint, false);
    return (entry != NULL);
}
```

## 5. Utility Functions

### 5.1 Uint to K-mer Conversion

**kmersearch_uint_to_kmer2():**
```c
void kmersearch_uint_to_kmer2(uint64 kmer2_as_uint, int kmer_size, unsigned char *kmer2_data) {
    int total_bits = kmer_size * 2;
    memset(kmer2_data, 0, (total_bits + 7) / 8);
    
    for (int i = 0; i < total_bits; i++) {
        if (kmer2_as_uint & (1ULL << (total_bits - 1 - i))) {
            int byte_idx = i / 8;
            int bit_idx = 7 - (i % 8);
            kmer2_data[byte_idx] |= (1 << bit_idx);
        }
    }
}
```

### 5.2 Occurrence Count Addition

## Understanding Occurrence Count

**Occurrence Count Definition**: The number of times the same kmer2 appears within a single row (sequence).

**Important**: Occurrence count is **0-based offset**, meaning:
- 1st occurrence: `00000000` (0)
- 2nd occurrence: `00000001` (1)  
- 3rd occurrence: `00000010` (2)
- etc.

**Example**: Sequence "ACGTACGT" with k=4:
- Extract k-mers: "ACGT", "CGTA", "GTAC", "TACG", "ACGT"
- Convert to kmer2 (2 bits per base): "00011011", "01101100", "10110001", "11000110", "00011011"
- Add occurrence count to create ngram_key2:
  - 1st "ACGT" (00011011): `0001101100000000` (occurrence count = 0)
  - 1st "CGTA" (01101100): `0110110000000000` (occurrence count = 0)
  - 1st "GTAC" (10110001): `1011000100000000` (occurrence count = 0)
  - 1st "TACG" (11000110): `1100011000000000` (occurrence count = 0)
  - 2nd "ACGT" (00011011): `0001101100000001` (occurrence count = 1)

**kmersearch_create_ngram_key2_from_kmer2_as_uint():**
```c
VarBit *kmersearch_create_ngram_key2_from_kmer2_as_uint(uint64 kmer2_as_uint, int kmer_size, int occurrence_count) {
    // Convert kmer2_as_uint to VarBit format
    VarBit *kmer2_varbit = kmersearch_uint_to_varbit(kmer2_as_uint, kmer_size * 2);
    
    // Add occurrence count to create ngram_key2 format for GIN index
    // NOTE: occurrence_count is 0-based (first occurrence = 0, second = 1, etc.)
    VarBit *ngram_key2 = kmersearch_create_ngram_key2_from_kmer2_and_count(kmer2_varbit, occurrence_count);
    
    pfree(kmer2_varbit);
    return ngram_key2;  // Returns ngram_key2 format for GIN index storage
}
```

**CRITICAL REFERENCE**: The occurrence count addition logic should follow the same pattern as implemented in:
- `kmersearch_extract_dna2_ngram_key2_direct()` (lines 787-846 in kmersearch.c)
- `kmersearch_extract_dna4_ngram_key2_with_expansion_direct()` (lines 931-989 in kmersearch.c)

These existing functions demonstrate the correct method for tracking k-mer occurrences within a single row and converting them to 0-based occurrence counts for ngram_key2 creation.

## 6. SIMD Optimization Updates

### 6.1 New SIMD Functions

**✅ ALREADY AVAILABLE**: SIMD-optimized versions of extraction functions are already implemented:

- `kmersearch_extract_dna2_kmer2_as_uint_direct_avx512()`
- `kmersearch_extract_dna2_kmer2_as_uint_direct_avx2()`
- `kmersearch_extract_dna2_kmer2_as_uint_direct_neon()`
- `kmersearch_extract_dna2_kmer2_as_uint_direct_sve()`

### 6.2 Dispatch Table Updates

**✅ ALREADY AVAILABLE**: Function dispatch tables are already configured to use existing SIMD functions:

```c
// In kmersearch.c - EXISTING dispatch table
static const SIMDCapabilities simd_capabilities[] = {
    {
        .extract_dna2_kmer2_as_uint_direct = kmersearch_extract_dna2_kmer2_as_uint_direct_avx512,
        .extract_dna4_kmer2_as_uint_with_expansion_direct = kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct_avx512,
        // ... other SIMD functions
    },
    // ... other capability entries
};
```

**STATUS**: **✅ NO CHANGES NEEDED** - SIMD dispatch is already implemented and functional.

## 7. Error Handling and Validation

### 7.1 Data Type Validation

**Add validation for k-mer size vs integer type:**
```c
void kmersearch_validate_kmer_size_for_uint_storage(int kmer_size) {
    if (kmer_size <= 8) {
        // Database: smallint, Memory: uint16 - OK
    } else if (kmer_size <= 16) {
        // Database: integer, Memory: uint32 - OK
    } else if (kmer_size <= 32) {
        // Database: bigint, Memory: uint64 - OK
    } else {
        elog(ERROR, "k-mer size %d exceeds maximum supported size for uint storage", kmer_size);
    }
}
```

**Data Type Mapping Validation:**
```c
typedef enum {
    KMER_UINT16 = 1,  // k ≤ 8: smallint ↔ uint16
    KMER_UINT32 = 2,  // k ≤ 16: integer ↔ uint32
    KMER_UINT64 = 3   // k ≤ 32: bigint ↔ uint64
} KmerUintType;

KmerUintType kmersearch_get_kmer_uint_type(int kmer_size) {
    if (kmer_size <= 8) return KMER_UINT16;
    if (kmer_size <= 16) return KMER_UINT32;
    if (kmer_size <= 32) return KMER_UINT64;
    elog(ERROR, "k-mer size %d exceeds maximum supported size", kmer_size);
}
```

### 7.2 Cache Consistency Validation

**Update cache validation to work with kmer2_as_uint:**
```c
bool kmersearch_validate_cache_consistency(HighfreqCacheKey *cache_key) {
    // Validate GUC parameters match cache metadata
    // Update to handle kmer2_as_uint format
    return (cache_key->kmer_size == kmersearch_kmer_size &&
            cache_key->occur_bitlen == kmersearch_occur_bitlen &&
            cache_key->max_appearance_rate == kmersearch_max_appearance_rate &&
            cache_key->max_appearance_nrow == kmersearch_max_appearance_nrow);
}
```

## 8. Testing and Validation

### 8.1 Unit Tests

**Add tests for new functions:**
- `test_kmersearch_uint_to_kmer2()`
- **✅ ALREADY TESTED**: `test_kmersearch_extract_dna2_kmer2_as_uint_direct()`
- **✅ ALREADY TESTED**: `test_kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()`
- `test_kmersearch_create_ngram_key2_from_kmer2_as_uint()`

### 8.2 Integration Tests

**Update existing tests:**
- Modify `09_highfreq_filter.sql` to work with new table schema
- Update `10_parallel_cache.sql` and `11_cache_hierarchy.sql` for new cache structure
- Verify GIN index creation and search functionality

### 8.3 Performance Testing

**Benchmark comparisons:**
- Memory usage: ngram_key2 vs kmer2_as_uint storage
- Cache loading time: VarBit vs uint64 processing
- GIN extraction performance: new vs old pipeline

## 9. Migration Considerations

### 9.1 Fresh Installation Only

**As per CLAUDE.md instructions:**
- No migration scripts needed (new installations only)
- Breaking changes to data types and SQL interfaces are acceptable
- No backward compatibility considerations

### 9.2 Database Schema Updates

**Update SQL extension definition:**
- Modify `pg_kmersearch--1.0.sql` to include new table schema
- Update system table definitions
- Add new function signatures

## 10. Implementation Priority

### Phase 1: Core Infrastructure
1. Database schema changes
2. **✅ ALREADY AVAILABLE**: Direct kmer2_as_uint extraction functions (kmersearch_extract_dna2_kmer2_as_uint_direct, kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct)
3. Cache structure updates

### Phase 2: Analysis Functions
1. Update `kmersearch_perform_highfreq_analysis()`
2. Update `kmersearch_undo_highfreq_analysis()`
3. Remove Phase 2 infrastructure

### Phase 3: Cache System
1. Update cache loading functions
2. Update cache lookup functions
3. Add new validation logic

### Phase 4: GIN Extraction
1. **✅ SKIP**: Extraction functions already implemented
2. Update GIN extract_value functions to use existing extraction functions
3. **✅ SKIP**: SIMD optimizations already included in existing functions

### Phase 5: Testing and Validation
1. Update regression tests
2. Performance benchmarking
3. Integration testing

## 11. Expected Benefits

### 11.1 Memory Efficiency
- **Reduced storage**: kmer2_as_uint (8 bytes max) vs ngram_key2 (variable, often 12+ bytes) in system tables and caches
- **Simplified caching**: Direct uint64 storage vs VarBit handling
- **Faster lookups**: Integer comparison vs bit pattern matching for high-frequency filtering

### 11.2 Performance Improvements
- **Eliminated Phase 2**: Single-pass analysis instead of two-pass
- **Reduced complexity**: Direct integer operations vs VarBit manipulation in cache operations
- **Better cache locality**: Compact uint64 arrays vs VarBit pointer arrays

### 11.3 Simplified Architecture
- **Fewer temporary tables**: Direct database storage from Phase 1
- **Cleaner code**: Removal of Phase 2 infrastructure
- **Maintainability**: Simpler data flow and fewer abstraction layers

### 11.4 Preserved Compatibility
- **GIN index format**: Continues to use ngram_key2 format (no breaking changes to index structure)
- **Search operations**: Existing search functionality remains unchanged
- **SQL interface**: No changes to user-facing query operations

## 12. Risk Assessment

### 12.1 Low Risk
- Database schema changes (new installations only)
- Cache structure updates (well-defined interfaces)
- Function signature changes (breaking changes acceptable)

### 12.2 Medium Risk
- SIMD optimization updates (requires careful testing)
- GIN extraction pipeline changes (complex logic)
- Cache consistency validation (depends on GUC parameter handling)

### 12.3 High Risk
- Performance regression (must benchmark against current implementation)
- Correctness of k-mer extraction (critical for search accuracy)
- Memory management (proper cleanup of new data structures)

## Conclusion

This optimization plan transforms the high-frequency k-mer analysis system from a complex two-phase approach with ngram_key2 storage to a streamlined single-phase approach with kmer2_as_uint storage in system tables and caches. The changes will significantly reduce memory usage, simplify the codebase, and improve performance while maintaining full functionality and search accuracy.

**Key Architecture Decision**: The optimization affects only the high-frequency k-mer analysis and caching layers. **The GIN index continues to store ngram_key2 format**, ensuring complete backward compatibility with existing search operations, query interfaces, and index structures.

The implementation should be done incrementally, with careful testing at each phase to ensure correctness and performance improvements are achieved while preserving the existing GIN index format and search functionality.