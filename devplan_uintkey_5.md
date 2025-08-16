# Development Plan: Complete Migration to Uintkey-based Implementation

## Overview
This document outlines the complete migration from VarBit-based ngram_key2 to direct uintkey (uint16/uint32/uint64) implementation for the pg_kmersearch extension. The goal is to eliminate unnecessary memory copying and data conversions by using native integer types throughout the processing pipeline.

## New and Modified Components

### New SQL Operator Classes
```sql
-- For DNA2 type
kmersearch_dna2_gin_ops_int2  -- STORAGE int2 (k <= 8)
kmersearch_dna2_gin_ops_int4  -- STORAGE int4 (k <= 16)
kmersearch_dna2_gin_ops_int8  -- STORAGE int8 (k <= 32)

-- For DNA4 type
kmersearch_dna4_gin_ops_int2  -- STORAGE int2 (k <= 8)
kmersearch_dna4_gin_ops_int4  -- STORAGE int4 (k <= 16)
kmersearch_dna4_gin_ops_int8  -- STORAGE int8 (k <= 32)
```

### New C Functions

#### Extract Value Functions (kmersearch_gin.c)
```c
Datum kmersearch_extract_value_dna2_int2(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna2_int4(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna2_int8(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna4_int2(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna4_int4(PG_FUNCTION_ARGS);
Datum kmersearch_extract_value_dna4_int8(PG_FUNCTION_ARGS);
```

#### Extract Query Functions (kmersearch_gin.c)
```c
Datum kmersearch_extract_query_int2(PG_FUNCTION_ARGS);
Datum kmersearch_extract_query_int4(PG_FUNCTION_ARGS);
Datum kmersearch_extract_query_int8(PG_FUNCTION_ARGS);
```

#### Consistent Functions (kmersearch_gin.c)
```c
Datum kmersearch_consistent_int2(PG_FUNCTION_ARGS);
Datum kmersearch_consistent_int4(PG_FUNCTION_ARGS);
Datum kmersearch_consistent_int8(PG_FUNCTION_ARGS);
```

#### Cache Functions (kmersearch_cache.c)
```c
void* get_cached_query_uintkey(const char *query_string, int k_size, int *nkeys);
int get_cached_actual_min_score_uintkey(void *uintkeys, int nkeys, int k_size);
int get_cached_actual_min_score_datum_int2(Datum *queryKeys, int nkeys);
int get_cached_actual_min_score_datum_int4(Datum *queryKeys, int nkeys);
int get_cached_actual_min_score_datum_int8(Datum *queryKeys, int nkeys);
```

#### Helper Functions (kmersearch_gin.c)
```c
void* filter_uintkey_and_set_actual_min_score(void *uintkeys, int *nkeys, 
                                               const char *query_string, int k_size);
bool kmersearch_is_uintkey_highfreq(uint64 uintkey, int k_size);
```

### Modified Existing Functions

#### Functions to Remove/Replace
```c
// Remove these VarBit-based functions
VarBit **get_cached_query_kmer()           -> get_cached_query_uintkey()
VarBit **filter_ngram_key2_and_set_actual_min_score() -> filter_uintkey_and_set_actual_min_score()
int get_cached_actual_min_score_uintarray() -> get_cached_actual_min_score_uintkey()
int get_cached_actual_min_score_or_error()  -> Remove (not needed)

// Keep but modify
int get_cached_actual_min_score_datum()     -> Split into int2/int4/int8 versions
```

### Modified Cache Structures (kmersearch.h)

```c
// query-kmer cache entry - store uintkeys instead of VarBit
typedef struct QueryPatternCacheEntry {
    uint64 hash_key;
    char *query_string;
    int k_size;
    void *extracted_uintkeys;  // Points to uint16/uint32/uint64 array
    int kmer_count;
} QueryPatternCacheEntry;
```

## Implementation Details

### 1. Uintkey Storage Format

The uintkey format combines k-mer bits with occurrence count:
- **k <= 8**: uint16 (k*2 bits for k-mer + remaining bits for count)
- **k <= 16**: uint32 (k*2 bits for k-mer + remaining bits for count)
- **k <= 32**: uint64 (k*2 bits for k-mer + remaining bits for count)

### 2. GIN Index Integration

Each operator class uses the appropriate integer storage type:

```sql
CREATE OPERATOR CLASS kmersearch_dna2_gin_ops_int2
    FOR TYPE DNA2 USING gin AS
        OPERATOR 1 =% (DNA2, text),
        FUNCTION 1 btint2cmp(int2, int2),
        FUNCTION 2 kmersearch_extract_value_dna2_int2(DNA2, internal),
        FUNCTION 3 kmersearch_extract_query_int2(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent_int2(internal, int2, text, int4, internal, internal),
        STORAGE int2;
```

### 3. Processing Flow

#### Extract Value (Index Build)
1. Extract DNA sequence from input
2. Use `kmersearch_extract_uintkey_from_dna2/4()` to get uintkey array
3. Filter high-frequency k-mers if enabled
4. Return Datum array of int2/int4/int8 values directly to GIN

#### Extract Query (Search)
1. Convert query text to uintkeys using `kmersearch_extract_uintkey_from_text()`
2. Filter high-frequency k-mers and cache actual_min_score
3. Return Datum array of int2/int4/int8 values directly to GIN

#### Consistent (Index Scan)
1. Receive Datum array from GIN (no conversion needed)
2. Pass directly to `get_cached_actual_min_score_datum_intX()`
3. Compare match count with cached actual_min_score

### 4. Cache Management

#### query-kmer cache
- Key: Query string hash
- Value: Uintkey array (void* pointing to uint16/uint32/uint64 array)
- No VarBit conversion needed

#### Actual Min Score Cache
- Key: Hash of filtered uintkey array
- Value: Calculated actual_min_score
- Direct processing of Datum arrays in consistent functions

### 5. Memory Optimization

Key optimizations:
- **No VarBit allocations**: Direct use of integer types
- **No data copying in consistent**: Datum arrays used as-is
- **Efficient caching**: Native integer arrays cached directly
- **SIMD compatibility**: Integer arrays work directly with SIMD instructions

### 6. Backward Compatibility

Since the extension assumes new installations only (as per CLAUDE.md):
- Existing indexes must be dropped and recreated
- Users select appropriate operator class based on their k value:
  - k <= 8: Use `*_int2` operator class
  - k <= 16: Use `*_int4` operator class
  - k <= 32: Use `*_int8` operator class

### 7. User Impact

```sql
-- Example: Creating index for k=6
SET kmersearch.kmer_size = 6;
CREATE INDEX idx_seq ON mytable USING gin(sequence kmersearch_dna2_gin_ops_int2);

-- Example: Creating index for k=12
SET kmersearch.kmer_size = 12;
CREATE INDEX idx_seq ON mytable USING gin(sequence kmersearch_dna2_gin_ops_int4);
```

## Implementation Order

1. **Phase 1**: Implement uintkey-based cache functions
   - `get_cached_query_uintkey()`
   - `get_cached_actual_min_score_uintkey()`
   - `get_cached_actual_min_score_datum_int2/4/8()`

2. **Phase 2**: Implement GIN support functions
   - Extract value functions for int2/int4/int8
   - Extract query functions for int2/int4/int8
   - Consistent functions for int2/int4/int8

3. **Phase 3**: Update SQL definitions
   - Create new operator classes
   - Register new functions

4. **Phase 4**: Testing and validation
   - Verify correct k-mer extraction
   - Validate high-frequency filtering
   - Performance benchmarking

## Benefits

1. **Performance**: Elimination of VarBit overhead and memory copying
2. **Memory**: Reduced memory footprint with native integer types
3. **Simplicity**: Direct integer comparisons instead of VarBit operations
4. **SIMD**: Better alignment with SIMD operations for future optimizations

## Notes

- The existing VarBit-based operator classes will be kept for backward compatibility
- Users must explicitly choose the appropriate integer-based operator class
- The k value must be set before creating the index and must match the operator class