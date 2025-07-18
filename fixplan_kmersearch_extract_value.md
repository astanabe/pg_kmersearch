# Fix Plan: kmersearch_extract_value Functions for kmer2_as_uint Integration

## Overview

The current implementation of `kmersearch_extract_value_dna2()` and `kmersearch_extract_value_dna4()` functions in `kmersearch_gin.c` uses the old approach of creating ngram_key2 first and then filtering high-frequency k-mers. This needs to be updated to align with the new kmer2_as_uint-based architecture.

## Current Implementation Issues

### 1. Outdated Workflow
**Current workflow:**
1. Extract ngram_key2 (kmer2 + occurrence count) from DNA data
2. Filter high-frequency k-mers using ngram_key2 comparison
3. Return filtered ngram_key2 array

**Required new workflow:**
1. Extract kmer2_as_uint from DNA data using existing functions
2. Filter high-frequency k-mers using kmer2_as_uint comparison with cache
3. Convert remaining kmer2_as_uint back to VarBit kmer2 format
4. Add occurrence count to create ngram_key2
5. Return ngram_key2 array for GIN index

### 2. Missing Utility Functions
The following utility functions are required but not implemented:
- `kmersearch_kmer2_as_uint_to_kmer2()` - Convert kmer2_as_uint to VarBit kmer2
- `kmersearch_create_ngram_key2_from_kmer2_as_uint()` - Create ngram_key2 from kmer2_as_uint
- `kmersearch_lookup_kmer2_as_uint_in_global_cache()` - Global cache lookup
- `kmersearch_lookup_kmer2_as_uint_in_parallel_cache()` - Parallel cache lookup

### 3. Existing Functions to Utilize
The following functions are already implemented and ready for use:
- `kmersearch_extract_dna2_kmer2_as_uint_direct()` - Extract kmer2_as_uint from DNA2
- `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()` - Extract kmer2_as_uint from DNA4

## Required Changes

### 1. Implement Utility Functions

#### 1.1 kmersearch_kmer2_as_uint_to_kmer2()
**Location**: `kmersearch_kmer.c`
**Purpose**: Convert kmer2_as_uint to VarBit kmer2 format (1 nucleotide = 2 bits)

```c
VarBit *kmersearch_kmer2_as_uint_to_kmer2(uint64 kmer2_as_uint, int kmer_size)
{
    VarBit *result;
    int total_bits = kmer_size * 2;  // 2 bits per nucleotide
    int byte_size = (total_bits + 7) / 8;
    int varbit_size = VARHDRSZ + byte_size;
    
    result = (VarBit *) palloc0(varbit_size);
    SET_VARSIZE(result, varbit_size);
    VARBITTOTALLEN(result) = total_bits;
    
    unsigned char *data = VARBITS(result);
    
    // Convert kmer2_as_uint to VarBit representation
    for (int i = 0; i < total_bits; i++) {
        if (kmer2_as_uint & (1ULL << (total_bits - 1 - i))) {
            int byte_idx = i / 8;
            int bit_idx = 7 - (i % 8);
            data[byte_idx] |= (1 << bit_idx);
        }
    }
    
    return result;
}
```

#### 1.2 kmersearch_create_ngram_key2_from_kmer2_as_uint()
**Location**: `kmersearch_kmer.c`
**Purpose**: Create ngram_key2 from kmer2_as_uint with occurrence count

```c
VarBit *kmersearch_create_ngram_key2_from_kmer2_as_uint(uint64 kmer2_as_uint, int kmer_size, int occurrence_count)
{
    // Convert kmer2_as_uint to VarBit kmer2 format
    VarBit *kmer2_varbit = kmersearch_kmer2_as_uint_to_kmer2(kmer2_as_uint, kmer_size);
    
    // Create ngram_key2 with occurrence count using existing function
    VarBit *ngram_key2 = kmersearch_create_ngram_key2(kmer2_varbit, occurrence_count);
    
    pfree(kmer2_varbit);
    return ngram_key2;
}
```

#### 1.3 Cache Lookup Functions
**Location**: `kmersearch_cache.c`
**Purpose**: Lookup kmer2_as_uint in cache structures

```c
bool kmersearch_lookup_kmer2_as_uint_in_global_cache(uint64 kmer2_as_uint)
{
    if (!global_highfreq_cache.is_valid || global_highfreq_cache.highfreq_count == 0)
        return false;
    
    bool found;
    HighfreqKmerHashEntry *entry = (HighfreqKmerHashEntry *) hash_search(
        global_highfreq_cache.highfreq_hash, 
        &kmer2_as_uint, 
        HASH_FIND, 
        &found
    );
    
    return found;
}

bool kmersearch_lookup_kmer2_as_uint_in_parallel_cache(uint64 kmer2_as_uint)
{
    if (!parallel_highfreq_cache || !parallel_highfreq_cache->is_initialized)
        return false;
    
    ParallelHighfreqKmerCacheEntry *entry = dshash_find(
        parallel_highfreq_cache->hash, 
        &kmer2_as_uint, 
        false
    );
    
    return (entry != NULL);
}
```

### 2. Update GIN Extract Value Functions

#### 2.1 kmersearch_extract_value_dna2()
**Location**: `kmersearch_gin.c`

**New implementation approach:**
```c
Datum kmersearch_extract_value_dna2(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *dna = (kmersearch_dna2 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    int k = kmersearch_kmer_size;
    
    if (k < 4 || k > 32)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 32")));
    
    // Step 1: Extract kmer2_as_uint values from DNA2 data using existing function
    uint64 *kmer2_as_uint_array;
    int total_kmers;
    void *output;
    
    kmersearch_extract_dna2_kmer2_as_uint_direct((VarBit *)dna, k, &output, &total_kmers);
    kmer2_as_uint_array = (uint64 *)output;
    
    if (!kmer2_as_uint_array || total_kmers == 0) {
        *nkeys = 0;
        PG_RETURN_POINTER(NULL);
    }
    
    // Step 2: Filter high-frequency k-mers using kmer2_as_uint comparison
    uint64 *filtered_kmer2_array;
    int filtered_count;
    
    if (kmersearch_preclude_highfreq_kmer) {
        filtered_kmer2_array = kmersearch_filter_highfreq_kmer2_as_uint(
            kmer2_as_uint_array, total_kmers, &filtered_count, k);
    } else {
        filtered_kmer2_array = kmer2_as_uint_array;
        filtered_count = total_kmers;
    }
    
    // Step 3: Convert filtered kmer2_as_uint back to ngram_key2 format
    keys = (Datum *) palloc(filtered_count * sizeof(Datum));
    
    for (int i = 0; i < filtered_count; i++) {
        VarBit *ngram_key2 = kmersearch_create_ngram_key2_from_kmer2_as_uint(
            filtered_kmer2_array[i], k, 0);  // occurrence_count = 0 for GIN index
        keys[i] = PointerGetDatum(ngram_key2);
    }
    
    *nkeys = filtered_count;
    
    // Cleanup
    if (filtered_kmer2_array != kmer2_as_uint_array)
        pfree(filtered_kmer2_array);
    pfree(kmer2_as_uint_array);
    
    PG_RETURN_POINTER(keys);
}
```

#### 2.2 kmersearch_extract_value_dna4()
**Location**: `kmersearch_gin.c`

**Similar implementation with DNA4-specific k-mer extraction:**
```c
Datum kmersearch_extract_value_dna4(PG_FUNCTION_ARGS)
{
    kmersearch_dna4 *dna = (kmersearch_dna4 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    int k = kmersearch_kmer_size;
    
    if (k < 4 || k > 32)
        ereport(ERROR, (errmsg("k-mer length must be between 4 and 32")));
    
    // Step 1: Extract kmer2_as_uint values from DNA4 data with degenerate expansion
    uint64 *kmer2_as_uint_array;
    int total_kmers;
    void *output;
    
    kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct((VarBit *)dna, k, &output, &total_kmers);
    kmer2_as_uint_array = (uint64 *)output;
    
    if (!kmer2_as_uint_array || total_kmers == 0) {
        *nkeys = 0;
        PG_RETURN_POINTER(NULL);
    }
    
    // Step 2: Filter high-frequency k-mers using kmer2_as_uint comparison
    uint64 *filtered_kmer2_array;
    int filtered_count;
    
    if (kmersearch_preclude_highfreq_kmer) {
        filtered_kmer2_array = kmersearch_filter_highfreq_kmer2_as_uint(
            kmer2_as_uint_array, total_kmers, &filtered_count, k);
    } else {
        filtered_kmer2_array = kmer2_as_uint_array;
        filtered_count = total_kmers;
    }
    
    // Step 3: Convert filtered kmer2_as_uint back to ngram_key2 format
    keys = (Datum *) palloc(filtered_count * sizeof(Datum));
    
    for (int i = 0; i < filtered_count; i++) {
        VarBit *ngram_key2 = kmersearch_create_ngram_key2_from_kmer2_as_uint(
            filtered_kmer2_array[i], k, 0);  // occurrence_count = 0 for GIN index
        keys[i] = PointerGetDatum(ngram_key2);
    }
    
    *nkeys = filtered_count;
    
    // Cleanup
    if (filtered_kmer2_array != kmer2_as_uint_array)
        pfree(filtered_kmer2_array);
    pfree(kmer2_as_uint_array);
    
    PG_RETURN_POINTER(keys);
}
```

### 3. Required Supporting Functions

#### 3.1 kmersearch_filter_highfreq_kmer2_as_uint()
**Location**: `kmersearch_gin.c`
**Purpose**: Filter high-frequency k-mers using kmer2_as_uint comparison

```c
uint64 *kmersearch_filter_highfreq_kmer2_as_uint(uint64 *original_array, int original_count, 
                                                int *filtered_count, int k)
{
    uint64 *filtered_array;
    int count = 0;
    
    filtered_array = (uint64 *) palloc(original_count * sizeof(uint64));
    
    for (int i = 0; i < original_count; i++) {
        bool is_highfreq = false;
        
        // Check against global cache
        if (global_highfreq_cache.is_valid) {
            is_highfreq = kmersearch_lookup_kmer2_as_uint_in_global_cache(original_array[i]);
        }
        
        // Check against parallel cache if needed
        if (!is_highfreq && (kmersearch_force_use_parallel_highfreq_kmer_cache || IsParallelWorker())) {
            is_highfreq = kmersearch_lookup_kmer2_as_uint_in_parallel_cache(original_array[i]);
        }
        
        if (!is_highfreq) {
            filtered_array[count++] = original_array[i];
        }
    }
    
    *filtered_count = count;
    return filtered_array;
}
```

## Implementation Priority

1. **High Priority**: Implement utility functions for kmer2_as_uint to VarBit conversion
2. **High Priority**: Update GIN extract_value functions to use new workflow
3. **Medium Priority**: Implement cache lookup functions for kmer2_as_uint
4. **Medium Priority**: Update filtering functions to work with kmer2_as_uint arrays

## Testing Requirements

1. Verify that GIN index creation works correctly with new extract_value functions
2. Test high-frequency k-mer filtering during index creation
3. Ensure compatibility with existing search operations
4. Validate that occurrence count handling is correct

## Notes

- The occurrence count for GIN index entries is typically 0 (first occurrence)
- Multiple occurrences of the same k-mer within a sequence should be handled appropriately
- Memory management must be careful to avoid leaks with the new VarBit creation functions
- All new functions should have proper error handling and validation