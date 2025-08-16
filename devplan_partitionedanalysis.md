# Development Plan: Partitioned Table Support for kmersearch_perform_highfreq_analysis()

## Executive Summary

Enable `kmersearch_perform_highfreq_analysis()` to work correctly with partitioned tables by implementing partition-aware parallel processing that produces **identical results** to analyzing an equivalent non-partitioned table with the same data.

## Current Implementation Issues

### Problem Analysis
1. **Physical Block Scanning Limitation**
   - Current implementation uses `RelationGetNumberOfBlocks()` on the parent table
   - Parent partitioned tables have 0 blocks (metadata only)
   - Workers scan non-existent blocks, resulting in no k-mer analysis

2. **Row Count Mismatch**
   - `SELECT COUNT(*)` correctly aggregates from child partitions
   - Physical scanning only targets parent table blocks
   - Creates inconsistency between reported rows and analyzed data

3. **Parallel Worker Distribution**
   - Workers are assigned to scan blocks of a single relation
   - Cannot leverage partition-level parallelism
   - Inefficient resource utilization for partitioned tables

## Core Requirement

**The implementation MUST produce identical results when analyzing a partitioned table as it would when analyzing a non-partitioned table with the same data.** This means:
- Same total row count
- Same k-mer frequency counts
- Same high-frequency k-mer identification
- Same threshold calculations
- Deterministic and reproducible results

## Proposed Architecture

### High-Level Design
```
kmersearch_perform_highfreq_analysis()
    ├── Detect if table is partitioned
    ├── If partitioned:
    │   ├── Get list of child partitions
    │   ├── Create unified shared k-mer hash table
    │   ├── Launch parallel workers
    │   ├── Workers process blocks from ALL partitions
    │   └── Aggregate in single shared hash table (same as regular tables)
    └── If regular table:
        └── Use existing block-based parallel processing
```

### Parallel Processing Strategy

#### Selected Approach: Unified Block Processing Across Partitions
- Treat all partition blocks as a single logical block sequence
- Workers process blocks from any partition using shared work queue
- Single unified k-mer hash table for all partitions
- Benefits:
  - **Guarantees identical results to non-partitioned tables**
  - Same k-mer counting logic
  - Same frequency threshold application
  - Maintains deterministic behavior

#### Why Not Partition-Level Parallelism?
- Would require complex merging of per-partition k-mer counts
- Risk of rounding errors or count discrepancies
- Difficult to ensure identical threshold calculations
- Cannot guarantee identical results

## Implementation Plan

### Phase 1: Infrastructure (Week 1)

#### 1.1 Partition Detection
```c
typedef enum {
    KMERSEARCH_TABLE_REGULAR,
    KMERSEARCH_TABLE_PARTITIONED,
    KMERSEARCH_TABLE_PARTITION_CHILD
} KmerSearchTableType;

KmerSearchTableType kmersearch_get_table_type(Oid table_oid);
List *kmersearch_get_partition_oids(Oid parent_oid);
```

#### 1.2 Shared State Extension
```c
typedef struct PartitionBlockInfo {
    Oid partition_oid;
    BlockNumber start_block;
    BlockNumber end_block;
} PartitionBlockInfo;

typedef struct KmerAnalysisSharedState {
    /* Existing fields */
    ...
    
    /* Partition-specific fields */
    bool is_partitioned;
    int num_partitions;
    PartitionBlockInfo *partition_blocks;  /* Array of partition block ranges */
    BlockNumber total_blocks_all_partitions;
    
    /* Unified work queue for all blocks across partitions */
    pg_atomic_uint32 next_global_block;    /* Global block counter */
    
    /* Single unified k-mer hash table (same as regular tables) */
    /* No per-partition statistics needed - maintain unified counts */
} KmerAnalysisSharedState;
```

### Phase 2: Core Implementation (Week 2)

#### 2.1 Modified Main Function Flow
```c
KmerAnalysisResult
kmersearch_perform_highfreq_analysis_parallel(Oid table_oid, ...) {
    KmerSearchTableType table_type = kmersearch_get_table_type(table_oid);
    
    if (table_type == KMERSEARCH_TABLE_PARTITIONED) {
        return kmersearch_analyze_partitioned_table(table_oid, ...);
    } else {
        return kmersearch_analyze_regular_table(table_oid, ...);
    }
}
```

#### 2.2 Partitioned Table Analysis
```c
static KmerAnalysisResult
kmersearch_analyze_partitioned_table(Oid parent_oid, ...) {
    List *partition_oids = kmersearch_get_partition_oids(parent_oid);
    BlockNumber total_blocks = 0;
    int partition_idx = 0;
    
    /* Build unified block map for all partitions */
    shared_state->is_partitioned = true;
    shared_state->num_partitions = list_length(partition_oids);
    shared_state->partition_blocks = kmersearch_allocate_partition_block_array();
    
    foreach(lc, partition_oids) {
        Oid part_oid = lfirst_oid(lc);
        Relation part_rel = table_open(part_oid, AccessShareLock);
        BlockNumber part_blocks = RelationGetNumberOfBlocks(part_rel);
        
        shared_state->partition_blocks[partition_idx].partition_oid = part_oid;
        shared_state->partition_blocks[partition_idx].start_block = total_blocks;
        shared_state->partition_blocks[partition_idx].end_block = total_blocks + part_blocks - 1;
        
        total_blocks += part_blocks;
        partition_idx++;
        table_close(part_rel, AccessShareLock);
    }
    
    shared_state->total_blocks_all_partitions = total_blocks;
    shared_state->next_global_block = 0;
    
    /* Use same k-mer hash table initialization as regular tables */
    /* Launch workers - they will process blocks from all partitions */
    LaunchParallelWorkers(pcxt);
    WaitForParallelWorkersToFinish(pcxt);
    
    /* Results are already in unified hash table - same as regular tables */
    return kmersearch_extract_analysis_results(shared_state);
}
```

#### 2.3 Worker Processing Logic
```c
static void
kmersearch_highfreq_analysis_worker(dsm_segment *seg, shm_toc *toc) {
    BlockNumber global_block;
    
    while ((global_block = kmersearch_claim_next_global_block(shared_state)) < 
           shared_state->total_blocks_all_partitions) {
        
        if (shared_state->is_partitioned) {
            /* Map global block to partition and local block */
            PartitionBlockMapping mapping = kmersearch_map_global_to_partition_block(
                global_block, shared_state);
            
            /* Process the block from the specific partition */
            kmersearch_process_partition_block(mapping.partition_oid, 
                                              mapping.local_block_number);
        } else {
            /* Regular table - process block directly */
            kmersearch_process_table_block(shared_state->table_oid, global_block);
        }
        
        /* K-mer counting uses same unified hash table for both cases */
        kmersearch_update_unified_kmer_counts(extracted_kmers);
    }
}

static PartitionBlockMapping
kmersearch_map_global_to_partition_block(BlockNumber global_block, 
                                        KmerAnalysisSharedState *state) {
    for (int i = 0; i < state->num_partitions; i++) {
        if (global_block >= state->partition_blocks[i].start_block &&
            global_block <= state->partition_blocks[i].end_block) {
            return (PartitionBlockMapping) {
                .partition_oid = state->partition_blocks[i].partition_oid,
                .local_block_number = global_block - 
                                    state->partition_blocks[i].start_block
            };
        }
    }
    elog(ERROR, "Invalid global block number");
}
```

### Phase 3: Result Validation & Optimization (Week 3)

#### 3.1 Result Consistency Validation
```c
/* Validation function to ensure identical results */
bool kmersearch_validate_partitioned_vs_regular_results(
    KmerAnalysisResult partitioned_result,
    KmerAnalysisResult regular_result) {
    
    return (partitioned_result.total_rows == regular_result.total_rows &&
            partitioned_result.highfreq_kmers_count == regular_result.highfreq_kmers_count &&
            fabs(partitioned_result.max_appearance_rate_used - 
                 regular_result.max_appearance_rate_used) < 0.0001 &&
            partitioned_result.max_appearance_nrow_used == 
                 regular_result.max_appearance_nrow_used);
}
```

#### 3.2 Performance Optimization
- Block prefetching for better cache utilization
- Batch processing of small partitions
- Optimize global-to-partition block mapping with binary search

#### 3.3 Memory Management
- Use same memory pools as regular table processing
- No additional memory overhead for partition metadata
- Maintain existing memory limits

### Phase 4: Testing & Validation (Week 4)

#### 4.1 Test Scenarios
1. **Basic Functionality**
   - Small partitioned table (2-4 partitions)
   - Large partitioned table (100+ partitions)
   - Mixed partition sizes

2. **Edge Cases**
   - Empty partitions
   - Single partition with all data
   - Partitions with different DNA types (DNA2/DNA4)

3. **Performance Tests**
   - Compare with non-partitioned equivalent
   - Parallel worker scaling (1, 2, 4, 8 workers)
   - Memory usage profiling

4. **Regression Tests**
   - Ensure regular tables still work
   - Verify backward compatibility
   - Check GUC parameter handling

#### 4.2 Test Implementation - Result Identity Verification
```sql
-- Create identical data in both partitioned and regular tables
CREATE TABLE test_regular (
    id serial PRIMARY KEY,
    sequence dna2 NOT NULL
);

CREATE TABLE test_partitioned (
    id serial,
    sequence dna2 NOT NULL
) PARTITION BY HASH (id);

-- Create 4 partitions
CREATE TABLE test_partitioned_0 PARTITION OF test_partitioned
    FOR VALUES WITH (modulus 4, remainder 0);
CREATE TABLE test_partitioned_1 PARTITION OF test_partitioned
    FOR VALUES WITH (modulus 4, remainder 1);
CREATE TABLE test_partitioned_2 PARTITION OF test_partitioned
    FOR VALUES WITH (modulus 4, remainder 2);
CREATE TABLE test_partitioned_3 PARTITION OF test_partitioned
    FOR VALUES WITH (modulus 4, remainder 3);

-- Insert IDENTICAL data into both tables
WITH test_data AS (
    SELECT i as id, generate_random_dna2(100) as sequence 
    FROM generate_series(1, 10000) i
)
INSERT INTO test_regular SELECT * FROM test_data;

WITH test_data AS (
    SELECT i as id, generate_random_dna2(100) as sequence 
    FROM generate_series(1, 10000) i
)
INSERT INTO test_partitioned SELECT * FROM test_data;

-- Run analysis on both tables
SELECT * FROM kmersearch_perform_highfreq_analysis('test_regular', 'sequence') 
INTO regular_result;

SELECT * FROM kmersearch_perform_highfreq_analysis('test_partitioned', 'sequence') 
INTO partitioned_result;

-- Verify IDENTICAL results
SELECT 
    CASE 
        WHEN r.total_rows = p.total_rows AND
             r.highfreq_kmers_count = p.highfreq_kmers_count AND
             r.max_appearance_rate_used = p.max_appearance_rate_used AND
             r.max_appearance_nrow_used = p.max_appearance_nrow_used
        THEN 'PASS: Results are identical'
        ELSE 'FAIL: Results differ'
    END as test_result,
    r.* as regular,
    p.* as partitioned
FROM regular_result r, partitioned_result p;

-- Verify k-mer data is identical
SELECT 
    CASE 
        WHEN COUNT(*) = 0 THEN 'PASS: K-mer data identical'
        ELSE 'FAIL: K-mer data differs'
    END as kmer_comparison
FROM (
    SELECT kmer_data FROM kmersearch_highfreq_kmer 
    WHERE table_oid = 'test_regular'::regclass
    EXCEPT
    SELECT kmer_data FROM kmersearch_highfreq_kmer 
    WHERE table_oid = 'test_partitioned'::regclass
) diff;
```

## Technical Considerations

### 1. Transaction Isolation
- Maintain consistent snapshot across all partitions
- Handle concurrent modifications during analysis

### 2. Lock Management
- Acquire appropriate locks on all partitions
- Prevent deadlocks with ordered locking

### 3. Error Handling
- Graceful handling of partition drops during analysis
- Worker failure recovery

### 4. Backward Compatibility
- Maintain existing API
- Transparent upgrade path for existing users

## Performance Impact

### Expected Benefits
- Transparent support for partitioned tables
- Same parallelism benefits as regular tables
- No additional memory overhead (uses same hash table structure)

### Potential Overhead
- Minimal: Global-to-partition block mapping computation
- Negligible: Partition metadata storage in shared state
- No impact on result accuracy or consistency

## Implementation Timeline

| Week | Phase | Deliverables |
|------|-------|--------------|
| 1 | Infrastructure | Partition detection, shared state design |
| 2 | Core Implementation | Partitioned table analysis, worker logic |
| 3 | Optimization | Load balancing, memory management |
| 4 | Testing | Test suite, performance validation |

## Risk Mitigation

| Risk | Mitigation Strategy |
|------|-------------------|
| Complex partition hierarchies | Start with simple hash/range partitions |
| Performance regression | Maintain separate code paths |
| Memory exhaustion | Implement memory limits and spilling |
| Uneven partition sizes | Dynamic work stealing between workers |

## Success Criteria

1. **PRIMARY REQUIREMENT**: Results from partitioned table analysis are 100% identical to non-partitioned table analysis with same data
2. `kmersearch_perform_highfreq_analysis()` successfully analyzes partitioned tables
3. All k-mer counts and frequency calculations match exactly
4. High-frequency k-mer identification is deterministic and consistent
5. No regression in regular table analysis performance
6. All existing tests pass
7. New partition-specific identity verification tests pass

## Future Enhancements

1. **Incremental Analysis**
   - Analyze only modified partitions
   - Maintain per-partition analysis timestamps

2. **Partition Pruning**
   - Skip partitions based on analysis predicates
   - Leverage partition constraints

3. **Distributed Processing**
   - Support for distributed PostgreSQL (Citus, etc.)
   - Cross-node parallel processing

4. **Adaptive Parallelism**
   - Dynamic worker allocation based on partition characteristics
   - Cost-based optimization for worker distribution