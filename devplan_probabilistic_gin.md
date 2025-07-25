# PGIN (Probabilistic GIN) Index Development Plan

## Overview

This document outlines the development plan for implementing PGIN as a new PostgreSQL access method, independent of traditional GIN indexes. PGIN integrates probabilistic data structures (Count-min sketch, Feature hashing, and Bloom filters) to achieve significant performance improvements and memory reduction while maintaining search accuracy through PostgreSQL's recheck mechanism.

## Motivation

The current GIN index implementation faces scalability challenges at TB-scale:
- **Index size**: 2-8x larger than raw data (500GB-2TB for 1 billion rows)
- **Memory constraints**: TB-scale indexes exceed available RAM
- **I/O bottleneck**: Frequent disk access for large posting lists
- **Search latency**: Minutes to hours for complex queries on large datasets
- **Storage overhead**: Combining probabilistic structures with traditional GIN would increase storage requirements

**Key Innovation**: PGIN eliminates the need for traditional GIN indexes entirely, using PostgreSQL's recheck mechanism for accuracy while achieving massive storage savings.

## Solution Architecture

### PGIN Independent Architecture

```
Query → Feature Hashing → Bloom Filter → Count-min Sketch → PGIN Index → Candidate TIDs → Recheck → Results
   ↓           ↓              ↓               ↓              ↓              ↓           ↓
4^16 keys → 10^6 keys → Existence → Frequency Est. → Reduced Index → TID List → Exact Score → Final Score
```

**No Traditional GIN Required**: PGIN operates entirely independently, using PostgreSQL's recheck mechanism for final accuracy verification.

### Core Components

#### 1. Feature Hashing Layer
- **Purpose**: Reduce k-mer key space from 4^16 to ~10^6 dimensions
- **Method**: Multiple hash functions with XOR combination
- **Memory**: ~8MB (10^6 × 8 bytes)
- **False positive rate**: Controlled through hash function selection

```c
typedef struct FeatureHashReducer {
    uint32_t target_dimensions;    // ~10^6
    uint32_t* hash_seeds;         // 2-4 hash functions
    int num_hash_functions;
} FeatureHashReducer;
```

#### 2. Bloom Filter Cascade
- **Coarse Filter**: High-speed, 10% false positive rate
- **Fine Filter**: Precise filtering, 1% false positive rate  
- **Quality Filter**: High-score candidate identification
- **Memory**: ~100MB total for all filters

```c
typedef struct BloomFilterCascade {
    BloomFilter* coarse_filter;    // Fast screening
    BloomFilter* fine_filter;      // Precision filtering
    BloomFilter* quality_filter;   // High-quality candidates
} BloomFilterCascade;
```

#### 3. Count-min Sketch
- **Purpose**: Multi-stage k-mer frequency management
- **Dimensions**: 64K width × 8 depth  
- **Memory**: ~16MB (64K × 8 × 4 bytes)
- **Property**: Over-estimates but never under-estimates

**Three Usage Scenarios**:
1. **Index Construction**: High-frequency k-mer filtering (replaces max_appearance_rate)
2. **Search Filtering**: Candidate pre-scoring for early elimination
3. **Recheck Optimization**: Priority ordering for exact score calculation

```c
typedef struct CountMinSketch {
    uint32_t width;        // 2^16
    uint32_t depth;        // 4-8 hash functions
    uint32_t** counters;   // [depth][width] counter array
    uint32_t* hash_seeds;  // Independent hash functions
} CountMinSketch;
```

## Implementation Phases

### Phase 1: PGIN Access Method Foundation
**Timeline**: 3-4 months

**Goals**:
- Implement new PostgreSQL access method "pgin"
- Create basic PGIN index structure with feature hashing
- Implement recheck mechanism for exact scoring
- Achieve basic functionality with 50% storage reduction

**Deliverables**:
- New access method handler in `kmersearch_pgin.c`
- Feature hashing implementation in `kmersearch_hash.c`
- PGIN index structure and basic search algorithm
- Recheck integration with existing scoring functions

**Key Functions**:
```c
Datum pgin_handler(PG_FUNCTION_ARGS);
Datum pgin_build(PG_FUNCTION_ARGS);
Datum pgin_gettuple(PG_FUNCTION_ARGS);
Datum pgin_recheck(PG_FUNCTION_ARGS);
```

### Phase 2: Bloom Filter Integration
**Timeline**: 2-3 months

**Goals**:
- Add Bloom filter existence checking to PGIN
- Implement candidate pre-filtering
- Achieve 5-10x search speed improvement

**Deliverables**:
- Bloom filter implementation in `kmersearch_bloom.c`
- Integration with PGIN search algorithm
- Performance benchmarks vs traditional GIN

**Key Functions**:
```c
BloomFilter* bloom_create(size_t bits, int num_hash_funcs);
bool bloom_test(BloomFilter* bf, uint32_t key);
void bloom_add(BloomFilter* bf, uint32_t key);
```

### Phase 3: Count-min Sketch Integration  
**Timeline**: 2-3 months

**Goals**:
- Add CMS for frequency-based filtering and optimization
- Implement intelligent recheck ordering
- Achieve maximum performance with intelligent candidate management

**Deliverables**:
- Count-min sketch implementation in `kmersearch_cms.c`
- Advanced candidate filtering and ordering
- Complete PGIN architecture

**Key Functions**:
```c
CountMinSketch* cms_create(uint32_t width, uint32_t depth);
void cms_increment(CountMinSketch* cms, uint32_t key);
uint32_t cms_estimate_frequency(CountMinSketch* cms, uint32_t key);
float cms_estimate_score(CountMinSketch* cms, uint32_t* query_keys, int num_keys);
```

## Performance Expectations

### Storage Requirements Comparison
| Component | Traditional GIN | PGIN Only | Reduction |
|-----------|-----------------|-----------|-----------|
| Main Index | 500GB-2TB | 50GB-100GB | 80-90% |
| Feature Hash Metadata | N/A | 8MB | N/A |
| Bloom Filter | N/A | 12.5MB | N/A |
| Count-min Sketch | N/A | 16MB | N/A |
| **Total Storage** | **500GB-2TB** | **50GB-100GB** | **80-90%** |

### Search Performance Improvement
| Stage | Speedup | Rationale |
|-------|---------|-----------|
| Feature Hash Lookup | 1000x | O(1) hash vs O(log N) tree |
| Bloom Filter | 100-500x | In-memory bit operations |
| CMS Pre-filtering | 10-50x | Early candidate elimination |
| Recheck Only High-Quality | 5-20x | Selective exact scoring |
| **Overall** | **50-200x** | **Compound improvement** |

### Accuracy Analysis with Recheck
| Component | False Positive Rate | Final Impact |
|-----------|-------------------|---------------|
| Feature Hash | ~0.1% | Eliminated by recheck |
| Bloom Filter | ~1% | Eliminated by recheck |
| CMS Over-estimation | ~5% | Eliminated by recheck |
| **Final Accuracy** | **100%** | **Recheck ensures correctness** |

## Technical Implementation Details

### PostgreSQL Access Method Integration
```sql
-- Create new access method
CREATE ACCESS METHOD pgin TYPE INDEX HANDLER pgin_handler;

-- PGIN operator class definition  
CREATE OPERATOR CLASS pgin_kmersearch_ops
FOR TYPE dna2 USING pgin AS
    OPERATOR 1 =%,
    FUNCTION 1 pgin_options,
    FUNCTION 2 pgin_build,
    FUNCTION 3 pgin_gettuple,
    FUNCTION 4 pgin_recheck;

-- Index creation with PGIN
CREATE INDEX ON sequences USING pgin (dna_sequence pgin_kmersearch_ops)
WITH (
    feature_hash_dimensions = 1000000,
    bloom_filter_bits = 100000000,
    cms_width = 65536,
    cms_depth = 8
);
```

### PGIN Data Structure
```c
typedef struct PGINIndex {
    // Feature hashing for key space reduction
    FeatureHashTable* feature_hash;
    uint32_t hash_dimensions;           // ~10^6
    
    // Bloom filter for existence checking
    BloomFilter* existence_filter;
    size_t bloom_bits;                  // ~100M bits
    
    // Count-min sketch for frequency estimation
    CountMinSketch* frequency_cms;
    uint32_t cms_width, cms_depth;      // 64K × 8
    
    // Reduced key index (core PGIN structure)
    PGINKeyIndex* reduced_index;        // Hash-reduced k-mer index
    
    // Metadata
    uint64_t total_kmers_processed;
    double collision_rate_estimate;
} PGINIndex;
```

### PGIN Search Algorithm
```c
Datum pgin_gettuple(PG_FUNCTION_ARGS) {
    IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
    PGINScanOpaque* so = (PGINScanOpaque*) scan->opaque;
    
    if (so->candidate_list == NIL) {
        // Initialize candidate search
        DNA2* query = so->query_sequence;
        
        // Step 1: Feature hash k-mers from query
        uint32_t* hashed_kmers = feature_hash_query_kmers(so->pgin->feature_hash, query);
        
        // Step 2: Bloom filter existence check
        List* bloom_candidates = bloom_filter_candidates(so->pgin->existence_filter, hashed_kmers);
        
        // Step 3: CMS frequency estimation and filtering
        so->candidate_list = cms_filter_candidates(so->pgin->frequency_cms, bloom_candidates);
        so->candidate_iterator = list_head(so->candidate_list);
    }
    
    // Return next candidate (requires recheck)
    if (so->candidate_iterator != NULL) {
        ItemPointer tid = (ItemPointer) lfirst(so->candidate_iterator);
        so->candidate_iterator = lnext(so->candidate_iterator);
        
        scan->xs_ctup.t_self = *tid;
        scan->xs_recheck = true;  // Always require recheck
        PG_RETURN_BOOL(true);
    }
    
    PG_RETURN_BOOL(false);
}

// Recheck for exact accuracy
Datum pgin_recheck(PG_FUNCTION_ARGS) {
    IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
    HeapTuple tuple = (HeapTuple) PG_GETARG_POINTER(1);
    
    PGINScanOpaque* so = (PGINScanOpaque*) scan->opaque;
    
    // Extract actual sequence and calculate exact score
    DNA2* target_sequence = extract_dna2_from_tuple(tuple);
    float exact_score = kmersearch_rawscore_direct(so->original_query, target_sequence);
    
    // Apply exact threshold check
    PG_RETURN_BOOL(exact_score >= kmersearch_min_score);
}
```

## Configuration Parameters

### New GUC Variables
```c
// PGIN configuration parameters
static int pgin_feature_hash_dimensions = 1000000;
static int pgin_bloom_filter_bits = 100000000;
static int pgin_cms_width = 65536;
static int pgin_cms_depth = 8;
static double pgin_bloom_false_positive_rate = 0.01;
static double pgin_cms_threshold_ratio = 0.8;      // fraction of min_score for CMS filtering
static bool pgin_enable_cms_reorder = true;        // reorder candidates by estimated score
```

### Memory Management
- **Feature Hash Table**: Index build context allocation
- **Bloom Filter**: Persistent index storage with memory mapping
- **Count-min Sketch**: Index-persistent with WAL logging support
- **PGIN Index**: Standard PostgreSQL index storage
- **Recheck Integration**: Utilizes existing kmersearch scoring functions

## Testing and Validation

### Accuracy Testing
1. **Recheck Verification**: Ensure 100% accuracy through recheck mechanism
2. **Candidate Quality**: Measure precision/recall of candidate selection
3. **Score Correlation**: Validate CMS estimated vs exact scores correlation
4. **Feature Hash Collision**: Analyze impact of hash collisions on result quality

### Performance Testing  
1. **Storage Efficiency**: Measure actual index size vs traditional GIN
2. **Search Latency**: Benchmark query response times across data scales
3. **Candidate Filtering**: Evaluate effectiveness of Bloom filter and CMS
4. **Recheck Overhead**: Measure cost of exact score computation
5. **Scalability**: Performance with billion-row datasets

### Test Datasets
- **Small Scale**: 1M sequences for development and correctness testing
- **Medium Scale**: 100M sequences for performance validation  
- **Large Scale**: 1B sequences for scalability and storage efficiency testing
- **Genomic Diversity**: Multiple organism types for biological relevance

## Risk Assessment and Mitigation

### Technical Risks
| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Feature hash collision degradation | Medium | Medium | Multiple hash functions and collision monitoring |
| PostgreSQL access method complexity | High | Medium | Follow existing access method patterns (GIN/GiST) |
| Recheck performance overhead | Medium | Medium | Intelligent candidate ordering with CMS |
| Index build time increase | Low | Low | Acceptable trade-off for storage savings |

### Performance Risks
| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| Storage savings less than projected | Low | Medium | Conservative estimates, even 50% savings valuable |
| Candidate filtering efficiency low | Medium | Medium | Tunable Bloom filter and CMS parameters |
| Recheck bottleneck at scale | Medium | High | CMS-based priority ordering for high-quality candidates |

## Success Criteria

### Phase 1 (PGIN Foundation)
- [ ] New access method "pgin" successfully implemented
- [ ] Basic feature hashing and recheck mechanism working
- [ ] 50% storage reduction compared to traditional GIN
- [ ] Functional compatibility with existing `=%` operator
- [ ] Index build and search operations stable

### Phase 2 (Bloom Filter Integration)
- [ ] Bloom filter candidate filtering implemented
- [ ] 5-10x search speed improvement achieved
- [ ] False positive rate <1% for Bloom filter stage
- [ ] Memory usage within projected bounds (<15MB overhead)

### Phase 3 (Count-min Sketch Completion)
- [ ] CMS frequency estimation and candidate reordering
- [ ] 50-200x overall performance improvement vs traditional GIN
- [ ] 80-90% total storage reduction achieved
- [ ] 100% accuracy maintained through recheck mechanism
- [ ] Production-ready for billion-row datasets

## Future Enhancements

### Advanced Probabilistic Techniques
- **Locality-Sensitive Hashing (LSH)**: Better similarity preservation in feature hashing
- **Cuckoo Filters**: Space-efficient alternative to Bloom filters with deletion support
- **HyperLogLog**: Cardinality estimation for PGIN index statistics

### PGIN + Vector Search Hybrid
- **Dual Index Architecture**: PGIN for exact matches + sparse vector search for similarity
- **Feature Hash to Vector**: Convert PGIN feature hash outputs to sparse vectors
- **Multi-modal Queries**: Support both exact k-mer and approximate similarity search

### Distributed PGIN
- **Probabilistic Structure Sharding**: Distribute feature hash, Bloom filter, and CMS across nodes
- **PGIN Replication**: Replicate probabilistic structures for fault tolerance
- **Incremental PGIN Updates**: Efficient maintenance without full rebuild

## Conclusion

PGIN represents a revolutionary approach to genomic sequence indexing by implementing a completely independent PostgreSQL access method that eliminates the storage overhead of traditional GIN indexes while maintaining perfect accuracy through recheck mechanisms.

### Key Innovations

- **Independent Access Method**: PGIN operates without traditional GIN indexes, achieving massive storage savings
- **Recheck-Based Accuracy**: Leverages PostgreSQL's recheck mechanism for 100% accuracy guarantee
- **Probabilistic Filtering Pipeline**: Feature hashing → Bloom filter → Count-min sketch for efficient candidate selection
- **Storage Efficiency**: 80-90% reduction in index storage requirements

### Expected Impact

- **80-90% storage reduction** compared to traditional GIN indexes (50-100GB vs 500GB-2TB)
- **50-200x performance improvement** through multi-stage probabilistic filtering
- **100% accuracy** maintained via PostgreSQL recheck mechanism
- **Billion-row scalability** with sub-minute query response times
- **Zero traditional GIN dependency** - completely self-contained solution

PGIN will enable pg_kmersearch to handle truly massive genomic datasets efficiently while remaining within practical storage and memory constraints, making it the ideal solution for next-generation genomic sequence search at planetary scale.