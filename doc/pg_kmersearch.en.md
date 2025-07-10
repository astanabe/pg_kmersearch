# pg_kmersearch

PostgreSQL extension for DNA sequence data types and k-mer search

## Overview

pg_kmersearch is a PostgreSQL extension that provides custom data types for efficiently storing and processing DNA sequence data with k-mer search capabilities. This extension implements two data types and high-performance search functionality:

### DNA2 Type
- **Purpose**: Store standard DNA sequences (A, C, G, T)
- **Encoding**: 2-bit compression (A=00, C=01, G=10, T=11)
- **Supported characters**: A, C, G, T, U (U is treated as T)
- **Storage efficiency**: 2 bits per character

### DNA4 Type
- **Purpose**: Store DNA sequences with degenerate codes
- **Encoding**: 4-bit compression
- **Supported characters**:
  - Standard bases: A(0001), C(0010), G(0100), T(1000)
  - Degenerate codes: M(0011), R(0101), W(1001), S(0110), Y(1010), K(1100), V(0111), H(1011), D(1101), B(1110), N(1111)
  - U (treated as T)
- **Storage efficiency**: 4 bits per character

### K-mer Search Features
- **K-mer length**: 4-64 bases (specified at index creation)
- **GIN indexing**: Fast search using n-gram keys
- **Degenerate code support**: MRWSYKVHDBN expansion for DNA4 type
- **Occurrence tracking**: Counts k-mer occurrences within rows (default 8-bit)
- **Scoring search**: Retrieve top matches by similarity, not just exact matches
- **High-frequency k-mer exclusion**: Automatically excludes overly common k-mers during index creation
- **Score-based filtering**: Minimum score thresholds with automatic adjustment for excluded k-mers
- **Score calculation functions**: `kmersearch_rawscore()` and `kmersearch_correctedscore()` for individual sequence scoring
- **High-frequency k-mer management**: `kmersearch_analyze_table()` for high-frequency k-mer analysis and `kmersearch_highfreq_kmers_cache_load()` and `kmersearch_highfreq_kmers_cache_free()` for cache management

## Installation

### Prerequisites
- PostgreSQL 16 or higher
- postgresql-server-dev-16
- make
- gcc or clang

### Installation Steps

1. Download the source code:
```bash
git clone <repository-url>
cd pg_kmersearch
```

2. Compile and install:
```bash
make
sudo make install
```

3. Install the extension in your database:
```sql
CREATE EXTENSION pg_kmersearch;
```

## Usage

### Basic Usage Examples

```sql
-- DNA2 type usage
CREATE TABLE sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    dna_seq DNA2
);

INSERT INTO sequences (name, dna_seq) VALUES 
    ('seq1', 'ATCGATCG'::DNA2),
    ('seq2', 'GCTAGCTA'::DNA2);

SELECT name, dna_seq FROM sequences;

-- DNA4 type usage (with degenerate codes)
CREATE TABLE degenerate_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    dna_seq DNA4
);

INSERT INTO degenerate_sequences (name, dna_seq) VALUES 
    ('seq1', 'ATCGMRWSYKN'::DNA4),
    ('seq2', 'VHDBNATCG'::DNA4);

SELECT name, dna_seq FROM degenerate_sequences;
```

### K-mer Search Usage Examples

```sql
-- Configure k-mer size (default 8-mer)
SET kmersearch.kmer_size = 8;

-- Create GIN index (uses current kmersearch.kmer_size setting)
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq);

-- K-mer search using =% operator
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore,
       kmersearch_correctedscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS correctedscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') DESC 
LIMIT 10;

-- Search with degenerate codes in query
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG') AS rawscore
FROM degenerate_sequences 
WHERE dna_seq =% 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'
ORDER BY rawscore DESC 
LIMIT 5;

-- Configure occurrence bit length (default 8-bit)
SET kmersearch.occur_bitlen = 12; -- Change to 12-bit (max 4095 occurrences)

-- Check current configuration
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
```

### Configuration Variables

pg_kmersearch provides several configuration variables that can be set using PostgreSQL's `SET` command:

| Variable | Default | Range | Description |
|----------|---------|-------|-------------|
| `kmersearch.kmer_size` | 8 | 4-64 | K-mer length for index creation and search |
| `kmersearch.occur_bitlen` | 8 | 0-16 | Bits for occurrence count storage |
| `kmersearch.max_appearance_rate` | 0.05 | 0.0-1.0 | Maximum k-mer appearance rate for indexing |
| `kmersearch.max_appearance_nrow` | 0 | 0-∞ | Maximum rows containing k-mer (0=unlimited) |
| `kmersearch.min_score` | 1 | 0-∞ | Minimum similarity score for search results |
| `kmersearch.min_shared_ngram_key_rate` | 0.9 | 0.0-1.0 | Minimum threshold for shared n-gram key rate |
| `kmersearch.rawscore_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for rawscore cache |
| `kmersearch.query_pattern_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for query pattern cache |
| `kmersearch.actual_min_score_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for actual min score cache |
| `kmersearch.preclude_highfreq_kmer` | false | true/false | Enable high-frequency k-mer exclusion during GIN index construction |
| `kmersearch.force_use_dshash` | false | true/false | Force use of dshash parallel cache for high-frequency k-mer lookups |

### High-Frequency K-mer Exclusion

The extension includes automatic exclusion of high-frequency k-mers to improve index performance:

```sql
-- Configure exclusion parameters (before index creation)
SET kmersearch.max_appearance_rate = 0.05;  -- Default: 5% max appearance rate
SET kmersearch.max_appearance_nrow = 1000;  -- Default: 0 (disabled)

-- Create index with frequency analysis
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq);

-- View excluded k-mers for an index
SELECT kmer_key, frequency_count, exclusion_reason 
FROM kmersearch_excluded_kmers 
WHERE index_oid = 'sequences_kmer_idx'::regclass;

-- View index statistics
SELECT total_rows, excluded_kmers_count, max_appearance_rate 
FROM kmersearch_index_info 
WHERE index_oid = 'sequences_kmer_idx'::regclass;
```

### Score-based Search Filtering

Control search quality with minimum score thresholds, automatically adjusted for excluded k-mers:

```sql
-- Set minimum score for GIN search results
SET kmersearch.min_score = 50;  -- Default: 1

-- Check current minimum score setting
SHOW kmersearch.min_score;

-- Search with automatic score adjustment
-- If query contains 3 excluded k-mers and min_score=50, 
-- actual threshold becomes 47 for that query
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY rawscore DESC LIMIT 10;
```

### Score Calculation Functions

Retrieve match scores for individual sequences:

```sql
-- Get raw scores for matched sequences
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY rawscore DESC 
LIMIT 10;

-- Get corrected scores (accounting for excluded k-mers)
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS raw_score,
       kmersearch_correctedscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS corrected_score
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY corrected_score DESC;

-- Example with both DNA2 and DNA4 types
SELECT 'DNA2' as type, id, kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score
FROM dna2_sequences WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
UNION ALL
SELECT 'DNA4' as type, id, kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score  
FROM dna4_sequences WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY score DESC;
```

## Length Functions

pg_kmersearch provides several length functions for DNA2 and DNA4 types that correctly handle padding and return accurate measurements:

### Available Length Functions

- **`bit_length(DNA2/DNA4)`**: Returns the actual bit length (excluding padding)
- **`nuc_length(DNA2/DNA4)`**: Returns the number of nucleotides
- **`char_length(DNA2/DNA4)`**: Same as `nuc_length()` (character count)
- **`length(DNA2/DNA4)`**: Same as `nuc_length()` (standard length function)

### Function Relationships

- **DNA2**: `nuc_length() = bit_length() / 2` (2 bits per nucleotide)
- **DNA4**: `nuc_length() = bit_length() / 4` (4 bits per nucleotide)
- **Consistency**: `char_length() = length() = nuc_length()`

### Usage Examples

```sql
-- Basic length measurements
SELECT 
    bit_length('ATCGA'::DNA2) AS bits,      -- Returns: 10
    nuc_length('ATCGA'::DNA2) AS nucs,      -- Returns: 5
    char_length('ATCGA'::DNA2) AS chars,    -- Returns: 5
    length('ATCGA'::DNA2) AS len;           -- Returns: 5

-- DNA4 with degenerate codes
SELECT 
    bit_length('ATCGMRWSYKN'::DNA4) AS bits,  -- Returns: 44
    nuc_length('ATCGMRWSYKN'::DNA4) AS nucs;  -- Returns: 11

-- Padding verification (non-multiples of 4/8)
SELECT 
    bit_length('ATCGATCGA'::DNA2) AS bits,    -- 9 nucleotides * 2 = 18 bits
    nuc_length('ATCGATCGA'::DNA2) AS nucs;    -- Returns: 9

-- Using length functions in queries
SELECT id, name, length(dna_seq) AS sequence_length
FROM sequences
WHERE length(dna_seq) >= 50
ORDER BY length(dna_seq) DESC;
```

### Degenerate Code Meanings

| Code | Meaning | Bit Representation |
|------|---------|-------------------|
| A | Adenine | 0001 |
| C | Cytosine | 0010 |
| G | Guanine | 0100 |
| T | Thymine | 1000 |
| M | A or C | 0011 |
| R | A or G | 0101 |
| W | A or T | 1001 |
| S | C or G | 0110 |
| Y | C or T | 1010 |
| K | G or T | 1100 |
| V | A or C or G | 0111 |
| H | A or C or T | 1011 |
| D | A or G or T | 1101 |
| B | C or G or T | 1110 |
| N | A or C or G or T | 1111 |

## Technical Details

### Architecture
- Built on PostgreSQL's varbit (bit varying) type
- Custom encoding/decoding functions
- Memory-efficient storage format
- GIN indexing support for k-mer search

### Internal Implementation
- **DNA2 type**: 2 bits per character encoding
- **DNA4 type**: 4 bits per character encoding
- **N-gram keys**: k-mer (2k bits) + occurrence count (8-16 bits)
- **Degenerate expansion**: Automatic expansion up to 10 combinations
- **Parallel index creation**: Supports max_parallel_maintenance_workers
- **High-frequency exclusion**: Full table scan before index creation
- **System tables**: Metadata storage for excluded k-mers and index statistics
- **Cache system**: TopMemoryContext-based high-performance caching
- Binary input/output support

### K-mer Search Mechanism
1. **Frequency analysis**: Full table scan to identify high-frequency k-mers
2. **K-mer extraction**: Sliding window with specified k-length
3. **High-frequency filtering**: Exclude k-mers exceeding appearance thresholds
4. **N-gram key generation**: Binary encoding of k-mer + occurrence count
5. **Degenerate processing**: Expansion of MRWSYKVHDBN to standard bases
6. **Scoring**: Similarity calculation based on shared n-gram key count

## Cache Management Features

pg_kmersearch provides three types of high-performance cache systems to improve search performance:

### Actual Min Score Cache
- **Purpose**: Optimization of search condition evaluation for the `=%` operator
- **Mechanism**: Pre-calculates and caches `actual_min_score = max(kmersearch_min_score, ceil(kmersearch_min_shared_ngram_key_rate × query_total_kmers))`
- **Use cases**: 
  - Matching condition evaluation in `=%` operator
  - Value judgment for rawscore cache storage
- **Memory management**: TopMemoryContext-based implementation

### Rawscore Cache
- **Purpose**: Fast retrieval of calculated rawscores
- **Mechanism**: Caches sequence and query combination results
- **Memory management**: PortalContext-based implementation

### Query Pattern Cache
- **Purpose**: Performance improvement through query pattern reuse
- **Memory management**: PortalContext-based implementation

### Cache Statistics and Management Functions

```sql
-- Check cache statistics
SELECT * FROM kmersearch_actual_min_score_cache_stats();
SELECT * FROM kmersearch_rawscore_cache_stats();
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- Clear caches
SELECT kmersearch_actual_min_score_cache_free();
SELECT kmersearch_query_pattern_cache_free();

-- Configure cache sizes
SET kmersearch.actual_min_score_cache_max_entries = 25000;
SET kmersearch.rawscore_cache_max_entries = 25000;
SET kmersearch.query_pattern_cache_max_entries = 25000;
```

## High-Frequency K-mer Hierarchical Cache System

pg_kmersearch provides a multi-layered cache architecture optimized for PostgreSQL 16/18, maximizing high-frequency k-mer search performance:

### Cache Hierarchy

High-frequency k-mer searches utilize caches in the following priority order:

1. **Global Cache** (`global_highfreq_cache`)
   - Fastest access
   - Managed in TopMemoryContext
   - Single-process reuse

2. **Parallel Cache** (`parallel_highfreq_cache`)
   - PostgreSQL dshash (dynamic shared hash tables) implementation
   - DSM (Dynamic Shared Memory) sharing across multiple processes
   - Future PostgreSQL 18 parallel GIN index support

3. **Table Reference Fallback** (`kmersearch_highfreq_kmers`)
   - Direct system table access
   - Final fallback when caches are unavailable

### GUC Validation Feature

The following GUC variables are automatically validated during cache loading:

- `kmersearch.max_appearance_rate`
- `kmersearch.max_appearance_nrow`
- `kmersearch.occur_bitlen`
- `kmersearch.kmer_size`

Detailed error messages and hints are provided when settings mismatch.

### Parallel Cache Technical Implementation

- **Memory Management**: DSM segment pinning prevents automatic cleanup
- **Process Identification**: Differentiated resource management for main process and workers
- **Error Handling**: Comprehensive PG_TRY/PG_CATCH for safe operations
- **Lock Management**: Proper concurrent control with dshash_release_lock()

### Cache Usage Examples

```sql
-- Execute high-frequency k-mer analysis (create metadata)
SELECT kmersearch_analyze_table(
    (SELECT oid FROM pg_class WHERE relname = 'sequences'), 
    'dna_seq', 
    8,    -- k-mer length
    100   -- max appearance row threshold
);

-- Set GUC variables to match metadata
SET kmersearch.max_appearance_rate = 0.05;
SET kmersearch.max_appearance_nrow = 100;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.kmer_size = 8;

-- Load global cache
SELECT kmersearch_highfreq_kmers_cache_load(
    (SELECT oid FROM pg_class WHERE relname = 'sequences'),
    'dna_seq', 8
);

-- Load parallel cache (optional)
SELECT kmersearch_parallel_highfreq_kmers_cache_load(
    (SELECT oid FROM pg_class WHERE relname = 'sequences'),
    'dna_seq', 8
);

-- High-speed search using cache hierarchy
SELECT id, name FROM sequences 
WHERE dna_seq =% 'ATCGATCG'
ORDER BY kmersearch_rawscore(dna_seq, 'ATCGATCG') DESC;

-- Free caches
SELECT kmersearch_highfreq_kmers_cache_free();
SELECT kmersearch_parallel_highfreq_kmers_cache_free();
```

### Parallel Cache Functions

- **`kmersearch_parallel_highfreq_kmers_cache_load(table_oid, column_name, k_value)`**: Load high-frequency k-mers into shared dshash cache
- **`kmersearch_parallel_highfreq_kmers_cache_free()`**: Free all entries from the parallel cache and destroy shared memory structures

### Usage Scenarios

The parallel cache system is particularly useful for:

1. **Large-scale genomic analysis** with multiple concurrent queries
2. **Batch processing** of DNA sequence searches
3. **Future PostgreSQL 18 parallel GIN index operations**
4. **High-throughput bioinformatics applications**

### Technical Implementation

- **Memory Context**: TopMemoryContext for proper cleanup
- **Shared Memory**: DSA (Dynamic Shared Areas) with automatic pinning
- **Error Handling**: Comprehensive PG_TRY/PG_CATCH blocks
- **Lock Management**: Proper dshash_release_lock() calls for all operations
- **Process Detection**: IsParallelWorker() for proper worker identification

## Limitations

- Query sequences must be at least 8 bases long
- Degenerate code expansion limited to 10 combinations (skipped if exceeded)
- Occurrence counts capped at maximum value for configured bit length
- Case-insensitive input, uppercase output
- High-frequency k-mer exclusion requires full table scan during index creation
- Tables with GIN k-mer indexes have restricted INSERT/UPDATE/DELETE operations

## License

This project is released under the MIT License.

## Contributing

Please report bugs and feature requests on the GitHub Issues page.