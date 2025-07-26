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
- **Score calculation functions**: `kmersearch_rawscore()` and `kmersearch_correctedscore()` for individual sequence scoring (both functions return identical values in current implementation)
- **High-frequency k-mer management**: `kmersearch_perform_highfreq_analysis()` for high-frequency k-mer analysis and `kmersearch_highfreq_kmer_cache_load()` and `kmersearch_highfreq_kmer_cache_free()` for cache management

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

3. Configure shared_preload_libraries (required for GUC variables):
```bash
# Add to postgresql.conf
shared_preload_libraries = 'pg_kmersearch'
```

4. Restart PostgreSQL and install the extension in your database:
```sql
CREATE EXTENSION pg_kmersearch;
```

**Important Note on GUC Variables:**
PostgreSQL sessions reset GUC variables to default values when a new session starts. Always set required GUC variables at the beginning of each session:

```sql
-- Example: Setting GUC variables
SET kmersearch.kmer_size = 4;
SET kmersearch.max_appearance_rate = 0.25;
-- Your analysis commands here
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
-- Configure k-mer size (default 16-mer)
SET kmersearch.kmer_size = 16;

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
| `kmersearch.kmer_size` | 16 | 4-64 | K-mer length for index creation and search |
| `kmersearch.occur_bitlen` | 8 | 0-16 | Bits for occurrence count storage |
| `kmersearch.max_appearance_rate` | 0.5 | 0.0-1.0 | Maximum k-mer appearance rate for indexing |
| `kmersearch.max_appearance_nrow` | 0 | 0-∞ | Maximum rows containing k-mer (0=unlimited) |
| `kmersearch.min_score` | 1 | 0-∞ | Minimum similarity score for search results |
| `kmersearch.min_shared_ngram_key_rate` | 0.9 | 0.0-1.0 | Minimum threshold for shared n-gram key rate |
| `kmersearch.rawscore_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for rawscore cache |
| `kmersearch.query_pattern_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for query pattern cache |
| `kmersearch.actual_min_score_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for actual min score cache |
| `kmersearch.preclude_highfreq_kmer` | false | true/false | Enable high-frequency k-mer exclusion during GIN index construction |
| `kmersearch.force_use_parallel_highfreq_kmer_cache` | false | true/false | Force use of dshash parallel cache for high-frequency k-mer lookups |
| `kmersearch.highfreq_kmer_cache_load_batch_size` | 10000 | 1000-1000000 | Batch size for loading high-frequency k-mers into cache |
| `kmersearch.highfreq_analysis_batch_size` | 10000 | 1000-1000000 | Batch size for high-frequency k-mer analysis |

### High-Frequency K-mer Exclusion

The extension includes automatic exclusion of high-frequency k-mers to improve index performance:

```sql
-- Configure exclusion parameters (before index creation)
SET kmersearch.max_appearance_rate = 0.5;  -- Default: 50% max appearance rate
SET kmersearch.max_appearance_nrow = 1000;  -- Default: 0 (disabled)

-- Create index with frequency analysis
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq);

-- View excluded k-mers for a table/column
SELECT ngram_key, detection_reason 
FROM kmersearch_highfreq_kmer 
WHERE table_oid = 'sequences'::regclass AND column_name = 'dna_seq';

-- View index statistics
SELECT table_oid, column_name, kmer_size, occur_bitlen, max_appearance_rate, max_appearance_nrow 
FROM kmersearch_highfreq_kmer_meta 
WHERE table_oid = 'sequences'::regclass;
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

### BYTEA Conversion Functions

- **`kmersearch_dna2_to_bytea(DNA2)`**: Converts DNA2 to BYTEA format suitable for hashing (includes bit length prefix)
- **`kmersearch_dna4_to_bytea(DNA4)`**: Converts DNA4 to BYTEA format suitable for hashing (includes bit length prefix)

These functions are optimized for use with cryptographic hash functions like pgcrypto's `digest()`. They include a 4-byte bit length prefix in network byte order, followed by the actual bit data, ensuring unique hash values even for sequences that differ only in padding.

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

-- BYTEA conversion for hashing (requires pgcrypto extension)
SELECT digest(kmersearch_dna2_to_bytea('ATCG'::DNA2), 'sha256');
SELECT digest(kmersearch_dna4_to_bytea('ATCGN'::DNA4), 'sha256');

-- Verifying padding collision avoidance
SELECT 
    digest(kmersearch_dna2_to_bytea('ATG'::DNA2), 'sha256') AS hash_atg,
    digest(kmersearch_dna2_to_bytea('ATGA'::DNA2), 'sha256') AS hash_atga;
-- These will produce different hash values despite similar padding

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

## System Tables and Views

pg_kmersearch creates several system tables and views for metadata management and monitoring.

### System Tables

#### kmersearch_highfreq_kmer
Stores high-frequency k-mers that are excluded from GIN indexes:

```sql
-- View excluded k-mers for a specific table/column
SELECT table_oid, column_name, ngram_key, detection_reason, created_at
FROM kmersearch_highfreq_kmer 
WHERE table_oid = 'sequences'::regclass AND column_name = 'dna_seq';

-- Count excluded k-mers per table/column
SELECT table_oid, column_name, COUNT(*) as excluded_count
FROM kmersearch_highfreq_kmer
GROUP BY table_oid, column_name;
```

#### kmersearch_highfreq_kmer_meta
Stores metadata about k-mer frequency analysis:

```sql
-- View analysis metadata for all tables
SELECT table_oid, column_name, kmer_size, occur_bitlen, 
       max_appearance_rate, max_appearance_nrow, analysis_timestamp
FROM kmersearch_highfreq_kmer_meta;

-- Check analysis parameters for specific table
SELECT * FROM kmersearch_highfreq_kmer_meta 
WHERE table_oid = 'sequences'::regclass AND column_name = 'dna_seq';
```

#### kmersearch_gin_index_meta
Stores GIN index metadata including high-frequency filtering information:

```sql
-- View GIN index metadata
SELECT index_oid, table_oid, column_name, highfreq_filtered, 
       kmer_size, occur_bitlen, created_at
FROM kmersearch_gin_index_meta;

-- Check if an index uses high-frequency filtering
SELECT highfreq_filtered FROM kmersearch_gin_index_meta 
WHERE index_oid = 'sequences_kmer_idx'::regclass;
```

#### kmersearch_index_info
Stores comprehensive index statistics and configuration:

```sql
-- View index statistics
SELECT index_oid, table_oid, column_name, kmer_size, total_nrow,
       highfreq_kmer_count, max_appearance_rate, created_at
FROM kmersearch_index_info;
```

### Management Views

#### kmersearch_cache_summary
Provides unified cache statistics across all cache types:

```sql
-- View cache performance summary
SELECT cache_type, total_entries, total_hits, total_misses, hit_rate
FROM kmersearch_cache_summary;

-- Monitor cache efficiency
SELECT cache_type, 
       ROUND(hit_rate * 100, 2) as hit_rate_percent,
       total_hits + total_misses as total_requests
FROM kmersearch_cache_summary
WHERE total_hits + total_misses > 0;
```

#### kmersearch_analysis_status
Shows high-frequency k-mer analysis status for all analyzed tables:

```sql
-- View analysis status for all tables
SELECT table_name, column_name, kmer_size, highfreq_kmer_count,
       analysis_timestamp
FROM kmersearch_analysis_status;

-- Check analysis parameters and results
SELECT table_name, column_name, kmer_size, occur_bitlen,
       max_appearance_rate, max_appearance_nrow,
       highfreq_kmer_count, analysis_timestamp
FROM kmersearch_analysis_status
WHERE table_name = 'sequences';
```

## Analysis and Management Functions

### K-mer Frequency Analysis

#### kmersearch_perform_highfreq_analysis()
Performs parallel k-mer frequency analysis on a table:

```sql
-- Basic analysis using table name and column name
SELECT kmersearch_perform_highfreq_analysis(
    'sequences',                   -- table name or OID
    'dna_seq'                     -- column name or attnum
);

-- Analysis using OID and attnum
SELECT kmersearch_perform_highfreq_analysis(
    '16384'::text,                -- table OID as text
    '3'::text                     -- column attnum as text
);

-- Example result interpretation
SELECT (result).total_rows,
       (result).highfreq_kmers_count,
       (result).parallel_workers_used,
       (result).max_appearance_rate_used
FROM (
    SELECT kmersearch_perform_highfreq_analysis('sequences', 'dna_seq') as result
) t;
```

This function uses PostgreSQL's standard parallel execution framework to distribute k-mer extraction and counting across multiple CPU cores.

#### kmersearch_undo_highfreq_analysis()
Removes analysis data and frees storage:

```sql
-- Undo analysis for specific table/column combination
SELECT kmersearch_undo_highfreq_analysis(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- Example result interpretation
SELECT (result).dropped_analyses,
       (result).dropped_highfreq_kmers,
       (result).freed_storage_bytes
FROM (
    SELECT kmersearch_undo_highfreq_analysis('sequences', 'dna_seq') as result
) t;
```

### High-Frequency K-mer Cache Management

#### Global Cache Functions

```sql
-- Load high-frequency k-mers into global cache
SELECT kmersearch_highfreq_kmer_cache_load(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- Free global cache for specific table/column
SELECT kmersearch_highfreq_kmer_cache_free(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- Free all entries from global cache
SELECT kmersearch_highfreq_kmer_cache_free_all();
```

#### Parallel Cache Functions

```sql
-- Load high-frequency k-mers into parallel cache (for multi-process sharing)
SELECT kmersearch_parallel_highfreq_kmer_cache_load(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- Free parallel cache for specific table/column
SELECT kmersearch_parallel_highfreq_kmer_cache_free(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- Free all entries from parallel cache
SELECT kmersearch_parallel_highfreq_kmer_cache_free_all();
```

### Cache Statistics and Management

#### Rawscore Cache

```sql
-- View rawscore cache statistics
SELECT * FROM kmersearch_rawscore_cache_stats();

-- Monitor cache performance
SELECT dna2_hits, dna2_misses, 
       CASE WHEN (dna2_hits + dna2_misses) > 0 
            THEN dna2_hits::float / (dna2_hits + dna2_misses)::float * 100
            ELSE 0 END as dna2_hit_rate_percent,
       dna4_hits, dna4_misses,
       CASE WHEN (dna4_hits + dna4_misses) > 0 
            THEN dna4_hits::float / (dna4_hits + dna4_misses)::float * 100
            ELSE 0 END as dna4_hit_rate_percent
FROM kmersearch_rawscore_cache_stats();

-- Clear rawscore cache
SELECT kmersearch_rawscore_cache_free();
```

#### Query Pattern Cache

```sql
-- View query pattern cache statistics
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- Monitor cache efficiency
SELECT hits, misses, current_entries, max_entries,
       CASE WHEN (hits + misses) > 0 
            THEN hits::float / (hits + misses)::float * 100
            ELSE 0 END as hit_rate_percent
FROM kmersearch_query_pattern_cache_stats();

-- Clear query pattern cache
SELECT kmersearch_query_pattern_cache_free();
```

#### Actual Min Score Cache

```sql
-- View actual min score cache statistics
SELECT * FROM kmersearch_actual_min_score_cache_stats();

-- Clear actual min score cache  
SELECT kmersearch_actual_min_score_cache_free();
```


## Complete Workflow Examples

### Setting Up K-mer Analysis and Indexing

```sql
-- 1. Configure parameters
SET kmersearch.kmer_size = 16;
SET kmersearch.max_appearance_rate = 0.5;
SET kmersearch.max_appearance_nrow = 1000;
SET kmersearch.occur_bitlen = 8;

-- 2. Perform frequency analysis
SELECT kmersearch_perform_highfreq_analysis(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- 3. Check analysis results
SELECT * FROM kmersearch_analysis_status WHERE table_name = 'sequences';

-- 4. Create GIN index (uses analysis results automatically)
CREATE INDEX sequences_kmer_idx ON sequences USING gin(dna_seq);

-- 5. Load caches for optimal performance
SELECT kmersearch_highfreq_kmer_cache_load('sequences', 'dna_seq');
SELECT kmersearch_parallel_highfreq_kmer_cache_load('sequences', 'dna_seq');

-- 6. Perform searches
SELECT id, name, kmersearch_rawscore(dna_seq, 'ATCGATCG') as score
FROM sequences 
WHERE dna_seq =% 'ATCGATCG'
ORDER BY score DESC;
```

### Cache Performance Monitoring

```sql
-- Monitor all cache types
SELECT cache_type, hit_rate, total_entries, 
       total_hits + total_misses as total_requests
FROM kmersearch_cache_summary;

-- Detailed rawscore cache analysis
WITH cache_stats AS (
    SELECT dna2_hits, dna2_misses, dna2_entries,
           dna4_hits, dna4_misses, dna4_entries
    FROM kmersearch_rawscore_cache_stats()
)
SELECT 
    'DNA2' as type,
    dna2_entries as entries,
    dna2_hits as hits,
    dna2_misses as misses,
    CASE WHEN (dna2_hits + dna2_misses) > 0 
         THEN dna2_hits::float / (dna2_hits + dna2_misses)::float 
         ELSE 0 END as hit_rate
FROM cache_stats
UNION ALL
SELECT 
    'DNA4' as type,
    dna4_entries as entries, 
    dna4_hits as hits,
    dna4_misses as misses,
    CASE WHEN (dna4_hits + dna4_misses) > 0 
         THEN dna4_hits::float / (dna4_hits + dna4_misses)::float 
         ELSE 0 END as hit_rate
FROM cache_stats;
```

## Complex Type Definitions

pg_kmersearch defines custom composite types for function return values:

### kmersearch_analysis_result
Returned by `kmersearch_analyze_table()`:

```sql
-- Type definition equivalent:
-- CREATE TYPE kmersearch_analysis_result AS (
--     total_rows bigint,
--     highfreq_kmers_count integer,
--     parallel_workers_used integer,
--     max_appearance_rate_used real,
--     max_appearance_nrow_used integer
-- );

-- Usage example:
SELECT (result).total_rows,
       (result).highfreq_kmers_count,
       (result).max_appearance_rate_used
FROM (
    SELECT kmersearch_analyze_table('sequences', 'dna_seq') as result
) t;
```

### kmersearch_drop_result
Returned by `kmersearch_undo_highfreq_analysis()`:

```sql
-- Type definition equivalent:
-- CREATE TYPE kmersearch_drop_result AS (
--     dropped_analyses integer,
--     dropped_highfreq_kmers integer,
--     freed_storage_bytes bigint
-- );

-- Usage example:
SELECT (result).dropped_analyses,
       (result).dropped_highfreq_kmers,
       pg_size_pretty((result).freed_storage_bytes) as freed_storage
FROM (
    SELECT kmersearch_undo_highfreq_analysis('sequences', 'dna_seq') as result
) t;
```

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
- **High-frequency exclusion**: Parallel table scan using multiple workers
- **Parallel k-mer analysis**: True parallel processing with PostgreSQL's ParallelContext
- **System tables**: Metadata storage for excluded k-mers and index statistics (`kmersearch_highfreq_kmer`, `kmersearch_highfreq_kmer_meta`)
- **Cache system**: TopMemoryContext-based high-performance caching
- Binary input/output support

### K-mer Search Mechanism
1. **Frequency analysis**: Parallel table scan using multiple workers to identify high-frequency k-mers
2. **K-mer extraction**: Sliding window with specified k-length (parallel processing)
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
- **Memory management**: TopMemoryContext-based implementation

### Query Pattern Cache
- **Purpose**: Performance improvement through query pattern reuse
- **Memory management**: TopMemoryContext-based implementation

### Cache Statistics and Management Functions

```sql
-- Check cache statistics
SELECT * FROM kmersearch_actual_min_score_cache_stats();
SELECT * FROM kmersearch_rawscore_cache_stats();
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- Clear caches
SELECT kmersearch_actual_min_score_cache_free();
SELECT kmersearch_query_pattern_cache_free();
SELECT kmersearch_highfreq_kmer_cache_free_all();

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

3. **Table Reference Fallback** (`kmersearch_highfreq_kmer`)
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
SELECT kmersearch_perform_highfreq_analysis(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- Set GUC variables to match metadata
SET kmersearch.max_appearance_rate = 0.5;
SET kmersearch.max_appearance_nrow = 100;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.kmer_size = 16;

-- Load global cache
SELECT kmersearch_highfreq_kmer_cache_load(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- Load parallel cache (optional)
SELECT kmersearch_parallel_highfreq_kmer_cache_load(
    'sequences',                   -- table name
    'dna_seq'                     -- column name
);

-- High-speed search using cache hierarchy
SELECT id, name FROM sequences 
WHERE dna_seq =% 'ATCGATCG'
ORDER BY kmersearch_rawscore(dna_seq, 'ATCGATCG') DESC;

-- Free caches
SELECT kmersearch_highfreq_kmer_cache_free('sequences', 'dna_seq');
SELECT kmersearch_parallel_highfreq_kmer_cache_free('sequences', 'dna_seq');
```

### Parallel Cache Functions

- **`kmersearch_parallel_highfreq_kmer_cache_load(table_name, column_name)`**: Load high-frequency k-mers into shared dshash cache
- **`kmersearch_parallel_highfreq_kmer_cache_free(table_name, column_name)`**: Free specific entries from the parallel cache
- **`kmersearch_parallel_highfreq_kmer_cache_free_all()`**: Free all entries from the parallel cache and destroy shared memory structures

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