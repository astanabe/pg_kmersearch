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
- **K-mer length**: 4-32 bases (configurable via GUC variable)
- **GIN indexing**: Multiple operator classes with different key storage strategies (int2/int4/int8)
- **Degenerate code support**: Full IUPAC code expansion for DNA4 type
- **Occurrence tracking**: Configurable bit length (0-16 bits, default 8-bit)
- **Scoring search**: Retrieve matches by similarity score
- **High-frequency k-mer filtering**: Optional exclusion of common k-mers
- **Score-based filtering**: Minimum score thresholds for search quality control
- **Score calculation functions**: `kmersearch_matchscore()` for sequence scoring
- **High-frequency k-mer management**: Analysis and cache management functions
- **Table partitioning**: Hash partitioning support for large databases

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
       kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS matchscore,
       kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS matchscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') DESC 
LIMIT 10;

-- Search with degenerate codes in query
SELECT id, name, dna_seq,
       kmersearch_matchscore(dna_seq, 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG') AS matchscore
FROM degenerate_sequences 
WHERE dna_seq =% 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'
ORDER BY matchscore DESC 
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
| `kmersearch.kmer_size` | 16 | 4-32 | K-mer length for index creation and search |
| `kmersearch.occur_bitlen` | 8 | 0-16 | Bits for occurrence count storage |
| `kmersearch.max_appearance_rate` | 0.5 | 0.0-1.0 | Maximum k-mer appearance rate for indexing |
| `kmersearch.max_appearance_nrow` | 0 | 0-∞ | Maximum rows containing k-mer (0=unlimited) |
| `kmersearch.min_score` | 1 | 0-∞ | Minimum similarity score for search results |
| `kmersearch.min_shared_kmer_rate` | 0.5 | 0.0-1.0 | Minimum threshold for shared k-mer rate |
| `kmersearch.query_kmer_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for query-kmer cache |
| `kmersearch.actual_min_score_cache_max_entries` | 50000 | 1000-10000000 | Maximum entries for actual min score cache |
| `kmersearch.preclude_highfreq_kmer` | false | true/false | Enable high-frequency k-mer exclusion during GIN index construction |
| `kmersearch.force_use_parallel_highfreq_kmer_cache` | false | true/false | Force use of dshash parallel cache for high-frequency k-mer lookups |
| `kmersearch.force_simd_capability` | -1 | -1-100 | Force SIMD capability level (-1 = auto-detect) |
| `kmersearch.highfreq_kmer_cache_load_batch_size` | 10000 | 1000-1000000 | Batch size for loading high-frequency k-mers into cache |
| `kmersearch.highfreq_analysis_hashtable_size` | 1000000 | 10000-100000000 | Initial hash table size for high-frequency k-mer analysis |

**Note:** High-frequency k-mer analysis batch size is automatically calculated from `maintenance_work_mem`, and ring buffer size is calculated from `shared_buffers`. No manual configuration is required for optimal I/O performance.

### High-Frequency K-mer Exclusion

The extension includes automatic exclusion of high-frequency k-mers to improve index performance:

```sql
-- Configure exclusion parameters (before index creation)
SET kmersearch.max_appearance_rate = 0.5;  -- Default: 50% max appearance rate
SET kmersearch.max_appearance_nrow = 1000;  -- Default: 0 (disabled)

-- Create index with appropriate operator class
-- For DNA2 with k-mer size <= 16:
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq kmersearch_dna2_gin_ops_int4);

-- For DNA4 with k-mer size <= 8:
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq kmersearch_dna4_gin_ops_int2);

-- For DNA4 with k-mer size <= 16:
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq kmersearch_dna4_gin_ops_int4);

-- View excluded k-mers for a table/column
SELECT uintkey, detection_reason 
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
       kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS matchscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY matchscore DESC LIMIT 10;
```

### Score Calculation Functions

Retrieve match scores for individual sequences:

```sql
-- Get raw scores for matched sequences
SELECT id, name, dna_seq,
       kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS matchscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY matchscore DESC 
LIMIT 10;

-- Get corrected scores (accounting for excluded k-mers)
SELECT id, name, dna_seq,
       kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS matchscore,
       kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS match_score
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY corrected_score DESC;

-- Example with both DNA2 and DNA4 types
SELECT 'DNA2' as type, id, kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score
FROM dna2_sequences WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
UNION ALL
SELECT 'DNA4' as type, id, kmersearch_matchscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score  
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

pg_kmersearch creates several system tables and views for metadata management and monitoring. These tables support multiple GIN indexes on the same table/column with different parameters.

### System Tables

#### kmersearch_highfreq_kmer
Stores high-frequency k-mers identified during analysis. Each k-mer is stored with its appearance count.

| Column | Type | Description |
|--------|------|-------------|
| table_oid | oid | Target table OID |
| column_name | name | Target column name |
| uintkey | bigint | K-mer encoded as integer |
| appearance_nrow | bigint | Number of rows containing this k-mer |
| detection_reason | text | Reason for detection (rate/nrow threshold) |
| created_at | timestamptz | Creation timestamp |

Primary Key: (table_oid, column_name, uintkey)

```sql
-- View excluded k-mers with appearance count
SELECT table_oid, column_name, uintkey, appearance_nrow, detection_reason, created_at
FROM kmersearch_highfreq_kmer
WHERE table_oid = 'sequences'::regclass AND column_name = 'dna_seq';

-- Find k-mers appearing in more than 1000 rows
SELECT uintkey, appearance_nrow
FROM kmersearch_highfreq_kmer
WHERE table_oid = 'sequences'::regclass
  AND column_name = 'dna_seq'
  AND appearance_nrow > 1000;
```

#### kmersearch_highfreq_kmer_meta
Stores metadata about k-mer frequency analysis. Supports multiple analyses with different k-mer sizes on the same table/column.

| Column | Type | Description |
|--------|------|-------------|
| table_oid | oid | Target table OID |
| column_name | name | Target column name |
| kmer_size | integer | K-mer size used for analysis |
| occur_bitlen | integer | Occurrence bit length setting |
| max_appearance_rate | real | Max appearance rate threshold |
| max_appearance_nrow | integer | Max appearance row count threshold |
| analysis_timestamp | timestamptz | When analysis was performed |

Primary Key: (table_oid, column_name, kmer_size)

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
Stores GIN index metadata including high-frequency filtering configuration. Tracks each GIN index separately.

| Column | Type | Description |
|--------|------|-------------|
| index_oid | oid | GIN index OID (Primary Key) |
| table_oid | oid | Target table OID |
| column_name | name | Target column name |
| highfreq_filtered | boolean | Whether high-frequency filtering is enabled |
| highfreq_source_table | name | Source table for high-frequency k-mers |
| kmer_size | integer | K-mer size for this index |
| occur_bitlen | integer | Occurrence bit length |
| max_appearance_rate | real | Max appearance rate setting |
| max_appearance_nrow | integer | Max appearance row count setting |
| created_at | timestamptz | Index creation timestamp |

```sql
-- View all GIN index metadata
SELECT index_oid, table_oid, column_name, highfreq_filtered,
       highfreq_source_table, kmer_size, occur_bitlen, created_at
FROM kmersearch_gin_index_meta;

-- Check if an index uses high-frequency filtering
SELECT highfreq_filtered, highfreq_source_table
FROM kmersearch_gin_index_meta
WHERE index_oid = 'sequences_kmer_idx'::regclass;
```

#### kmersearch_index_info
Stores comprehensive index statistics and GUC settings at index creation time. Automatically populated via event trigger when creating GIN indexes.

| Column | Type | Description |
|--------|------|-------------|
| index_oid | oid | GIN index OID (Primary Key) |
| table_oid | oid | Target table OID |
| column_name | name | Target column name |
| kmer_size | integer | K-mer size setting |
| occur_bitlen | integer | Occurrence bit length setting |
| total_nrow | bigint | Total rows in table at creation |
| highfreq_kmer_count | integer | Number of high-frequency k-mers |
| max_appearance_rate | real | Max appearance rate setting |
| max_appearance_nrow | integer | Max appearance row count setting |
| preclude_highfreq_kmer | boolean | Whether to exclude high-frequency k-mers |
| created_at | timestamptz | Index creation timestamp |

```sql
-- View index statistics with all settings
SELECT index_oid, table_oid, column_name, kmer_size, occur_bitlen,
       total_nrow, highfreq_kmer_count, max_appearance_rate,
       max_appearance_nrow, preclude_highfreq_kmer, created_at
FROM kmersearch_index_info;

-- Find indexes using high-frequency k-mer exclusion
SELECT index_oid, table_oid, column_name, kmer_size
FROM kmersearch_index_info
WHERE preclude_highfreq_kmer = true;
```

### Multiple GIN Index Support

pg_kmersearch supports creating multiple GIN indexes on the same table/column with different configurations. This is useful for:

- Testing different k-mer sizes for optimal performance
- Using different high-frequency filtering thresholds
- Comparing index performance with different settings

```sql
-- Create multiple indexes with different k-mer sizes
SET kmersearch.kmer_size = 8;
CREATE INDEX seq_kmer8_idx ON sequences USING gin(dna_seq kmersearch_dna4_gin_ops_int2);

SET kmersearch.kmer_size = 16;
CREATE INDEX seq_kmer16_idx ON sequences USING gin(dna_seq kmersearch_dna4_gin_ops_int4);

-- Each index is tracked separately in kmersearch_index_info
SELECT index_oid, kmer_size, preclude_highfreq_kmer
FROM kmersearch_index_info
WHERE table_oid = 'sequences'::regclass;
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

### Table Partitioning Functions

#### kmersearch_partition_table()
Converts a non-partitioned table to a hash-partitioned table based on DNA2/DNA4 column:

```sql
-- Convert a table to partitioned table with 4 partitions
SELECT kmersearch_partition_table(
    'sequences',        -- table name
    4                  -- number of partitions
);

-- Specify a tablespace for the partitioned table
SELECT kmersearch_partition_table(
    'sequences',        -- table name
    4,                 -- number of partitions
    'fast_ssd'         -- tablespace name (optional)
);

-- Example: Converting a large sequence table
-- Note: This preserves all data during conversion
SELECT kmersearch_partition_table('large_sequences', 16);

-- Use NULL to explicitly use the source table's tablespace
SELECT kmersearch_partition_table('large_sequences', 16, NULL);
```

Requirements:
- Table must have exactly one DNA2 or DNA4 column
- Table must not already be partitioned
- Sufficient disk space for temporary data during migration

Note: PostgreSQL does not allow explicitly specifying 'pg_default' tablespace for partitioned tables. Use NULL or omit the parameter to use the default tablespace.

#### kmersearch_unpartition_table()
Converts a hash-partitioned table back to a regular (non-partitioned) table:

```sql
-- Convert a partitioned table back to regular table
SELECT kmersearch_unpartition_table('sequences');

-- Specify a tablespace for the regular table
SELECT kmersearch_unpartition_table(
    'sequences',        -- table name
    'fast_ssd'         -- tablespace name (optional)
);

-- Use NULL to explicitly use the source table's tablespace
SELECT kmersearch_unpartition_table('sequences', NULL);
```

Requirements:
- Table must be a partitioned table (relkind = 'p')
- Table must have exactly one DNA2 or DNA4 column
- Table must have at least one partition
- Sufficient disk space for temporary data during migration

Features:
- Preserves all data during conversion
- Preserves high-frequency k-mer analysis metadata
- Progress reporting during migration
- Automatic cleanup on error

### Utility Functions

#### kmersearch_simd_capability()
Returns the detected SIMD capability of the current system:

```sql
-- Check SIMD support level
SELECT kmersearch_simd_capability();
-- Returns one of:
--   x86_64: 'None', 'AVX2', 'AVX2+BMI2', 'AVX512F', 'AVX512F+AVX512BW',
--           'AVX512F+AVX512BW+AVX512VBMI', 'AVX512F+AVX512BW+AVX512VBMI+AVX512VBMI2'
--   ARM64:  'None', 'NEON', 'NEON+SVE', 'NEON+SVE+SVE2'
```

#### kmersearch_show_buildno()
Returns the build version information:

```sql
-- Display build version
SELECT kmersearch_show_buildno();
-- Returns: '1.0.2025.12.13' (example)
```

#### kmersearch_delete_tempfiles()
Cleans up temporary files created during high-frequency k-mer analysis operations:

```sql
-- Clean up temporary files
SELECT * FROM kmersearch_delete_tempfiles();

-- Returns:
-- deleted_count: number of files deleted
-- deleted_size: total bytes freed
-- error_count: number of files that could not be deleted
```

This function removes temporary files from the `pgsql_tmp` directory that were created during `kmersearch_perform_highfreq_analysis()` operations. These files use a file-based hash table implementation for efficient k-mer counting during parallel analysis.

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

#### Query-kmer Cache

```sql
-- View query-kmer cache statistics
SELECT * FROM kmersearch_query_kmer_cache_stats();

-- Monitor cache efficiency
SELECT hits, misses, current_entries, max_entries,
       CASE WHEN (hits + misses) > 0 
            THEN hits::float / (hits + misses)::float * 100
            ELSE 0 END as hit_rate_percent
FROM kmersearch_query_kmer_cache_stats();

-- Clear query-kmer cache
SELECT kmersearch_query_kmer_cache_free();
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

-- 4. Create GIN index with appropriate operator class
-- Choose based on k-mer size and data type:
CREATE INDEX sequences_kmer_idx ON sequences USING gin(dna_seq kmersearch_dna4_gin_ops_int4);

-- 5. Load caches for optimal performance
SELECT kmersearch_highfreq_kmer_cache_load('sequences', 'dna_seq');
SELECT kmersearch_parallel_highfreq_kmer_cache_load('sequences', 'dna_seq');

-- 6. Perform searches
SELECT id, name, kmersearch_matchscore(dna_seq, 'ATCGATCG') as score
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

-- Detailed cache performance analysis
SELECT cache_type, total_entries, total_hits, total_misses,
       ROUND(hit_rate * 100, 2) as hit_rate_percent
FROM kmersearch_cache_summary
ORDER BY cache_type;
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
- **File-based hash table**: Efficient temporary storage for k-mer counting during analysis (supports uint16/uint32/uint64 keys)
- **System tables**: Metadata storage for excluded k-mers and index statistics (`kmersearch_highfreq_kmer`, `kmersearch_highfreq_kmer_meta`)
- **Cache system**: TopMemoryContext-based high-performance caching
- **SIMD optimization**: Platform-specific acceleration for encoding/decoding
  - x86_64: AVX2, BMI2, AVX512F, AVX512BW, AVX512VBMI, AVX512VBMI2
  - ARM64: NEON, SVE, SVE2
- Binary input/output support

### K-mer Search Mechanism
1. **Frequency analysis**: Parallel table scan using multiple workers to identify high-frequency k-mers
2. **K-mer extraction**: Sliding window with specified k-length (parallel processing)
3. **High-frequency filtering**: Exclude k-mers exceeding appearance thresholds
4. **N-gram key generation**: Binary encoding of k-mer + occurrence count
5. **Degenerate processing**: Expansion of MRWSYKVHDBN to standard bases
6. **Scoring**: Similarity calculation based on shared n-gram key count

## Cache Management Features

pg_kmersearch provides two types of high-performance cache systems to improve search performance:

### Actual Min Score Cache
- **Purpose**: Optimization of search condition evaluation for the `=%` operator
- **Mechanism**: Pre-calculates and caches `actual_min_score = max(kmersearch_min_score, ceil(kmersearch_min_shared_kmer_rate × query_total_kmers))`
- **Use cases**: 
  - Matching condition evaluation in `=%` operator
  - Optimization of search condition evaluation
- **Memory management**: TopMemoryContext-based implementation

### Query-kmer Cache
- **Purpose**: Performance improvement through query k-mer pattern reuse
- **Memory management**: TopMemoryContext-based implementation

### Cache Statistics and Management Functions

```sql
-- Check cache statistics
SELECT * FROM kmersearch_actual_min_score_cache_stats();
SELECT * FROM kmersearch_query_kmer_cache_stats();

-- Clear caches
SELECT kmersearch_actual_min_score_cache_free();
SELECT kmersearch_query_kmer_cache_free();
SELECT kmersearch_highfreq_kmer_cache_free_all();

-- Configure cache sizes
SET kmersearch.actual_min_score_cache_max_entries = 25000;
SET kmersearch.query_kmer_cache_max_entries = 25000;
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
ORDER BY kmersearch_matchscore(dna_seq, 'ATCGATCG') DESC;

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

## Performance Tuning

This section describes OS and PostgreSQL parameter settings recommended for optimizing pg_kmersearch performance.

### OS Parameter Settings

#### Shared Memory

```bash
# /etc/sysctl.conf or /etc/sysctl.d/99-postgresql.conf

# Maximum shared memory size (set larger than PostgreSQL's shared_buffers)
kernel.shmmax = 17179869184  # 16GB (adjust according to actual RAM)
kernel.shmall = 4194304      # shmmax / PAGE_SIZE

# Semaphore settings
kernel.sem = 250 32000 100 128
```

#### Huge Pages (Recommended)

Effective when using large amounts of memory. GIN indexes can grow large, making Huge Pages highly beneficial.

```bash
# Enable Huge Pages
vm.nr_hugepages = 5120  # shared_buffers / 2MB (Huge Page size)

# Verification command
grep Huge /proc/meminfo
```

#### Memory Management

```bash
# Suppress swap usage (when memory is sufficient)
vm.swappiness = 10

# Dirty page flush settings
vm.dirty_background_ratio = 5
vm.dirty_ratio = 10

# Or specify absolute values (recommended for large RAM systems)
vm.dirty_background_bytes = 268435456  # 256MB
vm.dirty_bytes = 1073741824            # 1GB

# Memory overcommit settings
vm.overcommit_memory = 2
vm.overcommit_ratio = 90
```

#### I/O Related

```bash
# File descriptor limit
fs.file-max = 2097152

# Maximum AIO (asynchronous I/O) requests
fs.aio-max-nr = 1048576
```

#### Disk Scheduler (for SSD)

```bash
# Set noop or mq scheduler for NVMe/SSD
echo none > /sys/block/nvme0n1/queue/scheduler
# or
echo mq-deadline > /sys/block/sda/queue/scheduler

# Adjust read-ahead (for sequential read-heavy workloads)
blockdev --setra 4096 /dev/nvme0n1
```

#### ulimit Settings (for PostgreSQL user)

```bash
# /etc/security/limits.d/postgresql.conf
postgres soft nofile 65536
postgres hard nofile 65536
postgres soft nproc 65536
postgres hard nproc 65536
postgres soft memlock unlimited
postgres hard memlock unlimited  # Required for Huge Pages
```

**Important:** When PostgreSQL is started via systemd, the settings in `/etc/security/limits.d/` are not applied. You must configure limits in systemd instead:

```bash
# Edit the systemd service override
sudo systemctl edit postgresql@16-main

# Add the following content:
[Service]
LimitNOFILE=65536

# Reload and restart to apply changes
sudo systemctl daemon-reload
sudo systemctl restart postgresql@16-main

# Verify the setting is applied
cat /proc/$(sudo head -1 /var/lib/postgresql/16/main/postmaster.pid)/limits | grep "open files"
```

#### NUMA Settings (for multi-socket CPU)

```bash
# Disable NUMA interleave to prioritize local memory
# Specify when starting PostgreSQL
numactl --interleave=all postgres -D /var/lib/pgsql/data
```

#### Applying Settings

```bash
# Apply settings immediately
sudo sysctl -p /etc/sysctl.d/99-postgresql.conf

# Restart for persistence after reboot
sudo systemctl restart postgresql
```

### PostgreSQL Parameter Settings

The following shows recommended settings for a system with 188GB RAM. Adjust values according to your environment.

#### Memory Configuration

```conf
# Memory settings (example for 188GB RAM system)
shared_buffers = 32GB                   # 17-25% of RAM (critical for GIN index caching)
temp_buffers = 1GB
effective_cache_size = 140GB            # Approximately RAM - shared_buffers
work_mem = 2GB                          # Consider parallel worker count
maintenance_work_mem = 8GB              # For GIN index creation
min_dynamic_shared_memory = 1GB
hash_mem_multiplier = 2.0               # Memory multiplier for hash joins
wal_buffers = 256MB                     # About 1/32 of shared_buffers
```

**Important: work_mem Considerations**

Maximum memory usage can be calculated as:
```
Maximum memory ≈ work_mem × max_connections × max_parallel_workers_per_gather
```

For example, with 8GB × 100 connections × 8 parallel = 6.4TB potential usage. Keep work_mem around 2GB and set per-session for specific queries:

```sql
SET work_mem = '8GB';
-- Execute heavy query
RESET work_mem;
```

#### WAL Configuration

```conf
max_wal_size = 16GB
min_wal_size = 4GB
checkpoint_timeout = 30min
checkpoint_completion_target = 0.9
```

#### I/O Cost Settings (for SSD)

```conf
random_page_cost = 1.1                  # Default 4.0 (keep 4.0 for HDD)
seq_page_cost = 1.0
effective_io_concurrency = 200
maintenance_io_concurrency = 200
```

#### Parallel Processing

```conf
max_worker_processes = 64
max_parallel_workers = 32
max_parallel_workers_per_gather = 8     # Default 2
max_parallel_maintenance_workers = 8    # Parallelism for index creation
parallel_setup_cost = 100               # Default 1000
parallel_tuple_cost = 0.001             # Default 0.01
```

#### GIN Index Settings

```conf
gin_pending_list_limit = 64MB           # Default 4MB (effective for batch inserts)
gin_fuzzy_search_limit = 0              # No limit
```

#### Planner Settings

```conf
enable_partitionwise_join = on          # For partitioned tables
enable_partitionwise_aggregate = on
default_statistics_target = 500         # 100→500 (for large tables)
```

#### JIT Settings

```conf
# Recommended to disable for k-mer search (overhead outweighs benefits)
jit = off
```

#### Huge Pages (PostgreSQL side)

```conf
huge_pages = try    # or 'on'
```

Calculate required Huge Pages:
```bash
# Check after starting PostgreSQL
grep ^VmPeak /proc/$(head -1 /var/lib/pgsql/data/postmaster.pid)/status
# Result / 2MB = required nr_hugepages
```

#### Background Writer

```conf
bgwriter_delay = 200ms
bgwriter_lru_maxpages = 100
bgwriter_lru_multiplier = 2.0
```

#### Locks and Limits

```conf
max_locks_per_transaction = 10000
max_files_per_process = 10000
```

### pg_kmersearch GUC Settings

Optimize pg_kmersearch-specific settings according to your workload:

```sql
-- Cache size settings (adjust based on query patterns)
ALTER SYSTEM SET kmersearch.query_kmer_cache_max_entries = 100000;
ALTER SYSTEM SET kmersearch.actual_min_score_cache_max_entries = 100000;

-- High-frequency k-mer cache loading batch size
ALTER SYSTEM SET kmersearch.highfreq_kmer_cache_load_batch_size = 50000;

-- High-frequency k-mer analysis hash table size (for large tables)
ALTER SYSTEM SET kmersearch.highfreq_analysis_hashtable_size = 10000000;

-- Apply settings
SELECT pg_reload_conf();
```

### Recommended Configuration Summary

| Parameter | Recommended Value | Notes |
|-----------|-------------------|-------|
| `shared_buffers` | 17-25% of RAM | Critical for GIN index performance |
| `effective_cache_size` | 70-80% of RAM | Used by planner for cost estimation |
| `work_mem` | 1-4GB | Careful with parallel queries |
| `maintenance_work_mem` | 4-8GB | For index creation |
| `random_page_cost` | 1.1 (SSD) / 4.0 (HDD) | I/O cost for planner |
| `jit` | off | Disable for k-mer search |
| `huge_pages` | try or on | Reduces TLB misses |
| `vm.swappiness` | 10 | Minimize swapping |

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