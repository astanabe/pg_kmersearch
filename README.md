# pg_kmersearch

A PostgreSQL extension for DNA sequence similarity search using k-mer indexing.

## Features

- **DNA2** and **DNA4** data types for storing DNA sequences
- **GIN index** support with multiple storage strategies (int2/int4/int8 keys)
- **Configurable parameters** via GUC variables for k-mer size and search thresholds
- **Similarity scoring** functions with degenerate base support
- **High-frequency k-mer filtering** for improved index performance
- **Parallel processing** support for large-scale analysis
- **SIMD optimizations** for encoding/decoding and k-mer extraction
- **Multi-tier caching** system for query patterns and analysis results

## Quick Start

```sql
CREATE EXTENSION pg_kmersearch;

-- Create a table with DNA sequences
CREATE TABLE sequences (id SERIAL, name TEXT, dna DNA4);

-- Create GIN index (default: int4 keys for k-mer size 16)
CREATE INDEX ON sequences USING gin(dna kmersearch_dna4_gin_ops_int4);

-- Insert sequences
INSERT INTO sequences (name, dna) VALUES ('seq1', 'ATCGATCG');

-- Search using k-mer similarity
SELECT * FROM sequences WHERE dna =% 'ATCGATCG';

-- Calculate similarity scores
SELECT name, kmersearch_correctedscore(dna, 'ATCGATCG') AS score 
FROM sequences ORDER BY score DESC;
```

## Key Features in Detail

### Data Types
- **DNA2**: 2-bit encoding for standard DNA sequences (ACGT, U treated as T)
- **DNA4**: 4-bit encoding for full DNA (ACGT) with IUPAC degenerate base support (MRWSYKVHDBN)

### Index Support
- **GIN operator classes** with different key storage strategies:
  - `kmersearch_dna2_gin_ops_int2`: 16-bit keys (k-mer size ≤ 8 for DNA2)
  - `kmersearch_dna2_gin_ops_int4`: 32-bit keys (k-mer size ≤ 16 for DNA2)
  - `kmersearch_dna2_gin_ops_int8`: 64-bit keys (k-mer size ≤ 32 for DNA2)
  - `kmersearch_dna4_gin_ops_int2`: 16-bit keys (k-mer size ≤ 8 for DNA4)
  - `kmersearch_dna4_gin_ops_int4`: 32-bit keys (k-mer size ≤ 16 for DNA4)
  - `kmersearch_dna4_gin_ops_int8`: 64-bit keys (k-mer size ≤ 32 for DNA4)
- **BTree and Hash** operator classes for standard operations

### Key Functions

#### Search and Scoring
- `=%` operator: K-mer based sequence search
- `kmersearch_rawscore()`: Calculate raw matching score
- `kmersearch_correctedscore()`: Calculate score with degenerate base adjustments

#### High-frequency K-mer Management
- `kmersearch_perform_highfreq_analysis()`: Analyze and identify high-frequency k-mers
- `kmersearch_undo_highfreq_analysis()`: Remove analysis results
- `kmersearch_highfreq_kmer_cache_load()`: Load high-frequency k-mers into cache
- `kmersearch_highfreq_kmer_cache_free()`: Free cache memory

#### Cache Management
- `kmersearch_query_pattern_cache_stats()`: View query pattern cache statistics
- `kmersearch_actual_min_score_cache_stats()`: View minimum score cache statistics
- `kmersearch_query_pattern_cache_free()`: Clear query pattern cache
- `kmersearch_actual_min_score_cache_free()`: Clear minimum score cache

#### Utility Functions
- `kmersearch_simd_capability()`: Check SIMD support level
- `kmersearch_partition_table()`: Convert table to hash partitions
- `bit_length()`, `nuc_length()`, `char_length()`, `length()`: Get sequence lengths

### Configuration Parameters (GUC Variables)

- `kmersearch.kmer_size` (default: 16): K-mer size for indexing (4-32)
- `kmersearch.occur_bitlen` (default: 8): Bits for occurrence count (0-16)
- `kmersearch.max_appearance_rate` (default: 0.5): Max appearance rate threshold (0.0-1.0)
- `kmersearch.max_appearance_nrow` (default: 0): Max rows threshold (0 = unlimited)
- `kmersearch.min_score` (default: 1): Minimum score for search results
- `kmersearch.min_shared_ngram_key_rate` (default: 0.9): Min shared n-gram rate (0.0-1.0)
- `kmersearch.preclude_highfreq_kmer` (default: false): Enable high-frequency filtering
- `kmersearch.query_pattern_cache_max_entries` (default: 50000): Query cache size
- `kmersearch.actual_min_score_cache_max_entries` (default: 50000): Score cache size

### System Tables and Views

#### System Tables
- `kmersearch_highfreq_kmer`: Stores identified high-frequency k-mers
- `kmersearch_highfreq_kmer_meta`: Metadata for high-frequency analysis
- `kmersearch_gin_index_meta`: GIN index metadata
- `kmersearch_index_info`: General index information

#### Management Views
- `kmersearch_cache_summary`: Overview of cache statistics
- `kmersearch_analysis_status`: Status of high-frequency k-mer analyses

For detailed documentation, see `doc/pg_kmersearch.en.md`.