# pg_kmersearch

A PostgreSQL extension for DNA sequence similarity search using k-mer indexing.

## Features

- **DNA2** and **DNA4** data types for storing DNA sequences
- **GIN index** support for fast k-mer based similarity search
- **Configurable parameters** for k-mer size and search thresholds
- **Similarity scoring** functions for sequence comparison
- **Table partitioning** support for large-scale sequence databases
- **Parallel high-frequency k-mer analysis** for optimal index performance
- **SIMD optimizations** for k-mer extraction and sequence operations
- **High-performance caching** for query patterns and search conditions

## Quick Start

```sql
CREATE EXTENSION pg_kmersearch;
CREATE TABLE sequences (id SERIAL, name TEXT, dna DNA2);
CREATE INDEX ON sequences USING gin(dna);
INSERT INTO sequences (name, dna) VALUES ('seq1', 'ATCGATCG...');
SELECT * FROM sequences WHERE dna =% 'ATCGATCG';
```

## Key Features in Detail

### Data Types
- **DNA2**: 2-bit encoding for standard DNA (ACGT)
- **DNA4**: 4-bit encoding with degenerate code support (MRWSYKVHDBN)

### Performance Features
- **GIN indexing** with k-mer based search (k=4-32)
- **High-frequency k-mer exclusion** for improved search performance
- **Parallel processing** for k-mer frequency analysis
- **SIMD optimizations** (AVX2/AVX512/NEON/SVE)
- **Multi-tier caching** system

### Analysis Functions
- `kmersearch_perform_highfreq_analysis()` - Parallel k-mer frequency analysis
- `kmersearch_partition_table()` - Convert tables to hash partitions
- `kmersearch_rawscore()`/`kmersearch_correctedscore()` - Similarity scoring

For detailed documentation, see `doc/pg_kmersearch.en.md`.