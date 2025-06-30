# pg_kmersearch

A PostgreSQL extension for DNA sequence similarity search using k-mer indexing.

## Features

- **DNA2** and **DNA4** data types for storing DNA sequences
- **GIN index** support for fast k-mer based similarity search
- **Configurable parameters** for k-mer size and search thresholds
- **Similarity scoring** functions for sequence comparison

## Quick Start

```sql
CREATE EXTENSION pg_kmersearch;
CREATE TABLE sequences (id SERIAL, name TEXT, dna dna2);
CREATE INDEX ON sequences USING gin(dna);
INSERT INTO sequences (name, dna) VALUES ('seq1', 'ATCGATCG...');
SELECT * FROM sequences WHERE dna =% 'ATCGATCG';
```

For detailed documentation, see `doc/pg_kmersearch.en.md`.