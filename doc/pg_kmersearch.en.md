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
- Binary input/output support

### K-mer Search Mechanism
1. **Frequency analysis**: Full table scan to identify high-frequency k-mers
2. **K-mer extraction**: Sliding window with specified k-length
3. **High-frequency filtering**: Exclude k-mers exceeding appearance thresholds
4. **N-gram key generation**: Binary encoding of k-mer + occurrence count
5. **Degenerate processing**: Expansion of MRWSYKVHDBN to standard bases
6. **Scoring**: Similarity calculation based on shared n-gram key count

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