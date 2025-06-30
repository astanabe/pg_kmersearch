# pg_kmersearch

PostgreSQL extension for DNA sequence data types

## Overview

pg_kmersearch is a PostgreSQL extension that provides custom data types for efficiently storing and processing DNA sequence data. This extension implements two data types:

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

### Internal Implementation
- DNA2 type: 2 bits per character encoding
- DNA4 type: 4 bits per character encoding
- Both types implemented as PostgreSQL variable-length data types
- Binary input/output support

## Limitations

- Currently implements only basic data type functionality
- Comparison operators and search functions are planned for future releases
- Case-insensitive input, uppercase output

## License

This project is released under the MIT License.

## Contributing

Please report bugs and feature requests on the GitHub Issues page.