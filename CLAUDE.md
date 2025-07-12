# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Development Commands

### Building the Extension
```bash
make                    # Build the extension
sudo make install       # Install to PostgreSQL
make clean             # Clean build artifacts
```

### Testing

**IMPORTANT NOTE FOR CLAUDE CODE USERS:**
Claude Code cannot execute `sudo` commands. When testing is required, Claude will provide the necessary commands for manual execution. The user must run these commands with appropriate privileges.

```bash
# IMPORTANT: Before running installcheck, you MUST rebuild and reinstall the extension
# especially after SQL function name changes or other modifications
sudo systemctl stop postgresql
make clean
make
sudo make install
sudo systemctl start postgresql
make installcheck      # Run regression tests using REGRESS variable
```

**Complete workflow for testing after modifications:**
```bash
# Drop existing extension and rebuild completely
psql -d postgres -c "DROP EXTENSION IF EXISTS pg_kmersearch CASCADE;"
sudo systemctl stop postgresql
make clean
make
sudo make install
sudo systemctl start postgresql
make installcheck
```

**Claude Code Workflow:**
When Claude Code needs to perform testing, it will:
1. Execute non-sudo commands (make clean, make, make installcheck)
2. Present required sudo commands for manual execution:
   - `sudo systemctl stop postgresql`
   - `sudo make install` 
   - `sudo systemctl start postgresql`
3. Wait for user confirmation before proceeding with testing

The test suite includes 12 test files covering:
- Basic DNA2/DNA4 data types (01_basic_types.sql)
- Configuration management (02_configuration.sql) 
- Table and index operations (03_tables_indexes.sql)
- Search operators (04_search_operators.sql)
- Scoring functions (05_scoring_functions.sql)
- Advanced search features (06_advanced_search.sql)
- Length functions (07_length_functions.sql)
- Cache management (08_cache_management.sql)
- High-frequency k-mer filtering (09_highfreq_filter.sql)
- Parallel cache operations (10_parallel_cache.sql)
- Cache hierarchy (11_cache_hierarchy.sql)
- Management views (12_management_views.sql)

Test results are compared between `results/` and `expected/` directories.

## Architecture Overview

### Core Components

**kmersearch.c** - Main module with initialization, configuration management, and SIMD dispatch setup
**kmersearch_datatype.c** - DNA2/DNA4 data type implementation with encoding/decoding
**kmersearch_gin.c** - GIN index support for k-mer indexing and search operators
**kmersearch_kmer.c** - K-mer extraction, encoding, and utility functions
**kmersearch_cache.c** - Multi-layered cache system (rawscore, query pattern, actual min score)
**kmersearch_freq.c** - High-frequency k-mer analysis and parallel processing

### Data Types

- **DNA2**: 2-bit encoding for standard bases (A,C,G,T) - 2 bits per nucleotide
- **DNA4**: 4-bit encoding supporting degenerate codes (MRWSYKVHDBN) - 4 bits per nucleotide

### Key Features

- **GIN indexing**: Fast k-mer based similarity search using n-gram keys
- **High-frequency k-mer exclusion**: Automatic filtering during index creation to improve performance
- **Multi-level caching**: Global cache, parallel cache (dshash), and table fallback for high-frequency k-mers
- **SIMD optimization**: Platform-specific optimizations for x86_64 (AVX2/AVX512) and ARM64 (NEON/SVE)
- **Parallel analysis**: Support for parallel workers in k-mer frequency analysis

### Configuration System

Important GUC variables:
- `kmersearch.kmer_size` (4-64): K-mer length for indexing
- `kmersearch.occur_bitlen` (0-16): Bits for occurrence count storage
- `kmersearch.max_appearance_rate` (0.0-1.0): Max k-mer appearance threshold
- `kmersearch.min_score`: Minimum similarity score for search results
- `kmersearch.preclude_highfreq_kmer`: Enable high-frequency k-mer exclusion

### Search Operations

- **=% operator**: Main similarity search operator
- **rawscore functions**: Calculate shared k-mer counts
- **correctedscore functions**: Adjust scores for excluded high-frequency k-mers

### Cache Architecture

Three-tier cache system for performance optimization:
1. **Global cache**: TopMemoryContext-based for single process
2. **Parallel cache**: DSM/dshash for cross-process sharing  
3. **Table fallback**: Direct system table access

All caches validate GUC parameters against stored metadata to ensure consistency.

## Development Notes

### File Organization
- SQL extension definition: `pg_kmersearch--1.0.sql`
- Test cases: `sql/*.sql` with expected outputs in `expected/`
- Documentation: `doc/pg_kmersearch.{en,ja}.md`

### Key Technical Details
- Built on PostgreSQL's varbit type for efficient bit storage
- N-gram keys combine k-mer data (2k bits) + occurrence count (8-16 bits)
- Degenerate expansion limited to 10 combinations to prevent explosion
- System tables `kmersearch_highfreq_kmer` and `kmersearch_highfreq_kmer_meta` store analysis metadata

### SIMD Implementation
The codebase includes platform-specific SIMD optimizations with function dispatch tables for encoding/decoding operations. SIMD capability is detected at runtime and appropriate implementation is selected.