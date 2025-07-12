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
1. Execute non-sudo commands (make clean, make)
2. Present required sudo commands for manual execution:
   - `sudo systemctl stop postgresql`
   - `sudo make install` 
   - `sudo systemctl start postgresql`
3. Wait for user confirmation before proceeding with testing

**CRITICAL TESTING NOTE:**
Running `make installcheck` after source code modifications without executing `sudo make install` first is meaningless and will not test the modified code. Always follow this exact sequence:
1. `make clean && make` (Build)
2. `sudo systemctl stop postgresql` (Stop PostgreSQL)
3. `sudo make install` (Install)
4. `sudo systemctl start postgresql` (Start PostgreSQL)
5. `make installcheck` (Run tests)

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
- **High-frequency k-mer management**: `kmersearch_analyze_table()` for high-frequency k-mer analysis and cache management functions

## Configuration System

Important GUC variables:
- `kmersearch.kmer_size` (4-64): K-mer length for indexing (default: 16)
- `kmersearch.occur_bitlen` (0-16): Bits for occurrence count storage (default: 8)
- `kmersearch.max_appearance_rate` (0.0-1.0): Max k-mer appearance threshold (default: 0.5)
- `kmersearch.max_appearance_nrow` (0-∞): Maximum rows containing k-mer (default: 0, unlimited)
- `kmersearch.min_score` (0-∞): Minimum similarity score for search results (default: 1)
- `kmersearch.min_shared_ngram_key_rate` (0.0-1.0): Minimum threshold for shared n-gram key rate (default: 0.9)
- `kmersearch.rawscore_cache_max_entries` (1000-10000000): Maximum entries for rawscore cache (default: 50000)
- `kmersearch.query_pattern_cache_max_entries` (1000-10000000): Maximum entries for query pattern cache (default: 50000)
- `kmersearch.actual_min_score_cache_max_entries` (1000-10000000): Maximum entries for actual min score cache (default: 50000)
- `kmersearch.preclude_highfreq_kmer` (true/false): Enable high-frequency k-mer exclusion during GIN index construction (default: false)
- `kmersearch.force_use_parallel_highfreq_kmer_cache` (true/false): Force use of dshash parallel cache for high-frequency k-mer lookups (default: false)

## Architecture Overview

### Core Components

**kmersearch.c** - Main module with initialization, configuration management, and SIMD dispatch setup
**kmersearch_datatype.c** - DNA2/DNA4 data type implementation with encoding/decoding
**kmersearch_gin.c** - GIN index support for k-mer indexing and search operators
**kmersearch_kmer.c** - K-mer extraction, encoding, and utility functions
**kmersearch_cache.c** - Multi-layered cache system (rawscore, query pattern, actual min score)
**kmersearch_freq.c** - High-frequency k-mer analysis and parallel processing

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

## System Tables and Views

### System Tables

#### kmersearch_highfreq_kmer
Stores high-frequency k-mers that are excluded from GIN indexes

#### kmersearch_highfreq_kmer_meta
Stores metadata about k-mer frequency analysis

#### kmersearch_gin_index_meta
Stores GIN index metadata including high-frequency filtering information

#### kmersearch_index_info
Stores comprehensive index statistics and configuration

### Management Views

#### kmersearch_cache_summary
Provides unified cache statistics across all cache types

#### kmersearch_analysis_status
Shows high-frequency k-mer analysis status for all analyzed tables

## Analysis and Management Functions

### K-mer Frequency Analysis

#### kmersearch_analyze_table()
Performs parallel k-mer frequency analysis on a table

#### kmersearch_drop_analysis()
Removes analysis data and frees storage

### High-Frequency K-mer Cache Management

#### Global Cache Functions
- `kmersearch_highfreq_kmer_cache_load(table_name, column_name)`: Load high-frequency k-mers into global cache
- `kmersearch_highfreq_kmer_cache_free(table_name, column_name)`: Free global cache for specific table/column
- `kmersearch_highfreq_kmer_cache_free_all()`: Free all entries from global cache

#### Parallel Cache Functions
- `kmersearch_parallel_highfreq_kmer_cache_load(table_name, column_name)`: Load high-frequency k-mers into parallel cache
- `kmersearch_parallel_highfreq_kmer_cache_free(table_name, column_name)`: Free parallel cache for specific table/column
- `kmersearch_parallel_highfreq_kmer_cache_free_all()`: Free all entries from parallel cache

### Cache Statistics and Management

#### Cache Statistics Functions
- `kmersearch_rawscore_cache_stats()`: View rawscore cache statistics
- `kmersearch_query_pattern_cache_stats()`: View query pattern cache statistics
- `kmersearch_actual_min_score_cache_stats()`: View actual min score cache statistics

#### Cache Management Functions
- `kmersearch_rawscore_cache_free()`: Clear rawscore cache
- `kmersearch_query_pattern_cache_free()`: Clear query pattern cache
- `kmersearch_actual_min_score_cache_free()`: Clear actual min score cache

## Complex Type Definitions

### kmersearch_analysis_result
Returned by `kmersearch_analyze_table()`:
- total_rows: bigint
- highfreq_kmers_count: integer
- parallel_workers_used: integer
- analysis_duration: real
- max_appearance_rate_used: real
- max_appearance_nrow_used: integer

### kmersearch_drop_result
Returned by `kmersearch_drop_analysis()`:
- dropped_analyses: integer
- dropped_highfreq_kmers: integer
- freed_storage_bytes: bigint

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

### Recent Code Improvements
- **Unified ngram_key2 creation**: Removed duplicate `kmersearch_create_ngram_key2_with_occurrence()` function and consolidated implementation in `kmersearch_create_ngram_key2()` for better maintainability and consistency
- **Enhanced error handling**: Improved negative occurrence value handling in ngram key creation functions
- **Removed testing artifacts**: Eliminated dummy data generation code from parallel analysis functions to ensure production data integrity
- **Fixed placeholder implementations**: Updated non-functional worker function with proper initialization and clear documentation of remaining work
- **Consolidated cache management**: Eliminated duplicate LRU eviction function and unified cache operation interfaces

## Limitations

- Query sequences must be at least 8 bases long
- Degenerate code expansion limited to 10 combinations (skipped if exceeded)
- Occurrence counts capped at maximum value for configured bit length
- Case-insensitive input, uppercase output
- High-frequency k-mer exclusion requires full table scan during index creation
- Tables with GIN k-mer indexes have restricted INSERT/UPDATE/DELETE operations

## License

This project is released under the MIT License.