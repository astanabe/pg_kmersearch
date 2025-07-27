# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Implementation Constraints

**IMPORTANT**: The following implementations are PROHIBITED:
- Stub implementations are forbidden
- Dummy data generation code is forbidden
- Incomplete function implementations are forbidden
- Unnecessary wrapper functions or helper functions are forbidden
- Workaround implementations are forbidden
- Implementing functions, structures, or variables with the same functionality as existing ones is forbidden
- When GUC variables exist, use them directly - assigning GUC variables to local variables is forbidden

## Important Development Assumptions

**Extension Version Management:**
This extension development assumes **NEW INSTALLATIONS ONLY**. Version upgrades, migration scripts, and backward compatibility are NOT considered. When making changes:
- Always assume users will perform fresh installations
- No need to create ALTER EXTENSION scripts
- No need to maintain compatibility with previous versions
- Breaking changes to function signatures, data types, or SQL interfaces are acceptable

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

The test suite includes 13 test files covering:
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
- Partition functions (13_partition_functions.sql)

Test results are compared between `results/` and `expected/` directories.

**IMPORTANT GUC Configuration Note:**
PostgreSQL sessions reset GUC variables to their default values when a new session starts. Therefore, when using `psql` commands, you must set the required GUC variables in each new session:

```bash
# Example: Setting GUC variables for each psql session
psql -d postgres -c "
SET kmersearch.kmer_size = 4;
SET kmersearch.max_appearance_rate = 0.25;
-- Your analysis commands here
"
```

This is particularly important when testing analysis functions and cache operations that depend on specific GUC variable values.

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
- **K-mer length**: 4-32 bases (specified at index creation)
- **GIN indexing**: Fast search using n-gram keys
- **Degenerate code support**: MRWSYKVHDBN expansion for DNA4 type
- **Occurrence tracking**: Counts k-mer occurrences within rows (default 8-bit)
- **Scoring search**: Retrieve top matches by similarity, not just exact matches
- **High-frequency k-mer exclusion**: Automatically excludes overly common k-mers during index creation
- **Score-based filtering**: Minimum score thresholds with automatic adjustment for excluded k-mers
- **Score calculation functions**: `kmersearch_rawscore()` and `kmersearch_correctedscore()` for individual sequence scoring
- **High-frequency k-mer management**: `kmersearch_perform_highfreq_analysis()` for high-frequency k-mer analysis and `kmersearch_undo_highfreq_analysis()` for analysis data removal, plus cache management functions

## Configuration System

Important GUC variables:
- `kmersearch.kmer_size` (4-32): K-mer length for indexing (default: 16)
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
- `kmersearch.highfreq_kmer_cache_load_batch_size` (1000-1000000): Batch size for loading high-frequency k-mers into cache (default: 10000)
- `kmersearch.highfreq_analysis_batch_size` (1000-1000000): Batch size for high-frequency k-mer analysis (default: 10000)

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

#### kmersearch_perform_highfreq_analysis()
Performs parallel k-mer frequency analysis on a table
- Supports table name or OID as first argument
- Supports column name or attnum as second argument
- Uses PostgreSQL's standard parallel execution framework
- Dynamically distributes work across multiple parallel workers

#### kmersearch_undo_highfreq_analysis()
Removes analysis data and frees storage

### Table Partitioning Functions

#### kmersearch_partition_table(table_name text, partition_count int, tablespace_name text DEFAULT NULL)
Converts a non-partitioned table to a hash-partitioned table based on DNA2/DNA4 column
- Requires exactly one DNA2 or DNA4 column in the table
- Creates hash partitions based on the DNA column
- Preserves all data during conversion
- Handles SERIAL column dependencies with CASCADE
- Optional tablespace_name parameter to specify target tablespace (NULL uses source table's tablespace)
- Note: Cannot explicitly specify 'pg_default' tablespace for partitioned tables (use NULL instead)


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


## Scoring Functions: rawscore vs correctedscore

**IMPORTANT**: The distinction between rawscore and correctedscore functions is critical for understanding high-frequency k-mer exclusion functionality:

### rawscore Functions
- **Purpose**: Calculate match count based on GIN index capabilities
- **Behavior**: 
  - When GIN index exists WITH high-frequency k-mer exclusion (`kmersearch.preclude_highfreq_kmer = true`): Returns match count with high-frequency k-mers excluded
  - When GIN index exists WITHOUT high-frequency k-mer exclusion: Returns full match count
  - When no GIN index exists: Returns full match count
- **Performance**: Fast (utilizes GIN index when available)
- **Use case**: Primary scoring function for search operations

### correctedscore Functions  
- **Purpose**: Calculate unfiltered match count through direct comparison
- **Behavior**: Always returns full match count by directly comparing ngram_key2 extracted from query sequence and DNA2/DNA4 data, regardless of GIN index configuration
- **Performance**: Slower (direct sequence comparison, no index utilization)
- **Use case**: Baseline comparison to understand the effect of high-frequency k-mer exclusion

### Expected Relationship
**IMPORTANT NOTE**: In the current implementation, `rawscore` and `correctedscore` functions always return the same values. The original design intention was for correctedscore to provide unfiltered results, but both functions now operate identically. This behavior is expected and not a bug.

## Complex Type Definitions

### kmersearch_analysis_result
Returned by `kmersearch_perform_highfreq_analysis()`:
- total_rows: bigint
- highfreq_kmers_count: integer
- parallel_workers_used: integer
- max_appearance_rate_used: real
- max_appearance_nrow_used: integer

### kmersearch_drop_result
Returned by `kmersearch_undo_highfreq_analysis()`:
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

### Length Functions

pg_kmersearch provides several length functions for DNA2 and DNA4 types that correctly handle padding and return accurate measurements:

#### Available Length Functions
- **`bit_length(DNA2/DNA4)`**: Returns the actual bit length (excluding padding)
- **`nuc_length(DNA2/DNA4)`**: Returns the number of nucleotides
- **`char_length(DNA2/DNA4)`**: Same as `nuc_length()` (character count)
- **`length(DNA2/DNA4)`**: Same as `nuc_length()` (standard length function)

#### Function Relationships
- **DNA2**: `nuc_length() = bit_length() / 2` (2 bits per nucleotide)
- **DNA4**: `nuc_length() = bit_length() / 4` (4 bits per nucleotide)
- **Consistency**: `char_length() = length() = nuc_length()`

### SIMD Implementation
The codebase includes platform-specific SIMD optimizations with function dispatch tables for encoding/decoding operations. SIMD capability is detected at runtime and appropriate implementation is selected.

### Parallel Cache Technical Implementation

- **Memory Management**: DSM segment pinning prevents automatic cleanup
- **Process Identification**: Differentiated resource management for main process and workers
- **Error Handling**: Comprehensive PG_TRY/PG_CATCH for safe operations
- **Lock Management**: Proper concurrent control with dshash_release_lock()

## Development Notes

### Repository File Management

**IMPORTANT: Files to EXCLUDE from Git Repository**

The following types of files should NEVER be committed to the git repository:

**Compiled Files and Build Artifacts:**
- `*.bc` - LLVM bitcode files
- `*.o` - Object files  
- `*.so` - Shared library files (e.g., `pg_kmersearch.so`)
- Any other compiled binaries or intermediate build files

**Regression Test Output Files:**
- `results/*.out` - Generated test result files
- `results/` directory itself (created during testing)

**Temporary and Debug Files:**
- `test_*.sql` - Temporary test files created for debugging
- `debug_*.sql` - Debug scripts
- `*.tmp` - Temporary files
- `*.log` - Log files
- `*~` - Backup files created by editors
- `.DS_Store` - macOS system files
- `Thumbs.db` - Windows thumbnail cache

**Note:** The `expected/` directory and its `*.out` files ARE required as they contain the expected test results for regression testing.

Use `make clean` to remove build artifacts before committing changes. Always verify with `git status` that no temporary or compiled files are being tracked.

### File Organization
- SQL extension definition: `pg_kmersearch--1.0.sql`
- Test cases: `sql/*.sql` with expected outputs in `expected/`
- Documentation: `doc/pg_kmersearch.{en,ja}.md`

### Key Technical Details
- Built on PostgreSQL's varbit type for efficient bit storage
- N-gram keys combine k-mer data (2k bits) + occurrence count (8-16 bits)
- Degenerate expansion limited to 10 combinations to prevent explosion
- System tables `kmersearch_highfreq_kmer` and `kmersearch_highfreq_kmer_meta` store analysis metadata

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

## PostgreSQL Parallel Processing Guidelines

### Critical Lessons Learned

Based on extensive debugging of parallel processing issues in `kmersearch_perform_highfreq_analysis()`, the following guidelines must be followed when implementing parallel processing in PostgreSQL extensions:

#### 1. Avoid PG_TRY/PG_CATCH in Parallel Workers
**Problem**: PG_TRY/PG_CATCH blocks in parallel workers can interfere with PostgreSQL's parallel execution context, causing unexpected behavior including attempts to execute SQL operations during cleanup.

**Solution**: Remove all PG_TRY/PG_CATCH blocks from parallel worker code. Use simple error checking with elog(ERROR) instead.

```c
/* BAD - Don't use in parallel workers */
PG_TRY();
{
    /* worker code */
}
PG_CATCH();
{
    /* cleanup that might trigger SQL */
    PG_RE_THROW();
}
PG_END_TRY();

/* GOOD - Simple error handling */
if (!condition)
    elog(ERROR, "Error message");
```

#### 2. Parallel Worker Restrictions
Parallel workers CANNOT:
- Execute SQL commands (INSERT, UPDATE, DELETE, SELECT)
- Use SPI (Server Programming Interface)
- Perform transaction control operations
- Access certain global variables safely

Always check IsParallelWorker() before any potentially restricted operation:
```c
if (IsParallelWorker())
    return;  /* or error out */
```

#### 3. TOAST Data Handling
**Problem**: Large data (>2KB) may be TOAST-compressed. Parallel workers crash when accessing TOAST pointers directly.

**Solution**: Always detoast data in parallel workers:
```c
/* Get properly detoasted sequence */
seq = DatumGetVarBitP(sequence_datum);

/* Use the detoasted data */
process_sequence(seq);

/* Free if a copy was made */
if ((void *)seq != DatumGetPointer(sequence_datum))
    pfree(seq);
```

#### 4. Global Variable Management
**Problem**: Static global variables can be inherited by parallel workers but should not be used.

**Solution**: 
- Mark globals as main-process only with comments
- Set to NULL at worker startup
- Add IsParallelWorker() checks before access

```c
/* IMPORTANT: These are only valid in the main process, NOT in parallel workers */
static dsm_segment *analysis_dsm_segment = NULL;

/* In worker function */
if (!IsParallelWorker())
{
    /* Access global variables */
}
```

#### 5. DSM/DSA/dshash Architecture
**Correct Pattern**:
1. Create DSM/DSA/dshash resources BEFORE EnterParallelMode()
2. Pass handles through shm_toc
3. Workers only attach and read, never create or modify structure

```c
/* Parent process - before EnterParallelMode() */
create_dshash_resources(&ctx);

/* After InitializeParallelDSM() */
handles->dsm_handle = ctx.dsm_handle;
shm_toc_insert(toc, KEY, handles);

/* Worker - read only */
handles = shm_toc_lookup(toc, KEY, false);
dsm_attach(handles->dsm_handle);
```

#### 6. Resource Cleanup
**Problem**: Cleanup functions may execute SQL or other restricted operations.

**Solution**: Guard all cleanup with IsParallelWorker() checks:
```c
static void cleanup_resources(void)
{
    /* Never perform cleanup in parallel workers */
    if (IsParallelWorker())
        return;
    
    /* Cleanup code */
}
```

#### 7. Error Handling Best Practices
- Use simple NULL checks and elog(ERROR) in workers
- Avoid complex error handling that might trigger SQL
- Let PostgreSQL handle worker errors through shared state

#### 8. Parallel Mode Exit Timing
**Critical Issue**: PostgreSQL remains in parallel mode until ExitParallelMode() is called, even after WaitForParallelWorkersToFinish().

**Problem**: Attempting SQL operations after workers finish but before ExitParallelMode() causes "cannot execute INSERT during a parallel operation" errors.

**Solution**: Always exit parallel mode BEFORE executing any SQL operations:
```c
/* BAD - SQL operations while still in parallel mode */
WaitForParallelWorkersToFinish(pcxt);
SPI_connect();
SPI_execute("INSERT ...", false, 0);  /* ERROR! */
DestroyParallelContext(pcxt);
ExitParallelMode();

/* GOOD - Exit parallel mode first */
WaitForParallelWorkersToFinish(pcxt);
DestroyParallelContext(pcxt);
ExitParallelMode();
/* Now safe to execute SQL */
SPI_connect();
SPI_execute("INSERT ...", false, 0);
```

This is a common pitfall because IsParallelWorker() returns false in the main process, but PostgreSQL still considers it to be in parallel mode until ExitParallelMode() is called.

#### 9. Common Parallel Processing Patterns

**Correct Parallel Execution Flow**:
```c
/* 1. Create shared resources (before parallel mode) */
create_shared_resources();

/* 2. Enter parallel mode */
EnterParallelMode();

/* 3. Create and configure parallel context */
pcxt = CreateParallelContext("worker_function", nworkers);
InitializeParallelDSM(pcxt);

/* 4. Launch workers */
LaunchParallelWorkers(pcxt);

/* 5. Wait for completion */
WaitForParallelWorkersToFinish(pcxt);

/* 6. Exit parallel mode BEFORE any SQL operations */
DestroyParallelContext(pcxt);
ExitParallelMode();

/* 7. Now safe to perform SQL operations */
SPI_connect();
/* SQL operations */
SPI_finish();

/* 8. Clean up shared resources */
cleanup_shared_resources();
```

#### 10. Testing Parallel Code
Always test with:
- Small data (non-TOAST)
- Large data (TOAST-compressed)
- Multiple workers
- Worker failures
- Check PostgreSQL logs for warnings about DSM leaks or SQL execution attempts
- Verify no "cannot execute ... during a parallel operation" errors

## License

This project is released under the MIT License.