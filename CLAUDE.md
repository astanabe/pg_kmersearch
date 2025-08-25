# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## CRITICAL DATABASE SAFETY RULES

**ABSOLUTELY FORBIDDEN - NEVER EXECUTE THESE COMMANDS:**
- **NEVER** use `DROP EXTENSION ... CASCADE` - this will destroy user data tables
- **NEVER** use `DROP TABLE` on any user data tables
- **NEVER** use `TRUNCATE` on any user data tables
- **NEVER** execute any destructive operations without explicit user confirmation
- **NEVER** assume data can be easily restored - always assume NO BACKUPS EXIST

**SAFE EXTENSION UPDATE PROCEDURES:**
- To update extension code, use one of these methods:
  1. Restart PostgreSQL to reload shared libraries: `sudo systemctl restart postgresql`
  2. Use `CREATE OR REPLACE FUNCTION` for individual function updates
  3. Drop and recreate ONLY the extension functions, never the data
- **ALWAYS** ask for explicit permission before any DROP operations
- **ALWAYS** warn about potential data loss before suggesting risky commands

## Language Usage Guidelines

- **Dialogue/Conversation**: Use Japanese (日本語) for all interactions with the user
- **Source Code**: Use English only, including all comments in source files
- **Documentation**: Use English only for:
  - Source code comments
  - Git commit messages
  - CLAUDE.md file
  - README.md file
  - Other technical documentation files
- **Note**: Japanese should NOT be used in any code, documentation, or version control messages

## Implementation Constraints

**IMPORTANT**: The following implementations are PROHIBITED:
- Stub implementations are forbidden
- Dummy data generation code is forbidden
- Incomplete function implementations are forbidden
- Unnecessary wrapper functions or helper functions are forbidden
- Workaround implementations are forbidden
- Implementing functions, structures, or variables with the same functionality as existing ones is forbidden
- When GUC variables exist, use them directly - assigning GUC variables to local variables is forbidden
- Unnecessary memory copying and data conversions are forbidden - always minimize memory allocations and data copies
- In GIN consistent functions, any data conversion or transformation is absolutely prohibited - use data as-is
- **CRITICAL**: Separate implementations for int2/int4/int8 and uint16/uint32/uint64 exist for performance and memory optimization reasons. These MUST NEVER be merged or consolidated - they are intentionally separate for speed and memory efficiency

## Code Analysis Requirements

**IMPORTANT**: When analyzing code:
- Statements about "possibilities" must always be verified and confirmed with certainty
- Never use uncertain phrases like "might be used", "possibly", "may be" without immediate verification
- Always check actual usage with grep, read, or other tools to provide definitive answers
- If something cannot be determined with certainty, explicitly state what needs to be checked
- **Never make descriptions or changes based on speculation or assumptions - always verify before describing or modifying**

## Coding Standards

### Comment Guidelines
**PROHIBITED**: Insertion of change history comments or modification tracking comments.
- Comments must ONLY describe the current source code content
- Do NOT add comments about what was changed, when it was changed, or why it was changed
- Do NOT add TODO comments about past modifications
- Do NOT add FIXME comments referencing previous versions
- Comments should explain what the code does NOW, not what it used to do

### C Function Naming Convention
**MANDATORY**: All C language functions in this extension MUST use the `kmersearch_` prefix.
- This applies to ALL functions: static, exported, utility, helper functions
- Examples: `kmersearch_init()`, `kmersearch_process_data()`, `kmersearch_validate_input()`
- This naming convention ensures namespace consistency and prevents conflicts with PostgreSQL internals
- NOTE: Structures, enums, and typedefs do NOT require the `kmersearch_` prefix

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

## Extension Architecture

### Data Types
- **DNA2**: Standard DNA sequences with 2-bit encoding (A, C, G, T, U)
- **DNA4**: Full 4-bit encoding (A, C, G, T) with support for IUPAC degenerate bases (M, R, W, S, Y, K, V, H, D, B, N)
- Both types use compact bit-packed storage with variable-length support

### GUC Variables (Configuration Parameters)
- `kmersearch.kmer_size` (default: 16): K-mer size for index operations (4-32)
- `kmersearch.occur_bitlen` (default: 8): Bit length for storing k-mer occurrence information (0-16)
- `kmersearch.max_appearance_rate` (default: 0.5): Maximum appearance rate for high-frequency k-mer detection (0.0-1.0)
- `kmersearch.max_appearance_nrow` (default: 0): Maximum number of rows for high-frequency k-mer detection (0 = unlimited)
- `kmersearch.min_score` (default: 1): Minimum score threshold for search operations
- `kmersearch.min_shared_kmer_rate` (default: 0.5): Minimum shared k-mer rate for scoring (0.0-1.0)
- `kmersearch.preclude_highfreq_kmer` (default: false): Enable high-frequency k-mer filtering
- `kmersearch.force_use_parallel_highfreq_kmer_cache` (default: false): Force parallel cache usage
- `kmersearch.force_simd_capability` (default: -1): Force specific SIMD capability level (-1 = auto-detect)
- `kmersearch.query_kmer_cache_max_entries` (default: 50000): Maximum cache entries for query patterns
- `kmersearch.actual_min_score_cache_max_entries` (default: 50000): Maximum cache entries for minimum scores
- `kmersearch.highfreq_kmer_cache_load_batch_size` (default: 10000): Batch size for loading high-frequency k-mers
- `kmersearch.highfreq_analysis_batch_size` (default: 10000): Batch size for high-frequency k-mer analysis
- `kmersearch.highfreq_analysis_hashtable_size` (default: 1000000): Initial hash table size for high-frequency k-mer analysis (10000-100000000)

### Core Functions
- **Search Operators**: `=%` operator for k-mer based sequence search
- **Scoring Functions**: 
  - `kmersearch_matchscore()`: Calculate similarity score by counting shared k-mers
- **Length Functions**: `bit_length()`, `nuc_length()`, `char_length()`, `length()`
- **High-frequency Analysis**:
  - `kmersearch_perform_highfreq_analysis()`: Analyze and identify high-frequency k-mers
  - `kmersearch_undo_highfreq_analysis()`: Remove high-frequency k-mer analysis
- **Cache Management**:
  - `kmersearch_query_kmer_cache_stats()`: Query-kmer cache statistics
  - `kmersearch_query_kmer_cache_free()`: Clear query-kmer cache
  - `kmersearch_actual_min_score_cache_stats()`: Minimum score cache statistics
  - `kmersearch_actual_min_score_cache_free()`: Clear minimum score cache
  - `kmersearch_highfreq_kmer_cache_load()`: Load high-frequency k-mers into cache
  - `kmersearch_highfreq_kmer_cache_free()`: Free high-frequency k-mer cache
  - `kmersearch_highfreq_kmer_cache_free_all()`: Free all high-frequency k-mer cache entries
  - `kmersearch_parallel_highfreq_kmer_cache_load()`: Load high-frequency k-mers into parallel dshash cache
  - `kmersearch_parallel_highfreq_kmer_cache_free()`: Free specific entries from parallel cache
  - `kmersearch_parallel_highfreq_kmer_cache_free_all()`: Free all entries from parallel cache

### System Tables
- `kmersearch_highfreq_kmer`: Stores identified high-frequency k-mers
- `kmersearch_highfreq_kmer_meta`: Metadata for high-frequency k-mer analysis
- `kmersearch_gin_index_meta`: GIN index metadata with high-frequency filtering information
- `kmersearch_index_info`: General index tracking information

### Management Views
- `kmersearch_cache_summary`: Overview of all cache statistics and hit rates
- `kmersearch_analysis_status`: Status of high-frequency k-mer analyses

### GIN Index Support
Multiple operator classes with different key storage strategies:
- `kmersearch_dna2_gin_ops_int2`: 16-bit integer keys for DNA2
- `kmersearch_dna2_gin_ops_int4`: 32-bit integer keys for DNA2
- `kmersearch_dna2_gin_ops_int8`: 64-bit integer keys for DNA2
- `kmersearch_dna4_gin_ops_int2`: 16-bit integer keys for DNA4
- `kmersearch_dna4_gin_ops_int4`: 32-bit integer keys for DNA4
- `kmersearch_dna4_gin_ops_int8`: 64-bit integer keys for DNA4

## SIMD Implementation

**SIMD (Single Instruction, Multiple Data) Implementation**:
- Platform-specific optimizations for encoding/decoding operations
- Runtime SIMD capability detection via `kmersearch_simd_capability()`
- Function dispatch tables for selecting optimal implementation
- Supports various instruction sets: SSE, AVX, AVX2, AVX-512, BMI2
- Automatic fallback to standard implementation if SIMD not available
- Optimizes critical path operations like:
  - DNA sequence encoding/decoding
  - K-mer extraction
  - N-gram key generation
  - High-frequency k-mer filtering
- Designed to maximize performance across different CPU architectures
- Careful boundary condition handling to prevent overflow
- Comprehensive runtime checks to ensure safe SIMD usage

Key aspects of SIMD implementation:
- Modular design allowing easy extension to new instruction sets
- Performance-critical sections identified and optimized
- Extensive testing to validate SIMD vs. standard implementation equivalence
- Minimal overhead for SIMD capability detection
- Portable across different PostgreSQL compilation environments

## Memory Tracking

**Historical Implementation Details**:
- Detailed tracking of SIMD memory management strategies
- Continuous optimization of memory allocation patterns
- Memory alignment considerations for SIMD operations
- Efficient memory pool management for SIMD contexts

## Error Recovery and Version Control

**CRITICAL**: When errors or incorrect deletions occur:
- **DO NOT** attempt to fix errors by creating workarounds or patches
- **DO NOT** try to cover up mistakes with compensatory changes
- **IMMEDIATELY** restore from backup or git repository to the previous working state
- **ALWAYS** create backups before making significant deletions or changes
- Use `git status`, `git diff`, and `git log` to track changes
- When a function is accidentally deleted, restore it from git history or backup, do not recreate it from scratch

Example recovery workflow:
```bash
# Check what was changed
git status
git diff

# If accidental deletion occurred, restore from git
git checkout -- <filename>  # Restore single file
git reset --hard HEAD       # Restore entire working directory

# If no git history, restore from backup
cp <backup_file> <original_file>
```

## Git Operations

**CRITICAL**: Git operations are strictly controlled:
- **NEVER** execute git commit without explicit user instruction
- **NEVER** use `git add -A` or `git add .` when temporary files exist in the repository
- **ALWAYS** wait for explicit user permission before any git operations
- **ALWAYS** update KMERSEARCH_BUILD_VERSION in kmersearch.h before committing changes
  - Format: `"1.0.YYYY.MM.DD"` where YYYY.MM.DD is the current date
  - Example: `#define KMERSEARCH_BUILD_VERSION "1.0.2025.08.19"`
  - Use `date +%Y.%m.%d` to get the correct date format
- When preparing to commit, list the specific files to be committed and wait for user approval
- Be aware of temporary files, build artifacts, and test outputs that should not be committed