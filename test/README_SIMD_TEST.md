# SIMD Test Scripts

This directory contains test scripts for verifying SIMD optimizations in pg_kmersearch.

## Scripts

### test_simd_amd64.pl
Tests SIMD capabilities on AMD64/x86-64 architecture in the following order:
- AVX512VBMI2 (6)
- AVX512VBMI (5)
- AVX512BW (4)
- AVX512F (3)
- BMI2 (2)
- AVX2 (1)
- None (0)

### test_simd_arm64.pl
Tests SIMD capabilities on ARM64/AArch64 architecture in the following order:
- SVE2 (23)
- SVE+NEON (22)
- NEON (21)
- None (0)

## Prerequisites

1. Perl with DBI and DBD::Pg modules:
   ```bash
   sudo apt-get install libdbi-perl libdbd-pg-perl
   # or
   cpan install DBI DBD::Pg
   ```

2. PostgreSQL with pg_kmersearch extension installed

3. Database connection environment variables (optional):
   - `PGDATABASE`: Database name (default: postgres)
   - `PGHOST`: Database host (default: localhost)
   - `PGPORT`: Database port (default: 5432)
   - `PGUSER`: Database user (default: current user)

## Usage

### AMD64 Systems
```bash
./test_simd_amd64.pl
```

### ARM64 Systems
```bash
./test_simd_arm64.pl
```

## What the Scripts Test

1. **SIMD Capability Switching**: Verifies that `kmersearch.force_simd_capability` correctly changes the active SIMD implementation

2. **DNA2 Encoding/Decoding**: Tests 2-bit DNA sequence encoding and decoding with various SIMD levels

3. **DNA4 Encoding/Decoding**: Tests 4-bit DNA sequence encoding (including degenerate bases) with various SIMD levels

4. **Data Integrity**: Ensures that sequences are correctly stored and retrieved regardless of SIMD level

5. **Performance**: Measures encoding/decoding throughput at each SIMD level

6. **K-mer Operations**: Tests GIN index creation and k-mer searches at different SIMD levels

## Expected Output

The scripts will:
1. Display the auto-detected SIMD capability
2. For each SIMD level:
   - Set the forced capability
   - Verify the setting took effect
   - Test DNA2 operations and show throughput
   - Test DNA4 operations and show throughput
   - Test k-mer index operations
3. Reset to auto-detected capability
4. Clean up test data

## Troubleshooting

1. **Permission Denied**: Make sure the scripts are executable:
   ```bash
   chmod +x test_simd_*.pl
   ```

2. **Database Connection Failed**: Check your PostgreSQL connection settings and ensure the extension is installed:
   ```sql
   CREATE EXTENSION pg_kmersearch;
   ```

3. **Invalid SIMD Capability Error**: The script may be trying to set a capability not supported by your CPU. This is expected behavior - the test should continue with supported levels.

4. **Performance Variations**: Performance differences between SIMD levels depend on:
   - Sequence length (longer sequences benefit more from SIMD)
   - CPU architecture and specific instruction support
   - System load and other factors

## Interpreting Results

- Higher SIMD levels should generally show better performance for longer sequences
- Short sequences may not show significant differences
- The "None" level provides a baseline for comparison
- Data integrity should be maintained at all SIMD levels