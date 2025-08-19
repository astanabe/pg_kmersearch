-- Test partitioning support functions

SET client_min_messages = WARNING;

-- Limit parallel workers to 2 for consistent test results
SET max_parallel_workers = 2;
SET max_parallel_maintenance_workers = 2;
SET max_parallel_workers_per_gather = 2;

CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Set maintenance_work_mem for batch processing
SET maintenance_work_mem = '64MB';

-- Test 1: Create test table with DNA2 column
CREATE TABLE test_sequences (
    id serial PRIMARY KEY,
    sequence dna2 NOT NULL,
    description text
);

-- Insert test data
INSERT INTO test_sequences (sequence, description) VALUES
    ('ATCGATCGATCGATCG', 'Test sequence 1'),
    ('GCTAGCTAGCTAGCTA', 'Test sequence 2'),
    ('TTTTTTTTTTTTTTTT', 'Test sequence 3'),
    ('AAAAAAAAAAAAAAAA', 'Test sequence 4'),
    ('ATCGATCGATCGATCG', 'Test sequence 5'),
    ('GCTAGCTAGCTAGCTA', 'Test sequence 6'),
    ('CCCCCCCCCCCCCCCC', 'Test sequence 7'),
    ('GGGGGGGGGGGGGGGG', 'Test sequence 8');

-- Test 2: Verify table is not partitioned
SELECT c.relname, c.relkind 
FROM pg_class c 
WHERE c.oid = 'test_sequences'::regclass;

-- Test 3: Convert table to partitioned table with 4 partitions
SET client_min_messages = NOTICE;
\set VERBOSITY terse
SELECT kmersearch_partition_table('test_sequences', 4);
\set VERBOSITY default
SET client_min_messages = WARNING;

-- Test 4: Verify table is now partitioned
SELECT c.relname, c.relkind 
FROM pg_class c 
WHERE c.oid = 'test_sequences'::regclass;

-- Test 5: List partitions
SELECT c.relname as partition_name
FROM pg_inherits i
JOIN pg_class c ON i.inhrelid = c.oid
WHERE i.inhparent = 'test_sequences'::regclass
ORDER BY c.relname;

-- Test 6: Verify data was migrated correctly
SELECT COUNT(*) as total_rows FROM test_sequences;
SELECT sequence, description FROM test_sequences ORDER BY id;

-- Test 7: Error case - try to partition already partitioned table
SELECT kmersearch_partition_table('test_sequences', 4);

-- Test 8: Test search functionality on partitioned table
SET kmersearch.kmer_size = 4;
SET kmersearch.preclude_highfreq_kmer = false;
SELECT COUNT(*) FROM test_sequences WHERE sequence =% 'ATCGATCG';

-- Test 9: Create test table with DNA4 column
CREATE TABLE test_sequences_dna4 (
    id serial PRIMARY KEY,
    sequence dna4 NOT NULL
);

INSERT INTO test_sequences_dna4 (sequence) VALUES
    ('ATCGATCGATCGATCG'),
    ('GCTAGCTAGCTAGCTA'),
    ('MMMMMMMMMMMMMMMM'),
    ('NNNNNNNNNNNNNNNN');

-- Test 10: Convert DNA4 table to partitioned
SET client_min_messages = NOTICE;
\set VERBOSITY terse
SELECT kmersearch_partition_table('test_sequences_dna4', 2);
\set VERBOSITY default
SET client_min_messages = WARNING;

-- Test 11: Error case - table without DNA column
CREATE TABLE test_no_dna_table (
    id serial PRIMARY KEY,
    name text
);

SELECT kmersearch_partition_table('test_no_dna_table', 2);

-- Test 12: Error case - table with multiple DNA columns
CREATE TABLE test_multi_dna_table (
    id serial PRIMARY KEY,
    seq1 dna2,
    seq2 dna4
);

SELECT kmersearch_partition_table('test_multi_dna_table', 2);

-- Test 13: Test partition_count validation
CREATE TABLE test_regular_table (
    id serial PRIMARY KEY,
    sequence dna2
);

SELECT kmersearch_partition_table('test_regular_table', 0);
SELECT kmersearch_partition_table('test_regular_table', -1);

-- Test 14: Test tablespace parameter (test with default behavior)
CREATE TABLE test_tablespace_table (
    id serial PRIMARY KEY,
    sequence dna2 NOT NULL
);

INSERT INTO test_tablespace_table (sequence) VALUES
    ('ATCGATCGATCGATCG'),
    ('GCTAGCTAGCTAGCTA');

-- Convert to partitioned table without specifying tablespace (uses default)
SET client_min_messages = NOTICE;
\set VERBOSITY terse
SELECT kmersearch_partition_table('test_tablespace_table', 2);
\set VERBOSITY default
SET client_min_messages = WARNING;

-- Verify partitions were created
SELECT c.relname as partition_name
FROM pg_inherits i
JOIN pg_class c ON i.inhrelid = c.oid
WHERE i.inhparent = 'test_tablespace_table'::regclass
ORDER BY c.relname;

-- Test 15: Test with NULL tablespace parameter (should work same as default)
CREATE TABLE test_null_tablespace (
    id serial PRIMARY KEY,
    sequence dna2 NOT NULL
);

INSERT INTO test_null_tablespace (sequence) VALUES
    ('ATCGATCGATCGATCG'),
    ('GCTAGCTAGCTAGCTA');

-- Convert to partitioned table with NULL tablespace
SELECT kmersearch_partition_table('test_null_tablespace', 2, NULL);

-- Verify it worked
SELECT c.relname as partition_name
FROM pg_inherits i
JOIN pg_class c ON i.inhrelid = c.oid
WHERE i.inhparent = 'test_null_tablespace'::regclass
ORDER BY c.relname;

-- Test 16: Test high-frequency k-mer analysis on partitioned table
-- Create a larger test table
CREATE TABLE test_highfreq_regular (
    id serial PRIMARY KEY,
    sequence dna2 NOT NULL
);

CREATE TABLE test_highfreq_partitioned (
    id serial PRIMARY KEY,
    sequence dna2 NOT NULL
);

-- Insert identical data into both tables
INSERT INTO test_highfreq_regular (sequence)
SELECT 
    CASE (i % 10)
        WHEN 0 THEN 'ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2  -- High frequency pattern
        WHEN 1 THEN 'ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2  -- Same high frequency
        WHEN 2 THEN 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'::dna2
        WHEN 3 THEN 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'::dna2
        WHEN 4 THEN 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'::dna2
        WHEN 5 THEN 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'::dna2
        WHEN 6 THEN 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'::dna2
        WHEN 7 THEN 'ATGCATGCATGCATGCATGCATGCATGCATGC'::dna2
        WHEN 8 THEN 'CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG'::dna2
        ELSE 'TATATATATATATATATATATATATATATATA'::dna2
    END
FROM generate_series(1, 100) i;

-- Copy same data to partitioned table
INSERT INTO test_highfreq_partitioned (sequence)
SELECT sequence FROM test_highfreq_regular ORDER BY id;

-- Convert the second table to partitioned with 4 partitions
SELECT kmersearch_partition_table('test_highfreq_partitioned', 4);

-- Configure analysis parameters
SET kmersearch.kmer_size = 8;
SET kmersearch.max_appearance_rate = 0.15;  -- 15% threshold
SET kmersearch.max_appearance_nrow = 0;     -- Use only rate-based threshold

-- Test 17: Perform high-frequency analysis on regular table
SELECT 
    total_rows,
    highfreq_kmers_count,
    max_appearance_rate_used,
    max_appearance_nrow_used,
    parallel_workers_used
FROM kmersearch_perform_highfreq_analysis('test_highfreq_regular', 'sequence');

-- Test 18: Perform high-frequency analysis on partitioned table
SELECT 
    total_rows,
    highfreq_kmers_count,
    max_appearance_rate_used,
    max_appearance_nrow_used,
    parallel_workers_used
FROM kmersearch_perform_highfreq_analysis('test_highfreq_partitioned', 'sequence');

-- Test 19: Verify that both analyses produced identical results
-- Compare high-frequency k-mers from both tables
SELECT 
    CASE 
        WHEN (
            SELECT COUNT(*) FROM kmersearch_highfreq_kmer 
            WHERE table_oid = 'test_highfreq_regular'::regclass
        ) = (
            SELECT COUNT(*) FROM kmersearch_highfreq_kmer 
            WHERE table_oid = 'test_highfreq_partitioned'::regclass
        )
        THEN 'Same count of high-frequency k-mers'
        ELSE 'Different count of high-frequency k-mers'
    END as comparison_result;

-- Test 20: Clean up analysis results
SELECT * FROM kmersearch_undo_highfreq_analysis('test_highfreq_regular', 'sequence');
SELECT * FROM kmersearch_undo_highfreq_analysis('test_highfreq_partitioned', 'sequence');

-- Test 21: Test with different k-mer size and parallel workers
SET kmersearch.kmer_size = 4;
SET max_parallel_workers_per_gather = 2;

-- Analyze regular table with 2 parallel workers
SELECT 
    total_rows,
    highfreq_kmers_count,
    parallel_workers_used
FROM kmersearch_perform_highfreq_analysis('test_highfreq_regular', 'sequence');

-- Analyze partitioned table with 2 parallel workers
SELECT 
    total_rows,
    highfreq_kmers_count,
    parallel_workers_used
FROM kmersearch_perform_highfreq_analysis('test_highfreq_partitioned', 'sequence');

-- Clean up analysis results
SELECT * FROM kmersearch_undo_highfreq_analysis('test_highfreq_regular', 'sequence');
SELECT * FROM kmersearch_undo_highfreq_analysis('test_highfreq_partitioned', 'sequence');

-- Reset parameters
RESET kmersearch.kmer_size;
RESET kmersearch.max_appearance_rate;
RESET kmersearch.max_appearance_nrow;
RESET max_parallel_workers_per_gather;

-- Cleanup
DROP TABLE IF EXISTS test_sequences CASCADE;
DROP TABLE IF EXISTS test_sequences_dna4 CASCADE;
DROP TABLE IF EXISTS test_regular_table CASCADE;
DROP TABLE IF EXISTS test_no_dna_table CASCADE;
DROP TABLE IF EXISTS test_multi_dna_table CASCADE;
DROP TABLE IF EXISTS test_tablespace_table CASCADE;
DROP TABLE IF EXISTS test_null_tablespace CASCADE;
DROP TABLE IF EXISTS test_highfreq_regular CASCADE;
DROP TABLE IF EXISTS test_highfreq_partitioned CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;