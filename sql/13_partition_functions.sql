-- Test partitioning support functions

SET client_min_messages = WARNING;
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

-- Cleanup
DROP TABLE IF EXISTS test_sequences CASCADE;
DROP TABLE IF EXISTS test_sequences_dna4 CASCADE;
DROP TABLE IF EXISTS test_regular_table CASCADE;
DROP TABLE IF EXISTS test_no_dna_table CASCADE;
DROP TABLE IF EXISTS test_multi_dna_table CASCADE;
DROP TABLE IF EXISTS test_tablespace_table CASCADE;
DROP TABLE IF EXISTS test_null_tablespace CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;