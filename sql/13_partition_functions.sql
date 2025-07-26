-- Test partitioning support functions

SET client_min_messages = WARNING;
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

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
SELECT kmersearch_partition_table('test_sequences', 4);
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

-- Test 8: Create GIN indexes on all partitions
SET kmersearch.kmer_size = 4;
SET kmersearch.preclude_highfreq_kmer = false;
SET client_min_messages = NOTICE;
SELECT * FROM kmersearch_parallel_create_index('test_sequences', 'sequence');
SET client_min_messages = WARNING;

-- Test 9: Verify indexes were created
SELECT i.indexrelid::regclass as index_name, 
       c.relname as table_name
FROM pg_index i
JOIN pg_class c ON i.indrelid = c.oid
WHERE c.relname LIKE 'test_sequences_%'
  AND i.indrelid IN (
    SELECT inhrelid FROM pg_inherits WHERE inhparent = 'test_sequences'::regclass
  )
ORDER BY c.relname;

-- Test 10: Test search functionality on partitioned table
SELECT COUNT(*) FROM test_sequences WHERE sequence =% 'ATCGATCG';

-- Test 11: Create test table with DNA4 column
CREATE TABLE test_sequences_dna4 (
    id serial PRIMARY KEY,
    sequence dna4 NOT NULL
);

INSERT INTO test_sequences_dna4 (sequence) VALUES
    ('ATCGATCGATCGATCG'),
    ('GCTAGCTAGCTAGCTA'),
    ('MMMMMMMMMMMMMMMM'),
    ('NNNNNNNNNNNNNNNN');

-- Test 12: Convert DNA4 table to partitioned
SET client_min_messages = NOTICE;
SELECT kmersearch_partition_table('test_sequences_dna4', 2);
SET client_min_messages = WARNING;

-- Test 13: Create indexes on DNA4 partitioned table
SET client_min_messages = NOTICE;
SELECT * FROM kmersearch_parallel_create_index('test_sequences_dna4', 'sequence');
SET client_min_messages = WARNING;

-- Test 14: Error case - non-partitioned table
CREATE TABLE test_regular_table (
    id serial PRIMARY KEY,
    sequence dna2
);

SELECT * FROM kmersearch_parallel_create_index('test_regular_table', 'sequence');

-- Test 15: Error case - table without DNA column
CREATE TABLE test_no_dna_table (
    id serial PRIMARY KEY,
    name text
);

SELECT kmersearch_partition_table('test_no_dna_table', 2);

-- Test 16: Error case - table with multiple DNA columns
CREATE TABLE test_multi_dna_table (
    id serial PRIMARY KEY,
    seq1 dna2,
    seq2 dna4
);

SELECT kmersearch_partition_table('test_multi_dna_table', 2);

-- Test 17: Test with high-frequency k-mer exclusion
SET kmersearch.preclude_highfreq_kmer = true;
SET kmersearch.force_use_parallel_highfreq_kmer_cache = false;

-- This should fail because force_use_parallel_highfreq_kmer_cache must be true
SELECT * FROM kmersearch_parallel_create_index('test_sequences', 'sequence');

-- Test 18: Test partition_count validation
SELECT kmersearch_partition_table('test_regular_table', 0);
SELECT kmersearch_partition_table('test_regular_table', -1);

-- Cleanup
DROP TABLE IF EXISTS test_sequences CASCADE;
DROP TABLE IF EXISTS test_sequences_dna4 CASCADE;
DROP TABLE IF EXISTS test_regular_table CASCADE;
DROP TABLE IF EXISTS test_no_dna_table CASCADE;
DROP TABLE IF EXISTS test_multi_dna_table CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;