SET client_min_messages = WARNING;
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test advanced search functionality and configuration
-- This test covers different k-mer sizes and configuration options

-- Test with different k-mer sizes
SET kmersearch.kmer_size = 6;

-- Create a new index with k=6
CREATE TABLE test_k6_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

INSERT INTO test_k6_sequences (name, sequence) VALUES
    ('k6_test1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('k6_test2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA');

CREATE INDEX idx_k6_gin ON test_k6_sequences USING gin (sequence);

-- Test search with k=6
SELECT id, name FROM test_k6_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;

-- Test scoring with k=6
SELECT id, name,
       kmersearch_rawscore(sequence, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM test_k6_sequences ORDER BY id;

-- Test with k=6
SET kmersearch.kmer_size = 6;

CREATE TABLE test_k6_sequences_2 (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

INSERT INTO test_k6_sequences_2 (name, sequence) VALUES
    ('k6_test1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('k6_test2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA');

CREATE INDEX idx_k6_gin_2 ON test_k6_sequences_2 USING gin (sequence);

-- Test search with k=6
SELECT id, name FROM test_k6_sequences_2 WHERE sequence =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA' ORDER BY id;

-- Test complete search workflow with ORDER BY score
SELECT id, name, 
       kmersearch_rawscore(sequence, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score
FROM test_k6_sequences_2 
WHERE sequence =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY score DESC, id;

-- Reset k-mer size for other tests
SET kmersearch.kmer_size = 4;

-- Clean up test tables
DROP TABLE IF EXISTS test_k6_sequences CASCADE;
DROP TABLE IF EXISTS test_k6_sequences_2 CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;