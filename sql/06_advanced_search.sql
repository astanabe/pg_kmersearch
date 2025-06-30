CREATE EXTENSION pg_kmersearch;

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

-- Test with k=10
SET kmersearch.kmer_size = 10;

CREATE TABLE test_k10_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

INSERT INTO test_k10_sequences (name, sequence) VALUES
    ('k10_test1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('k10_test2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA');

CREATE INDEX idx_k10_gin ON test_k10_sequences USING gin (sequence);

-- Test search with k=10
SELECT id, name FROM test_k10_sequences WHERE sequence =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA' ORDER BY id;

-- Test complete search workflow with ORDER BY score
SELECT id, name, 
       kmersearch_rawscore(sequence, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score
FROM test_k10_sequences 
WHERE sequence =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY score DESC, id;

-- Reset k-mer size for other tests
SET kmersearch.kmer_size = 8;

DROP EXTENSION pg_kmersearch CASCADE;