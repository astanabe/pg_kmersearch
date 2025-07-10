CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test scoring functions
-- This test covers kmersearch_rawscore and kmersearch_correctedscore functions

-- Clean up any existing tables
DROP TABLE IF EXISTS test_DNA2_sequences, test_DNA4_sequences CASCADE;

-- Create test tables for this test
CREATE TABLE test_DNA2_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

CREATE TABLE test_DNA4_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA4
);

-- Insert test data
INSERT INTO test_DNA2_sequences (name, sequence) VALUES
    ('seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('seq2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'),
    ('seq3', 'ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT');

INSERT INTO test_DNA4_sequences (name, sequence) VALUES
    ('seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('seq2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'),
    ('seq3', 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC');

-- Test rawscore function with DNA2
SELECT id, name, 
       kmersearch_rawscore(sequence, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM test_DNA2_sequences ORDER BY id;

-- Test rawscore function with DNA4
SELECT id, name,
       kmersearch_rawscore(sequence, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM test_DNA4_sequences ORDER BY id;

-- Test correctedscore function with DNA2
SELECT id, name,
       kmersearch_correctedscore(sequence, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS correctedscore
FROM test_DNA2_sequences ORDER BY id;

-- Test correctedscore function with DNA4
SELECT id, name,
       kmersearch_correctedscore(sequence, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS correctedscore
FROM test_DNA4_sequences ORDER BY id;

-- Test with shorter sequences for exact match
SELECT kmersearch_rawscore('ATCGATCG'::DNA2, 'ATCGATCG') AS exact_match_score;
SELECT kmersearch_rawscore('ATCGATCG'::DNA4, 'ATCGATCG') AS exact_match_score;

-- Test with completely different sequences
SELECT kmersearch_rawscore('ATCGATCG'::DNA2, 'GCTAGCTA') AS no_match_score;
SELECT kmersearch_rawscore('ATCGATCG'::DNA4, 'GCTAGCTA') AS no_match_score;

-- Test with partial matches
SELECT kmersearch_rawscore('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA2, 'ATCGATCG') AS partial_match;
SELECT kmersearch_rawscore('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA4, 'ATCGATCG') AS partial_match;

DROP EXTENSION pg_kmersearch CASCADE;