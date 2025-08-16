SET client_min_messages = WARNING;
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Set consistent GUC variables for this test
SET kmersearch.kmer_size = 4;
SHOW kmersearch.kmer_size;
SET kmersearch.occur_bitlen = 8;
SHOW kmersearch.occur_bitlen;
SET kmersearch.max_appearance_rate = 0.5;
SHOW kmersearch.max_appearance_rate;
SET kmersearch.max_appearance_nrow = 0;
SHOW kmersearch.max_appearance_nrow;
SET kmersearch.min_shared_kmer_rate = 0.2;  -- Allow matches with 20% shared k-mers
SHOW kmersearch.min_shared_kmer_rate;

-- Test k-mer search operators and functionality
-- This test covers the =% operator and search functionality

-- Clean up any existing tables
DROP TABLE IF EXISTS test_dna2_sequences, test_dna4_sequences CASCADE;

-- Create test tables for this test
CREATE TABLE test_dna2_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

CREATE TABLE test_dna4_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA4
);

-- Insert test data
INSERT INTO test_dna2_sequences (name, sequence) VALUES
    ('seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('seq2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'),
    ('seq3', 'ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT');

INSERT INTO test_dna4_sequences (name, sequence) VALUES
    ('seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('seq2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'),
    ('seq3', 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC');

-- Create GIN indexes
-- With k=4 and occur_bitlen=8, total bits = 4*2+8 = 16, so we can use int2 operator class
CREATE INDEX idx_DNA2_gin ON test_dna2_sequences USING gin (sequence kmersearch_dna2_gin_ops_int2);
CREATE INDEX idx_DNA4_gin ON test_dna4_sequences USING gin (sequence kmersearch_dna4_gin_ops_int2);

-- Test =% operator with DNA2
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'GCTAG' ORDER BY id;
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'TTTTTTTT' ORDER BY id;

-- Test =% operator with DNA4
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'GCTAG' ORDER BY id;
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'NNATCG' ORDER BY id;

-- Test query length validation (should work with 4+ characters when k=4)
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;

-- Test too short query (should error)
\set ON_ERROR_STOP off
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'ATC';  -- 3 chars, k=4, should error
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'AT';   -- 2 chars, k=4, should error  
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'A';    -- 1 char, k=4, should error
SELECT id, name FROM test_dna2_sequences WHERE sequence =% '';     -- empty, k=4, should error
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'ATC';  -- 3 chars, k=4, should error
\set ON_ERROR_STOP on

-- Test exact k-mer size query (should work)
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'ATCG' ORDER BY id;  -- exactly 4 chars

-- Test case insensitivity in queries
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'atcgatcg' ORDER BY id;
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'atcgatcg' ORDER BY id;

-- Test different query patterns
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'GCTAGCTA' ORDER BY id; -- should match seq2
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'TTTTTTTT' ORDER BY id; -- should match seq3
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'AAAAAAAA' ORDER BY id; -- should not match any

-- Test DNA4 with degenerate bases
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'NNATCG' ORDER BY id;   -- contains N
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'ATCGNN' ORDER BY id;   -- contains N

-- Test scoring functions with different patterns
SELECT id, name, kmersearch_matchscore(sequence, 'ATCGATCG') as score 
FROM test_dna2_sequences ORDER BY score DESC, id;

-- Clean up test tables
DROP TABLE IF EXISTS test_dna2_sequences CASCADE;
DROP TABLE IF EXISTS test_dna4_sequences CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;