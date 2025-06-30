CREATE EXTENSION pg_kmersearch;

-- Test k-mer search operators and functionality
-- This test covers the =% operator and search functionality

-- Clean up any existing tables
DROP TABLE IF EXISTS test_dna2_sequences, test_dna4_sequences CASCADE;

-- Create test tables for this test
CREATE TABLE test_dna2_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence dna2
);

CREATE TABLE test_dna4_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence dna4
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
CREATE INDEX idx_dna2_gin ON test_dna2_sequences USING gin (sequence);
CREATE INDEX idx_dna4_gin ON test_dna4_sequences USING gin (sequence);

-- Test =% operator with DNA2
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'GCTAG' ORDER BY id;
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'TTTTTTTT' ORDER BY id;

-- Test =% operator with DNA4
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'GCTAG' ORDER BY id;
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'NNATCG' ORDER BY id;

-- Test query length validation (should work with 8+ characters)
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;

-- Test too short query (should error)
\set ON_ERROR_STOP off
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'ATCG';
\set ON_ERROR_STOP on

-- Test case insensitivity in queries
SELECT id, name FROM test_dna2_sequences WHERE sequence =% 'atcgatcg' ORDER BY id;
SELECT id, name FROM test_dna4_sequences WHERE sequence =% 'atcgatcg' ORDER BY id;

DROP EXTENSION pg_kmersearch CASCADE;