CREATE EXTENSION pg_kmersearch;

-- Test k-mer search operators and functionality
-- This test covers the =% operator and search functionality

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

-- Create GIN indexes
CREATE INDEX idx_DNA2_gin ON test_DNA2_sequences USING gin (sequence);
CREATE INDEX idx_DNA4_gin ON test_DNA4_sequences USING gin (sequence);

-- Test =% operator with DNA2
SELECT id, name FROM test_DNA2_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;
SELECT id, name FROM test_DNA2_sequences WHERE sequence =% 'GCTAG' ORDER BY id;
SELECT id, name FROM test_DNA2_sequences WHERE sequence =% 'TTTTTTTT' ORDER BY id;

-- Test =% operator with DNA4
SELECT id, name FROM test_DNA4_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;
SELECT id, name FROM test_DNA4_sequences WHERE sequence =% 'GCTAG' ORDER BY id;
SELECT id, name FROM test_DNA4_sequences WHERE sequence =% 'NNATCG' ORDER BY id;

-- Test query length validation (should work with 8+ characters)
SELECT id, name FROM test_DNA2_sequences WHERE sequence =% 'ATCGATCG' ORDER BY id;

-- Test too short query (should error)
\set ON_ERROR_STOP off
SELECT id, name FROM test_DNA2_sequences WHERE sequence =% 'ATCG';
\set ON_ERROR_STOP on

-- Test case insensitivity in queries
SELECT id, name FROM test_DNA2_sequences WHERE sequence =% 'atcgatcg' ORDER BY id;
SELECT id, name FROM test_DNA4_sequences WHERE sequence =% 'atcgatcg' ORDER BY id;

DROP EXTENSION pg_kmersearch CASCADE;