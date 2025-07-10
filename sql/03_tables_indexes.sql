CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test table creation and GIN index functionality
-- This test covers DDL operations and index creation

-- Create test tables
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

-- Test data retrieval
SELECT id, name, sequence FROM test_DNA2_sequences ORDER BY id;
SELECT id, name, sequence FROM test_DNA4_sequences ORDER BY id;

-- Create GIN indexes
CREATE INDEX idx_DNA2_gin ON test_DNA2_sequences USING gin (sequence);
CREATE INDEX idx_DNA4_gin ON test_DNA4_sequences USING gin (sequence);

-- Verify indexes exist
SELECT indexname, tablename FROM pg_indexes 
WHERE indexname IN ('idx_DNA2_gin', 'idx_DNA4_gin')
ORDER BY indexname;

DROP EXTENSION pg_kmersearch CASCADE;