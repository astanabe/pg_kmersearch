CREATE EXTENSION pg_kmersearch;

-- Test table creation and GIN index functionality
-- This test covers DDL operations and index creation

-- Create test tables
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

-- Test data retrieval
SELECT id, name, sequence FROM test_dna2_sequences ORDER BY id;
SELECT id, name, sequence FROM test_dna4_sequences ORDER BY id;

-- Create GIN indexes
CREATE INDEX idx_dna2_gin ON test_dna2_sequences USING gin (sequence);
CREATE INDEX idx_dna4_gin ON test_dna4_sequences USING gin (sequence);

-- Verify indexes exist
SELECT indexname, tablename FROM pg_indexes 
WHERE indexname IN ('idx_dna2_gin', 'idx_dna4_gin')
ORDER BY indexname;

DROP EXTENSION pg_kmersearch CASCADE;