-- Test high-frequency k-mer filtering functionality
-- This test verifies the high-frequency k-mer cache system and filtering

-- Create test table
CREATE TABLE test_highfreq_filter (
    id SERIAL PRIMARY KEY,
    dna2_seq DNA2,
    dna4_seq DNA4,
    description TEXT
);

-- Insert test data with some high-frequency patterns
INSERT INTO test_highfreq_filter (dna2_seq, dna4_seq, description) VALUES
    ('ATCGATCGATCGATCG'::DNA2, 'ATCGATCGATCGATCG'::DNA4, 'Repetitive ATCG pattern'),
    ('GGGGGGGGGGGGGGGG'::DNA2, 'GGGGGGGGGGGGGGGG'::DNA4, 'High-frequency G pattern'),
    ('AAAAAAAAAAAAAAAA'::DNA2, 'AAAAAAAAAAAAAAAA'::DNA4, 'High-frequency A pattern'),
    ('ACGTACGTACGTACGT'::DNA2, 'ACGTACGTACGTACGT'::DNA4, 'Repetitive ACGT pattern'),
    ('TGCATGCATGCATGCA'::DNA2, 'TGCATGCATGCATGCA'::DNA4, 'Repetitive TGCA pattern'),
    ('CCCCCCCCCCCCCCC'::DNA2, 'CCCCCCCCCCCCCCCC'::DNA4, 'High-frequency C pattern'),
    ('TTTTTTTTTTTTTTTT'::DNA2, 'TTTTTTTTTTTTTTTT'::DNA4, 'High-frequency T pattern'),
    ('ATGCATGCATGCATGC'::DNA2, 'ATGCATGCATGCATGC'::DNA4, 'Repetitive ATGC pattern');

-- Test metadata tables exist
SELECT COUNT(*) as meta_table_count FROM information_schema.tables 
WHERE table_name IN ('kmersearch_highfreq_kmers_meta', 'kmersearch_gin_index_meta');

-- Test cache management functions exist
SELECT proname FROM pg_proc WHERE proname LIKE '%cache%' AND proname LIKE 'kmersearch%';

-- Create GIN indexes to test filtering
CREATE INDEX idx_test_dna2_gin ON test_highfreq_filter USING gin (dna2_seq);
CREATE INDEX idx_test_dna4_gin ON test_highfreq_filter USING gin (dna4_seq);

-- Verify indexes were created
SELECT indexname, tablename FROM pg_indexes WHERE tablename = 'test_highfreq_filter';

-- Test basic search functionality before high-frequency analysis
SELECT id, description FROM test_highfreq_filter WHERE dna2_seq =% 'ATCG';
SELECT id, description FROM test_highfreq_filter WHERE dna4_seq =% 'GGGG';

-- Simulate high-frequency k-mer analysis (manual metadata insertion for testing)
-- Insert sample metadata to test cache functionality
INSERT INTO kmersearch_highfreq_kmers_meta (table_oid, column_name, k_value, max_appearance_rate, max_appearance_nrow)
SELECT 
    'test_highfreq_filter'::regclass::oid,
    'dna2_seq',
    8,
    0.05,
    1000
WHERE NOT EXISTS (
    SELECT 1 FROM kmersearch_highfreq_kmers_meta 
    WHERE table_oid = 'test_highfreq_filter'::regclass::oid 
    AND column_name = 'dna2_seq' 
    AND k_value = 8
);

INSERT INTO kmersearch_highfreq_kmers_meta (table_oid, column_name, k_value, max_appearance_rate, max_appearance_nrow)
SELECT 
    'test_highfreq_filter'::regclass::oid,
    'dna4_seq',
    8,
    0.05,
    1000
WHERE NOT EXISTS (
    SELECT 1 FROM kmersearch_highfreq_kmers_meta 
    WHERE table_oid = 'test_highfreq_filter'::regclass::oid 
    AND column_name = 'dna4_seq' 
    AND k_value = 8
);

-- Verify metadata was inserted
SELECT table_oid::regclass, column_name, k_value, max_appearance_rate 
FROM kmersearch_highfreq_kmers_meta 
WHERE table_oid = 'test_highfreq_filter'::regclass::oid;

-- Test search after metadata insertion (should trigger cache loading)
SELECT id, description FROM test_highfreq_filter WHERE dna2_seq =% 'ATCG' LIMIT 3;
SELECT id, description FROM test_highfreq_filter WHERE dna4_seq =% 'GGGG' LIMIT 3;

-- Test with different k-mer patterns
SELECT id, description FROM test_highfreq_filter WHERE dna2_seq =% 'ACGT';
SELECT id, description FROM test_highfreq_filter WHERE dna4_seq =% 'TGCA';

-- Test scoring functions
SELECT id, description, kmersearch_rawscore(dna2_seq, 'ATCGATCG') as raw_score 
FROM test_highfreq_filter 
WHERE dna2_seq =% 'ATCG'
ORDER BY raw_score DESC;

SELECT id, description, kmersearch_rawscore(dna4_seq, 'GGGGGGGG') as raw_score 
FROM test_highfreq_filter 
WHERE dna4_seq =% 'GGGG'
ORDER BY raw_score DESC;

-- Test edge cases
-- Empty search pattern (should handle gracefully)
SELECT COUNT(*) FROM test_highfreq_filter WHERE dna2_seq =% '';

-- Very short pattern
SELECT COUNT(*) FROM test_highfreq_filter WHERE dna2_seq =% 'A';

-- Very long pattern
SELECT COUNT(*) FROM test_highfreq_filter WHERE dna2_seq =% 'ATCGATCGATCGATCGATCGATCG';

-- Cleanup test data
DROP TABLE test_highfreq_filter CASCADE;

-- Clean up metadata (optional - keep for other tests)
-- DELETE FROM kmersearch_highfreq_kmers_meta WHERE table_oid = 'test_highfreq_filter'::regclass::oid;