CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test parallel high-frequency k-mer cache functionality
-- This test verifies dshash-based parallel cache operations

-- Create test table with DNA sequences
CREATE TABLE test_dna_parallel (
    id SERIAL PRIMARY KEY,
    seq dna2
);

-- Insert test data with repeated k-mers to create high-frequency patterns
INSERT INTO test_dna_parallel (seq) VALUES 
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),  -- Contains repeated ATCG
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),   -- Contains repeated GCTA  
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),  -- Duplicate
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),   -- Duplicate
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),  -- More duplicates to ensure high frequency
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),
    ('TTTTAAAACCCCGGGGTTTTAAAACCCCGGGG'::dna2),  -- Different pattern
    ('GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT'::dna2);   -- Different pattern

-- Create GIN index to trigger high-frequency k-mer detection
CREATE INDEX test_dna_parallel_gin_idx ON test_dna_parallel USING gin (seq);

-- Insert test metadata for cache loading
INSERT INTO kmersearch_highfreq_kmers_meta (table_oid, column_name, k_value, occur_bitlen, max_appearance_rate, max_appearance_nrow)
VALUES ((SELECT oid FROM pg_class WHERE relname = 'test_dna_parallel'), 'seq', 8, 8, 0.05, 0);

-- Load high-frequency k-mers into global cache (for comparison)
SELECT kmersearch_highfreq_kmers_cache_load('test_dna_parallel'::regclass, 'seq', 8);

-- Show current GUC settings
SHOW kmersearch.preclude_highfreq_kmer;
SHOW kmersearch.force_use_dshash;

-- Test 1: Load parallel cache with default settings
SELECT kmersearch_parallel_highfreq_kmers_cache_load('test_dna_parallel'::regclass, 'seq', 8) AS parallel_cache_loaded;

-- Test 2: Enable force_use_dshash and test extraction
SET kmersearch.force_use_dshash = true;
SHOW kmersearch.force_use_dshash;

-- Test 3: Enable high-frequency k-mer exclusion
SET kmersearch.preclude_highfreq_kmer = true;
SHOW kmersearch.preclude_highfreq_kmer;

-- Test 4: Search with parallel cache (should use dshash)
-- This should use parallel cache due to force_use_dshash = true
SELECT COUNT(*) AS results_with_parallel_cache
FROM test_dna_parallel 
WHERE seq =% 'ATCGATCG';

-- Test 5: Compare with global cache behavior
SET kmersearch.force_use_dshash = false;
SHOW kmersearch.force_use_dshash;

SELECT COUNT(*) AS results_with_global_cache
FROM test_dna_parallel 
WHERE seq =% 'ATCGATCG';

-- Test 6: Test cache cleanup
SELECT kmersearch_parallel_highfreq_kmers_cache_free() AS parallel_cache_freed_entries;

-- Test 7: Test error handling when parallel cache is not loaded
SET kmersearch.force_use_dshash = true;
-- This should cause an error since parallel cache was freed
SELECT COUNT(*) AS error_test
FROM test_dna_parallel 
WHERE seq =% 'ATCGATCG'
LIMIT 1;

-- Reset settings
SET kmersearch.preclude_highfreq_kmer = false;
SET kmersearch.force_use_dshash = false;

-- Cleanup
DROP TABLE test_dna_parallel;
SELECT kmersearch_highfreq_kmers_cache_free();

DROP EXTENSION pg_kmersearch CASCADE;