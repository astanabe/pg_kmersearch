CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test high-frequency k-mer cache management functionality
-- This test covers manual cache loading and clearing

-- Clean up any existing tables
DROP TABLE IF EXISTS test_highfreq_dna2 CASCADE;

-- Create test table for high-frequency k-mer analysis
CREATE TABLE test_highfreq_dna2 (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

-- Insert test data with 10 sequences of 100 bp each
INSERT INTO test_highfreq_dna2 (name, sequence) VALUES
    ('seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'),
    ('seq2', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'),
    ('seq3', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'),
    ('seq4', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
    ('seq5', 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'),
    ('seq6', 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'),
    ('seq7', 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'),
    ('seq8', 'TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA'),
    ('seq9', 'CAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGT'),
    ('seq10', 'GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC');

-- Create GIN index for high-frequency k-mer analysis
CREATE INDEX idx_test_highfreq_dna2_gin ON test_highfreq_dna2 USING gin (sequence);

-- Test kmersearch_analyze_table functionality
SELECT 'Testing kmersearch_analyze_table...' as test_phase;

-- Analyze the table to identify high-frequency k-mers
-- This will create entries in kmersearch_highfreq_kmer and kmersearch_highfreq_kmer_meta tables
WITH analysis_result AS (
    SELECT kmersearch_analyze_table(
        (SELECT oid FROM pg_class WHERE relname = 'test_highfreq_dna2'), 
        'sequence', 
        8, 
        2  -- parallel_workers (default 0)
    ) AS result
)
SELECT (result).* FROM analysis_result;

-- No need to manually insert data since analyze_table should handle it

-- Check if analysis created metadata
SELECT 'Checking analysis metadata...' as test_phase;
SELECT COUNT(*) as meta_count FROM kmersearch_highfreq_kmer_meta 
WHERE table_oid = (SELECT oid FROM pg_class WHERE relname = 'test_highfreq_dna2');

-- Check if analysis created high-frequency k-mer data
SELECT 'Checking high-frequency k-mers...' as test_phase;
SELECT COUNT(*) as highfreq_kmer_count FROM kmersearch_highfreq_kmer
WHERE index_oid IN (
    SELECT indexrelid FROM pg_stat_user_indexes 
    WHERE schemaname = 'public' AND relname = 'test_highfreq_dna2'
);

-- Test high-frequency k-mer cache functions with analysis data
SELECT 'Testing cache loading with analysis data...' as test_phase;

-- Test high-frequency k-mer cache loading
SELECT 'Testing cache loading...' as test_phase;

-- Set GUC variables to match metadata
SET kmersearch.max_appearance_rate = 0.05;
SET kmersearch.max_appearance_nrow = 0;

SELECT kmersearch_highfreq_kmer_cache_load(test_highfreq_dna2.tableoid, 'sequence', 8) as cache_loaded
FROM test_highfreq_dna2 LIMIT 1;

-- Test high-frequency k-mer cache clearing
SELECT 'Testing cache clearing...' as test_phase;
SELECT kmersearch_highfreq_kmers_cache_free() as freed_entries;

-- Test cache loading with non-existent data
SELECT 'Testing cache loading with invalid parameters...' as test_phase;
SELECT kmersearch_highfreq_kmer_cache_load(0, 'nonexistent', 8) as cache_loaded;

-- Test cache clearing when no cache exists
SELECT 'Testing cache clearing when no cache exists...' as test_phase;
SELECT kmersearch_highfreq_kmers_cache_free() as freed_entries;

-- Test GUC validation and cache hierarchy
SELECT 'Testing GUC validation and cache hierarchy...' as test_phase;

-- Reset GUC settings to match metadata for testing
SET kmersearch.max_appearance_rate = 0.05;
SET kmersearch.max_appearance_nrow = 0;
SET kmersearch.occur_bitlen = 8;

-- Insert some test high-frequency k-mer data for cache hierarchy testing
INSERT INTO kmersearch_highfreq_kmer (index_oid, ngram_key, detection_reason)
VALUES (
  (SELECT indexrelid FROM pg_stat_user_indexes WHERE relname = 'test_highfreq_dna2' AND indexrelname = 'idx_test_highfreq_dna2_gin'),
  '01010101'::bit(16),  -- Sample n-gram key  
  'regression_test'
);

-- Test 1: GUC validation error handling
SELECT 'Testing GUC validation errors...' as test_substep;

-- Change GUC to cause mismatch and test cache loading (should fail)
SET kmersearch.max_appearance_rate = 0.9;  -- Different from metadata (0.05)
SELECT kmersearch_highfreq_kmer_cache_load(
  (SELECT oid FROM pg_class WHERE relname = 'test_highfreq_dna2'),
  'sequence', 8
) as cache_load_with_mismatch;

-- Reset to correct values
SET kmersearch.max_appearance_rate = 0.05;

-- Test 2: Cache hierarchy - Global cache
SELECT 'Testing global cache hierarchy...' as test_substep;
SELECT kmersearch_highfreq_kmer_cache_load(
  (SELECT oid FROM pg_class WHERE relname = 'test_highfreq_dna2'),
  'sequence', 8
) as global_cache_loaded;

-- Test query with global cache
SELECT sequence =% 'ATCGATCG' as global_cache_query FROM test_highfreq_dna2 LIMIT 1;

-- Test 3: Cache hierarchy - Parallel cache  
SELECT 'Testing parallel cache hierarchy...' as test_substep;
SELECT kmersearch_parallel_highfreq_kmer_cache_load(
  (SELECT oid FROM pg_class WHERE relname = 'test_highfreq_dna2'),
  'sequence', 8
) as parallel_cache_loaded;

-- Test 4: Cache hierarchy fallback
SELECT 'Testing cache hierarchy fallback...' as test_substep;

-- Clear global cache
SELECT kmersearch_highfreq_kmers_cache_free() as global_freed;

-- Test with only parallel cache
SET kmersearch.force_use_parallel_highfreq_kmer_cache = true;
SELECT sequence =% 'ATCGATCG' as parallel_only_query FROM test_highfreq_dna2 LIMIT 1;

-- Clear parallel cache  
SELECT kmersearch_parallel_highfreq_kmer_cache_free() as parallel_freed;

-- Test with no cache (table lookup fallback)
SET kmersearch.force_use_parallel_highfreq_kmer_cache = false;
SELECT sequence =% 'ATCGATCG' as table_fallback_query FROM test_highfreq_dna2 LIMIT 1;

-- Clean up test data
DELETE FROM kmersearch_highfreq_kmer WHERE detection_reason = 'regression_test';

-- Clean up analysis data
SELECT 'Cleaning up analysis data...' as test_phase;
DELETE FROM kmersearch_highfreq_kmer
WHERE index_oid IN (
    SELECT indexrelid FROM pg_stat_user_indexes 
    WHERE schemaname = 'public' AND relname = 'test_highfreq_dna2'
);
DELETE FROM kmersearch_highfreq_kmer_meta 
WHERE table_oid = (SELECT oid FROM pg_class WHERE relname = 'test_highfreq_dna2');

-- Clean up
DROP TABLE test_highfreq_dna2 CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;