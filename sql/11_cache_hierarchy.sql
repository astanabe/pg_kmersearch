SET client_min_messages = WARNING;
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test comprehensive cache hierarchy and GUC validation functionality
-- This test covers the enhanced kmersearch_is_kmer_highfreq implementation

-- Set GUC variables for index creation (will be changed later for testing)
SET kmersearch.kmer_size = 4;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.max_appearance_rate = 0.2;
SET kmersearch.max_appearance_nrow = 3;

-- Clean up any existing test data
DROP TABLE IF EXISTS test_cache_hierarchy CASCADE;
DELETE FROM kmersearch_highfreq_kmer WHERE detection_reason LIKE 'cache_hierarchy_%';
DELETE FROM kmersearch_highfreq_kmer_meta WHERE column_name = 'test_seq';

-- Create test table for cache hierarchy testing
CREATE TABLE test_cache_hierarchy (
    id SERIAL PRIMARY KEY,
    test_seq dna2
);

-- Insert test data with repeated patterns to ensure high-frequency k-mers
INSERT INTO test_cache_hierarchy (test_seq) VALUES 
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),
    ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2),
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),
    ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2),
    ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'::dna2),
    ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'::dna2),
    ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'::dna2),
    ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'::dna2);

-- Create GIN index (using kmer_size=4)
-- With k=4 and occur_bitlen=8, total bits = 4*2+8 = 16, so we can use int2 operator class
CREATE INDEX test_cache_hierarchy_gin_idx ON test_cache_hierarchy USING gin (test_seq kmersearch_dna2_gin_ops_int2);

-- Test Phase 1: No metadata scenario (should work with table fallback)
SELECT 'Phase 1: Testing without metadata table...' as test_phase;
SELECT test_seq =% 'ATCGATCG' as no_metadata_query FROM test_cache_hierarchy LIMIT 1;

-- Test Phase 2: Perform high-frequency k-mer analysis
SELECT 'Phase 2: Performing high-frequency k-mer analysis...' as test_phase;

-- Perform actual high-frequency k-mer analysis
SELECT kmersearch_perform_highfreq_analysis(
    'test_cache_hierarchy'::text,
    'test_seq'::text
) AS analysis_result;

-- Verify analysis results
SELECT 'Verifying analysis results:' AS status;
SELECT highfreq_kmer_count, max_appearance_rate, max_appearance_nrow 
FROM kmersearch_analysis_status 
WHERE table_name = 'test_cache_hierarchy' AND column_name = 'test_seq';

-- Test 2.1: Matching GUC settings (should work)
SET kmersearch.kmer_size = 4;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.max_appearance_rate = 0.2;
SET kmersearch.max_appearance_nrow = 3;
SET kmersearch.min_shared_ngram_key_rate = 0.2;  -- Allow matches with 20% shared k-mers
SELECT test_seq =% 'ATCGATCG' as matching_guc_query FROM test_cache_hierarchy LIMIT 1;

-- Test 2.2: Mismatched GUC - occur_bitlen (should fail during cache load)
SET kmersearch.occur_bitlen = 16;
SELECT kmersearch_highfreq_kmer_cache_load(
  'test_cache_hierarchy',
  'test_seq'
) as mismatched_occur_bitlen;

-- Test 2.3: Mismatched GUC - max_appearance_rate (should fail during cache load)
SET kmersearch.kmer_size = 4;  -- Reset
SET kmersearch.occur_bitlen = 8;  -- Reset
SET kmersearch.max_appearance_rate = 0.8;
SELECT kmersearch_highfreq_kmer_cache_load(
  'test_cache_hierarchy',
  'test_seq'
) as mismatched_rate;

-- Test 2.4: Mismatched GUC - max_appearance_nrow (should fail during cache load)
SET kmersearch.kmer_size = 4;  -- Reset
SET kmersearch.max_appearance_rate = 0.4;  -- Reset
SET kmersearch.max_appearance_nrow = 10;
SELECT kmersearch_highfreq_kmer_cache_load(
  'test_cache_hierarchy',
  'test_seq'
) as mismatched_nrow;

-- Reset to correct values for subsequent tests
SET kmersearch.kmer_size = 4;
SET kmersearch.max_appearance_nrow = 1;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.max_appearance_rate = 0.01;

-- Test Phase 3: Cache hierarchy testing
SELECT 'Phase 3: Testing cache hierarchy...' as test_phase;

-- High-frequency k-mer data should already exist from the analysis above

-- Test 3.1: Global cache loading and usage
SELECT 'Testing global cache...' as test_substep;
SELECT kmersearch_highfreq_kmer_cache_load(
  'test_cache_hierarchy',
  'test_seq'
) as global_cache_load_result;

-- Query with global cache loaded
SELECT test_seq =% 'ATCGATCG' as global_cache_query FROM test_cache_hierarchy LIMIT 1;

-- Test 3.2: Parallel cache loading and usage
SELECT 'Testing parallel cache...' as test_substep;
SELECT kmersearch_parallel_highfreq_kmer_cache_load(
  'test_cache_hierarchy',
  'test_seq'
) as parallel_cache_load_result;

-- Test 3.3: Cache hierarchy fallback sequence
SELECT 'Testing cache hierarchy fallback...' as test_substep;

-- Clear global cache first
SELECT kmersearch_highfreq_kmer_cache_free_all() as global_cache_cleared;

-- Test with parallel cache only (should work)
SET kmersearch.force_use_parallel_highfreq_kmer_cache = true;
SELECT test_seq =% 'ATCGATCG' as parallel_only_query FROM test_cache_hierarchy LIMIT 1;

-- Clear parallel cache too
SELECT kmersearch_parallel_highfreq_kmer_cache_free_all() as parallel_cache_cleared;

-- Test with no cache (should fall back to table lookup)
SET kmersearch.force_use_parallel_highfreq_kmer_cache = false;
SELECT test_seq =% 'ATCGATCG' as table_fallback_query FROM test_cache_hierarchy LIMIT 1;

-- Test Phase 4: Table lookup without high-frequency data
SELECT 'Phase 4: Testing table lookup without high-frequency data...' as test_phase;

-- Remove all high-frequency k-mer data
DELETE FROM kmersearch_highfreq_kmer WHERE table_oid = (SELECT oid FROM pg_class WHERE relname = 'test_cache_hierarchy') AND column_name = 'test_seq';

-- Test query (should work but find no high-frequency k-mers)
SELECT test_seq =% 'ATCGATCG' as no_highfreq_data_query FROM test_cache_hierarchy LIMIT 1;

-- Test Phase 5: Non-existent table scenario
SELECT 'Phase 5: Testing non-existent table scenario...' as test_phase;

-- Remove metadata to simulate non-existent table scenario
DELETE FROM kmersearch_highfreq_kmer_meta WHERE column_name = 'test_seq';

-- Test query (should work with no validation)
SELECT test_seq =% 'ATCGATCG' as no_table_query FROM test_cache_hierarchy LIMIT 1;

-- Clean up test data
SELECT 'Cleaning up test data...' as cleanup_phase;
DROP TABLE IF EXISTS test_cache_hierarchy CASCADE;
DELETE FROM kmersearch_highfreq_kmer WHERE table_oid = (SELECT oid FROM pg_class WHERE relname = 'test_cache_hierarchy') AND column_name = 'test_seq';
DELETE FROM kmersearch_highfreq_kmer_meta WHERE column_name = 'test_seq';

SELECT 'Cache hierarchy test completed successfully!' as final_result;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;