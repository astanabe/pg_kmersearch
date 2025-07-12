-- Enhanced parallel cache test with actual high-frequency k-mers

CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Create test table with many repeated sequences to ensure high-frequency k-mers
CREATE TABLE test_parallel_enhanced (
    id SERIAL PRIMARY KEY,
    seq dna2
);

-- Insert many sequences with the same k-mers to create high frequency
DO $$
BEGIN
    FOR i IN 1..50 LOOP
        INSERT INTO test_parallel_enhanced (seq) VALUES ('ATCGATCGATCGATCGATCGATCGATCGATCG'::dna2);
        INSERT INTO test_parallel_enhanced (seq) VALUES ('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT'::dna2);
    END LOOP;
END $$;

-- Set parameters to detect high-frequency k-mers
SET kmersearch.kmer_size = 4;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.max_appearance_rate = 0.1;  -- Allow up to 10% appearance rate
SET kmersearch.max_appearance_nrow = 10;   -- k-mers appearing in > 10 rows are high-frequency

-- Perform analysis
SELECT 'Performing k-mer frequency analysis...' AS status;
SELECT kmersearch_analyze_table(
    'test_parallel_enhanced'::text,
    'seq'::text
) AS analysis_result;

-- Check analysis results
SELECT 'Analysis results:' AS status;
SELECT highfreq_kmer_count FROM kmersearch_analysis_status WHERE table_name = 'test_parallel_enhanced';

-- Insert some manual high-frequency k-mer entries to ensure we have data
INSERT INTO kmersearch_highfreq_kmer (table_oid, column_name, ngram_key, detection_reason)
SELECT 
    (SELECT oid FROM pg_class WHERE relname = 'test_parallel_enhanced'),
    'seq'::text,
    substring(('01010101010101010101010101010101'::bit(32))::text::bit varying, 1, 24),  -- Sample 24-bit n-gram key
    'manual_test_entry'
WHERE EXISTS (SELECT 1 FROM pg_stat_user_indexes WHERE relname = 'test_parallel_enhanced');

-- Test Phase 1: Load global cache
SELECT 'Phase 1: Loading global cache...' AS test_phase;
SELECT kmersearch_highfreq_kmer_cache_load(
    'test_parallel_enhanced'::text,
    'seq'::text
) AS global_cache_loaded;

-- Test Phase 2: Load parallel cache
SELECT 'Phase 2: Loading parallel cache...' AS test_phase;
SELECT kmersearch_parallel_highfreq_kmer_cache_load(
    'test_parallel_enhanced'::text,
    'seq'::text
) AS parallel_cache_loaded;

-- Test Phase 3: Functional testing
SELECT 'Phase 3: Testing search functionality...' AS test_phase;

-- Enable high-frequency k-mer exclusion and test different cache modes
SET kmersearch.preclude_highfreq_kmer = true;

-- Test with global cache only
SET kmersearch.force_use_parallel_highfreq_kmer_cache = false;
SELECT 'Testing with global cache:' AS cache_type;
SELECT COUNT(*) AS results FROM test_parallel_enhanced WHERE seq =% 'ATCGATCG';

-- Test with parallel cache (force dshash)
SET kmersearch.force_use_parallel_highfreq_kmer_cache = true;
SELECT 'Testing with parallel cache:' AS cache_type;
SELECT COUNT(*) AS results FROM test_parallel_enhanced WHERE seq =% 'ATCGATCG';

-- Test Phase 4: Cache isolation testing
SELECT 'Phase 4: Testing cache isolation...' AS test_phase;

-- Clear global cache and test parallel cache only
SELECT kmersearch_highfreq_kmer_cache_free_all() AS global_cache_freed;
SELECT 'Testing parallel cache after global cache freed:' AS test_description;
SELECT COUNT(*) AS results FROM test_parallel_enhanced WHERE seq =% 'ATCGATCG';

-- Test Phase 5: Raw score verification
SELECT 'Phase 5: Testing rawscore calculations...' AS test_phase;
SELECT id, 
       kmersearch_rawscore(seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM test_parallel_enhanced 
WHERE seq =% 'ATCGATCG' 
LIMIT 5;

-- Test Phase 6: Cache cleanup and fallback
SELECT 'Phase 6: Testing cache cleanup and fallback...' AS test_phase;

-- Free parallel cache
SELECT kmersearch_parallel_highfreq_kmer_cache_free_all() AS parallel_cache_freed;

-- Test fallback to table lookup
SET kmersearch.force_use_parallel_highfreq_kmer_cache = false;
SELECT 'Testing table lookup fallback:' AS test_description;
SELECT COUNT(*) AS results FROM test_parallel_enhanced WHERE seq =% 'ATCGATCG';

-- Final cleanup
SELECT 'Cleanup: Removing test data...' AS cleanup_phase;
DROP TABLE test_parallel_enhanced CASCADE;
SELECT kmersearch_highfreq_kmer_cache_free_all() AS final_cleanup;

SELECT 'Enhanced parallel cache test completed!' AS final_status;

DROP EXTENSION pg_kmersearch CASCADE;