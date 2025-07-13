CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Set k-mer size to 4 for efficient testing (must be after CREATE EXTENSION)
SET kmersearch.kmer_size = 4;

-- Test management views and new functions
-- This test covers the newly added views and complex return types

-- Test cache summary view
SELECT 'Testing kmersearch_cache_summary view...' as test_phase;
SELECT * FROM kmersearch_cache_summary;

-- Test analysis status view (should be empty initially)
SELECT 'Testing kmersearch_analysis_status view...' as test_phase;
SELECT COUNT(*) as analysis_count FROM kmersearch_analysis_status;

-- Create test table for analysis functions
CREATE TABLE test_analysis_dna2 (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

-- Insert test data with patterns that create high-frequency k-mers
INSERT INTO test_analysis_dna2 (name, sequence) VALUES
    ('seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'),
    ('seq2', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'),
    ('seq3', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG');

-- Create GIN index
CREATE INDEX idx_test_analysis_gin ON test_analysis_dna2 USING gin (sequence);

-- Test kmersearch_perform_highfreq_analysis with complex return type
SELECT 'Testing kmersearch_perform_highfreq_analysis return type...' as test_phase;
SELECT 
    total_rows,
    highfreq_kmers_count,
    parallel_workers_used,
    analysis_duration >= 0 as duration_valid,
    max_appearance_rate_used,
    max_appearance_nrow_used
FROM kmersearch_perform_highfreq_analysis(
    'test_analysis_dna2', 
    'sequence'
);

-- Check if the analysis is reflected in the view
SELECT 'Checking analysis status after analyze...' as test_phase;
SELECT 
    table_name,
    column_name,
    kmer_size,
    occur_bitlen,
    max_appearance_rate,
    max_appearance_nrow,
    highfreq_kmer_count >= 0 as has_kmers
FROM kmersearch_analysis_status
WHERE table_name = 'test_analysis_dna2';

-- Test kmersearch_undo_highfreq_analysis with complex return type
SELECT 'Testing kmersearch_undo_highfreq_analysis return type...' as test_phase;
SELECT 
    dropped_analyses,
    dropped_highfreq_kmers,
    freed_storage_bytes >= 0 as storage_freed
FROM kmersearch_undo_highfreq_analysis(
    'test_analysis_dna2', 
    'sequence'
);

-- Verify analysis was dropped
SELECT 'Checking analysis status after drop...' as test_phase;
SELECT COUNT(*) as remaining_analyses 
FROM kmersearch_analysis_status
WHERE table_name = 'test_analysis_dna2';

-- Test dropping non-existent analysis
SELECT 'Testing drop of non-existent analysis...' as test_phase;
SELECT 
    dropped_analyses,
    dropped_highfreq_kmers,
    freed_storage_bytes
FROM kmersearch_undo_highfreq_analysis(
    'test_analysis_dna2', 
    'sequence'
);

-- Test the cache summary view again after operations
SELECT 'Testing cache summary after operations...' as test_phase;
SELECT 
    cache_type,
    total_entries >= 0 as has_entries,
    total_hits >= 0 as has_hits,
    total_misses >= 0 as has_misses,
    hit_rate >= 0 AND hit_rate <= 1 as hit_rate_valid
FROM kmersearch_cache_summary
ORDER BY cache_type;

-- Clean up
DROP TABLE test_analysis_dna2 CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;