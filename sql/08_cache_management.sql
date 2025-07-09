CREATE EXTENSION pg_kmersearch;

-- Test cache management functionality
-- This test covers actual min score cache, rawscore cache, and query pattern cache

-- Clean up any existing tables
DROP TABLE IF EXISTS test_cache_DNA2, test_cache_DNA4 CASCADE;

-- Create test tables for cache testing
CREATE TABLE test_cache_DNA2 (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA2
);

CREATE TABLE test_cache_DNA4 (
    id SERIAL PRIMARY KEY,
    name TEXT,
    sequence DNA4
);

-- Insert test data for cache testing
INSERT INTO test_cache_DNA2 (name, sequence) VALUES
    ('cache_seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('cache_seq2', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'),
    ('cache_seq3', 'ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'),
    ('cache_seq4', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'),
    ('cache_seq5', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA');

INSERT INTO test_cache_DNA4 (name, sequence) VALUES
    ('cache_seq1', 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'),
    ('cache_seq2', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'),
    ('cache_seq3', 'ATCGATCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'),
    ('cache_seq4', 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA'),
    ('cache_seq5', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA');

-- Create GIN indexes for testing
CREATE INDEX test_cache_dna2_gin_idx ON test_cache_DNA2 USING gin (sequence);
CREATE INDEX test_cache_dna4_gin_idx ON test_cache_DNA4 USING gin (sequence);

-- Test GUC configuration for cache settings
SHOW kmersearch.actual_min_score_cache_max_entries;
SET kmersearch.actual_min_score_cache_max_entries = 1000;
SHOW kmersearch.actual_min_score_cache_max_entries;

SHOW kmersearch.rawscore_cache_max_entries;
SET kmersearch.rawscore_cache_max_entries = 2000;
SHOW kmersearch.rawscore_cache_max_entries;

SHOW kmersearch.query_pattern_cache_max_entries;
SET kmersearch.query_pattern_cache_max_entries = 3000;
SHOW kmersearch.query_pattern_cache_max_entries;

-- Reset to defaults for consistent testing
SET kmersearch.actual_min_score_cache_max_entries = 50000;
SET kmersearch.rawscore_cache_max_entries = 50000;
SET kmersearch.query_pattern_cache_max_entries = 50000;

-- Test initial cache statistics (should be empty)
SELECT 'Initial actual min score cache:' as test_phase;
SELECT * FROM kmersearch_actual_min_score_cache_stats();

SELECT 'Initial rawscore cache:' as test_phase;
SELECT * FROM kmersearch_rawscore_cache_stats();

SELECT 'Initial query pattern cache:' as test_phase;
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- Execute some queries to populate caches
SELECT 'Populating caches with queries...' as test_phase;

-- Same query multiple times to test actual min score caching
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'ATCGATCG';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'ATCGATCG';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'ATCGATCG';

-- Different queries to test cache diversity
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'TTTTTTTT';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'TTTTTTTT';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'GCTAGCTA';

-- DNA4 queries
SELECT COUNT(*) FROM test_cache_DNA4 WHERE sequence =% 'ATCGATCG';
SELECT COUNT(*) FROM test_cache_DNA4 WHERE sequence =% 'NNNNNNNN';

-- Test scoring functions to populate rawscore cache
SELECT kmersearch_rawscore(sequence, 'ATCGATCG') FROM test_cache_DNA2 WHERE id <= 2;
SELECT kmersearch_rawscore(sequence, 'TTTTTTTT') FROM test_cache_DNA2 WHERE id <= 2;

-- Check cache statistics after usage
SELECT 'After query execution - actual min score cache:' as test_phase;
SELECT * FROM kmersearch_actual_min_score_cache_stats();

SELECT 'After query execution - rawscore cache:' as test_phase;
SELECT * FROM kmersearch_rawscore_cache_stats();

SELECT 'After query execution - query pattern cache:' as test_phase;
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- Test cache clearing functions
SELECT 'Testing cache clear functions...' as test_phase;

SELECT 'Clearing actual min score cache:' as action, kmersearch_actual_min_score_cache_free() as freed_entries;
SELECT 'Clearing rawscore cache:' as action, kmersearch_rawscore_cache_free() as freed_entries;
SELECT 'Clearing query pattern cache:' as action, kmersearch_query_pattern_cache_free() as freed_entries;

-- Verify caches are cleared
SELECT 'After clearing - actual min score cache:' as test_phase;
SELECT * FROM kmersearch_actual_min_score_cache_stats();

SELECT 'After clearing - rawscore cache:' as test_phase;
SELECT * FROM kmersearch_rawscore_cache_stats();

SELECT 'After clearing - query pattern cache:' as test_phase;
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- Test cache behavior with different min_score settings
SELECT 'Testing with different min_score settings...' as test_phase;

SET kmersearch.min_score = 5;
SET kmersearch.min_shared_ngram_key_rate = 0.8;

-- Execute queries with new settings
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'ATCGATCG';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'ATCGATCG';

-- Check cache after configuration change
SELECT 'With changed settings - actual min score cache:' as test_phase;
SELECT * FROM kmersearch_actual_min_score_cache_stats();

-- Reset to defaults
SET kmersearch.min_score = 1;
SET kmersearch.min_shared_ngram_key_rate = 0.9;

-- Test edge cases for cache limits
SELECT 'Testing cache limit behavior...' as test_phase;

-- Set very small cache limit (minimum allowed is 1000)
SET kmersearch.actual_min_score_cache_max_entries = 1000;

-- Clear cache first
SELECT kmersearch_actual_min_score_cache_free();

-- Execute multiple different queries to test eviction
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'AAAAAAAA';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'TTTTTTTT';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'CCCCCCCC';
SELECT COUNT(*) FROM test_cache_DNA2 WHERE sequence =% 'GGGGGGGG';

-- Check final cache state
SELECT 'Final cache state:' as test_phase;
SELECT * FROM kmersearch_actual_min_score_cache_stats();

-- Clean up
DROP TABLE test_cache_DNA2, test_cache_DNA4 CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;