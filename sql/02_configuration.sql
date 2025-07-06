CREATE EXTENSION pg_kmersearch;

-- Test GUC configuration variables
-- This test covers all configuration parameters

-- Show default values
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
SHOW kmersearch.max_appearance_rate;
SHOW kmersearch.max_appearance_nrow;
SHOW kmersearch.min_score;
SHOW kmersearch.min_shared_ngram_key_rate;
SHOW kmersearch.rawscore_cache_max_entries;

SHOW kmersearch.query_pattern_cache_max_entries;

-- Test setting valid values
SET kmersearch.kmer_size = 12;
SHOW kmersearch.kmer_size;

SET kmersearch.occur_bitlen = 10;
SHOW kmersearch.occur_bitlen;

SET kmersearch.max_appearance_rate = 0.1;
SHOW kmersearch.max_appearance_rate;

SET kmersearch.max_appearance_nrow = 1000;
SHOW kmersearch.max_appearance_nrow;

SET kmersearch.min_score = 5;
SHOW kmersearch.min_score;

SET kmersearch.min_shared_ngram_key_rate = 0.8;
SHOW kmersearch.min_shared_ngram_key_rate;

SET kmersearch.rawscore_cache_max_entries = 25000;
SHOW kmersearch.rawscore_cache_max_entries;

SET kmersearch.query_pattern_cache_max_entries = 25000;
SHOW kmersearch.query_pattern_cache_max_entries;

-- Test boundary values
SET kmersearch.kmer_size = 4;  -- minimum
SHOW kmersearch.kmer_size;

SET kmersearch.kmer_size = 64; -- maximum
SHOW kmersearch.kmer_size;

SET kmersearch.occur_bitlen = 0;  -- minimum
SHOW kmersearch.occur_bitlen;

SET kmersearch.occur_bitlen = 16; -- maximum
SHOW kmersearch.occur_bitlen;

SET kmersearch.min_shared_ngram_key_rate = 0.0;  -- minimum
SHOW kmersearch.min_shared_ngram_key_rate;

SET kmersearch.min_shared_ngram_key_rate = 1.0;  -- maximum
SHOW kmersearch.min_shared_ngram_key_rate;

SET kmersearch.rawscore_cache_max_entries = 1000;  -- minimum
SHOW kmersearch.rawscore_cache_max_entries;

SET kmersearch.rawscore_cache_max_entries = 10000000;  -- maximum
SHOW kmersearch.rawscore_cache_max_entries;

SET kmersearch.query_pattern_cache_max_entries = 1000;  -- minimum
SHOW kmersearch.query_pattern_cache_max_entries;

SET kmersearch.query_pattern_cache_max_entries = 10000000;  -- maximum
SHOW kmersearch.query_pattern_cache_max_entries;

-- Test invalid values (should error)
\set ON_ERROR_STOP off
SET kmersearch.kmer_size = 3;   -- below minimum
SET kmersearch.kmer_size = 65;  -- above maximum
SET kmersearch.occur_bitlen = 17; -- above maximum
SET kmersearch.max_appearance_rate = 1.1; -- above maximum
SET kmersearch.min_shared_ngram_key_rate = -0.1;  -- below minimum
SET kmersearch.min_shared_ngram_key_rate = 1.1;   -- above maximum
SET kmersearch.rawscore_cache_max_entries = 999;  -- below minimum
SET kmersearch.rawscore_cache_max_entries = 10000001;  -- above maximum
SET kmersearch.query_pattern_cache_max_entries = 999;  -- below minimum
SET kmersearch.query_pattern_cache_max_entries = 10000001;  -- above maximum
\set ON_ERROR_STOP on

-- Reset to defaults for other tests
SET kmersearch.kmer_size = 8;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.max_appearance_rate = 0.05;
SET kmersearch.max_appearance_nrow = 0;
SET kmersearch.min_score = 1;
SET kmersearch.min_shared_ngram_key_rate = 0.9;
SET kmersearch.rawscore_cache_max_entries = 50000;
SET kmersearch.query_pattern_cache_max_entries = 50000;

-- Test cache management functions
SELECT kmersearch_rawscore_cache_stats(); -- Should show all zeros initially
SELECT kmersearch_query_pattern_cache_stats(); -- Should show all zeros initially
SELECT kmersearch_query_pattern_cache_free();  -- Should return 0 (no entries to free)

DROP EXTENSION pg_kmersearch CASCADE;