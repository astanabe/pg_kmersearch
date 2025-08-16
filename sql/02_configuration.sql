SET client_min_messages = WARNING;
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test GUC configuration variables
-- This test covers all configuration parameters

-- Show default values
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
SHOW kmersearch.max_appearance_rate;
SHOW kmersearch.max_appearance_nrow;
SHOW kmersearch.min_score;
SHOW kmersearch.min_shared_kmer_rate;
SHOW kmersearch.query_kmer_cache_max_entries;

-- Test setting valid values
SET kmersearch.kmer_size = 6;
SHOW kmersearch.kmer_size;

SET kmersearch.occur_bitlen = 10;
SHOW kmersearch.occur_bitlen;

SET kmersearch.max_appearance_rate = 0.1;
SHOW kmersearch.max_appearance_rate;

SET kmersearch.max_appearance_nrow = 1000;
SHOW kmersearch.max_appearance_nrow;

SET kmersearch.min_score = 5;
SHOW kmersearch.min_score;

SET kmersearch.min_shared_kmer_rate = 0.8;
SHOW kmersearch.min_shared_kmer_rate;

SET kmersearch.query_kmer_cache_max_entries = 25000;
SHOW kmersearch.query_kmer_cache_max_entries;

-- Test boundary values
SET kmersearch.kmer_size = 4;  -- minimum
SHOW kmersearch.kmer_size;

-- Test maximum k-mer size with occur_bitlen=10 (should error: 32*2+10=74 > 64)
\set ON_ERROR_STOP off
SET kmersearch.kmer_size = 32; -- should fail with occur_bitlen=10
\set ON_ERROR_STOP on
SHOW kmersearch.kmer_size; -- should still be 4

-- Test maximum k-mer size with occur_bitlen=0 (should succeed: 32*2+0=64)
SET kmersearch.occur_bitlen = 0;
SET kmersearch.kmer_size = 32; -- should succeed now
SHOW kmersearch.kmer_size;

-- Reset occur_bitlen (should fail because kmer_size is now 32)
\set ON_ERROR_STOP off
SET kmersearch.occur_bitlen = 10;
\set ON_ERROR_STOP on

SET kmersearch.occur_bitlen = 0;  -- minimum
SHOW kmersearch.occur_bitlen;

\set ON_ERROR_STOP off
SET kmersearch.occur_bitlen = 16; -- maximum (should fail with kmer_size=32)
\set ON_ERROR_STOP on

SET kmersearch.min_shared_kmer_rate = 0.0;  -- minimum
SHOW kmersearch.min_shared_kmer_rate;

SET kmersearch.min_shared_kmer_rate = 1.0;  -- maximum
SHOW kmersearch.min_shared_kmer_rate;

SET kmersearch.query_kmer_cache_max_entries = 1000;  -- minimum
SHOW kmersearch.query_kmer_cache_max_entries;

SET kmersearch.query_kmer_cache_max_entries = 10000000;  -- maximum
SHOW kmersearch.query_kmer_cache_max_entries;

-- Test invalid values (should error)
\set ON_ERROR_STOP off
SET kmersearch.kmer_size = 3;   -- below minimum
SET kmersearch.kmer_size = 33;  -- above maximum
SET kmersearch.occur_bitlen = -1; -- below minimum  
SET kmersearch.occur_bitlen = 17; -- above maximum
SET kmersearch.max_appearance_rate = -0.1; -- below minimum
SET kmersearch.max_appearance_rate = 1.1; -- above maximum
SET kmersearch.min_shared_kmer_rate = -0.1;  -- below minimum
SET kmersearch.min_shared_kmer_rate = 1.1;   -- above maximum
SET kmersearch.query_kmer_cache_max_entries = 999;  -- below minimum
SET kmersearch.query_kmer_cache_max_entries = 10000001;  -- above maximum
SET kmersearch.actual_min_score_cache_max_entries = 999;  -- below minimum
SET kmersearch.actual_min_score_cache_max_entries = 10000001;  -- above maximum
\set ON_ERROR_STOP on

-- Test that valid values still work after errors
SHOW kmersearch.kmer_size; -- should still be 32
SHOW kmersearch.occur_bitlen; -- should still be 0

-- Test combinations of kmer_size and occur_bitlen
SET kmersearch.kmer_size = 4;
SET kmersearch.occur_bitlen = 8;
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
-- Total: 4*2+8=16 bits (valid)

SET kmersearch.kmer_size = 8;
SET kmersearch.occur_bitlen = 8;  
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
-- Total: 8*2+8=24 bits (valid)

SET kmersearch.kmer_size = 16;
SET kmersearch.occur_bitlen = 8;
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
-- Total: 16*2+8=40 bits (valid)

SET kmersearch.kmer_size = 27;
SET kmersearch.occur_bitlen = 10;
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
-- Total: 27*2+10=64 bits (valid, maximum)

\set ON_ERROR_STOP off
SET kmersearch.kmer_size = 28; -- 28*2+10=66 > 64 (should fail)
\set ON_ERROR_STOP on
SHOW kmersearch.kmer_size; -- should still be 27

-- Reset to defaults for other tests (ensure consistency across tests)
SET kmersearch.kmer_size = 16;
SHOW kmersearch.kmer_size;
SET kmersearch.occur_bitlen = 8;
SHOW kmersearch.occur_bitlen;
SET kmersearch.max_appearance_rate = 0.5;
SHOW kmersearch.max_appearance_rate;
SET kmersearch.max_appearance_nrow = 0;
SHOW kmersearch.max_appearance_nrow;
SET kmersearch.min_score = 1;
SHOW kmersearch.min_score;
SET kmersearch.min_shared_kmer_rate = 0.2;
SHOW kmersearch.min_shared_kmer_rate;
SET kmersearch.query_kmer_cache_max_entries = 50000;
SHOW kmersearch.query_kmer_cache_max_entries;
SET kmersearch.actual_min_score_cache_max_entries = 50000;
SHOW kmersearch.actual_min_score_cache_max_entries;

-- Test cache management functions
SELECT kmersearch_query_kmer_cache_stats(); -- Should show all zeros initially
SELECT kmersearch_actual_min_score_cache_stats(); -- Should show all zeros initially
SELECT kmersearch_query_kmer_cache_free();  -- Should return 0 (no entries to free)

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;