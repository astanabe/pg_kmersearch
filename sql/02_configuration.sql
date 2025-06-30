CREATE EXTENSION pg_kmersearch;

-- Test GUC configuration variables
-- This test covers all configuration parameters

-- Show default values
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
SHOW kmersearch.max_appearance_rate;
SHOW kmersearch.max_appearance_nrow;
SHOW kmersearch.min_score;

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

-- Test boundary values
SET kmersearch.kmer_size = 4;  -- minimum
SHOW kmersearch.kmer_size;

SET kmersearch.kmer_size = 64; -- maximum
SHOW kmersearch.kmer_size;

SET kmersearch.occur_bitlen = 0;  -- minimum
SHOW kmersearch.occur_bitlen;

SET kmersearch.occur_bitlen = 16; -- maximum
SHOW kmersearch.occur_bitlen;

-- Test invalid values (should error)
\set ON_ERROR_STOP off
SET kmersearch.kmer_size = 3;   -- below minimum
SET kmersearch.kmer_size = 65;  -- above maximum
SET kmersearch.occur_bitlen = 17; -- above maximum
SET kmersearch.max_appearance_rate = 1.1; -- above maximum
\set ON_ERROR_STOP on

-- Reset to defaults for other tests
SET kmersearch.kmer_size = 8;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.max_appearance_rate = 0.05;
SET kmersearch.max_appearance_nrow = 0;
SET kmersearch.min_score = 1;

DROP EXTENSION pg_kmersearch CASCADE;