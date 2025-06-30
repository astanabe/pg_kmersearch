CREATE EXTENSION pg_kmersearch;

-- Test basic DNA2 and DNA4 type functionality
-- This test covers input/output, basic operations

-- Test DNA2 type
SELECT 'ATCG'::dna2 AS dna2_simple;
SELECT 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::dna2 AS dna2_long;

-- Test DNA4 type with degenerate codes
SELECT 'ATCGMRWSYKN'::dna4 AS dna4_with_degenerate;
SELECT 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::dna4 AS dna4_long;

-- Test invalid characters (should error)
\set ON_ERROR_STOP off
SELECT 'ATCGX'::dna2;
SELECT 'ATCGZ'::dna4;
\set ON_ERROR_STOP on

-- Test case insensitivity
SELECT 'atcg'::dna2 AS dna2_lowercase;
SELECT 'atcgmrwsykn'::dna4 AS dna4_lowercase;

-- Test with U (treated as T)
SELECT 'AUCG'::dna2 AS dna2_with_u;
SELECT 'AUCG'::dna4 AS dna4_with_u;

DROP EXTENSION pg_kmersearch CASCADE;