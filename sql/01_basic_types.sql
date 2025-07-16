SET client_min_messages = WARNING;
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test basic DNA2 and DNA4 type functionality
-- This test covers input/output, basic operations

-- Test DNA2 type
SELECT 'ATCG'::DNA2 AS DNA2_simple;
SELECT 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA2 AS DNA2_long;

-- Test DNA4 type with degenerate codes
SELECT 'ATCGMRWSYKN'::DNA4 AS DNA4_with_degenerate;
SELECT 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA4 AS DNA4_long;

-- Test invalid characters (should error)
\set ON_ERROR_STOP off
SELECT 'ATCGX'::DNA2;
SELECT 'ATCGZ'::DNA4;
\set ON_ERROR_STOP on

-- Test case insensitivity
SELECT 'atcg'::DNA2 AS DNA2_lowercase;
SELECT 'atcgmrwsykn'::DNA4 AS DNA4_lowercase;

-- Test with U (treated as T)
SELECT 'AUCG'::DNA2 AS DNA2_with_u;
SELECT 'AUCG'::DNA4 AS DNA4_with_u;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;