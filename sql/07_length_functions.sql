CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- Test length functions for DNA2 and DNA4 types
-- This test covers bit_length, nuc_length, char_length, and length functions

-- Test bit_length function with non-multiples to test padding
SELECT bit_length('ATCGA'::DNA2) AS dna2_5bp_bits;  -- 4+1
SELECT bit_length('ATCGATCGA'::DNA2) AS dna2_9bp_bits;  -- 8+1
SELECT bit_length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA2) AS dna2_61bp_bits;  -- 60+1

SELECT bit_length('ATCGA'::DNA4) AS dna4_5bp_bits;  -- 4+1
SELECT bit_length('ATCGATCGA'::DNA4) AS dna4_9bp_bits;  -- 8+1
SELECT bit_length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA4) AS dna4_61bp_bits;  -- 60+1

-- Test nuc_length function with non-multiples to test padding
SELECT nuc_length('ATCGA'::DNA2) AS dna2_5bp_nucs;  -- 4+1
SELECT nuc_length('ATCGATCGA'::DNA2) AS dna2_9bp_nucs;  -- 8+1
SELECT nuc_length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA2) AS dna2_61bp_nucs;  -- 60+1

SELECT nuc_length('ATCGA'::DNA4) AS dna4_5bp_nucs;  -- 4+1
SELECT nuc_length('ATCGATCGA'::DNA4) AS dna4_9bp_nucs;  -- 8+1
SELECT nuc_length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA4) AS dna4_61bp_nucs;  -- 60+1

-- Test char_length function (should be same as nuc_length) with non-multiples
SELECT char_length('ATCGA'::DNA2) AS dna2_5bp_chars;  -- 4+1
SELECT char_length('ATCGATCGA'::DNA2) AS dna2_9bp_chars;  -- 8+1
SELECT char_length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA2) AS dna2_61bp_chars;  -- 60+1

SELECT char_length('ATCGA'::DNA4) AS dna4_5bp_chars;  -- 4+1
SELECT char_length('ATCGATCGA'::DNA4) AS dna4_9bp_chars;  -- 8+1
SELECT char_length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA4) AS dna4_61bp_chars;  -- 60+1

-- Test length function (should be same as nuc_length) with non-multiples
SELECT length('ATCGA'::DNA2) AS dna2_5bp_length;  -- 4+1
SELECT length('ATCGATCGA'::DNA2) AS dna2_9bp_length;  -- 8+1
SELECT length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA2) AS dna2_61bp_length;  -- 60+1

SELECT length('ATCGA'::DNA4) AS dna4_5bp_length;  -- 4+1
SELECT length('ATCGATCGA'::DNA4) AS dna4_9bp_length;  -- 8+1
SELECT length('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'::DNA4) AS dna4_61bp_length;  -- 60+1

-- Test relationship: nuc_length = bit_length / 2 for DNA2 (with padding test)
SELECT 
    bit_length('ATCGATCGA'::DNA2) / 2 = nuc_length('ATCGATCGA'::DNA2) AS dna2_bit_nuc_relation_9bp;

-- Test relationship: nuc_length = bit_length / 4 for DNA4 (with padding test)
SELECT 
    bit_length('ATCGATCGA'::DNA4) / 4 = nuc_length('ATCGATCGA'::DNA4) AS dna4_bit_nuc_relation_9bp;

-- Test that char_length = nuc_length = length (with padding test)
SELECT 
    nuc_length('ATCGATCGA'::DNA2) = char_length('ATCGATCGA'::DNA2) AND
    char_length('ATCGATCGA'::DNA2) = length('ATCGATCGA'::DNA2) AS dna2_length_consistency_9bp;

SELECT 
    nuc_length('ATCGATCGA'::DNA4) = char_length('ATCGATCGA'::DNA4) AND
    char_length('ATCGATCGA'::DNA4) = length('ATCGATCGA'::DNA4) AS dna4_length_consistency_9bp;

DROP EXTENSION pg_kmersearch CASCADE;