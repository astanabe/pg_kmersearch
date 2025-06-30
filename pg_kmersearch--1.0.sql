-- pg_kmersearch extension SQL definitions

-- Create shell types first
CREATE TYPE dna2;
CREATE TYPE dna4;

-- DNA2 input/output functions
CREATE FUNCTION kmersearch_dna2_in(cstring) RETURNS dna2
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_in'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_out(dna2) RETURNS cstring
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_out'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_recv(internal) RETURNS dna2
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_recv'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_send(dna2) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_send'
    LANGUAGE C IMMUTABLE STRICT;

-- DNA4 input/output functions
CREATE FUNCTION kmersearch_dna4_in(cstring) RETURNS dna4
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_in'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_out(dna4) RETURNS cstring
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_out'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_recv(internal) RETURNS dna4
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_recv'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_send(dna4) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_send'
    LANGUAGE C IMMUTABLE STRICT;

-- Complete DNA2 type definition
CREATE TYPE dna2 (
    INPUT = kmersearch_dna2_in,
    OUTPUT = kmersearch_dna2_out,
    RECEIVE = kmersearch_dna2_recv,
    SEND = kmersearch_dna2_send,
    STORAGE = extended,
    ALIGNMENT = int4
);

-- Complete DNA4 type definition
CREATE TYPE dna4 (
    INPUT = kmersearch_dna4_in,
    OUTPUT = kmersearch_dna4_out,
    RECEIVE = kmersearch_dna4_recv,
    SEND = kmersearch_dna4_send,
    STORAGE = extended,
    ALIGNMENT = int4
);

-- Equality operators
CREATE FUNCTION kmersearch_dna2_eq(dna2, dna2) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_eq'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR = (
    LEFTARG = dna2,
    RIGHTARG = dna2,
    FUNCTION = kmersearch_dna2_eq,
    COMMUTATOR = =,
    HASHES
);

CREATE FUNCTION kmersearch_dna4_eq(dna4, dna4) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_eq'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR = (
    LEFTARG = dna4,
    RIGHTARG = dna4,
    FUNCTION = kmersearch_dna4_eq,
    COMMUTATOR = =,
    HASHES
);

-- =% operators for k-mer search
CREATE FUNCTION kmersearch_dna2_match(dna2, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_match'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR =% (
    LEFTARG = dna2,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna2_match
);

CREATE FUNCTION kmersearch_dna4_match(dna4, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_match'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR =% (
    LEFTARG = dna4,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna4_match
);

-- Note: GIN operator classes commented out for initial testing
-- Will be added back after basic functionality is confirmed
-- 
-- GIN operator class support functions
-- CREATE FUNCTION kmersearch_extract_value(dna2, internal)
--     RETURNS internal
--     AS 'MODULE_PATHNAME', 'kmersearch_extract_value'
--     LANGUAGE C IMMUTABLE STRICT;
-- 
-- CREATE FUNCTION kmersearch_extract_value_dna4(dna4, internal)
--     RETURNS internal
--     AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna4'
--     LANGUAGE C IMMUTABLE STRICT;
-- 
-- CREATE FUNCTION kmersearch_extract_query(text, internal, int2, internal, internal)
--     RETURNS internal
--     AS 'MODULE_PATHNAME', 'kmersearch_extract_query'
--     LANGUAGE C IMMUTABLE STRICT;
-- 
-- CREATE FUNCTION kmersearch_consistent(internal, int2, text, int4, internal, internal)
--     RETURNS boolean
--     AS 'MODULE_PATHNAME', 'kmersearch_consistent'
--     LANGUAGE C IMMUTABLE STRICT;
-- 
-- -- GIN operator classes for DNA2 and DNA4
-- CREATE OPERATOR CLASS kmersearch_dna2_gin_ops
--     DEFAULT FOR TYPE dna2 USING gin AS
--         OPERATOR 1 LIKE,
--         FUNCTION 1 varbit_cmp(varbit, varbit),
--         FUNCTION 2 kmersearch_extract_value(dna2, internal),
--         FUNCTION 3 kmersearch_extract_query(text, internal, int2, internal, internal),
--         FUNCTION 4 kmersearch_consistent(internal, int2, text, int4, internal, internal),
--         STORAGE varbit;
-- 
-- CREATE OPERATOR CLASS kmersearch_dna4_gin_ops
--     DEFAULT FOR TYPE dna4 USING gin AS
--         OPERATOR 1 LIKE,
--         FUNCTION 1 varbit_cmp(varbit, varbit),
--         FUNCTION 2 kmersearch_extract_value_dna4(dna4, internal),
--         FUNCTION 3 kmersearch_extract_query(text, internal, int2, internal, internal),
--         FUNCTION 4 kmersearch_consistent(internal, int2, text, int4, internal, internal),
--         STORAGE varbit;

-- Configuration function
CREATE FUNCTION set_kmersearch_occur_bitlen(integer) RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_set_occur_bitlen'
    LANGUAGE C IMMUTABLE STRICT;

-- System table for storing excluded k-mers
CREATE TABLE kmersearch_excluded_kmers (
    index_oid oid NOT NULL,
    kmer_key varbit NOT NULL,
    frequency_count integer NOT NULL,
    exclusion_reason text,
    created_at timestamp with time zone DEFAULT now(),
    PRIMARY KEY (index_oid, kmer_key)
);

-- Index tracking table
CREATE TABLE kmersearch_index_info (
    index_oid oid PRIMARY KEY,
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    k_value integer NOT NULL,
    total_rows bigint NOT NULL,
    excluded_kmers_count integer NOT NULL,
    max_appearance_rate real NOT NULL,
    max_appearance_nrow integer NOT NULL,
    created_at timestamp with time zone DEFAULT now()
);

-- K-mer frequency analysis functions
CREATE FUNCTION kmersearch_analyze_table_frequency(table_oid oid, column_name text, k integer, index_oid oid) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_analyze_table_frequency'
    LANGUAGE C VOLATILE STRICT;

CREATE FUNCTION kmersearch_get_excluded_kmers(index_oid oid) 
    RETURNS varbit[]
    AS 'MODULE_PATHNAME', 'kmersearch_get_excluded_kmers'
    LANGUAGE C STABLE STRICT;

-- Show minimum score function
CREATE FUNCTION show_kmersearch_min_score() 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'show_kmersearch_min_score'
    LANGUAGE C STABLE;

-- Score calculation functions
CREATE FUNCTION kmersearch_rawscore(dna2, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_rawscore'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_rawscore(dna4, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_rawscore'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_correctedscore(dna2, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_correctedscore'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_correctedscore(dna4, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_correctedscore'
    LANGUAGE C IMMUTABLE STRICT;