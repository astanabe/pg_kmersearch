-- pg_kmersearch extension SQL definitions

-- Create shell types first
CREATE TYPE DNA2;
CREATE TYPE DNA4;

-- DNA2 input/output functions
CREATE FUNCTION kmersearch_dna2_in(cstring) RETURNS DNA2
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_in'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_out(DNA2) RETURNS cstring
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_out'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_recv(internal) RETURNS DNA2
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_recv'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_send(DNA2) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_send'
    LANGUAGE C IMMUTABLE STRICT;

-- DNA4 input/output functions
CREATE FUNCTION kmersearch_dna4_in(cstring) RETURNS DNA4
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_in'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_out(DNA4) RETURNS cstring
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_out'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_recv(internal) RETURNS DNA4
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_recv'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_send(DNA4) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_send'
    LANGUAGE C IMMUTABLE STRICT;

-- Complete DNA2 type definition
CREATE TYPE DNA2 (
    INPUT = kmersearch_dna2_in,
    OUTPUT = kmersearch_dna2_out,
    RECEIVE = kmersearch_dna2_recv,
    SEND = kmersearch_dna2_send,
    STORAGE = extended,
    ALIGNMENT = int4
);

-- Complete DNA4 type definition
CREATE TYPE DNA4 (
    INPUT = kmersearch_dna4_in,
    OUTPUT = kmersearch_dna4_out,
    RECEIVE = kmersearch_dna4_recv,
    SEND = kmersearch_dna4_send,
    STORAGE = extended,
    ALIGNMENT = int4
);

-- Equality operators
CREATE FUNCTION kmersearch_dna2_eq(DNA2, DNA2) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_eq'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR = (
    LEFTARG = DNA2,
    RIGHTARG = DNA2,
    FUNCTION = kmersearch_dna2_eq,
    COMMUTATOR = =,
    HASHES
);

CREATE FUNCTION kmersearch_dna4_eq(DNA4, DNA4) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_eq'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR = (
    LEFTARG = DNA4,
    RIGHTARG = DNA4,
    FUNCTION = kmersearch_dna4_eq,
    COMMUTATOR = =,
    HASHES
);

-- =% operators for k-mer search
CREATE FUNCTION kmersearch_dna2_match(DNA2, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_match'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR =% (
    LEFTARG = DNA2,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna2_match
);

CREATE FUNCTION kmersearch_dna4_match(DNA4, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_match'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR =% (
    LEFTARG = DNA4,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna4_match
);

-- GIN operator class support functions
CREATE FUNCTION kmersearch_extract_value_dna2(DNA2, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna2'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_value_dna4(DNA4, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna4'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_query(text, internal, int2, internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_query'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_consistent(internal, int2, text, int4, internal, internal)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_consistent'
    LANGUAGE C IMMUTABLE STRICT;

-- GIN operator classes for DNA2 and DNA4
CREATE OPERATOR CLASS kmersearch_dna2_gin_ops
    DEFAULT FOR TYPE DNA2 USING gin AS
        OPERATOR 1 =% (DNA2, text),
        FUNCTION 1 varbitcmp(varbit, varbit),
        FUNCTION 2 kmersearch_extract_value_dna2(DNA2, internal),
        FUNCTION 3 kmersearch_extract_query(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent(internal, int2, text, int4, internal, internal),
        STORAGE varbit;

CREATE OPERATOR CLASS kmersearch_dna4_gin_ops
    DEFAULT FOR TYPE DNA4 USING gin AS
        OPERATOR 1 =% (DNA4, text),
        FUNCTION 1 varbitcmp(varbit, varbit),
        FUNCTION 2 kmersearch_extract_value_dna4(DNA4, internal),
        FUNCTION 3 kmersearch_extract_query(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent(internal, int2, text, int4, internal, internal),
        STORAGE varbit;


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


-- Score calculation functions
CREATE FUNCTION kmersearch_rawscore(DNA2, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_rawscore'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_rawscore(DNA4, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_rawscore'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_correctedscore(DNA2, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_correctedscore'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_correctedscore(DNA4, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_correctedscore'
    LANGUAGE C IMMUTABLE STRICT;

-- Length functions for DNA2 and DNA4 types

-- bit_length functions
CREATE FUNCTION bit_length(DNA2) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_bit_length'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION bit_length(DNA4) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_bit_length'
    LANGUAGE C IMMUTABLE STRICT;

-- nuc_length functions
CREATE FUNCTION nuc_length(DNA2) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_nuc_length'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION nuc_length(DNA4) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_nuc_length'
    LANGUAGE C IMMUTABLE STRICT;

-- char_length functions (same as nuc_length)
CREATE FUNCTION char_length(DNA2) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_char_length'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION char_length(DNA4) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_char_length'
    LANGUAGE C IMMUTABLE STRICT;

-- length functions (same as nuc_length)
CREATE FUNCTION length(DNA2) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_char_length'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION length(DNA4) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_char_length'
    LANGUAGE C IMMUTABLE STRICT;