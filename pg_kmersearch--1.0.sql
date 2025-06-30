-- pg_kmersearch extension SQL definitions

-- DNA2 type (2-bit encoding for ACGT)
CREATE TYPE dna2 (
    INPUT = kmersearch_dna2_in,
    OUTPUT = kmersearch_dna2_out,
    RECEIVE = kmersearch_dna2_recv,
    SEND = kmersearch_dna2_send,
    STORAGE = extended,
    ALIGNMENT = int4
);

-- DNA4 type (4-bit encoding with degenerate codes)
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

-- LIKE operators for k-mer search
CREATE FUNCTION kmersearch_dna2_like(dna2, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_like'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR ~~ (
    LEFTARG = dna2,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna2_like
);

CREATE FUNCTION kmersearch_dna4_like(dna4, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_like'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR ~~ (
    LEFTARG = dna4,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna4_like
);

-- GIN operator class support functions
CREATE FUNCTION kmersearch_extract_value(anyarray, internal, internal, internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_query(anyarray, internal, int2, internal, internal, internal, internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_query'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_consistent(internal, int2, anyarray, int4, internal, internal, internal, internal)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_consistent'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_compare_partial(anyarray, anyarray, int2, internal)
    RETURNS int4
    AS 'MODULE_PATHNAME', 'kmersearch_compare_partial'
    LANGUAGE C IMMUTABLE STRICT;

-- GIN operator classes for DNA2 and DNA4
CREATE OPERATOR CLASS kmersearch_dna2_gin_ops
    DEFAULT FOR TYPE dna2 USING gin AS
        OPERATOR 1 ~~,
        FUNCTION 1 varbit_cmp(varbit, varbit),
        FUNCTION 2 kmersearch_extract_value(anyarray, internal, internal, internal, internal),
        FUNCTION 3 kmersearch_extract_query(anyarray, internal, int2, internal, internal, internal, internal, internal),
        FUNCTION 4 kmersearch_consistent(internal, int2, anyarray, int4, internal, internal, internal, internal),
        FUNCTION 6 kmersearch_compare_partial(anyarray, anyarray, int2, internal),
        STORAGE varbit;

CREATE OPERATOR CLASS kmersearch_dna4_gin_ops
    DEFAULT FOR TYPE dna4 USING gin AS
        OPERATOR 1 ~~,
        FUNCTION 1 varbit_cmp(varbit, varbit),
        FUNCTION 2 kmersearch_extract_value(anyarray, internal, internal, internal, internal),
        FUNCTION 3 kmersearch_extract_query(anyarray, internal, int2, internal, internal, internal, internal, internal),
        FUNCTION 4 kmersearch_consistent(internal, int2, anyarray, int4, internal, internal, internal, internal),
        FUNCTION 6 kmersearch_compare_partial(anyarray, anyarray, int2, internal),
        STORAGE varbit;

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