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

-- BTree comparison functions
CREATE FUNCTION kmersearch_dna2_cmp(DNA2, DNA2) RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_cmp'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_cmp(DNA4, DNA4) RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_cmp'
    LANGUAGE C IMMUTABLE STRICT;

-- DNA2 comparison operators
CREATE FUNCTION kmersearch_dna2_lt(DNA2, DNA2) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_lt'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_le(DNA2, DNA2) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_le'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_gt(DNA2, DNA2) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_gt'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_ge(DNA2, DNA2) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_ge'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna2_ne(DNA2, DNA2) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_ne'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR < (
    LEFTARG = DNA2,
    RIGHTARG = DNA2,
    FUNCTION = kmersearch_dna2_lt,
    COMMUTATOR = >,
    NEGATOR = >=
);

CREATE OPERATOR <= (
    LEFTARG = DNA2,
    RIGHTARG = DNA2,
    FUNCTION = kmersearch_dna2_le,
    COMMUTATOR = >=,
    NEGATOR = >
);

CREATE OPERATOR > (
    LEFTARG = DNA2,
    RIGHTARG = DNA2,
    FUNCTION = kmersearch_dna2_gt,
    COMMUTATOR = <,
    NEGATOR = <=
);

CREATE OPERATOR >= (
    LEFTARG = DNA2,
    RIGHTARG = DNA2,
    FUNCTION = kmersearch_dna2_ge,
    COMMUTATOR = <=,
    NEGATOR = <
);

CREATE OPERATOR <> (
    LEFTARG = DNA2,
    RIGHTARG = DNA2,
    FUNCTION = kmersearch_dna2_ne,
    COMMUTATOR = <>,
    NEGATOR = =
);

-- DNA4 comparison operators
CREATE FUNCTION kmersearch_dna4_lt(DNA4, DNA4) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_lt'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_le(DNA4, DNA4) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_le'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_gt(DNA4, DNA4) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_gt'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_ge(DNA4, DNA4) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_ge'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_ne(DNA4, DNA4) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_ne'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR < (
    LEFTARG = DNA4,
    RIGHTARG = DNA4,
    FUNCTION = kmersearch_dna4_lt,
    COMMUTATOR = >,
    NEGATOR = >=
);

CREATE OPERATOR <= (
    LEFTARG = DNA4,
    RIGHTARG = DNA4,
    FUNCTION = kmersearch_dna4_le,
    COMMUTATOR = >=,
    NEGATOR = >
);

CREATE OPERATOR > (
    LEFTARG = DNA4,
    RIGHTARG = DNA4,
    FUNCTION = kmersearch_dna4_gt,
    COMMUTATOR = <,
    NEGATOR = <=
);

CREATE OPERATOR >= (
    LEFTARG = DNA4,
    RIGHTARG = DNA4,
    FUNCTION = kmersearch_dna4_ge,
    COMMUTATOR = <=,
    NEGATOR = <
);

CREATE OPERATOR <> (
    LEFTARG = DNA4,
    RIGHTARG = DNA4,
    FUNCTION = kmersearch_dna4_ne,
    COMMUTATOR = <>,
    NEGATOR = =
);

-- =% operators for k-mer search
CREATE FUNCTION kmersearch_dna2_match(DNA2, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_match'
    LANGUAGE C IMMUTABLE STRICT
    COST 1000;

CREATE OPERATOR =% (
    LEFTARG = DNA2,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna2_match
);

CREATE FUNCTION kmersearch_dna4_match(DNA4, text) RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_match'
    LANGUAGE C IMMUTABLE STRICT
    COST 1000;

CREATE OPERATOR =% (
    LEFTARG = DNA4,
    RIGHTARG = text,
    FUNCTION = kmersearch_dna4_match
);


-- New uintkey-based GIN functions for DNA2
CREATE FUNCTION kmersearch_extract_value_dna2_int2(DNA2, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna2_int2'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_value_dna2_int4(DNA2, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna2_int4'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_value_dna2_int8(DNA2, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna2_int8'
    LANGUAGE C IMMUTABLE STRICT;

-- New uintkey-based GIN functions for DNA4
CREATE FUNCTION kmersearch_extract_value_dna4_int2(DNA4, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna4_int2'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_value_dna4_int4(DNA4, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna4_int4'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_value_dna4_int8(DNA4, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_value_dna4_int8'
    LANGUAGE C IMMUTABLE STRICT;

-- New uintkey-based extract_query functions
CREATE FUNCTION kmersearch_extract_query_int2(text, internal, int2, internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_query_int2'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_query_int4(text, internal, int2, internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_query_int4'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_extract_query_int8(text, internal, int2, internal, internal)
    RETURNS internal
    AS 'MODULE_PATHNAME', 'kmersearch_extract_query_int8'
    LANGUAGE C IMMUTABLE STRICT;

-- New uintkey-based consistent functions
CREATE FUNCTION kmersearch_consistent_int2(internal, int2, text, int4, internal, internal)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_consistent_int2'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_consistent_int4(internal, int2, text, int4, internal, internal)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_consistent_int4'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_consistent_int8(internal, int2, text, int4, internal, internal)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_consistent_int8'
    LANGUAGE C IMMUTABLE STRICT;

-- New uintkey-based GIN operator classes for DNA2
CREATE OPERATOR CLASS kmersearch_dna2_gin_ops_int2
    FOR TYPE DNA2 USING gin AS
        OPERATOR 1 =% (DNA2, text),
        FUNCTION 1 btint2cmp(int2, int2),
        FUNCTION 2 kmersearch_extract_value_dna2_int2(DNA2, internal),
        FUNCTION 3 kmersearch_extract_query_int2(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent_int2(internal, int2, text, int4, internal, internal),
        STORAGE int2;

CREATE OPERATOR CLASS kmersearch_dna2_gin_ops_int4
    FOR TYPE DNA2 USING gin AS
        OPERATOR 1 =% (DNA2, text),
        FUNCTION 1 btint4cmp(int4, int4),
        FUNCTION 2 kmersearch_extract_value_dna2_int4(DNA2, internal),
        FUNCTION 3 kmersearch_extract_query_int4(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent_int4(internal, int2, text, int4, internal, internal),
        STORAGE int4;

CREATE OPERATOR CLASS kmersearch_dna2_gin_ops_int8
    FOR TYPE DNA2 USING gin AS
        OPERATOR 1 =% (DNA2, text),
        FUNCTION 1 btint8cmp(int8, int8),
        FUNCTION 2 kmersearch_extract_value_dna2_int8(DNA2, internal),
        FUNCTION 3 kmersearch_extract_query_int8(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent_int8(internal, int2, text, int4, internal, internal),
        STORAGE int8;

-- New uintkey-based GIN operator classes for DNA4
CREATE OPERATOR CLASS kmersearch_dna4_gin_ops_int2
    FOR TYPE DNA4 USING gin AS
        OPERATOR 1 =% (DNA4, text),
        FUNCTION 1 btint2cmp(int2, int2),
        FUNCTION 2 kmersearch_extract_value_dna4_int2(DNA4, internal),
        FUNCTION 3 kmersearch_extract_query_int2(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent_int2(internal, int2, text, int4, internal, internal),
        STORAGE int2;

CREATE OPERATOR CLASS kmersearch_dna4_gin_ops_int4
    FOR TYPE DNA4 USING gin AS
        OPERATOR 1 =% (DNA4, text),
        FUNCTION 1 btint4cmp(int4, int4),
        FUNCTION 2 kmersearch_extract_value_dna4_int4(DNA4, internal),
        FUNCTION 3 kmersearch_extract_query_int4(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent_int4(internal, int2, text, int4, internal, internal),
        STORAGE int4;

CREATE OPERATOR CLASS kmersearch_dna4_gin_ops_int8
    FOR TYPE DNA4 USING gin AS
        OPERATOR 1 =% (DNA4, text),
        FUNCTION 1 btint8cmp(int8, int8),
        FUNCTION 2 kmersearch_extract_value_dna4_int8(DNA4, internal),
        FUNCTION 3 kmersearch_extract_query_int8(text, internal, int2, internal, internal),
        FUNCTION 4 kmersearch_consistent_int8(internal, int2, text, int4, internal, internal),
        STORAGE int8;

-- BTree operator classes for DNA2 and DNA4
CREATE OPERATOR CLASS kmersearch_dna2_btree_ops
    DEFAULT FOR TYPE DNA2 USING btree AS
        OPERATOR 1 <,
        OPERATOR 2 <=,
        OPERATOR 3 =,
        OPERATOR 4 >=,
        OPERATOR 5 >,
        FUNCTION 1 kmersearch_dna2_cmp(DNA2, DNA2);

CREATE OPERATOR CLASS kmersearch_dna4_btree_ops
    DEFAULT FOR TYPE DNA4 USING btree AS
        OPERATOR 1 <,
        OPERATOR 2 <=,
        OPERATOR 3 =,
        OPERATOR 4 >=,
        OPERATOR 5 >,
        FUNCTION 1 kmersearch_dna4_cmp(DNA4, DNA4);

-- Hash functions for GROUP BY and hash-based operations
CREATE FUNCTION kmersearch_dna2_hash(DNA2) RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_hash'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_hash(DNA4) RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_hash'
    LANGUAGE C IMMUTABLE STRICT;

-- Extended hash functions for improved collision resistance
CREATE FUNCTION kmersearch_dna2_hash_extended(DNA2, bigint) RETURNS bigint
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_hash_extended'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_hash_extended(DNA4, bigint) RETURNS bigint
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_hash_extended'
    LANGUAGE C IMMUTABLE STRICT;

-- Hash operator classes for DNA2 and DNA4
CREATE OPERATOR CLASS kmersearch_dna2_hash_ops
    DEFAULT FOR TYPE DNA2 USING hash AS
        OPERATOR 1 =,
        FUNCTION 1 kmersearch_dna2_hash(DNA2),
        FUNCTION 2 kmersearch_dna2_hash_extended(DNA2, bigint);

CREATE OPERATOR CLASS kmersearch_dna4_hash_ops
    DEFAULT FOR TYPE DNA4 USING hash AS
        OPERATOR 1 =,
        FUNCTION 1 kmersearch_dna4_hash(DNA4),
        FUNCTION 2 kmersearch_dna4_hash_extended(DNA4, bigint);

-- System table for storing highly frequent k-mers
CREATE TABLE kmersearch_highfreq_kmer (
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    uintkey bigint NOT NULL,
    detection_reason text,
    created_at timestamp with time zone DEFAULT now(),
    PRIMARY KEY (table_oid, column_name, uintkey)
);

-- High-frequency k-mers metadata table
CREATE TABLE kmersearch_highfreq_kmer_meta (
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    kmer_size integer NOT NULL,
    occur_bitlen integer NOT NULL,
    max_appearance_rate real NOT NULL,
    max_appearance_nrow integer NOT NULL,
    analysis_timestamp timestamp with time zone DEFAULT now(),
    PRIMARY KEY (table_oid, column_name, kmer_size)
);

-- GIN index metadata table
CREATE TABLE kmersearch_gin_index_meta (
    index_oid oid PRIMARY KEY,
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    highfreq_filtered boolean NOT NULL,
    highfreq_source_table name NOT NULL,
    kmer_size integer NOT NULL,
    occur_bitlen integer NOT NULL,
    max_appearance_rate real NOT NULL,
    max_appearance_nrow integer NOT NULL,
    created_at timestamp with time zone DEFAULT now()
);

-- Index tracking table
CREATE TABLE kmersearch_index_info (
    index_oid oid PRIMARY KEY,
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    kmer_size integer NOT NULL,
    occur_bitlen integer NOT NULL,
    total_nrow bigint NOT NULL,
    highfreq_kmer_count integer NOT NULL,
    max_appearance_rate real NOT NULL,
    max_appearance_nrow integer NOT NULL,
    preclude_highfreq_kmer boolean NOT NULL DEFAULT false,
    created_at timestamp with time zone DEFAULT now()
);


-- Score calculation functions
-- Type-specific matchscore functions
CREATE FUNCTION kmersearch_matchscore_dna2(DNA2, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_matchscore_dna2'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_matchscore_dna4(DNA4, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_matchscore_dna4'
    LANGUAGE C IMMUTABLE STRICT;

-- Overloaded matchscore functions for convenience
CREATE FUNCTION kmersearch_matchscore(DNA2, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_matchscore_dna2'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_matchscore(DNA4, text) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_matchscore_dna4'
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
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_nuc_length'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION length(DNA4) 
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_nuc_length'
    LANGUAGE C IMMUTABLE STRICT;

-- SIMD capability detection function
CREATE FUNCTION kmersearch_simd_capability() 
    RETURNS text
    AS 'MODULE_PATHNAME', 'kmersearch_simd_capability'
    LANGUAGE C IMMUTABLE;

-- BYTEA conversion functions for hash compatibility
CREATE FUNCTION kmersearch_dna2_to_bytea(DNA2) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_to_bytea'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmersearch_dna4_to_bytea(DNA4) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_to_bytea'
    LANGUAGE C IMMUTABLE STRICT;

-- Complex types for function return values
CREATE TYPE kmersearch_analysis_result AS (
    total_rows bigint,
    highfreq_kmers_count integer,
    parallel_workers_used integer,
    max_appearance_rate_used real,
    max_appearance_nrow_used integer
);

CREATE TYPE kmersearch_drop_result AS (
    dropped_analyses integer,
    dropped_highfreq_kmers integer,
    freed_storage_bytes bigint
);

-- Build information function
CREATE FUNCTION kmersearch_show_buildno() 
    RETURNS text
    AS 'MODULE_PATHNAME', 'kmersearch_show_buildno'
    LANGUAGE C IMMUTABLE STRICT;

-- Parallel k-mer analysis functions
CREATE FUNCTION kmersearch_perform_highfreq_analysis(table_name text, column_name text) 
    RETURNS kmersearch_analysis_result
    AS 'MODULE_PATHNAME', 'kmersearch_perform_highfreq_analysis'
    LANGUAGE C VOLATILE STRICT;

CREATE FUNCTION kmersearch_undo_highfreq_analysis(table_name text, column_name text) 
    RETURNS kmersearch_drop_result
    AS 'MODULE_PATHNAME', 'kmersearch_undo_highfreq_analysis'
    LANGUAGE C VOLATILE STRICT;


-- Query-kmer cache statistics function
CREATE FUNCTION kmersearch_query_kmer_cache_stats()
    RETURNS TABLE (
        hits bigint,
        misses bigint,
        current_entries integer,
        max_entries integer
    )
    AS 'MODULE_PATHNAME', 'kmersearch_query_kmer_cache_stats'
    LANGUAGE C STABLE;

-- Query-kmer cache management function
CREATE FUNCTION kmersearch_query_kmer_cache_free()
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_query_kmer_cache_free'
    LANGUAGE C VOLATILE;

-- Actual min score cache statistics function
CREATE FUNCTION kmersearch_actual_min_score_cache_stats()
    RETURNS TABLE (
        hits bigint,
        misses bigint,
        current_entries integer,
        max_entries integer
    )
    AS 'MODULE_PATHNAME', 'kmersearch_actual_min_score_cache_stats'
    LANGUAGE C STABLE;

-- Actual min score cache management function
CREATE FUNCTION kmersearch_actual_min_score_cache_free()
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_actual_min_score_cache_free'
    LANGUAGE C VOLATILE;

-- High-frequency k-mer cache management functions
CREATE FUNCTION kmersearch_highfreq_kmer_cache_load(table_name text, column_name text)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_highfreq_kmer_cache_load'
    LANGUAGE C VOLATILE STRICT;

CREATE FUNCTION kmersearch_highfreq_kmer_cache_free(table_name text, column_name text)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_highfreq_kmer_cache_free'
    LANGUAGE C VOLATILE STRICT;

-- Parallel high-frequency k-mer cache management functions
CREATE FUNCTION kmersearch_parallel_highfreq_kmer_cache_load(table_name text, column_name text)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmersearch_parallel_highfreq_kmer_cache_load'
    LANGUAGE C VOLATILE STRICT;

CREATE FUNCTION kmersearch_parallel_highfreq_kmer_cache_free(table_name text, column_name text)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_parallel_highfreq_kmer_cache_free'
    LANGUAGE C VOLATILE STRICT;

-- Cache free functions without parameters (for backwards compatibility)
CREATE FUNCTION kmersearch_highfreq_kmer_cache_free_all()
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_highfreq_kmer_cache_free_all'
    LANGUAGE C VOLATILE;

CREATE FUNCTION kmersearch_parallel_highfreq_kmer_cache_free_all()
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmersearch_parallel_highfreq_kmer_cache_free_all'
    LANGUAGE C VOLATILE;

-- Performance optimization: Add indexes on system tables
CREATE INDEX kmersearch_highfreq_kmer_idx 
    ON kmersearch_highfreq_kmer(table_oid, column_name, uintkey);

CREATE INDEX kmersearch_highfreq_kmer_meta_idx 
    ON kmersearch_highfreq_kmer_meta(table_oid, column_name);

CREATE INDEX kmersearch_gin_index_meta_idx
    ON kmersearch_gin_index_meta(table_oid, column_name);

CREATE INDEX kmersearch_index_info_idx
    ON kmersearch_index_info(table_oid, column_name);

-- Management views for easier administration
CREATE VIEW kmersearch_cache_summary AS
SELECT 
    'query-kmer' as cache_type,
    current_entries as total_entries,
    hits as total_hits,
    misses as total_misses,
    CASE WHEN (hits + misses) > 0 
         THEN hits::float / (hits + misses)::float 
         ELSE 0 END as hit_rate
FROM kmersearch_query_kmer_cache_stats()
UNION ALL
SELECT 
    'actual_min_score' as cache_type,
    current_entries as total_entries,
    hits as total_hits,
    misses as total_misses,
    CASE WHEN (hits + misses) > 0 
         THEN hits::float / (hits + misses)::float 
         ELSE 0 END as hit_rate
FROM kmersearch_actual_min_score_cache_stats();

-- View for high-frequency k-mer analysis status
CREATE VIEW kmersearch_analysis_status AS
SELECT 
    m.table_oid,
    pg_class.relname as table_name,
    m.column_name,
    m.kmer_size,
    m.occur_bitlen,
    m.max_appearance_rate,
    m.max_appearance_nrow,
    m.analysis_timestamp,
    COUNT(DISTINCT h.uintkey) as highfreq_kmer_count
FROM kmersearch_highfreq_kmer_meta m
LEFT JOIN pg_class ON pg_class.oid = m.table_oid
LEFT JOIN kmersearch_highfreq_kmer h ON (
    h.table_oid = m.table_oid 
    AND h.column_name = m.column_name
)
GROUP BY m.table_oid, pg_class.relname, m.column_name, m.kmer_size, 
         m.occur_bitlen, m.max_appearance_rate, m.max_appearance_nrow, 
         m.analysis_timestamp;

-- Partitioning support functions
CREATE FUNCTION kmersearch_partition_table(table_name text, partition_count int, tablespace_name text DEFAULT NULL)
RETURNS void
AS 'MODULE_PATHNAME', 'kmersearch_partition_table'
LANGUAGE C;

-- Temporary file cleanup function
CREATE FUNCTION kmersearch_delete_tempfiles()
RETURNS TABLE(
    deleted_count integer,
    deleted_size bigint,
    error_count integer
)
AS 'MODULE_PATHNAME', 'kmersearch_delete_tempfiles'
LANGUAGE C STRICT;

-- Event trigger for capturing GIN index creation and saving settings
CREATE FUNCTION kmersearch_on_index_create()
RETURNS event_trigger
LANGUAGE plpgsql
AS $$
DECLARE
    obj record;
    idx_oid oid;
    tbl_oid oid;
    col_name name;
    opclass_name name;
BEGIN
    FOR obj IN SELECT * FROM pg_event_trigger_ddl_commands()
        WHERE object_type = 'index' AND command_tag = 'CREATE INDEX'
    LOOP
        idx_oid := obj.objid;

        SELECT
            i.indrelid,
            a.attname,
            oc.opcname
        INTO tbl_oid, col_name, opclass_name
        FROM pg_index i
        JOIN pg_class ic ON ic.oid = i.indexrelid
        JOIN pg_am am ON am.oid = ic.relam
        JOIN pg_attribute a ON a.attrelid = i.indrelid
            AND a.attnum = ANY(i.indkey)
        JOIN pg_opclass oc ON oc.oid = i.indclass[0]
        WHERE i.indexrelid = idx_oid
          AND am.amname = 'gin'
          AND oc.opcname LIKE 'kmersearch_%';

        IF FOUND THEN
            INSERT INTO kmersearch_index_info (
                index_oid, table_oid, column_name,
                kmer_size, occur_bitlen,
                max_appearance_rate, max_appearance_nrow,
                preclude_highfreq_kmer,
                total_nrow, highfreq_kmer_count
            ) VALUES (
                idx_oid, tbl_oid, col_name,
                current_setting('kmersearch.kmer_size')::integer,
                current_setting('kmersearch.occur_bitlen')::integer,
                current_setting('kmersearch.max_appearance_rate')::real,
                current_setting('kmersearch.max_appearance_nrow')::integer,
                current_setting('kmersearch.preclude_highfreq_kmer')::boolean,
                0, 0
            )
            ON CONFLICT (index_oid) DO UPDATE SET
                kmer_size = EXCLUDED.kmer_size,
                occur_bitlen = EXCLUDED.occur_bitlen,
                max_appearance_rate = EXCLUDED.max_appearance_rate,
                max_appearance_nrow = EXCLUDED.max_appearance_nrow,
                preclude_highfreq_kmer = EXCLUDED.preclude_highfreq_kmer,
                created_at = now();
        END IF;
    END LOOP;
END;
$$;

CREATE EVENT TRIGGER kmersearch_index_create_trigger
ON ddl_command_end
WHEN TAG IN ('CREATE INDEX')
EXECUTE FUNCTION kmersearch_on_index_create();

-- Event trigger for cleaning up metadata when index is dropped
CREATE FUNCTION kmersearch_on_index_drop()
RETURNS event_trigger
LANGUAGE plpgsql
AS $$
DECLARE
    obj record;
BEGIN
    FOR obj IN SELECT * FROM pg_event_trigger_dropped_objects()
        WHERE object_type = 'index'
    LOOP
        DELETE FROM kmersearch_index_info WHERE index_oid = obj.objid;
        DELETE FROM kmersearch_gin_index_meta WHERE index_oid = obj.objid;
    END LOOP;
END;
$$;

CREATE EVENT TRIGGER kmersearch_index_drop_trigger
ON sql_drop
WHEN TAG IN ('DROP INDEX', 'DROP TABLE', 'DROP SCHEMA', 'DROP EXTENSION')
EXECUTE FUNCTION kmersearch_on_index_drop();