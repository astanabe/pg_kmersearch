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