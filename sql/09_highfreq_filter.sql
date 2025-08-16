SET client_min_messages = WARNING;
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- High-frequency k-mer exclusion functionality comprehensive test
-- kmersearch.kmer_size=4, kmersearch.occur_bitlen=4, kmersearch.max_appearance_rate=0.25,
-- kmersearch.max_appearance_nrow=0, kmersearch.min_score=1, kmersearch.min_shared_kmer_rate=0.9,
-- kmersearch.preclude_highfreq_kmer=true, kmersearch.force_use_parallel_highfreq_kmer_cache=false

-- Apply specified GUC settings
SET kmersearch.kmer_size = 4;
SET kmersearch.occur_bitlen = 4;
SET kmersearch.max_appearance_rate = 0.25;
SET kmersearch.max_appearance_nrow = 0;
SET kmersearch.min_score = 1;
SET kmersearch.min_shared_kmer_rate = 0.9;
SET kmersearch.preclude_highfreq_kmer = true;
SET kmersearch.force_use_parallel_highfreq_kmer_cache = false;

-- Configuration verification
SELECT 'GUC settings verification:' as info;
SHOW kmersearch.kmer_size;
SHOW kmersearch.max_appearance_rate;
SHOW kmersearch.preclude_highfreq_kmer;

-- Create test table
DROP TABLE IF EXISTS test_dna_highfreq CASCADE;
CREATE TABLE test_dna_highfreq (
    id SERIAL PRIMARY KEY,
    seq DNA2
);

-- Insert test data containing high-frequency k-mers for k=4
INSERT INTO test_dna_highfreq (seq) VALUES 
    ('AAAACTGTACGT'::DNA2),    -- Contains AAAA, AAAa, AACa, ACTa
    ('AAAAGCATGCAT'::DNA2),    -- Contains AAAA, AAAa, AAGa, AGCa  
    ('AAAATCGATCGA'::DNA2),    -- Contains AAAA, AAAa, AATa, ATCa
    ('AAAACCCCGGGG'::DNA2),    -- Contains AAAA, AAAa, AACc, CCCC
    ('AAAATTTTAAAA'::DNA2),    -- Contains AAAA, AAAa, AATt, TTTT
    ('TCGTAAAACGTA'::DNA2),    -- Contains AAAA
    ('GCATAAAATCGA'::DNA2),    -- Contains AAAA
    ('ATCGAAAACCCC'::DNA2),    -- Contains AAAA, CCCC
    ('CCCCAAAATTTT'::DNA2),    -- Contains CCCC, AAAA, TTTT
    ('TTTTAAAACCCC'::DNA2),    -- Contains TTTT, AAAA, CCCC
    ('CCCCCCCCCCCC'::DNA2),    -- Contains many CCCC
    ('TTTTTTTTTTTT'::DNA2),    -- Contains many TTTT
    ('ACGTACGTACGT'::DNA2),    -- Different pattern
    ('TGCATGCATGCA'::DNA2);    -- Different pattern

SELECT 'Test data creation completed:' as info, COUNT(*) as total_rows FROM test_dna_highfreq;

-- 1. Execute kmersearch_perform_highfreq_analysis()
SELECT '1. High-frequency k-mer analysis execution:' as step;
SELECT kmersearch_perform_highfreq_analysis('test_dna_highfreq', 'seq');

-- 2. Execute kmersearch_highfreq_kmer_cache_load()
SELECT '2. High-frequency k-mer cache load:' as step;
SELECT kmersearch_highfreq_kmer_cache_load('test_dna_highfreq', 'seq');

-- 3. CREATE INDEX (high-frequency ngram_key2 are excluded from GIN index)
SELECT '3. GIN index creation:' as step;
-- With k=4 and occur_bitlen=4, total bits = 4*2+4 = 12, so we can use int2 operator class
CREATE INDEX idx_test_dna_highfreq_seq ON test_dna_highfreq USING gin (seq kmersearch_dna2_gin_ops_int2);

-- Index creation verification
SELECT '   Index creation verification:' as info;
SELECT schemaname, tablename, indexname 
FROM pg_indexes 
WHERE tablename = 'test_dna_highfreq' AND indexname = 'idx_test_dna_highfreq_seq';

-- 4. GIN index-based search
SELECT '4. GIN index-based search:' as step;
SELECT id, seq, kmersearch_matchscore(seq, 'AAAACTGT') as matchscore
FROM test_dna_highfreq 
WHERE seq =% 'AAAACTGT'
ORDER BY matchscore DESC, id;

-- 5. Execute kmersearch_highfreq_kmer_cache_free()
SELECT '5. High-frequency k-mer cache release:' as step;
SELECT kmersearch_highfreq_kmer_cache_free('test_dna_highfreq', 'seq');

-- 6. GIN index-based search after cache release
SELECT '6. GIN index search after cache release:' as step;
SELECT id, seq, kmersearch_matchscore(seq, 'AAAACTGT') as matchscore
FROM test_dna_highfreq 
WHERE seq =% 'AAAACTGT'
ORDER BY matchscore DESC, id;

-- Additional test: High-frequency k-mer exclusion effect verification
SELECT '7. High-frequency k-mer exclusion effect verification:' as step;

-- High-frequency k-mer status verification
SELECT 'System table content verification:' as info;
SELECT COUNT(*) as highfreq_kmer_count 
FROM kmersearch_highfreq_kmer 
WHERE table_oid = (SELECT oid FROM pg_class WHERE relname = 'test_dna_highfreq');

-- Multiple queries search result verification
SELECT 'Multiple queries search results:' as info;
SELECT 'AAAACTGT' as query, COUNT(*) as match_count 
FROM test_dna_highfreq 
WHERE seq =% 'AAAACTGT';

SELECT 'CCCCTTTT' as query, COUNT(*) as match_count 
FROM test_dna_highfreq 
WHERE seq =% 'CCCCTTTT';

-- Test completion message
SELECT 'High-frequency k-mer exclusion functionality test completed' as test_result;

-- Cleanup
DROP TABLE test_dna_highfreq CASCADE;

DROP EXTENSION pg_kmersearch CASCADE;
SET client_min_messages = NOTICE;