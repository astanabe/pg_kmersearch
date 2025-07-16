-- 高頻度k-mer除外機能の並列キャッシュテスト
-- kmersearch.kmer_size=4、kmersearch.occur_bitlen=4、kmersearch.max_appearance_rate=0.25、
-- kmersearch.max_appearance_nrow=0、kmersearch.min_score=1、kmersearch.min_shared_ngram_key_rate=0.9、
-- kmersearch.preclude_highfreq_kmer=true、kmersearch.force_use_parallel_highfreq_kmer_cache=true

-- 拡張を作成
CREATE EXTENSION IF NOT EXISTS pg_kmersearch;

-- 指定されたGUC設定を適用
SET kmersearch.kmer_size = 4;
SET kmersearch.occur_bitlen = 4;
SET kmersearch.max_appearance_rate = 0.25;
SET kmersearch.max_appearance_nrow = 0;
SET kmersearch.min_score = 1;
SET kmersearch.min_shared_ngram_key_rate = 0.9;
SET kmersearch.preclude_highfreq_kmer = true;
SET kmersearch.force_use_parallel_highfreq_kmer_cache = true;

-- 設定確認
SELECT 'GUC設定確認:' as info;
SHOW kmersearch.kmer_size;
SHOW kmersearch.max_appearance_rate;
SHOW kmersearch.preclude_highfreq_kmer;

-- テスト用テーブルの作成
DROP TABLE IF EXISTS test_dna_highfreq CASCADE;
CREATE TABLE test_dna_highfreq (
    id SERIAL PRIMARY KEY,
    seq DNA2
);

-- k=4で高頻度k-merが含まれるテストデータを挿入
INSERT INTO test_dna_highfreq (seq) VALUES 
    ('AAAACTGTACGT'::DNA2),    -- AAAA, AAAa, AACa, ACTa が含まれる
    ('AAAAGCATGCAT'::DNA2),    -- AAAA, AAAa, AAGa, AGCa が含まれる  
    ('AAAATCGATCGA'::DNA2),    -- AAAA, AAAa, AATa, ATCa が含まれる
    ('AAAACCCCGGGG'::DNA2),    -- AAAA, AAAa, AACc, CCCCが含まれる
    ('AAAATTTTAAAA'::DNA2),    -- AAAA, AAAa, AATt, TTTTが含まれる
    ('TCGTAAAACGTA'::DNA2),    -- AAAA が含まれる
    ('GCATAAAATCGA'::DNA2),    -- AAAA が含まれる
    ('ATCGAAAACCCC'::DNA2),    -- AAAA, CCCC が含まれる
    ('CCCCAAAATTTT'::DNA2),    -- CCCC, AAAA, TTTT が含まれる
    ('TTTTAAAACCCC'::DNA2),    -- TTTT, AAAA, CCCC が含まれる
    ('CCCCCCCCCCCC'::DNA2),    -- CCCCが多数含まれる
    ('TTTTTTTTTTTT'::DNA2),    -- TTTTが多数含まれる
    ('ACGTACGTACGT'::DNA2),    -- 違うパターン
    ('TGCATGCATGCA'::DNA2);    -- 違うパターン

SELECT 'テストデータ作成完了:' as info, COUNT(*) as total_rows FROM test_dna_highfreq;

-- 1. kmersearch_perform_highfreq_analysis()を実行
SELECT '1. 高頻度k-mer分析実行:' as step;
SELECT kmersearch_perform_highfreq_analysis('test_dna_highfreq', 'seq');

-- 2. kmersearch_parallel_highfreq_kmer_cache_load()を実行
SELECT '2. 高頻度k-mer並列キャッシュロード:' as step;
SELECT kmersearch_parallel_highfreq_kmer_cache_load('test_dna_highfreq', 'seq');

-- 3. CREATE INDEX (高頻度出現ngram_key2がGINインデックスから除外される)
SELECT '3. GINインデックス作成:' as step;
CREATE INDEX idx_test_dna_highfreq_seq ON test_dna_highfreq USING gin (seq);

-- インデックス作成確認
SELECT '   インデックス作成確認:' as info;
SELECT schemaname, tablename, indexname 
FROM pg_indexes 
WHERE tablename = 'test_dna_highfreq' AND indexname = 'idx_test_dna_highfreq_seq';

-- 4. GINインデックスを使用した検索
SELECT '4. GINインデックスを使用した検索:' as step;
SELECT id, seq, kmersearch_rawscore_dna2(seq, 'AAAACTGT') as rawscore
FROM test_dna_highfreq 
WHERE seq =% 'AAAACTGT'
ORDER BY rawscore DESC, id;

-- 5. kmersearch_parallel_highfreq_kmer_cache_free()を実行
SELECT '5. 高頻度k-mer並列キャッシュ解放:' as step;
SELECT kmersearch_parallel_highfreq_kmer_cache_free('test_dna_highfreq', 'seq');

-- 6. 再度GINインデックスを使用した検索
SELECT '6. キャッシュ解放後の再度GINインデックス検索:' as step;
SELECT id, seq, kmersearch_rawscore_dna2(seq, 'AAAACTGT') as rawscore
FROM test_dna_highfreq 
WHERE seq =% 'AAAACTGT'
ORDER BY rawscore DESC, id;

-- 追加テスト: 高頻度k-mer除外の効果確認
SELECT '7. 高頻度k-mer除外効果の確認:' as step;

-- 高頻度k-merの状況確認
SELECT 'システムテーブルの内容確認:' as info;
SELECT COUNT(*) as highfreq_kmer_count 
FROM kmersearch_highfreq_kmer 
WHERE table_oid = (SELECT oid FROM pg_class WHERE relname = 'test_dna_highfreq');

-- 複数のクエリで検索結果を確認
SELECT '複数クエリでの検索結果:' as info;
SELECT 'AAAACTGT' as query, COUNT(*) as match_count 
FROM test_dna_highfreq 
WHERE seq =% 'AAAACTGT';

SELECT 'CCCCTTTT' as query, COUNT(*) as match_count 
FROM test_dna_highfreq 
WHERE seq =% 'CCCCTTTT';

-- テスト完了メッセージ
SELECT '高頻度k-mer除外機能のテスト完了' as test_result;

-- クリーンアップ
DROP TABLE test_dna_highfreq CASCADE;

-- 拡張を削除してリグレッションテストをクリーンアップ
DROP EXTENSION pg_kmersearch CASCADE;