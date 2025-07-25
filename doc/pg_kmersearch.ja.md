# pg_kmersearch

PostgreSQL用のDNA配列型拡張とk-mer検索機能

## 概要

pg_kmersearchは、PostgreSQL用のDNA配列データを効率的に格納・処理するためのカスタムデータ型とk-mer検索機能を提供する拡張です。この拡張は以下の2つのデータ型と高速検索機能を実装しています：

### DNA2型
- **用途**: 標準的なDNA配列（A、C、G、T）を格納
- **エンコーディング**: 2ビット圧縮（A=00, C=01, G=10, T=11）
- **対応文字**: A, C, G, T, U（UはTとして扱われます）
- **ストレージ効率**: 1文字あたり2ビット

### DNA4型
- **用途**: 縮重コードを含むDNA配列を格納
- **エンコーディング**: 4ビット圧縮
- **対応文字**:
  - 標準塩基: A(0001), C(0010), G(0100), T(1000)
  - 縮重コード: M(0011), R(0101), W(1001), S(0110), Y(1010), K(1100), V(0111), H(1011), D(1101), B(1110), N(1111)
  - U（UはTとして扱われます）
- **ストレージ効率**: 1文字あたり4ビット

### k-mer検索機能
- **k-mer長**: 4～64塩基（インデックス作成時に指定）
- **GINインデックス**: n-gramキーによる高速検索
- **縮重コード対応**: DNA4型でのMRWSYKVHDBN展開
- **出現回数追跡**: 同一行内でのk-mer出現回数を考慮（デフォルト8ビット）
- **スコアリング検索**: 完全一致だけでなく、類似度による上位結果取得
- **高頻出k-mer除外**: インデックス作成時に過度に頻出するk-merを自動除外
- **スコアベースフィルタリング**: 除外k-merに応じて自動調整される最小スコア閾値
- **スコア計算関数**: 個別配列のスコア算出用の`kmersearch_rawscore()`と`kmersearch_correctedscore()`関数（現在の実装では両関数は同一の値を返します）
- **高頻出k-mer管理**: `kmersearch_perform_highfreq_analysis()`による高頻出k-mer解析と`kmersearch_undo_highfreq_analysis()`による解析データ削除、`kmersearch_highfreq_kmer_cache_load()`および`kmersearch_highfreq_kmer_cache_free()`によるキャッシュ管理

## インストール

### 前提条件
- PostgreSQL 16以上
- postgresql-server-dev-16
- make
- gcc または clang

### インストール手順

1. ソースコードをダウンロード：
```bash
git clone <repository-url>
cd pg_kmersearch
```

2. コンパイルとインストール：
```bash
make
sudo make install
```

3. shared_preload_librariesの設定（GUC変数に必要）：
```bash
# postgresql.confに追加
shared_preload_libraries = 'pg_kmersearch'
```

4. PostgreSQLを再起動し、データベースに拡張をインストール：
```sql
CREATE EXTENSION pg_kmersearch;
```

**重要：GUC変数に関する注意事項**
PostgreSQLは新しいセッションが開始されるたびにGUC変数をデフォルト値にリセットします。各セッションの開始時に必要なGUC変数を設定してください：

```sql
-- 例：GUC変数の設定
SET kmersearch.kmer_size = 4;
SET kmersearch.max_appearance_rate = 0.25;
-- 解析コマンドをここに記述
```

## 使用方法

### 基本的な使用例

```sql
-- DNA2型の使用例
CREATE TABLE sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    dna_seq DNA2
);

INSERT INTO sequences (name, dna_seq) VALUES 
    ('seq1', 'ATCGATCG'::DNA2),
    ('seq2', 'GCTAGCTA'::DNA2);

SELECT name, dna_seq FROM sequences;

-- DNA4型の使用例（縮重コード対応）
CREATE TABLE degenerate_sequences (
    id SERIAL PRIMARY KEY,
    name TEXT,
    dna_seq DNA4
);

INSERT INTO degenerate_sequences (name, dna_seq) VALUES 
    ('seq1', 'ATCGMRWSYKN'::DNA4),
    ('seq2', 'VHDBNATCG'::DNA4);

SELECT name, dna_seq FROM degenerate_sequences;
```

### k-mer検索機能の使用例

```sql
-- k-merサイズを設定（デフォルト16-mer）
SET kmersearch.kmer_size = 16;

-- GINインデックスを作成（現在のkmersearch.kmer_size設定を使用）
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq);

-- k-mer検索（=%演算子を使用）
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore,
       kmersearch_correctedscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS correctedscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') DESC 
LIMIT 10;

-- 縮重コードを含むクエリでの検索
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG') AS rawscore
FROM degenerate_sequences 
WHERE dna_seq =% 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'
ORDER BY rawscore DESC 
LIMIT 5;

-- 出現回数ビット長の設定（デフォルト8ビット）
SET kmersearch.occur_bitlen = 12; -- 12ビットに変更（最大4095回の出現をカウント）

-- 現在の設定を確認
SHOW kmersearch.kmer_size;
SHOW kmersearch.occur_bitlen;
```

### 設定変数

pg_kmersearchは、PostgreSQLの`SET`コマンドで設定可能な複数の設定変数を提供します：

| 変数名 | デフォルト値 | 範囲 | 説明 |
|--------|-------------|------|------|
| `kmersearch.kmer_size` | 16 | 4-64 | インデックス作成と検索のk-mer長 |
| `kmersearch.occur_bitlen` | 8 | 0-16 | 出現回数格納のビット数 |
| `kmersearch.max_appearance_rate` | 0.5 | 0.0-1.0 | インデックス化するk-merの最大出現率 |
| `kmersearch.max_appearance_nrow` | 0 | 0-∞ | k-merが含まれる最大行数（0=無制限） |
| `kmersearch.min_score` | 1 | 0-∞ | 検索結果の最小類似度スコア |
| `kmersearch.min_shared_ngram_key_rate` | 0.9 | 0.0-1.0 | 共有n-gramキー率の最小閾値 |
| `kmersearch.rawscore_cache_max_entries` | 50000 | 1000-10000000 | rawscoreキャッシュの最大エントリ数 |
| `kmersearch.query_pattern_cache_max_entries` | 50000 | 1000-10000000 | クエリパターンキャッシュの最大エントリ数 |
| `kmersearch.actual_min_score_cache_max_entries` | 50000 | 1000-10000000 | actual min scoreキャッシュの最大エントリ数 |
| `kmersearch.preclude_highfreq_kmer` | false | true/false | GINインデックス構築時の高頻出k-mer除外の有効化 |
| `kmersearch.force_use_parallel_highfreq_kmer_cache` | false | true/false | 高頻出k-mer検索での並列dshashキャッシュの強制使用 |
| `kmersearch.highfreq_kmer_cache_load_batch_size` | 10000 | 1000-1000000 | 高頻出k-merをキャッシュに読み込む際のバッチサイズ |

### 高頻出k-mer除外機能

インデックスのパフォーマンス向上のため、高頻出k-merの自動除外機能を搭載：

```sql
-- 除外パラメータの設定（インデックス作成前）
SET kmersearch.max_appearance_rate = 0.5;  -- デフォルト: 50%の最大出現率
SET kmersearch.max_appearance_nrow = 1000;  -- デフォルト: 0（無効）

-- 頻度解析付きインデックス作成
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq);

-- テーブル/カラムの除外k-mer確認
SELECT ngram_key, detection_reason 
FROM kmersearch_highfreq_kmer 
WHERE table_oid = 'sequences'::regclass AND column_name = 'dna_seq';

-- インデックス統計情報
SELECT table_oid, column_name, kmer_size, occur_bitlen, max_appearance_rate, max_appearance_nrow 
FROM kmersearch_highfreq_kmer_meta 
WHERE table_oid = 'sequences'::regclass;
```

### スコアベース検索フィルタリング

除外k-merに応じて自動調整される最小スコア閾値で検索品質を制御：

```sql
-- GIN検索結果の最小スコア設定
SET kmersearch.min_score = 50;  -- デフォルト: 1

-- 現在の最小スコア設定を確認
SHOW kmersearch.min_score;

-- 自動スコア調整による検索
-- クエリに除外k-merが3個含まれ、min_score=50の場合、
-- そのクエリでは実際の閾値は47に調整される
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY rawscore DESC LIMIT 10;
```

### スコア計算関数

マッチした各配列のスコアを個別に取得：

```sql
-- マッチした配列の生スコアを取得
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS rawscore
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY rawscore DESC 
LIMIT 10;

-- 修正スコア（除外k-merを考慮）を取得
SELECT id, name, dna_seq,
       kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS raw_score,
       kmersearch_correctedscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS corrected_score
FROM sequences 
WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY corrected_score DESC;

-- DNA2型とDNA4型両方の例
SELECT 'DNA2' as type, id, kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score
FROM dna2_sequences WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
UNION ALL
SELECT 'DNA4' as type, id, kmersearch_rawscore(dna_seq, 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA') AS score  
FROM dna4_sequences WHERE dna_seq =% 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY score DESC;
```

## 長さ関数

pg_kmersearchは、DNA2型およびDNA4型に対してパディングを正しく処理し、正確な測定値を返す複数の長さ関数を提供します：

### 利用可能な長さ関数

- **`bit_length(DNA2/DNA4)`**: 実際のビット長を返す（パディングを除く）
- **`nuc_length(DNA2/DNA4)`**: 塩基数を返す
- **`char_length(DNA2/DNA4)`**: `nuc_length()`と同じ（文字数）
- **`length(DNA2/DNA4)`**: `nuc_length()`と同じ（標準の長さ関数）

### BYTEA変換関数

- **`kmersearch_dna2_to_bytea(DNA2)`**: DNA2をハッシュ化に適したBYTEA形式に変換（ビット長プレフィックス付き）
- **`kmersearch_dna4_to_bytea(DNA4)`**: DNA4をハッシュ化に適したBYTEA形式に変換（ビット長プレフィックス付き）

これらの関数は、pgcryptoの`digest()`などの暗号学的ハッシュ関数での使用に最適化されています。ネットワークバイトオーダーで4バイトのビット長プレフィックスを含み、その後に実際のビットデータが続くため、パディングのみが異なる配列でもユニークなハッシュ値を保証します。

### 関数の関係性

- **DNA2**: `nuc_length() = bit_length() / 2` (1塩基あたり2ビット)
- **DNA4**: `nuc_length() = bit_length() / 4` (1塩基あたり4ビット)
- **一貫性**: `char_length() = length() = nuc_length()`

### 使用例

```sql
-- 基本的な長さ測定
SELECT 
    bit_length('ATCGA'::DNA2) AS bits,      -- 戻り値: 10
    nuc_length('ATCGA'::DNA2) AS nucs,      -- 戻り値: 5
    char_length('ATCGA'::DNA2) AS chars,    -- 戻り値: 5
    length('ATCGA'::DNA2) AS len;           -- 戻り値: 5

-- 縮重コードを含むDNA4
SELECT 
    bit_length('ATCGMRWSYKN'::DNA4) AS bits,  -- 戻り値: 44
    nuc_length('ATCGMRWSYKN'::DNA4) AS nucs;  -- 戻り値: 11

-- ハッシュ化のためのBYTEA変換（pgcrypto拡張が必要）
SELECT digest(kmersearch_dna2_to_bytea('ATCG'::DNA2), 'sha256');
SELECT digest(kmersearch_dna4_to_bytea('ATCGN'::DNA4), 'sha256');

-- パディング衝突回避の確認
SELECT 
    digest(kmersearch_dna2_to_bytea('ATG'::DNA2), 'sha256') AS hash_atg,
    digest(kmersearch_dna2_to_bytea('ATGA'::DNA2), 'sha256') AS hash_atga;
-- 類似したパディングにも関わらず異なるハッシュ値が生成される

-- パディング検証（4/8の倍数でない場合）
SELECT 
    bit_length('ATCGATCGA'::DNA2) AS bits,    -- 9塩基 * 2 = 18ビット
    nuc_length('ATCGATCGA'::DNA2) AS nucs;    -- 戻り値: 9

-- クエリでの長さ関数の使用
SELECT id, name, length(dna_seq) AS sequence_length
FROM sequences
WHERE length(dna_seq) >= 50
ORDER BY length(dna_seq) DESC;
```

### 縮重コードの意味

| コード | 意味 | ビット表現 |
|--------|------|-----------|
| A | アデニン | 0001 |
| C | シトシン | 0010 |
| G | グアニン | 0100 |
| T | チミン | 1000 |
| M | AまたはC | 0011 |
| R | AまたはG | 0101 |
| W | AまたはT | 1001 |
| S | CまたはG | 0110 |
| Y | CまたはT | 1010 |
| K | GまたはT | 1100 |
| V | AまたはCまたはG | 0111 |
| H | AまたはCまたはT | 1011 |
| D | AまたはGまたはT | 1101 |
| B | CまたはGまたはT | 1110 |
| N | AまたはCまたはGまたはT | 1111 |

## システムテーブルとビュー

pg_kmersearchは、メタデータ管理と監視のための複数のシステムテーブルとビューを作成します。

### システムテーブル

#### kmersearch_highfreq_kmer
GINインデックスから除外される高頻出k-merを格納：

```sql
-- 特定テーブル/カラムの除外k-merを表示
SELECT table_oid, column_name, ngram_key, detection_reason, created_at
FROM kmersearch_highfreq_kmer 
WHERE table_oid = 'sequences'::regclass AND column_name = 'dna_seq';

-- テーブル/カラムごとの除外k-mer数をカウント
SELECT table_oid, column_name, COUNT(*) as excluded_count
FROM kmersearch_highfreq_kmer
GROUP BY table_oid, column_name;
```

#### kmersearch_highfreq_kmer_meta
k-mer頻度解析のメタデータを格納：

```sql
-- 全テーブルの解析メタデータを表示
SELECT table_oid, column_name, kmer_size, occur_bitlen, 
       max_appearance_rate, max_appearance_nrow, analysis_timestamp
FROM kmersearch_highfreq_kmer_meta;

-- 特定テーブルの解析パラメータを確認
SELECT * FROM kmersearch_highfreq_kmer_meta 
WHERE table_oid = 'sequences'::regclass AND column_name = 'dna_seq';
```

#### kmersearch_gin_index_meta
高頻出フィルタリング情報を含むGINインデックスメタデータを格納：

```sql
-- GINインデックスメタデータを表示
SELECT index_oid, table_oid, column_name, highfreq_filtered, 
       kmer_size, occur_bitlen, created_at
FROM kmersearch_gin_index_meta;

-- インデックスが高頻出フィルタリングを使用しているか確認
SELECT highfreq_filtered FROM kmersearch_gin_index_meta 
WHERE index_oid = 'sequences_kmer_idx'::regclass;
```

#### kmersearch_index_info
包括的なインデックス統計と設定を格納：

```sql
-- インデックス統計を表示
SELECT index_oid, table_oid, column_name, kmer_size, total_nrow,
       highfreq_kmer_count, max_appearance_rate, created_at
FROM kmersearch_index_info;
```

### 管理ビュー

#### kmersearch_cache_summary
全キャッシュタイプの統合統計を提供：

```sql
-- キャッシュパフォーマンス概要を表示
SELECT cache_type, total_entries, total_hits, total_misses, hit_rate
FROM kmersearch_cache_summary;

-- キャッシュ効率を監視
SELECT cache_type, 
       ROUND(hit_rate * 100, 2) as hit_rate_percent,
       total_hits + total_misses as total_requests
FROM kmersearch_cache_summary
WHERE total_hits + total_misses > 0;
```

#### kmersearch_analysis_status
解析済み全テーブルの高頻出k-mer解析状況を表示：

```sql
-- 全テーブルの解析状況を表示
SELECT table_name, column_name, kmer_size, highfreq_kmer_count,
       analysis_timestamp
FROM kmersearch_analysis_status;

-- 解析パラメータと結果を確認
SELECT table_name, column_name, kmer_size, occur_bitlen,
       max_appearance_rate, max_appearance_nrow,
       highfreq_kmer_count, analysis_timestamp
FROM kmersearch_analysis_status
WHERE table_name = 'sequences';
```

## 解析・管理関数

### k-mer頻度解析

#### kmersearch_perform_highfreq_analysis()
テーブルに対する並列k-mer頻度解析を実行：

```sql
-- テーブル名とカラム名を使用した基本解析
SELECT kmersearch_perform_highfreq_analysis(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- 結果の解釈例
SELECT (result).total_rows,
       (result).highfreq_kmers_count,
       (result).parallel_workers_used,
       (result).max_appearance_rate_used
FROM (
    SELECT kmersearch_perform_highfreq_analysis('sequences', 'dna_seq') as result
) t;
```

#### kmersearch_undo_highfreq_analysis()
解析データを削除してストレージを解放：

```sql
-- 特定のテーブル/カラム組み合わせの解析を削除
SELECT kmersearch_undo_highfreq_analysis(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- 結果の解釈例
SELECT (result).dropped_analyses,
       (result).dropped_highfreq_kmers,
       (result).freed_storage_bytes
FROM (
    SELECT kmersearch_undo_highfreq_analysis('sequences', 'dna_seq') as result
) t;
```

### 高頻出k-merキャッシュ管理

#### グローバルキャッシュ関数

```sql
-- 高頻出k-merをグローバルキャッシュに読み込み
SELECT kmersearch_highfreq_kmer_cache_load(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- 特定テーブル/カラムのグローバルキャッシュを解放
SELECT kmersearch_highfreq_kmer_cache_free(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- グローバルキャッシュの全エントリを解放
SELECT kmersearch_highfreq_kmer_cache_free_all();
```

#### 並列キャッシュ関数

```sql
-- 高頻出k-merを並列キャッシュに読み込み（マルチプロセス共有用）
SELECT kmersearch_parallel_highfreq_kmer_cache_load(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- 特定テーブル/カラムの並列キャッシュを解放
SELECT kmersearch_parallel_highfreq_kmer_cache_free(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- 並列キャッシュの全エントリを解放
SELECT kmersearch_parallel_highfreq_kmer_cache_free_all();
```

### キャッシュ統計・管理

#### Rawscoreキャッシュ

```sql
-- rawscoreキャッシュ統計を表示
SELECT * FROM kmersearch_rawscore_cache_stats();

-- キャッシュパフォーマンスを監視
SELECT dna2_hits, dna2_misses, 
       CASE WHEN (dna2_hits + dna2_misses) > 0 
            THEN dna2_hits::float / (dna2_hits + dna2_misses)::float * 100
            ELSE 0 END as dna2_hit_rate_percent,
       dna4_hits, dna4_misses,
       CASE WHEN (dna4_hits + dna4_misses) > 0 
            THEN dna4_hits::float / (dna4_hits + dna4_misses)::float * 100
            ELSE 0 END as dna4_hit_rate_percent
FROM kmersearch_rawscore_cache_stats();

-- rawscoreキャッシュをクリア
SELECT kmersearch_rawscore_cache_free();
```

#### クエリパターンキャッシュ

```sql
-- クエリパターンキャッシュ統計を表示
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- キャッシュ効率を監視
SELECT hits, misses, current_entries, max_entries,
       CASE WHEN (hits + misses) > 0 
            THEN hits::float / (hits + misses)::float * 100
            ELSE 0 END as hit_rate_percent
FROM kmersearch_query_pattern_cache_stats();

-- クエリパターンキャッシュをクリア
SELECT kmersearch_query_pattern_cache_free();
```

#### Actual Min Scoreキャッシュ

```sql
-- actual min scoreキャッシュ統計を表示
SELECT * FROM kmersearch_actual_min_score_cache_stats();

-- actual min scoreキャッシュをクリア
SELECT kmersearch_actual_min_score_cache_free();
```


## 完全なワークフロー例

### k-mer解析とインデックス構築のセットアップ

```sql
-- 1. パラメータを設定
SET kmersearch.kmer_size = 16;
SET kmersearch.max_appearance_rate = 0.5;
SET kmersearch.max_appearance_nrow = 1000;
SET kmersearch.occur_bitlen = 8;

-- 2. 頻度解析を実行
SELECT kmersearch_perform_highfreq_analysis(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- 3. 解析結果を確認
SELECT * FROM kmersearch_analysis_status WHERE table_name = 'sequences';

-- 4. GINインデックスを作成（解析結果を自動使用）
CREATE INDEX sequences_kmer_idx ON sequences USING gin(dna_seq);

-- 5. 最適なパフォーマンスのためキャッシュを読み込み
SELECT kmersearch_highfreq_kmer_cache_load('sequences', 'dna_seq');
SELECT kmersearch_parallel_highfreq_kmer_cache_load('sequences', 'dna_seq');

-- 6. 検索を実行
SELECT id, name, kmersearch_rawscore(dna_seq, 'ATCGATCG') as score
FROM sequences 
WHERE dna_seq =% 'ATCGATCG'
ORDER BY score DESC;
```

### キャッシュパフォーマンス監視

```sql
-- 全キャッシュタイプを監視
SELECT cache_type, hit_rate, total_entries, 
       total_hits + total_misses as total_requests
FROM kmersearch_cache_summary;

-- rawscoreキャッシュの詳細解析
WITH cache_stats AS (
    SELECT dna2_hits, dna2_misses, dna2_entries,
           dna4_hits, dna4_misses, dna4_entries
    FROM kmersearch_rawscore_cache_stats()
)
SELECT 
    'DNA2' as type,
    dna2_entries as entries,
    dna2_hits as hits,
    dna2_misses as misses,
    CASE WHEN (dna2_hits + dna2_misses) > 0 
         THEN dna2_hits::float / (dna2_hits + dna2_misses)::float 
         ELSE 0 END as hit_rate
FROM cache_stats
UNION ALL
SELECT 
    'DNA4' as type,
    dna4_entries as entries, 
    dna4_hits as hits,
    dna4_misses as misses,
    CASE WHEN (dna4_hits + dna4_misses) > 0 
         THEN dna4_hits::float / (dna4_hits + dna4_misses)::float 
         ELSE 0 END as hit_rate
FROM cache_stats;
```

## 複合型定義

pg_kmersearchは、関数の戻り値用にカスタム複合型を定義します：

### kmersearch_analysis_result
`kmersearch_perform_highfreq_analysis()`によって返される：

```sql
-- 型定義相当:
-- CREATE TYPE kmersearch_analysis_result AS (
--     total_rows bigint,
--     highfreq_kmers_count integer,
--     parallel_workers_used integer,
--     max_appearance_rate_used real,
--     max_appearance_nrow_used integer
-- );

-- 使用例:
SELECT (result).total_rows,
       (result).highfreq_kmers_count,
       (result).max_appearance_rate_used
FROM (
    SELECT kmersearch_perform_highfreq_analysis('sequences', 'dna_seq') as result
) t;
```

### kmersearch_drop_result
`kmersearch_undo_highfreq_analysis()`によって返される：

```sql
-- 型定義相当:
-- CREATE TYPE kmersearch_drop_result AS (
--     dropped_analyses integer,
--     dropped_highfreq_kmers integer,
--     freed_storage_bytes bigint
-- );

-- 使用例:
SELECT (result).dropped_analyses,
       (result).dropped_highfreq_kmers,
       pg_size_pretty((result).freed_storage_bytes) as freed_storage
FROM (
    SELECT kmersearch_undo_highfreq_analysis('sequences', 'dna_seq') as result
) t;
```

## 技術的詳細

### アーキテクチャ
- PostgreSQLのvarbit型（bit varying型）をベースとして実装
- カスタムのエンコーディング・デコーディング関数
- メモリ効率的なストレージ形式
- GINインデックスによるk-mer検索サポート

### 内部実装
- **DNA2型**: 2ビット/文字でエンコード
- **DNA4型**: 4ビット/文字でエンコード
- **n-gramキー**: k-mer（2k bit）+ 出現回数（8-16 bit）
- **縮重コード展開**: 最大10組み合わせまで自動展開
- **並列インデックス作成**: max_parallel_maintenance_workersに対応
- **高頻出除外**: インデックス作成前の全テーブルスキャン
- **システムテーブル**: 除外k-merとインデックス統計のメタデータ格納（`kmersearch_highfreq_kmer`, `kmersearch_highfreq_kmer_meta`）
- **キャッシュシステム**: TopMemoryContext-based高速キャッシュ
- バイナリ入出力サポート

### k-mer検索の仕組み
1. **頻度解析**: 高頻出k-mer特定のための全テーブルスキャン
2. **k-mer抽出**: 指定されたk長でスライディングウィンドウ
3. **高頻出フィルタリング**: 出現閾値を超えるk-merの除外
4. **n-gramキー生成**: k-mer + 出現回数をバイナリエンコード
5. **縮重コード処理**: MRWSYKVHDBN を標準塩基に展開
6. **スコアリング**: 共有n-gramキー数による類似度計算

## キャッシュ管理機能

pg_kmersearchは、検索性能を向上させるための3種類の高速キャッシュシステムを提供します：

### Actual Min Score Cache
- **目的**: `=%`演算子での検索条件評価の最適化
- **仕組み**: `actual_min_score = max(kmersearch_min_score, ceil(kmersearch_min_shared_ngram_key_rate × query_total_kmers))`を事前計算してキャッシュ
- **使用場面**: 
  - `=%`演算子でのマッチング条件判定
  - rawscore cacheへの格納価値判定
- **メモリ管理**: TopMemoryContext-based実装

### Rawscore Cache
- **目的**: 計算済みrawscoreの高速取得
- **仕組み**: 配列とクエリの組み合わせ結果をキャッシュ
- **メモリ管理**: TopMemoryContext-based実装

### Query Pattern Cache
- **目的**: クエリパターンの再利用による高速化
- **メモリ管理**: TopMemoryContext-based実装

### キャッシュ統計・管理関数

```sql
-- キャッシュ統計情報の確認
SELECT * FROM kmersearch_actual_min_score_cache_stats();
SELECT * FROM kmersearch_rawscore_cache_stats();
SELECT * FROM kmersearch_query_pattern_cache_stats();

-- キャッシュクリア
SELECT kmersearch_actual_min_score_cache_free();
SELECT kmersearch_query_pattern_cache_free();
SELECT kmersearch_highfreq_kmer_cache_free_all();

-- キャッシュサイズ設定
SET kmersearch.actual_min_score_cache_max_entries = 25000;
SET kmersearch.rawscore_cache_max_entries = 25000;
SET kmersearch.query_pattern_cache_max_entries = 25000;
```

## 高頻出k-mer階層キャッシュシステム

pg_kmersearchは、PostgreSQL 16/18に最適化された多層キャッシュアーキテクチャを提供し、高頻出k-mer検索のパフォーマンスを最大化します：

### キャッシュ階層

高頻出k-mer検索では、以下の優先順位でキャッシュが使用されます：

1. **グローバルキャッシュ** (`global_highfreq_cache`)
   - 最も高速なアクセス
   - TopMemoryContextで管理
   - 単一プロセス内での再利用

2. **並列キャッシュ** (`parallel_highfreq_cache`) 
   - PostgreSQL dshash（動的共有ハッシュテーブル）実装
   - DSM（Dynamic Shared Memory）による複数プロセス間共有
   - 将来のPostgreSQL 18並列GINインデックス対応

3. **テーブル参照フォールバック** (`kmersearch_highfreq_kmer`)
   - システムテーブルへの直接アクセス
   - キャッシュ不在時の最終手段

### GUC設定検証機能

キャッシュ読み込み時に以下のGUC変数が自動検証されます：

- `kmersearch.max_appearance_rate`
- `kmersearch.max_appearance_nrow` 
- `kmersearch.occur_bitlen`
- `kmersearch.kmer_size`

設定不一致時は詳細なエラーメッセージとヒントが提供されます。

### 並列キャッシュの技術実装

- **メモリ管理**: DSMセグメントピン機能による自動クリーンアップ防止
- **プロセス識別**: メインプロセスとワーカーでの差別化されたリソース管理  
- **エラーハンドリング**: 包括的なPG_TRY/PG_CATCHによる安全な操作
- **ロック管理**: dshash_release_lock()による適切な並行制御

### キャッシュ使用例

```sql
-- 高頻出k-mer解析の実行（メタデータ作成）
SELECT kmersearch_perform_highfreq_analysis(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- GUC設定をメタデータに合わせる
SET kmersearch.max_appearance_rate = 0.5;
SET kmersearch.max_appearance_nrow = 100;
SET kmersearch.occur_bitlen = 8;
SET kmersearch.kmer_size = 16;

-- グローバルキャッシュの読み込み
SELECT kmersearch_highfreq_kmer_cache_load(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- 並列キャッシュの読み込み（オプション）
SELECT kmersearch_parallel_highfreq_kmer_cache_load(
    'sequences',                   -- テーブル名
    'dna_seq'                     -- カラム名
);

-- キャッシュ階層を使用した高速検索
SELECT id, name FROM sequences 
WHERE dna_seq =% 'ATCGATCG'
ORDER BY kmersearch_rawscore(dna_seq, 'ATCGATCG') DESC;

-- キャッシュの解放
SELECT kmersearch_highfreq_kmer_cache_free('sequences', 'dna_seq');
SELECT kmersearch_parallel_highfreq_kmer_cache_free('sequences', 'dna_seq');
```

### 並列キャッシュ関数

- **`kmersearch_parallel_highfreq_kmer_cache_load(table_name, column_name)`**: 高頻出k-merを共有dshashキャッシュに読み込み
- **`kmersearch_parallel_highfreq_kmer_cache_free(table_name, column_name)`**: 並列キャッシュから特定エントリを解放
- **`kmersearch_parallel_highfreq_kmer_cache_free_all()`**: 並列キャッシュからすべてのエントリを解放し、共有メモリ構造を破棄

### 使用シナリオ

並列キャッシュシステムは以下の用途に特に有用です：

1. **大規模ゲノム解析**での複数の並行クエリ
2. **DNA配列検索のバッチ処理**
3. **将来のPostgreSQL 18並列GINインデックス操作**
4. **高スループットバイオインフォマティクスアプリケーション**

### 技術的実装

- **メモリコンテキスト**: 適切なクリーンアップのためのTopMemoryContext
- **共有メモリ**: 自動ピン機能付きDSA（Dynamic Shared Areas）
- **エラーハンドリング**: 包括的なPG_TRY/PG_CATCHブロック
- **ロック管理**: すべての操作での適切なdshash_release_lock()呼び出し
- **プロセス検出**: 適切なワーカー識別のためのIsParallelWorker()

## 制限事項

- クエリ配列は最小8塩基が必要
- 縮重コード展開は10組み合わせまで（超過時はスキップ）
- 出現回数は設定可能ビット長の最大値でキャップ
- 大文字小文字は区別されず、出力時は大文字で統一
- 高頻出k-mer除外にはインデックス作成時の全テーブルスキャンが必要
- GINのk-merインデックスが存在するテーブルはINSERT/UPDATE/DELETE操作が制限される

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。

## 貢献

バグレポートや機能要求は、GitHubのIssueページでお知らせください。