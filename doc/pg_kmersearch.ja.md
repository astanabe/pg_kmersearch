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
- **スコア計算関数**: 個別配列のスコア算出用の`kmersearch_rawscore()`と`kmersearch_correctedscore()`関数

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

3. データベースに拡張をインストール：
```sql
CREATE EXTENSION pg_kmersearch;
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
-- k-merサイズを設定（デフォルト8-mer）
SET kmersearch.kmer_size = 8;

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
| `kmersearch.kmer_size` | 8 | 4-64 | インデックス作成と検索のk-mer長 |
| `kmersearch.occur_bitlen` | 8 | 0-16 | 出現回数格納のビット数 |
| `kmersearch.max_appearance_rate` | 0.05 | 0.0-1.0 | インデックス化するk-merの最大出現率 |
| `kmersearch.max_appearance_nrow` | 0 | 0-∞ | k-merが含まれる最大行数（0=無制限） |
| `kmersearch.min_score` | 1 | 0-∞ | 検索結果の最小類似度スコア |

### 高頻出k-mer除外機能

インデックスのパフォーマンス向上のため、高頻出k-merの自動除外機能を搭載：

```sql
-- 除外パラメータの設定（インデックス作成前）
SET kmersearch.max_appearance_rate = 0.05;  -- デフォルト: 5%の最大出現率
SET kmersearch.max_appearance_nrow = 1000;  -- デフォルト: 0（無効）

-- 頻度解析付きインデックス作成
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq);

-- インデックスの除外k-mer確認
SELECT kmer_key, frequency_count, exclusion_reason 
FROM kmersearch_excluded_kmers 
WHERE index_oid = 'sequences_kmer_idx'::regclass;

-- インデックス統計情報
SELECT total_rows, excluded_kmers_count, max_appearance_rate 
FROM kmersearch_index_info 
WHERE index_oid = 'sequences_kmer_idx'::regclass;
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
- **システムテーブル**: 除外k-merとインデックス統計のメタデータ格納
- バイナリ入出力サポート

### k-mer検索の仕組み
1. **頻度解析**: 高頻出k-mer特定のための全テーブルスキャン
2. **k-mer抽出**: 指定されたk長でスライディングウィンドウ
3. **高頻出フィルタリング**: 出現閾値を超えるk-merの除外
4. **n-gramキー生成**: k-mer + 出現回数をバイナリエンコード
5. **縮重コード処理**: MRWSYKVHDBN を標準塩基に展開
6. **スコアリング**: 共有n-gramキー数による類似度計算

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