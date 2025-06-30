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
-- k=8でGINインデックスを作成（8-merを使用）
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq) WITH (k = 8);

-- k-mer検索（最小64塩基のクエリが必要）
SELECT * FROM sequences 
WHERE dna_seq LIKE 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA' 
LIMIT 10;

-- 出現回数ビット長の設定（デフォルト8ビット）
SELECT set_kmersearch_occur_bitlen(12); -- 12ビットに変更（最大4095回の出現をカウント）

-- 縮重コードを含むクエリでの検索
SELECT * FROM degenerate_sequences 
WHERE dna_seq LIKE 'ATCGATCGNNATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG' 
LIMIT 5;
```

### 高頻出k-mer除外機能

インデックスのパフォーマンス向上のため、高頻出k-merの自動除外機能を搭載：

```sql
-- 除外パラメータの設定（インデックス作成前）
SET kmersearch.max_appearance_rate = 0.05;  -- デフォルト: 5%の最大出現率
SET kmersearch.max_appearance_nrow = 1000;  -- デフォルト: 0（無効）

-- 頻度解析付きインデックス作成
CREATE INDEX sequences_kmer_idx ON sequences USING gin (dna_seq) WITH (k = 8);

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
SELECT show_kmersearch_min_score();

-- 自動スコア調整による検索
-- クエリに除外k-merが3個含まれ、min_score=50の場合、
-- そのクエリでは実際の閾値は47に調整される
SELECT * FROM sequences 
WHERE dna_seq LIKE 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA'
ORDER BY score DESC LIMIT 10;
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

- クエリ配列は最小64塩基が必要
- 縮重コード展開は10組み合わせまで（超過時はスキップ）
- 出現回数は設定可能ビット長の最大値でキャップ
- 大文字小文字は区別されず、出力時は大文字で統一
- 高頻出k-mer除外にはインデックス作成時の全テーブルスキャンが必要
- GINのk-merインデックスが存在するテーブルはINSERT/UPDATE/DELETE操作が制限される

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。

## 貢献

バグレポートや機能要求は、GitHubのIssueページでお知らせください。