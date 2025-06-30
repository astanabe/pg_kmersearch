# pg_kmersearch

PostgreSQL用のDNA配列型拡張

## 概要

pg_kmersearchは、PostgreSQL用のDNA配列データを効率的に格納・処理するためのカスタムデータ型を提供する拡張です。この拡張は以下の2つのデータ型を実装しています：

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

### 内部実装
- DNA2型：2ビット/文字でエンコード
- DNA4型：4ビット/文字でエンコード
- 両型ともにPostgreSQLの可変長データ型として実装
- バイナリ入出力サポート

## 制限事項

- 現在は基本的なデータ型機能のみ実装
- 比較演算子や検索機能は今後の実装予定
- 大文字小文字は区別されず、出力時は大文字で統一

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。

## 貢献

バグレポートや機能要求は、GitHubのIssueページでお知らせください。