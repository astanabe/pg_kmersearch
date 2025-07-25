# DNA2/DNA4型のBYTEA変換関数開発計画

## 概要

DNA2型およびDNA4型をpgcryptoの`digest()`関数でSHA256ハッシュ化するために、効率的なBYTEA変換関数を作成する。

## 現状分析

### 既存の関数

- `kmersearch_dna2_send(DNA2)`: BYTEA型を返すが、PostgreSQLバイナリプロトコル形式（ビット長メタデータ付き）
- `kmersearch_dna4_send(DNA4)`: BYTEA型を返すが、PostgreSQLバイナリプロトコル形式（ビット長メタデータ付き）

### 内部データ構造

- DNA2/DNA4型は`VarBit`型をベースとして実装
- `VARBITLEN(dna)`: ビット長を取得
- `VARBITS(dna)`: 実際のビットデータへのポインタ
- DNA2: 2ビット/塩基（A=00, C=01, G=10, T=11）
- DNA4: 4ビット/塩基（縮退塩基コード対応）

## 開発対象関数

### 1. kmersearch_dna2_to_bytea(DNA2)

**目的**: DNA2型をハッシュ化に適したBYTEA型として返す（パディング衝突回避）

**実装仕様**:
- 入力: DNA2型
- 出力: BYTEA型（ビット長[4バイト] + ビットデータ）
- 処理: ネットワークバイトオーダーでビット長を格納後、`VARBITS(dna)`からバイトデータを抽出
- パディング衝突対策: ビット長情報により一意性を保証

### 2. kmersearch_dna4_to_bytea(DNA4)

**目的**: DNA4型をハッシュ化に適したBYTEA型として返す

**実装仕様**:
- 入力: DNA4型
- 出力: BYTEA型（ビット長[4バイト] + ビットデータ）
- 処理: ネットワークバイトオーダーでビット長を格納後、`VARBITS(dna)`からバイトデータを抽出
- 一貫性: DNA2型との形式統一

## 実装計画

### ステップ1: C関数の実装

**ファイル**: `kmersearch_datatype.c`

```c
/*
 * DNA2 to BYTEA conversion (bit length + bit data for hash uniqueness)
 */
Datum
kmersearch_dna2_to_bytea(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int byte_len = (bit_len + 7) / 8;
    bytea *result;
    int32 net_bit_len;
    
    /* Allocate space for bit length (4 bytes) + bit data */
    result = (bytea *) palloc(VARHDRSZ + 4 + byte_len);
    SET_VARSIZE(result, VARHDRSZ + 4 + byte_len);
    
    /* Store bit length in network byte order for consistency */
    net_bit_len = htonl(bit_len);
    memcpy(VARDATA(result), &net_bit_len, 4);
    
    /* Store bit data */
    if (byte_len > 0) {
        memcpy(VARDATA(result) + 4, VARBITS(dna), byte_len);
    }
    
    PG_RETURN_BYTEA_P(result);
}

/*
 * DNA4 to BYTEA conversion (bit length + bit data for consistency)
 */
Datum
kmersearch_dna4_to_bytea(PG_FUNCTION_ARGS)
{
    VarBit *dna = PG_GETARG_VARBIT_P(0);
    int32 bit_len = VARBITLEN(dna);
    int byte_len = (bit_len + 7) / 8;
    bytea *result;
    int32 net_bit_len;
    
    /* Allocate space for bit length (4 bytes) + bit data */
    result = (bytea *) palloc(VARHDRSZ + 4 + byte_len);
    SET_VARSIZE(result, VARHDRSZ + 4 + byte_len);
    
    /* Store bit length in network byte order for consistency */
    net_bit_len = htonl(bit_len);
    memcpy(VARDATA(result), &net_bit_len, 4);
    
    /* Store bit data */
    if (byte_len > 0) {
        memcpy(VARDATA(result) + 4, VARBITS(dna), byte_len);
    }
    
    PG_RETURN_BYTEA_P(result);
}
```

### ステップ2: ヘッダファイルの更新

**ファイル**: `kmersearch.h`

関数宣言を追加:
```c
PG_FUNCTION_INFO_V1(kmersearch_dna2_to_bytea);
PG_FUNCTION_INFO_V1(kmersearch_dna4_to_bytea);
```

### ステップ3: SQL関数の定義

**ファイル**: `pg_kmersearch--1.0.sql`

```sql
-- DNA2 to BYTEA conversion function
CREATE FUNCTION kmersearch_dna2_to_bytea(DNA2) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna2_to_bytea'
    LANGUAGE C IMMUTABLE STRICT;

-- DNA4 to BYTEA conversion function  
CREATE FUNCTION kmersearch_dna4_to_bytea(DNA4) RETURNS bytea
    AS 'MODULE_PATHNAME', 'kmersearch_dna4_to_bytea'
    LANGUAGE C IMMUTABLE STRICT;
```

### ステップ4: テストケースの作成

**新規ファイル**: `sql/13_bytea_conversion.sql`

```sql
-- Test BYTEA conversion functions
SELECT kmersearch_dna2_to_bytea('ATCG'::DNA2);
SELECT kmersearch_dna4_to_bytea('ATCGMRWSYKVHDBN'::DNA4);

-- Test with pgcrypto digest
SELECT digest(kmersearch_dna2_to_bytea('ATCG'::DNA2), 'sha256');
SELECT digest(kmersearch_dna4_to_bytea('ATCGN'::DNA4), 'sha256');

-- Critical padding collision test for DNA2
SELECT digest(kmersearch_dna2_to_bytea('ATG'::DNA2), 'sha256') as hash_atg;
SELECT digest(kmersearch_dna2_to_bytea('ATGA'::DNA2), 'sha256') as hash_atga;

-- Verify different hashes for padding boundary cases
SELECT 
    digest(kmersearch_dna2_to_bytea('A'::DNA2), 'sha256') as hash_1char,
    digest(kmersearch_dna2_to_bytea('AA'::DNA2), 'sha256') as hash_2char,
    digest(kmersearch_dna2_to_bytea('AAA'::DNA2), 'sha256') as hash_3char,
    digest(kmersearch_dna2_to_bytea('AAAA'::DNA2), 'sha256') as hash_4char;

-- Compare with existing send functions (should be different due to protocol headers)
SELECT digest(kmersearch_dna2_send('ATCG'::DNA2), 'sha256') as send_hash;
SELECT digest(kmersearch_dna2_to_bytea('ATCG'::DNA2), 'sha256') as tobytea_hash;
```

**対応する期待結果ファイル**: `expected/13_bytea_conversion.out`

### ステップ5: ドキュメントの更新

**ファイル**: `doc/pg_kmersearch.en.md`, `doc/pg_kmersearch.ja.md`

新しい関数の説明を追加:
- 関数の目的と使用方法
- pgcryptoとの連携例
- 既存のsend関数との違い

## パフォーマンス考慮事項

### 利点
1. **一意性保証**: ビット長情報によりパディング衝突を完全回避
2. **効率的**: 直接メモリコピーによる高速処理  
3. **ハッシュ最適化**: digest()関数との組み合わせに最適
4. **プラットフォーム独立**: ネットワークバイトオーダーによる一貫性

### 既存send関数との比較
- `send`関数: 4バイトのビット長メタデータ + データ + プロトコルヘッダ
- `to_bytea`関数: 4バイトのビット長メタデータ + データ

`to_bytea`関数は`send`関数とほぼ同じデータ構造ですが、PostgreSQLプロトコルヘッダが不要なため、よりクリーンな形式でハッシュ化できます。

### パディング衝突問題の解決
**DNS2型での問題**:
- `ATG` (6ビット): `001011xx` (x=パディング0)
- `ATGA` (8ビット): `00101100`

**解決後**:
- `ATG`: `[6バイト長][001011xx]` → 異なるハッシュ
- `ATGA`: `[8バイト長][00101100]` → 異なるハッシュ

## 実装上の注意点

1. **メモリ管理**: `palloc()`による適切なメモリ割り当て（4バイト分の追加領域確保）
2. **エラーハンドリング**: NULL入力および不正データの処理
3. **バイトオーダー**: `htonl()`によるネットワークバイトオーダー変換で一貫性確保
4. **パディング衝突回避**: ビット長情報により一意性を完全保証
5. **テスト**: パディング境界での動作確認（例：6ビット vs 8ビットシーケンス）

## 使用例

```sql
-- 基本的なハッシュ化
SELECT digest(kmersearch_dna2_to_bytea('ATCGATCG'::DNA2), 'sha256');

-- 複数シーケンスの比較
WITH sequences AS (
    SELECT 'ATCG'::DNA2 as seq1, 'GCTA'::DNA2 as seq2
)
SELECT 
    digest(kmersearch_dna2_to_bytea(seq1), 'sha256') as hash1,
    digest(kmersearch_dna2_to_bytea(seq2), 'sha256') as hash2
FROM sequences;
```

この開発計画により、効率的でハッシュ化に最適化されたBYTEA変換関数を提供できる。