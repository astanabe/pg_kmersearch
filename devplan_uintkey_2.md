# kmersearch_perform_highfreq_analysis() uintkey対応改訂計画

## 概要
`kmersearch_perform_highfreq_analysis()`関数を、オカレンスカウント付きuintkey形式に対応させる改訂計画です。

## 現在の実装

### 処理フロー
1. DNA2/DNA4データからkmer2（オカレンスカウントなし）を抽出してuint化（kmer2_as_uint）
   - `kmersearch_extract_dna2_kmer2_as_uint_direct()`
   - `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()`
2. 各kmer2_as_uintの出現行数を並列処理でカウント
3. 閾値を超えた高頻度k-merを`kmersearch_highfreq_kmer`テーブルの`kmer2_as_uint`カラムに保存

### データベーススキーマ
```sql
CREATE TABLE kmersearch_highfreq_kmer (
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    kmer2_as_uint bigint NOT NULL,  -- オカレンスカウントなしのk-mer
    detection_reason text,
    created_at timestamp with time zone DEFAULT now(),
    PRIMARY KEY (table_oid, column_name, kmer2_as_uint)
);
```

## 改訂後の実装

### 処理フロー
1. DNA2/DNA4データからkmer2を抽出してオカレンスカウントを追加してuint化（uintkey）
   - `kmersearch_extract_uintkey_from_dna2()`
   - `kmersearch_extract_uintkey_from_dna4()`
   - **注意**: これらの関数が返す配列には同一uintkeyの重複はないため、重複除外処理は不要
2. 各uintkeyの出現行数を並列処理でカウント
3. 閾値を超えた高頻度uintkeyを`kmersearch_highfreq_kmer`テーブルの`uintkey`カラムに保存

### データベーススキーマ変更
```sql
CREATE TABLE kmersearch_highfreq_kmer (
    table_oid oid NOT NULL,
    column_name name NOT NULL,
    uintkey bigint NOT NULL,  -- オカレンスカウント付きk-mer
    detection_reason text,
    created_at timestamp with time zone DEFAULT now(),
    PRIMARY KEY (table_oid, column_name, uintkey)
);
```

## 実装詳細

### 1. SQLテーブル定義の変更
- `pg_kmersearch--1.0.sql`内の`kmersearch_highfreq_kmer`テーブル定義を修正
  - カラム名: `kmer2_as_uint` → `uintkey`

### 2. kmersearch_freq.c の修正

#### 2.1 kmersearch_update_kmer_counts_in_dshash() の変更
```c
static void
kmersearch_update_kmer_counts_in_dshash(Datum sequence_datum, int kmer_size, 
                                       dshash_table *hash, Oid column_type_oid)
{
    // 現在の実装:
    // - kmersearch_extract_dna2_kmer2_as_uint_direct()
    // - kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()
    
    // 改訂後:
    // - kmersearch_extract_uintkey_from_dna2()
    // - kmersearch_extract_uintkey_from_dna4()
    
    // ローカルのk-mer重複チェック処理は削除可能
    // （extract_uintkey関数が既に重複を除外しているため）
}
```

#### 2.2 kmersearch_insert_kmer2_as_uint_from_dshash() の変更
```c
static void
kmersearch_insert_uintkey_from_dshash(Oid table_oid, const char *column_name, 
                                     int k_size, uint64 threshold_rows)
{
    // 関数名変更: kmer2_as_uint → uintkey
    // SQLクエリのカラム名変更: kmer2_as_uint → uintkey
    
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer "
        "(table_oid, column_name, uintkey, detection_reason) VALUES ");
}
```

#### 2.3 関数呼び出しの更新
- `kmersearch_perform_highfreq_analysis_parallel()`内での関数呼び出しを更新

### 3. dshashエントリ構造の考慮事項

現在のdshashエントリ構造（KmerEntry16/32/64）では、k-merのビット表現のみを保存しています。uintkey形式では、オカレンスカウントがビット表現に含まれるため：

- **同一のk-mer（オカレンスカウントを除く）でも、異なるuintkeyとして扱われる**
- 一部のuintkeyは高頻度と判定され、一部は判定されない可能性がある
- これは設計上の意図された動作

### 4. 並列処理への影響

並列ワーカー（`kmersearch_analysis_worker`）での処理：
- 各ワーカーは独立してuintkeyを抽出
- dshashテーブルでuintkey単位でカウント
- 重複除外のためのローカルハッシュテーブルは不要になる

### 5. テストと検証

#### 5.1 単体テスト
- uintkey抽出関数の動作確認
- 高頻度判定ロジックの検証
- オカレンスカウント別の判定結果確認

#### 5.2 統合テスト
- 既存の`make installcheck`テストの更新
- 新しいテストケースの追加：
  - 同一k-merで異なるオカレンスカウントの場合
  - 閾値境界でのuintkey判定

#### 5.3 パフォーマンステスト
- uintkey形式での並列処理性能測定
- メモリ使用量の確認

## 実装順序

1. **フェーズ1: 基盤整備**
   - SQLスキーマ変更（`pg_kmersearch--1.0.sql`）
   - 関数名とデータ構造の更新

2. **フェーズ2: コア機能実装**
   - `kmersearch_update_kmer_counts_in_dshash()`の修正
   - `kmersearch_insert_uintkey_from_dshash()`への変更
   - 並列ワーカー関数の更新

3. **フェーズ3: テストと最適化**
   - テストケースの更新
   - パフォーマンス測定と最適化
   - ドキュメント更新

## 注意事項

1. **後方互換性**: 新規インストールのみを想定しているため、既存データとの互換性は考慮不要

2. **メモリ効率**: uintkey形式では、同一k-merでも異なるオカレンスカウントを持つ複数のエントリが存在する可能性があるため、メモリ使用量が増加する可能性がある

3. **検索性能への影響**: 高頻度k-merフィルタリングの精度が変わる可能性があるため、検索性能への影響を評価する必要がある

## 期待される効果

- より精密な高頻度k-merフィルタリング
- オカレンスカウント情報を活用した高度な分析が可能
- GINインデックスとの整合性向上（uintkey形式の統一）