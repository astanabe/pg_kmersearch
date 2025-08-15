# 高頻度k-merキャッシュ関数のuintkey対応改訂計画

## 概要
高頻度k-merキャッシュ関連関数を、`devplan_uintkey_2.md`で実装したuintkey形式（オカレンスカウント付き）に対応させる改訂計画です。

## 対象関数
- `kmersearch_highfreq_kmer_cache_load()`
- `kmersearch_highfreq_kmer_cache_free()`
- `kmersearch_highfreq_kmer_cache_free_all()`
- `kmersearch_parallel_highfreq_kmer_cache_load()`
- `kmersearch_parallel_highfreq_kmer_cache_free()`
- `kmersearch_parallel_highfreq_kmer_cache_free_all()`

## 現在の実装の問題点

### 1. データ取得時の不整合
現在の実装では、`kmersearch_highfreq_kmer`テーブルから`kmer2_as_uint`カラムを読み取っていますが、このカラムは`uintkey`に変更されています。

### 2. キャッシュ構造の不適合
- 現在: `kmer2_as_uint`（オカレンスカウントなし）をキーとしてハッシュテーブルに格納
- 改訂後: `uintkey`（オカレンスカウント付き）をキーとして格納する必要がある

### 3. 検索ロジックの不整合
高頻度k-mer判定時に`kmer2_as_uint`形式で検索しているが、`uintkey`形式での検索に変更が必要。

## 改訂内容

### 1. kmersearch_cache.c の修正

#### 1.1 kmersearch_highfreq_kmer_cache_load_internal()
```c
// 現在のSQL文（857-865行目）:
"SELECT COUNT(DISTINCT hkm.kmer2_as_uint) FROM kmersearch_highfreq_kmer hkm"
// ↓
"SELECT COUNT(DISTINCT hkm.uintkey) FROM kmersearch_highfreq_kmer hkm"

// 現在のSQL文（991-1000行目）:
"SELECT DISTINCT hkm.kmer2_as_uint FROM kmersearch_highfreq_kmer hkm"
// ↓
"SELECT DISTINCT hkm.uintkey FROM kmersearch_highfreq_kmer hkm"

// 変数名の変更:
uint64 kmer2_as_uint; → uint64 uintkey;

// ハッシュテーブルへの挿入（1066-1074行目）:
entry->hash_value = kmer2_as_uint; → entry->hash_value = uintkey;
```

#### 1.2 kmersearch_parallel_highfreq_kmer_cache_load_internal()
```c
// 同様の変更を並列キャッシュ版にも適用
// SQL文のカラム名変更: kmer2_as_uint → uintkey
// 変数名の変更: kmer2_as_uint → uintkey
```

#### 1.3 キャッシュ検索関数の更新
```c
// kmersearch_lookup_kmer2_as_uint_in_global_cache() → kmersearch_lookup_uintkey_in_global_cache()
// kmersearch_lookup_kmer2_as_uint_in_parallel_cache() → kmersearch_lookup_uintkey_in_parallel_cache()
// （既に実装済み）
```

### 2. キャッシュ構造体の考慮事項

#### 2.1 HighfreqKmerHashEntry構造体
```c
typedef struct HighfreqKmerHashEntry {
    uint64 hash_value;  // uintkeyの値を格納（オカレンスカウント付き）
    VarBit *kmer_key;   // 未使用（NULLを設定）
} HighfreqKmerHashEntry;
```

現在の構造体はそのまま使用可能ですが、`hash_value`フィールドには`uintkey`値を格納します。

#### 2.2 メモリ使用量への影響
- **増加要因**: 同一k-merでも異なるオカレンスカウントを持つ複数のuintkeyがキャッシュされる
- **対策**: バッチサイズ（`kmersearch_highfreq_kmer_cache_load_batch_size`）の調整が必要な場合がある

### 3. 高頻度判定ロジックの変更

#### 3.1 現在の判定方法
```c
// k-merのビット表現をuint化してキャッシュ内を検索
uint64 kmer2_as_uint = extract_kmer_as_uint(kmer);
bool is_highfreq = lookup_in_cache(kmer2_as_uint);
```

#### 3.2 改訂後の判定方法
```c
// uintkey（オカレンスカウント付き）でキャッシュ内を検索
uint64 uintkey = extract_uintkey(kmer, occurrence_count);
bool is_highfreq = lookup_in_cache(uintkey);
```

**重要な変更点**:
- 同一k-merでも、オカレンスカウントが異なれば別のエントリとして扱われる
- 一部のオカレンスカウント値のみが高頻度と判定される可能性がある

### 4. GINインデックス構築・検索時の影響

#### 4.1 構築時（kmersearch_gin_extract_value）
- 現在: `kmer2_as_uint`で高頻度判定し、該当するk-merを除外
- 改訂後: `uintkey`で高頻度判定し、該当するuintkeyを除外
- **注意**: オカレンスカウント別に除外判定が行われる

#### 4.2 検索時（kmersearch_gin_extract_query）
- 同様にuintkey形式での高頻度判定に変更が必要

### 5. 実装手順

#### フェーズ1: 基本的な修正
1. SQL文のカラム名変更（`kmer2_as_uint` → `uintkey`）
2. 変数名の更新
3. コメントの更新

#### フェーズ2: 機能検証
1. キャッシュのロード/フリー動作の確認
2. 高頻度k-mer判定の正確性検証
3. メモリ使用量の測定

#### フェーズ3: 最適化
1. バッチサイズの調整
2. ハッシュテーブルサイズの最適化
3. パフォーマンステスト

## テスト計画

### 1. 単体テスト
```sql
-- キャッシュロードテスト
SELECT kmersearch_highfreq_kmer_cache_load('test_table', 'seq_column');

-- キャッシュ内容の確認
SELECT COUNT(*) FROM kmersearch_highfreq_kmer 
WHERE table_oid = 'test_table'::regclass::oid;

-- 高頻度判定テスト（異なるオカレンスカウント）
-- TODO: 具体的なテストケースを作成
```

### 2. 統合テスト
- GINインデックス構築時の高頻度k-mer除外動作
- 検索時のフィルタリング動作
- 並列キャッシュの動作確認

### 3. パフォーマンステスト
- キャッシュロード時間の測定
- メモリ使用量の監視
- 検索性能への影響評価

## リスクと対策

### 1. メモリ使用量の増加
**リスク**: uintkey形式では同一k-merの複数バリエーションが格納される
**対策**: 
- キャッシュサイズの上限設定
- 最も頻度の高いuintkeyのみをキャッシュする選択的ロード

### 2. 検索精度の変化
**リスク**: オカレンスカウント別の除外により、検索結果が変わる可能性
**対策**:
- 閾値の調整機能を提供
- 検索結果の比較検証を実施

### 3. 後方互換性
**リスク**: 既存のキャッシュデータとの非互換
**対策**: 
- 新規インストールのみ対応（既存データの移行は不要）
- キャッシュの再構築を必須とする

## 期待される効果

1. **より精密な高頻度フィルタリング**: オカレンスカウント情報を活用した細かい制御
2. **GINインデックスとの整合性向上**: uintkey形式の統一による一貫性
3. **将来の拡張性**: オカレンスカウント情報を活用した高度な最適化が可能

## 実装優先順位

1. **必須**: kmersearch_highfreq_kmer_cache_load_internal()の修正
2. **必須**: kmersearch_parallel_highfreq_kmer_cache_load_internal()の修正
3. **推奨**: バッチサイズとメモリ管理の最適化
4. **オプション**: 選択的キャッシュロード機能の追加

## 補足事項

- デバッグログ（DEBUG1レベル）は変更後も維持し、uintkey値を適切に出力
- エラーハンドリングは既存のロジックを維持
- トランザクション管理やSPI接続は既存の方式を踏襲