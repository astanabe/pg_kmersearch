# GUC変数使用法の問題と修正計画

## 概要
pg_kmersearch拡張におけるGUC変数の使用法を調査し、不正確な使用パターンや潜在的な問題を特定した。以下に発見された問題と修正計画を記載する。

## 発見された問題

### 1. 【高優先度】 kmersearch_freq.c の型変換エラー (line 474)
**問題の内容:**
```c
rate_threshold = (int)(actual_row_count * kmersearch_max_appearance_rate);
```
`kmersearch_max_appearance_rate` は `double` 型（0.0-1.0の範囲）なので、この計算は意図された動作だが、0.0-1.0の範囲外の値が設定された場合の処理が不適切。

**修正計画:**
- **ファイル:** `kmersearch_freq.c`
- **行:** 474
- **修正内容:** 
  - `kmersearch_max_appearance_rate` の範囲検証を強化
  - 計算結果が0になった場合の適切なエラーハンドリングを追加
  - 型変換の安全性を向上

### 2. 【高優先度】 occurrence count の上限チェック不正確 (kmersearch.c)
**問題の内容:**
```c
if (current_count > (1 << kmersearch_occur_bitlen))
```
これは `2^kmersearch_occur_bitlen` と比較しているが、実際には `2^kmersearch_occur_bitlen - 1` が最大値のはず。

**修正計画:**
- **ファイル:** `kmersearch.c`
- **行:** 815, 959
- **修正内容:**
  ```c
  // 修正前
  if (current_count > (1 << kmersearch_occur_bitlen))
  
  // 修正後
  if (current_count >= (1 << kmersearch_occur_bitlen))
  ```
  または
  ```c
  if (current_count > ((1 << kmersearch_occur_bitlen) - 1))
  ```

### 3. 【中優先度】 GUC変数の assign hook の不統一
**問題の内容:**
一部のGUC変数にはassign hookが定義されていないため、値の変更時に適切な検証やキャッシュ更新が行われない。

**修正計画:**
- **ファイル:** `kmersearch.c`
- **対象変数:**
  - `kmersearch_preclude_highfreq_kmer` (line 527-528)
  - `kmersearch_force_use_parallel_highfreq_kmer_cache` (line 538-539)
  - `kmersearch_actual_min_score_cache_max_entries` (line 564-565)
  - `kmersearch_highfreq_kmer_cache_load_batch_size` (line 577-578)
- **修正内容:**
  - 各変数に対応するassign hookを実装
  - 値の変更時にキャッシュの整合性を確保

### 4. 【低優先度】 外部宣言の不足
**問題の内容:**
`kmersearch_force_use_parallel_highfreq_kmer_cache` の外部宣言が `kmersearch_cache.c` にもあるが、ヘッダーファイルに統一されていない。

**修正計画:**
- **ファイル:** `kmersearch_cache.c`
- **行:** 31
- **修正内容:**
  ```c
  // 削除対象
  bool kmersearch_force_use_parallel_highfreq_kmer_cache = false;
  ```
  この行を削除し、ヘッダーファイルの宣言のみを使用する

## 修正の優先順位と実装順序

### Phase 1: 高優先度修正（即座に実装）
1. `kmersearch_freq.c` line 474 の型変換とエラーハンドリング修正
2. `kmersearch.c` line 815, 959 の occurrence count 上限チェック修正

### Phase 2: 中優先度修正（次回リリース）
1. 欠落しているassign hookの実装
2. キャッシュ整合性の確保

### Phase 3: 低優先度修正（コードクリーンアップ）
1. 外部宣言の統一
2. コード品質の向上

## 実装時の注意点

### 一般的な注意事項
- 全てのGUC変数の変更はPostgreSQLセッションの再起動により初期化される
- assign hookの実装時は、既存のキャッシュとの整合性を確保する
- 型変換時は常に範囲チェックを行う

### テスト要件
- 各修正後は `make installcheck` で回帰テストを実行
- 特に occurrence count の上限チェック修正後は、境界値テストを追加
- GUC変数の動的変更のテストを追加

## コードレビューのポイント
- 型の安全性（特に整数オーバーフロー）
- エラーハンドリングの妥当性
- キャッシュ整合性の確保
- PostgreSQL標準に準拠したassign hook実装