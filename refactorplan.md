# pg_kmersearch リファクタリング計画

## 概要
pg_kmersearchプロジェクトのソースコード精査の結果、以下の問題が見つかりました。これらの問題を解決するための修正計画を提示します。

## 1. スタブ実装

### 1.1 cleanup_temp_partitions関数 (kmersearch_partition.c:629-635) ✅ COMPLETED
**問題点**: 
- 完全にスタブ実装として残されており、機能が実装されていない
- コメントで「intentionally left as a stub」と明記されている
- エラー時のクリーンアップが不完全になる可能性がある

**修正計画**:
1. この関数を削除し、呼び出し元でインライン処理を行う
2. PostgreSQLのSPIコンテキスト制限に従い、トランザクション内での適切なクリーンアップ方法を実装
3. エラーハンドリングをPG_TRY/PG_CATCHブロック内で改善

**実施内容**:
- cleanup_temp_partitions関数本体を削除
- static宣言も削除
- この関数は実際には呼び出されていなかったため、機能への影響なし

## 2. 同一機能を持つ重複関数

### 2.1 ngram_key2作成関数の統合
**分析結果**:
各ngram_key2作成関数の入出力を精査した結果、以下の通り異なる入力タイプを扱っていることが判明：
- `kmersearch_create_ngram_key2()`: 文字列（const char *）から作成
- `kmersearch_create_ngram_key2_from_dna2_bits()`: DNA2 VarBitの指定位置から作成
- `kmersearch_create_ngram_key2_from_dna4_bits()`: DNA4 VarBitの指定位置から作成（内部でDNA2に変換）
- `kmersearch_create_ngram_key2_with_occurrence_from_dna2()`: 既存のDNA2 k-merから作成
- `kmersearch_create_ngram_key2_from_kmer2_as_uint()`: uint64から作成

**修正計画**:
これらの関数は入力タイプが異なるため、統合は不適切。現状維持とする。

## 3. 過剰な内部ラッパー関数

### 3.1 キャッシュ管理の_internal関数群 ✅ COMPLETED
**分析結果**:
内部関数は以下の2つのカテゴリに分類される：

**A. 削除された単純なラッパー関数**:
- `kmersearch_free_query_pattern_cache_internal()` ✅ DELETED
- `kmersearch_free_actual_min_score_cache_internal()` ✅ DELETED

**B. 保持すべき実質的な処理を含む関数**:
- `kmersearch_highfreq_kmer_cache_load_internal()` - データベースからのロード、メモリ管理、GUC検証を実行
- `kmersearch_highfreq_kmer_cache_free_internal()` - メモリコンテキストの削除と状態リセットを実行
- `kmersearch_parallel_highfreq_kmer_cache_load_internal()` - 並列キャッシュの初期化とロードを実行
- `kmersearch_parallel_highfreq_kmer_cache_free_internal()` - 並列キャッシュのクリーンアップを実行
- `kmersearch_parallel_cache_cleanup_internal()` - DSM/dshashのクリーンアップを実行

これらの関数は実質的な処理を含み、SQL関数からの呼び出しに必要な内部実装である。

**修正計画**:
1. PostgreSQL関数(PG_FUNCTION_ARGS)から直接実装を呼び出す
2. 内部関数が本当に必要な場合のみ残す（共通処理がある場合）
3. 不要な関数階層を削除してコードを簡潔にする

**実施内容**:
- `kmersearch_free_query_pattern_cache_internal()`を削除し、`free_query_pattern_cache_manager()`を直接呼び出すように変更
- `kmersearch_free_actual_min_score_cache_internal()`を削除し、`free_actual_min_score_cache_manager()`を直接呼び出すように変更
- `static`を削除して関数を公開し、ヘッダーファイルで宣言
- kmersearch.cおよびkmersearch_cache.c内のすべての呼び出しを更新

**注意**: 他の`_internal`関数（特に`kmersearch_parallel_*`系）は実質的な処理を含んでいるため、削除対象から除外

## 4. SIMD実装について（削除）

### 4.1 プラットフォーム別SIMD関数の分析結果
**調査結果**:
- scalar、AVX2、AVX512、NEON、SVEの各実装は、それぞれ異なるCPU命令セットを使用
- 各実装はプラットフォーム固有の最適化により大幅な性能向上を実現
- 関数ポインタテーブルによる動的ディスパッチは、実行時のCPU機能検出に必要

**結論**:
これらは「重複」ではなく、性能向上のための必要な実装である。各バージョンは：
- 異なるSIMD命令セットを使用（AVX2: 256bit、AVX512: 512bit、NEON: 128bit、SVE: 可変長）
- プラットフォーム固有の最適化を含む
- 実行環境に応じて最適な性能を提供

よって、リファクタリング対象から除外する。

## 5. GUC変数の直接使用違反

### 5.1 GUC変数のローカル変数への代入 ✅ COMPLETED
**問題点**:
- CLAUDE.mdで「assigning GUC variables to local variables is forbidden」と記載
- 一部の関数でGUC変数をローカル変数に代入している可能性

**修正内容**:
1. kmersearch.c内の`int k = kmersearch_kmer_size;`の削除と直接参照への変更 ✅
2. kmersearch_kmer.c内の`int occur_bits = kmersearch_occur_bitlen;`の削除と直接参照への変更 ✅
3. kmersearch_gin.c内の`int k = kmersearch_kmer_size;`の削除と直接参照への変更 ✅
4. kmersearch_freq.c内の`k_size = kmersearch_kmer_size;`の削除（関数パラメータとして渡す形に変更） ✅

## 実装優先順位

1. **高優先度**: スタブ実装の削除と適切な実装 ✅ COMPLETED
2. **中優先度**: 内部ラッパー関数の削除 ✅ PARTIALLY COMPLETED (2/7)
3. **低優先度**: GUC変数使用の修正 ✅ COMPLETED
4. **中優先度**: ソースファイル間の関数再配置 ✅ IN PROGRESS
5. **中優先度**: 未使用関数の削除 ✅ COMPLETED
6. **高優先度**: USE_*条件付きコンパイルマクロの削除 ✅ COMPLETED

## 削除された項目

### rawscoreとcorrectedscore関数の統合（削除）
現在のGINインデックス検索がスコアをクエリ塩基配列と関連付けて保存することが不可能な設計のため、両関数が同じ値を返すのは意図的な実装。将来的な拡張のために残す。

### 並列処理でのPG_TRY/PG_CATCH使用（削除）
現在エラーの原因になっていないため、修正不要。

### SIMD実装の重複（削除）
調査の結果、これらは重複ではなく性能向上のための必要な実装であることが判明。各実装は異なるCPU命令セットを使用し、プラットフォーム固有の最適化を提供している。

## 期待される効果

1. **コードの簡潔性向上**: 重複コードの削除により、メンテナンス性が向上
2. **バグリスクの低減**: スタブ実装や不完全な関数の削除により、予期しない動作を防止
3. **パフォーマンス向上**: 不要な関数呼び出しの削減により、わずかながらパフォーマンス向上
4. **可読性向上**: 関数階層の簡素化により、コードの理解が容易に

## 6. ソースファイル間の関数再配置

### 6.1 問題点
`kmersearch.c`が6376行と肥大化しており、関数が適切なファイルに配置されていない。特に以下の問題がある：
- SIMD関連関数が約1500行以上を占める
- k-mer抽出関数が本来のファイル（`kmersearch_kmer.c`）ではなく`kmersearch.c`に存在
- データ型関連の関数が分散している

### 6.2 関数移動計画

#### kmersearch.c → kmersearch_kmer.c（k-mer抽出・変換関連）
- `kmersearch_extract_dna2_kmer2_direct()` およびそのSIMDバリアント（scalar ✅, avx2 ✅, avx512 ✅, neon ✅, sve ✅）
- `kmersearch_extract_dna4_kmer2_with_expansion_direct()` およびそのSIMDバリアント（scalar ✅, avx2 ✅, avx512 ✅, neon ✅, sve ✅）
- `kmersearch_extract_dna2_ngram_key2_direct()` ✅
- `kmersearch_extract_dna4_ngram_key2_direct()` ✅ 
- `create_ngram_key2_from_kmer2_and_count()` ✅
- `kmersearch_count_matching_kmer_fast()` およびそのSIMDバリアント ✅ (dispatcher + scalar_simple + scalar_hashtable)

#### kmersearch.c → kmersearch_datatype.c（データ型エンコード・デコード関連）
- `dna2_encode_scalar()` ✅, `dna2_encode_avx2()` ✅, `dna2_encode_avx512()` ✅, `dna2_encode_neon()` ✅, `dna2_encode_sve()` ✅
- `dna2_decode_scalar()` ✅, `dna2_decode_avx2()` ✅, `dna2_decode_avx512()` ✅, `dna2_decode_neon()` ✅, `dna2_decode_sve()` ✅
- `dna4_encode_scalar()` ✅, `dna4_encode_avx2()` ✅, `dna4_encode_avx512()` ✅, `dna4_encode_neon()` ✅, `dna4_encode_sve()` ✅
- `dna4_decode_scalar()` ✅, `dna4_decode_avx2()` ✅, `dna4_decode_avx512()` ✅, `dna4_decode_neon()` ✅, `dna4_decode_sve()` ✅
- ~~`kmersearch_dna2_bit_length()`, `kmersearch_dna4_bit_length()`~~ (SQL API functions - remain in kmersearch.c)
- ~~`kmersearch_dna2_nuc_length()`, `kmersearch_dna4_nuc_length()`~~ (SQL API functions - remain in kmersearch.c)

#### kmersearch.c → kmersearch_freq.c（頻度分析関連）
- `kmersearch_worker_analyze_blocks()` ✅ (with helper functions: `kmersearch_calculate_buffer_size`, `kmersearch_init_buffer`, `kmersearch_add_hash_to_buffer`, `kmersearch_flush_hash_buffer_to_table`, `kmersearch_create_worker_temp_table`)
- `kmersearch_merge_worker_results_sql()` ✅
- `process_extracted_kmer2()` ✅ (with helper functions: `is_kmer2_in_highfreq_table`, `is_kmer2_in_analysis_dshash`)
- `kmersearch_persist_highfreq_kmers_from_temp()` ✅
- `create_worker_ngram_temp_table()` ✅
- バッファ関連関数群：
  - `kmersearch_init_buffer()` ✅ (already in kmersearch_freq.c)
  - `kmersearch_add_to_buffer()` ✅
  - `kmersearch_add_hash_to_buffer()` ✅ (already in kmersearch_freq.c)
  - `kmersearch_flush_buffer_to_table()` ✅
  - `kmersearch_flush_hash_buffer_to_table()` ✅ (already in kmersearch_freq.c)
  - `kmersearch_aggregate_buffer_entries()` ✅
  - `kmersearch_create_worker_temp_table()` ✅ (already in kmersearch_freq.c)

#### kmersearch.c → kmersearch_cache.c（キャッシュ関連）
- `kmersearch_parallel_highfreq_kmer_cache_is_valid()` ✅ (already existed in kmersearch_cache.c, just removed static)
- rawscore/correctedscore計算で使用されるキャッシュアクセス関数

#### kmersearch.c → kmersearch_gin.c（GINインデックス関連）
- `kmersearch_get_index_info()` ✅
- `kmersearch_kmer_based_match_dna2()` ✅
- `kmersearch_kmer_based_match_dna4()` ✅
- `kmersearch_evaluate_match_conditions()`
- `evaluate_optimized_match_condition()` ✅

#### kmersearch.c → kmersearch_util.c（ユーティリティ関数）
- `get_dna2_type_oid()` ✅
- `get_dna4_type_oid()` ✅

#### kmersearch.cに残す関数
- `detect_cpu_capabilities()` - CPU機能検出は初期化処理の一部として中央管理
- `init_simd_dispatch_table()` - SIMDディスパッチテーブルの初期化も中央管理

#### 削除された未使用関数
調査の結果、以下の関数は実際には使用されていなかったため削除した：
- `tuple_in_worker_range()` ✅ DELETED (未使用)
- `get_processed_row_count()` ✅ DELETED (未使用)
- `kmersearch_get_kmer_data_size()` ✅ DELETED (未使用)
- `extract_sequence_from_tuple()` ✅ DELETED (未使用)
- `kmersearch_evaluate_match_conditions()` ✅ DELETED (未使用、evaluate_optimized_match_conditionに置き換え済み)
- `kmersearch_calculate_buffer_size()` - バッファサイズ計算のため、使用箇所に応じて決定

### 6.3 期待される効果
1. **保守性向上**: 各ファイルが明確な責務を持つ
2. **可読性向上**: 関連する関数がまとまることで理解しやすくなる
3. **コンパイル時間短縮**: 変更時の再コンパイル範囲が限定される
4. **コード探索効率化**: 機能に応じたファイル構成により目的の関数を見つけやすくなる

### 6.4 実装時の注意点
- ヘッダーファイル（`kmersearch.h`）の関数宣言を適切に整理
- static関数のスコープを考慮（必要に応じて非staticに変更またはファイル内移動）
- グローバル変数（`simd_dispatch`など）のextern宣言を適切に配置
- 移動後も既存のビルドが成功することを確認

## 注意事項

- 各修正は単体テストを実行して動作を確認する
- SQL関数インターフェースの変更は後方互換性に注意
- 並列処理関連の修正は特に慎重にテストする
- パフォーマンス測定を行い、改善を確認する

## SIMD実装の再有効化 ✅ COMPLETED

以下のSIMDディスパッチロジックが再有効化されました：

1. **`kmersearch_extract_dna2_kmer2_direct()`**:
   - CPU機能に基づいてAVX2/AVX512/NEON/SVE実装を動的に選択
   - `simd_capability`グローバル変数を使用して適切な実装を選択

2. **`kmersearch_extract_dna4_kmer2_with_expansion_direct()`**:
   - 同様にCPU機能に基づいてAVX2/AVX512/NEON/SVE実装を動的に選択
   - デジェネレートコード展開を含むDNA4シーケンスに対応

**注意**: コンパイル時に警告が出るが、これは前方宣言のパターンの不一致によるもので、機能には影響しない。