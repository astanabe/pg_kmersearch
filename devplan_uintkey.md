# kmersearch_kmer.c への uintkey 抽出関数追加計画

## 1. 概要

`kmersearch_kmer.c` に、k-mer を uint 型のキーとして抽出する3つの新しい関数を追加する。
これらの関数は、k-mer とオカレンスカウント（同一 k-mer の出現順序）を組み合わせた整数キーを生成する。

**重要**: これらの関数は他のCソースファイル（kmersearch_gin.c、kmersearch_freq.c など）からも利用されることを前提に設計する。
そのため、関数プロトタイプは `kmersearch.h` に宣言し、外部から呼び出し可能な公開APIとして実装する。

## 2. 追加する関数

### 2.1 kmersearch_extract_uintkey_from_dna2()

#### 機能
- DNA2 型データ（VarBit）から k-mer（kmer2）を抽出
- 各 kmer2 のオカレンスカウント（何回目の出現か）を計算
- kmer2 とオカレンスカウントをビット連結して uint に変換

#### 実装方針
- `kmersearch_extract_dna2_kmer2_as_uint_direct()` を参考に実装
- SIMD 化を見据えた構造とする
  - `kmersearch_extract_uintkey_from_dna2_scalar()` を作成
  - `kmersearch_extract_uintkey_from_dna2()` はディスパッチ関数として実装
- k 値は `kmersearch_kmer_size` GUC 変数を使用
- オカレンスカウントのビット長は `kmersearch_occur_bitlen` GUC 変数を使用

#### オカレンスカウントのエンコーディング
- `kmersearch_occur_bitlen` が 8 の場合：
  - オカレンスカウント 1 → 0b00000000
  - オカレンスカウント 256 → 0b11111111
  - 一般式：`encoded = (occurrence - 1) & ((1 << occur_bitlen) - 1)`

#### 出力型の決定
- kmer2 ビット長 + オカレンスカウントビット長の合計により決定：
  - 16 ビット以下 → uint16
  - 32 ビット以下 → uint32
  - 64 ビット以下 → uint64

### 2.2 kmersearch_extract_uintkey_from_dna4()

#### 機能
- DNA4 型データ（VarBit）から k-mer を抽出
- `kmersearch_expand_dna4_kmer2_to_dna2_direct()` で kmer2 に変換
- 各 kmer2 のオカレンスカウントを計算
- kmer2 とオカレンスカウントをビット連結して uint に変換

#### 実装方針
- `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()` を参考に実装
- SIMD 化を見据えた構造とする
  - `kmersearch_extract_uintkey_from_dna4_scalar()` を作成
  - `kmersearch_extract_uintkey_from_dna4()` はディスパッチ関数として実装
- デジェネレート塩基の展開処理を含む
- k 値は `kmersearch_kmer_size` GUC 変数を使用
- オカレンスカウントのビット長は `kmersearch_occur_bitlen` GUC 変数を使用

#### 処理フロー
1. DNA4 データから k-mer を抽出
2. デジェネレート塩基を含む場合は展開
3. 各展開された kmer2 に対してオカレンスカウントを計算
4. kmer2 + オカレンスカウントを uint に変換

### 2.3 kmersearch_extract_uintkey_from_text()

#### 機能
- テキスト文字列を DNA4 型に変換
- `kmersearch_extract_uintkey_from_dna4()` を呼び出して uint キーを取得

#### 実装方針
- `kmersearch_dna4_in()` の文字列→DNA4 変換ロジックを利用
- 変換後は `kmersearch_extract_uintkey_from_dna4()` に委譲
- エラーハンドリングを適切に実装

#### 処理フロー
1. 入力文字列の検証
2. DNA4 VarBit への変換
3. `kmersearch_extract_uintkey_from_dna4()` 呼び出し
4. 結果の返却

## 3. データ構造とメモリ管理

### オカレンスカウント管理
- ハッシュテーブルを使用して各 k-mer の出現回数を追跡
- PostgreSQL の HTAB (dynahash) を使用
- キー：kmer2 のビット列
- 値：現在のオカレンスカウント

### メモリ割り当て
- 結果配列は palloc() で割り当て
- 要素サイズは `kmersearch_get_kmer_uint_size()` で決定
- 最大 k-mer 数を事前計算してメモリを確保

## 4. エラー処理

### 検証項目
- k 値の範囲（4 ≤ k ≤ 32）
- オカレンスカウントビット長の範囲（1 ≤ occur_bitlen ≤ 16）
- 合計ビット長が 64 ビット以下であること
- デジェネレート塩基の展開限界チェック

### エラー処理方針
- 不正な入力に対しては ereport(ERROR) で適切なメッセージを返す
- デバッグ情報は elog(DEBUG2/DEBUG3) で出力

## 5. パフォーマンス考慮事項

### 最適化ポイント
- ビット操作の効率化（BMI2 命令の活用準備）
- ハッシュテーブルの初期サイズ最適化
- メモリアロケーションの最小化
- SIMD 化への準備（データアライメント等）

### ベンチマーク対象
- 異なる k 値での処理速度
- 異なるシーケンス長での処理速度
- デジェネレート塩基を含む/含まない場合の比較

## 6. テスト計画

### ユニットテスト
- 各関数の基本動作確認
- 境界値テスト（k=4, k=32, occur_bitlen=1, occur_bitlen=16）
- オカレンスカウントの正確性確認
- デジェネレート塩基展開の検証

### 統合テスト
- 既存の regression test への追加
- パフォーマンステストの実施
- メモリリークのチェック

## 7. 今後の拡張計画

### SIMD 実装
- AVX2/AVX512 実装（x86_64）
- NEON/SVE 実装（ARM）
- ランタイムディスパッチの実装

### 機能拡張
- バッチ処理の最適化
- キャッシュ機構との統合
- GIN インデックスとの連携

## 8. 実装順序

1. ヘッダファイルへの宣言追加
   - `kmersearch.h` に関数プロトタイプを追加
   - 他のCファイルから利用可能にする

2. スカラー実装の完成
   - `kmersearch_extract_uintkey_from_dna2_scalar()`（static関数）
   - `kmersearch_extract_uintkey_from_dna4_scalar()`（static関数）
   - `kmersearch_extract_uintkey_from_text()`（公開関数）

3. ディスパッチ関数の実装
   - `kmersearch_extract_uintkey_from_dna2()`（公開関数）
   - `kmersearch_extract_uintkey_from_dna4()`（公開関数）

4. テストコードの作成と実行

5. ドキュメント更新

## 9. API設計と他モジュールとの連携

### 公開API（kmersearch.hで宣言）
```c
/* 公開関数 - 他のCファイルから呼び出し可能 */
void kmersearch_extract_uintkey_from_dna2(VarBit *seq, void **output, int *nkeys);
void kmersearch_extract_uintkey_from_dna4(VarBit *seq, void **output, int *nkeys);
void kmersearch_extract_uintkey_from_text(const char *text, void **output, int *nkeys);
```

### 内部実装（static関数）
```c
/* スカラー実装 - kmersearch_kmer.c内部でのみ使用 */
static void kmersearch_extract_uintkey_from_dna2_scalar(VarBit *seq, void **output, int *nkeys);
static void kmersearch_extract_uintkey_from_dna4_scalar(VarBit *seq, void **output, int *nkeys);
```

### 利用想定モジュール
- **kmersearch_gin.c**: GINインデックスのキー生成で利用
- **kmersearch_freq.c**: 頻度解析でのk-mer識別で利用
- **kmersearch_cache.c**: キャッシュキーの生成で利用
- **kmersearch_partition.c**: パーティショニングキーの生成で利用

## 10. 注意事項

- GUC 変数は直接使用し、ローカル変数への代入は避ける
- スタブ実装や不完全な実装は禁止
- 既存の関数と重複する機能を実装しない
- エラーメッセージは明確で具体的にする
- 他モジュールからの呼び出しを考慮し、適切なエラーハンドリングを実装
- 関数の入出力仕様を明確にドキュメント化する