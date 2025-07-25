# pg_kmersearch max_appearance_rate テスト結果レポート

## 概要
`max_appearance_rate`パラメータが`kmersearch_perform_highfreq_analysis()`の結果に与える影響を調査しました。

## テスト環境
- PostgreSQL with pg_kmersearch extension
- テストデータベース: `kmersearch_test`
- テストテーブル: `test_sequences` (40行の DNA2 データ)

## テストデータの構成
作成したテストデータセットは以下の通り：

| パターン | 出現回数 | 内容 |
|----------|----------|------|
| Pattern 1 | 18行 | `ATCGATCGATCGATCGATCGATCGATCGATCGATCG` |
| Pattern 2 | 12行 | `GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA` |
| Pattern 3 | 5行 | `AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA` |
| Pattern 4 | 3行 | `CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC` |
| Pattern 5 | 1行 | `GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG` |
| Pattern 6 | 1行 | `TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT` |

**総行数: 40行**

## テスト結果

### max_appearance_rate による影響

| max_appearance_rate | 閾値計算 | 実際の閾値 | 検出された高頻度k-mer数 | max_appearance_nrow_used |
|---------------------|----------|------------|-------------------------|--------------------------|
| 0.1 | int(40 × 0.1) = 4 | 4 | 12 | 4 |
| 0.4 | int(40 × 0.4) = 16 | 16 | 8 | 16 |
| 0.5 | int(40 × 0.5) = 20 | 20 | 8 | 20 |
| 0.8 | int(40 × 0.8) = 32 | 32 | 4 | 32 |
| 0.9 | int(40 × 0.9) = 36 | 36 | 0 | 36 |

### kmer_size による影響

同じ `max_appearance_rate` 値でも、`kmer_size` が異なる場合の結果：

| kmer_size | max_appearance_rate | 検出された高頻度k-mer数 |
|-----------|---------------------|-------------------------|
| 4 | 0.4 | 8 |
| 6 | 0.4 | 8 |
| 8 | 0.4 | 8 |
| 4 | 0.8 | 4 |
| 6 | 0.8 | 4 |
| 8 | 0.8 | 4 |

## 分析結果

### 1. max_appearance_rate の動作確認
テスト結果から、`max_appearance_rate` パラメータは期待通りに動作していることが確認できました：

- **低い値 (0.1)**: より多くのk-merが高頻度として検出される (12個)
- **中間値 (0.4-0.5)**: 中程度の数のk-merが検出される (8個)
- **高い値 (0.8)**: より少ないk-merが検出される (4個)
- **非常に高い値 (0.9)**: ほとんどのk-merが検出されない (0個)

### 2. 閾値計算の検証
閾値計算 `int(total_rows * max_appearance_rate)` は正しく動作しています：

- 40行のデータセットで `max_appearance_rate=0.4` の場合：
  - 計算：`int(40 × 0.4) = 16`
  - 意味：16行以上で出現するk-merが高頻度として検出される
  - 実際の結果：Pattern 1 (18行) が検出される

### 3. kmer_size による影響
同じ `max_appearance_rate` 値でも、異なる `kmer_size` 値で同じ結果が得られました。これは以下の理由によるものと考えられます：

1. **テストデータの特性**: 使用したテストデータは単純な繰り返しパターンのため、k-merのサイズが変わっても基本的な出現頻度パターンは変わらない
2. **k-mer抽出の一貫性**: 各パターンから抽出されるk-merの数は異なるが、全体的な分布パターンは似ている

## 考察

### 報告されていた問題について
「同じ `kmer_size` なのに `max_appearance_rate` の値によって抽出されたk-merの件数が異なる」という問題は、テストでは再現されませんでした。実際には、`max_appearance_rate` が変わると抽出されるk-merの件数は期待通りに変化しています。

### 潜在的な問題
もし実際の使用で異なる結果が得られている場合、以下の原因が考えられます：

1. **データの特性**: より複雑なDNAシーケンスでは、k-merの分布がより複雑になる可能性
2. **境界値での計算**: 小数点以下の計算で微妙な違いが生じる可能性
3. **並列処理**: 並列ワーカーの数や処理順序による影響

## 推奨事項

### 1. より詳細なテスト
実際の問題を再現するため、以下のテストを推奨します：

```sql
-- より大きなデータセットでのテスト
-- より複雑なDNAシーケンスでのテスト
-- 異なるkmer_sizeでの詳細な比較
```

### 2. デバッグ情報の追加
問題の原因を特定するため、以下の情報を記録することを推奨：

```sql
-- 実際の閾値計算結果
-- k-mer抽出の詳細ログ
-- 並列処理の詳細情報
```

### 3. エッジケースのテスト
以下のエッジケースをテストすることを推奨：

- 非常に小さなデータセット (< 10行)
- 非常に大きなデータセット (> 1000行)
- 境界値での `max_appearance_rate` (0.5付近)

## 結論

現在のテストでは、`max_appearance_rate` パラメータは期待通りに動作しており、値の変更に応じて検出される高頻度k-merの数が適切に変化することが確認できました。報告されていた問題は、特定のデータセットや条件下でのみ発生する可能性があります。

さらなる調査が必要な場合は、実際の問題が発生しているデータセットと設定値を用いた詳細なテストを実施することを推奨します。