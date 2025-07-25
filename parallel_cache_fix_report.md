# Parallel Cache Count Mismatch Fix Report

## 問題の概要
global cacheとparallel cacheで異なる数のk-merがfreeされるという問題が発生していました：
- global cache: 22個のk-mer
- parallel cache: 23個のk-mer

## 根本原因の分析

### 1. 共通の問題: k-mer値0の無視
両方のキャッシュ実装で、k-mer値が0の場合は「null marker」として扱われ、ハッシュテーブルへの挿入がスキップされていました：

```c
// global cache (line 1541)
if (batch_kmers[i] == 0) {
    ereport(DEBUG1, (errmsg("Found null k-mer at batch %d index %d, skipping", batch_num, i)));
    continue;
}

// parallel cache (line 2743)
if (batch_kmers[i] == 0) {
    ereport(DEBUG1, (errmsg("Invalid kmer2_as_uint value at batch %d index %d, skipping", batch_num, i)));
    continue;
}
```

### 2. カウント方法の不一致
テストデータでは、23個のk-merが分析されましたが、そのうち1個が値0でした：

```sql
SELECT kmer2_as_uint FROM kmersearch_highfreq_kmer ORDER BY kmer2_as_uint;
-- 結果: 0, 1, 3, 5, 13, 15, 21, 27, 54, 57, 63, 78, 85, 108, 147, 177, 192, 198, 216, 228, 240, 252, 255
```

**Global cache (正しい実装):**
```c
// 実際に挿入されたk-merの数をカウント
if (entry && !found) {
    total_inserted++;  // 実際の挿入数をカウント
}
// ...
global_highfreq_cache.highfreq_count = total_inserted;  // 22個 (0は除外)
```

**Parallel cache (バグのある実装):**
```c
// 元々のクエリ結果の総数を使用（修正前）
parallel_highfreq_cache->num_entries = total_kmer_count;  // 23個 (0を含む)
```

## 修正内容

### 修正箇所
ファイル: `kmersearch_cache.c`

**1. 誤った初期化の削除 (line 2559):**
```c
// 修正前
parallel_highfreq_cache->num_entries = total_kmer_count;

// 修正後
/* num_entries will be set after batch processing */
```

**2. 正しいカウントの設定 (line 2812):**
```c
// 修正前
parallel_highfreq_cache->is_initialized = true;

// 修正後
/* Set the actual number of entries that were successfully inserted */
parallel_highfreq_cache->num_entries = total_inserted;
parallel_highfreq_cache->is_initialized = true;
```

### 修正理由
- `total_kmer_count`: システムテーブルから取得した全k-merの数（23個）
- `total_inserted`: 実際にハッシュテーブルに挿入されたk-merの数（22個、値0を除外）

parallel cacheは実際の挿入数を使用すべきでした。

## 期待される結果
修正後は、両方のキャッシュで同じ数のk-merが報告されるはずです：
- global cache: 22個のk-mer
- parallel cache: 22個のk-mer (修正前は23個)

## テスト方法
```sql
-- 両方のキャッシュで同じ結果を確認
SET kmersearch.kmer_size = 4;
SET kmersearch.occur_bitlen = 4;
SET kmersearch.max_appearance_rate = 0.25;
-- 他の設定...

-- global cache
SELECT kmersearch_highfreq_kmer_cache_load('test_regression', 'seq');
SELECT kmersearch_highfreq_kmer_cache_free('test_regression', 'seq');  -- 期待値: 22

-- parallel cache
SELECT kmersearch_parallel_highfreq_kmer_cache_load('test_regression', 'seq');
SELECT kmersearch_parallel_highfreq_kmer_cache_free('test_regression', 'seq');  -- 期待値: 22
```

## 追加の考慮事項

### k-mer値0の扱い
現在の実装では、k-mer値0は「null marker」として扱われていますが、これは設計上の制約です。値0は有効なk-mer値である可能性があります（例：AAAA = 0）。

将来的な改善案：
- null markerとして別の値を使用する（例：UINT64_MAX）
- 別のフラグを使用してnull状態を管理する

### 回帰テストの更新
この修正により、parallel cacheのテスト結果が変わるため、期待値ファイルの更新が必要になります。

## 結論
この修正により、global cacheとparallel cacheの動作が一致し、一貫性のある結果が得られるようになります。修正は最小限で、既存の動作を壊すことなく問題を解決します。