# GINインデックスと検索処理のuintkey完全移行計画

## 概要
ngram_key2とkmer2_as_uintを完全に廃止し、既存の`kmersearch_extract_uintkey_from_dna2/4()`および`kmersearch_extract_uintkey_from_text()`関数を活用してuintkey形式に統一する計画です。クエリパターンキャッシュ機能は維持し、uintkey形式に対応させます。

## 命名規則
- 関数や変数名でuintkeyを使用する場合、複数形にせず単数形`uintkey`のままで使用
- 既存の構造体名は変更せず、内部実装のみをuintkey対応に改訂

## 基本方針

### 廃止するデータ形式
1. **ngram_key2 (VarBit)**: k-mer2 + オカレンスカウントのビット列表現
2. **kmer2_as_uint**: オカレンスカウントなしのk-mer値

### 既に実装済みの機能を活用
- **kmersearch_extract_uintkey_from_dna2()**: DNA2からuintkey配列を抽出（実装済み）
- **kmersearch_extract_uintkey_from_dna4()**: DNA4からuintkey配列を抽出（実装済み）
- **kmersearch_extract_uintkey_from_text()**: テキスト文字列からuintkey配列を抽出（実装済み）
- これらの関数は既にオカレンスカウントを計算し、適切なデータ型で出力している

### 維持・改訂する機能
- **クエリパターンキャッシュ**: VarBit形式からuintkey形式に改訂
- **actual_min_scoreキャッシュ**: uintkey形式に対応

## 改訂対象と修正内容

### 1. GINインデックス構築関数の簡略化

#### kmersearch_extract_value_dna2() / kmersearch_extract_value_dna4()
**ファイル**: kmersearch_gin.c

**現在の複雑な処理**:
1. kmer2_as_uint形式で抽出
2. KmerOccurrence構造体でオカレンスカウント追跡
3. 高頻度判定（間違ったキーで検索）
4. ngram_key2形式のVarBitを生成
5. Datum配列として返す

**新しいシンプルな処理**:
```c
Datum
kmersearch_extract_value_dna2(PG_FUNCTION_ARGS)
{
    kmersearch_dna2 *dna = (kmersearch_dna2 *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    
    Datum *keys;
    void *uintkey_array;
    int uintkey_count;
    
    // 既存の関数でuintkey形式で抽出（オカレンスカウント込み）
    kmersearch_extract_uintkey_from_dna2((VarBit *)dna, &uintkey_array, &uintkey_count);
    
    if (kmersearch_preclude_highfreq_kmer) {
        // 高頻度k-merをフィルタリング
        void *filtered_array;
        int filtered_count;
        filter_highfreq_uintkeys(uintkey_array, uintkey_count, 
                                 &filtered_array, &filtered_count);
        pfree(uintkey_array);
        uintkey_array = filtered_array;
        uintkey_count = filtered_count;
    }
    
    // uintkeyを直接Datumとして返す
    keys = convert_uintkeys_to_datums(uintkey_array, uintkey_count);
    
    *nkeys = uintkey_count;
    pfree(uintkey_array);
    PG_RETURN_POINTER(keys);
}
```

### 2. クエリパターンキャッシュのuintkey対応

#### QueryPatternCacheEntry構造体の改訂
**ファイル**: kmersearch.h

```c
// 改訂後の構造体（VarBit形式からuintkey形式へ）
typedef struct QueryPatternCacheEntry {
    uint64 hash_key;
    char *query_pattern;
    int k_size;
    void *uintkey;       // uint16/32/64の配列（VarBit **kmersから変更）
    int nkeys;
    int query_len;
    uint64 last_access;
    uint64 access_count;
} QueryPatternCacheEntry;
```

#### get_cached_query_uintkey()の実装
**ファイル**: kmersearch_cache.c

```c
// 新しい関数: get_cached_query_uintkey()
void *get_cached_query_uintkey(const char *query_string, int k_size, int *nkeys)
{
    QueryPatternCacheEntry *cache_entry;
    void *result_uintkey = NULL;
    uint64 hash_key;
    MemoryContext old_context;
    
    *nkeys = 0;
    
    // キャッシュマネージャーの初期化（必要に応じて）
    if (query_pattern_cache_manager == NULL) {
        old_context = MemoryContextSwitchTo(TopMemoryContext);
        init_query_pattern_cache_manager(&query_pattern_cache_manager);
        MemoryContextSwitchTo(old_context);
    }
    
    // ハッシュキーの計算
    hash_key = hash_any((unsigned char *)query_string, strlen(query_string));
    
    // キャッシュから検索
    cache_entry = lookup_query_pattern_cache(hash_key, query_string, k_size);
    
    if (cache_entry != NULL) {
        // キャッシュヒット
        result_uintkey = copy_uintkey_array(cache_entry->uintkey, cache_entry->nkeys);
        *nkeys = cache_entry->nkeys;
        cache_entry->access_count++;
        cache_entry->last_access = GetCurrentTimestamp();
    } else {
        // キャッシュミス - 新規抽出
        kmersearch_extract_uintkey_from_text(query_string, &result_uintkey, nkeys);
        
        // キャッシュに追加
        if (result_uintkey != NULL && *nkeys > 0) {
            add_to_query_pattern_cache(hash_key, query_string, k_size, 
                                       result_uintkey, *nkeys);
        }
    }
    
    return result_uintkey;
}
```

### 3. GIN検索関数の修正

#### kmersearch_extract_query()
**ファイル**: kmersearch_gin.c

```c
Datum
kmersearch_extract_query(PG_FUNCTION_ARGS)
{
    Datum query = PG_GETARG_DATUM(0);
    int32 *nkeys = (int32 *) PG_GETARG_POINTER(1);
    StrategyNumber strategy = PG_GETARG_UINT16(2);
    bool **pmatch = (bool **) PG_GETARG_POINTER(3);
    Pointer **extra_data = (Pointer **) PG_GETARG_POINTER(4);
    bool **nullFlags = (bool **) PG_GETARG_POINTER(5);
    int32 *searchMode = (int32 *) PG_GETARG_POINTER(6);
    
    text *query_text = DatumGetTextP(query);
    char *query_string = text_to_cstring(query_text);
    Datum *keys;
    void *uintkey_array;
    int uintkey_count;
    
    // キャッシュを使用してuintkey抽出
    uintkey_array = get_cached_query_uintkey(query_string, kmersearch_kmer_size, &uintkey_count);
    
    // 高頻度フィルタリングとactual_min_scoreの計算
    if (uintkey_array != NULL && uintkey_count > 0) {
        uintkey_array = filter_uintkey_and_set_actual_min_score(uintkey_array, &uintkey_count, 
                                                                query_string);
    }
    
    // uintkeyをDatum配列に変換
    if (uintkey_array != NULL && uintkey_count > 0) {
        keys = convert_uintkey_to_datum(uintkey_array, uintkey_count);
    } else {
        keys = NULL;
    }
    
    *pmatch = NULL;
    *extra_data = NULL;
    *nullFlags = NULL;
    *searchMode = GIN_SEARCH_MODE_DEFAULT;
    
    pfree(query_string);
    // NOTE: uintkey_arrayはキャッシュ管理されているため、解放しない
    *nkeys = uintkey_count;
    
    if (*nkeys == 0)
        PG_RETURN_POINTER(NULL);
    
    PG_RETURN_POINTER(keys);
}
```

### 4. actual_min_scoreキャッシュのuintkey対応

#### filter_uintkey_and_set_actual_min_score()
```c
void *filter_uintkey_and_set_actual_min_score(void *uintkey_array, int *nkeys, 
                                              const char *query_string)
{
    void *filtered_keys;
    int filtered_count = 0;
    int original_nkeys = *nkeys;
    int actual_min_score;
    
    if (!uintkey_array || *nkeys <= 0)
        return uintkey_array;
    
    // Step 1: フィルタリング前のuintkey配列でactual_min_scoreを計算してキャッシュ
    actual_min_score = get_cached_actual_min_score(uintkey_array, *nkeys);
    
    // Step 2: 高頻度k-merのフィルタリング（有効な場合）
    if (!kmersearch_preclude_highfreq_kmer) {
        return uintkey_array;
    }
    
    // フィルタリング実行
    filter_highfreq_uintkey(uintkey_array, *nkeys, &filtered_keys, &filtered_count);
    
    *nkeys = filtered_count;
    
    elog(DEBUG1, "filter_uintkey_and_set_actual_min_score: filtered %d high-freq k-mers from %d total, "
                 "cached actual_min_score=%d", 
         original_nkeys - filtered_count, original_nkeys, actual_min_score);
    
    if (filtered_count == 0) {
        if (filtered_keys)
            pfree(filtered_keys);
        return NULL;
    }
    
    return filtered_keys;
}
```

### 5. 新規実装が必要な関数

#### filter_highfreq_uintkey()
```c
void filter_highfreq_uintkey(void *uintkey, int count, void **filtered, int *filtered_count)
{
    int total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;
    int filtered_idx = 0;
    void *result;
    
    // 結果配列を確保
    if (total_bits <= 16) {
        uint16 *src = (uint16 *)uintkey;
        uint16 *dst = (uint16 *)palloc(count * sizeof(uint16));
        
        for (int i = 0; i < count; i++) {
            if (!kmersearch_is_uintkey_highfreq(src[i])) {
                dst[filtered_idx++] = src[i];
            }
        }
        result = dst;
    } else if (total_bits <= 32) {
        uint32 *src = (uint32 *)uintkey;
        uint32 *dst = (uint32 *)palloc(count * sizeof(uint32));
        
        for (int i = 0; i < count; i++) {
            if (!kmersearch_is_uintkey_highfreq(src[i])) {
                dst[filtered_idx++] = src[i];
            }
        }
        result = dst;
    } else {
        uint64 *src = (uint64 *)uintkey;
        uint64 *dst = (uint64 *)palloc(count * sizeof(uint64));
        
        for (int i = 0; i < count; i++) {
            if (!kmersearch_is_uintkey_highfreq(src[i])) {
                dst[filtered_idx++] = src[i];
            }
        }
        result = dst;
    }
    
    *filtered = result;
    *filtered_count = filtered_idx;
}
```

#### kmersearch_is_uintkey_highfreq()
```c
bool kmersearch_is_uintkey_highfreq(uint64 uintkey)
{
    // 直接uintkeyでキャッシュを検索
    if (global_highfreq_cache.is_valid) {
        return kmersearch_lookup_uintkey_in_global_cache(uintkey, NULL, NULL);
    }
    if (kmersearch_is_parallel_highfreq_cache_loaded()) {
        return kmersearch_lookup_uintkey_in_parallel_cache(uintkey, NULL, NULL);
    }
    return false;
}
```

#### convert_uintkey_to_datum()
```c
Datum *convert_uintkey_to_datum(void *uintkey_array, int count)
{
    int total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;
    Datum *keys = (Datum *) palloc(count * sizeof(Datum));
    
    if (total_bits <= 16) {
        uint16 *arr = (uint16 *)uintkey_array;
        for (int i = 0; i < count; i++)
            keys[i] = Int16GetDatum(arr[i]);
    } else if (total_bits <= 32) {
        uint32 *arr = (uint32 *)uintkey_array;
        for (int i = 0; i < count; i++)
            keys[i] = Int32GetDatum(arr[i]);
    } else {
        uint64 *arr = (uint64 *)uintkey_array;
        for (int i = 0; i < count; i++)
            keys[i] = Int64GetDatum(arr[i]);
    }
    
    return keys;
}
```

#### copy_uintkey_array()
```c
void *copy_uintkey_array(void *src_array, int count)
{
    int total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;
    void *result;
    
    if (total_bits <= 16) {
        size_t size = count * sizeof(uint16);
        result = palloc(size);
        memcpy(result, src_array, size);
    } else if (total_bits <= 32) {
        size_t size = count * sizeof(uint32);
        result = palloc(size);
        memcpy(result, src_array, size);
    } else {
        size_t size = count * sizeof(uint64);
        result = palloc(size);
        memcpy(result, src_array, size);
    }
    
    return result;
}
```

### 6. その他のGIN関数の修正

#### kmersearch_consistent()
- VarBit配列の処理をuintkey配列の処理に変更
- actual_min_scoreの取得方法を更新

#### kmersearch_compare_partial()
```c
Datum
kmersearch_compare_partial(PG_FUNCTION_ARGS)
{
    Datum a = PG_GETARG_DATUM(0);
    Datum b = PG_GETARG_DATUM(1);
    
    // uintkey同士を直接比較
    int total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;
    
    if (total_bits <= 16) {
        int16 va = DatumGetInt16(a);
        int16 vb = DatumGetInt16(b);
        PG_RETURN_INT32(va - vb);
    } else if (total_bits <= 32) {
        int32 va = DatumGetInt32(a);
        int32 vb = DatumGetInt32(b);
        PG_RETURN_INT32(va - vb);
    } else {
        int64 va = DatumGetInt64(a);
        int64 vb = DatumGetInt64(b);
        if (va < vb) PG_RETURN_INT32(-1);
        if (va > vb) PG_RETURN_INT32(1);
        PG_RETURN_INT32(0);
    }
}
```

### 7. 削除・改訂が必要な関数

#### 削除する関数（VarBit専用）
- `kmersearch_extract_dna2_kmer2_as_uint_direct()`
- `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()`
- `kmersearch_extract_dna2_ngram_key2_direct()`
- `kmersearch_extract_dna4_ngram_key2_with_expansion_direct()`
- `kmersearch_create_ngram_key2_from_kmer2_as_uint()`
- `kmersearch_kmer2_as_uint_to_kmer2()`
- `kmersearch_remove_occurrence_from_ngram_key2()`
- `kmersearch_find_or_add_kmer_occurrence()`
- `kmersearch_is_kmer_highfreq()` (VarBit版)
- `filter_ngram_key2_and_set_actual_min_score()` (VarBit版)
- `get_cached_query_kmer()` (VarBit版)

#### 新規実装する関数
- `get_cached_query_uintkey()` - uintkey形式でクエリパターンをキャッシュ

#### 改訂する関数
- `get_cached_actual_min_score()` - uintkey形式に対応
- `init_query_pattern_cache_manager()` - uintkey形式に対応
- `lookup_query_pattern_cache()` - uintkey形式に対応
- `add_to_query_pattern_cache()` - uintkey形式に対応

#### 改訂する構造体
- `QueryPatternCacheEntry` - VarBit **kmersをvoid *uintkeyに変更
- `ActualMinScoreCacheEntry` - uintkey形式に対応

### 8. GINオペレータクラスの変更

```sql
-- STORAGEをint8に変更（最大64ビットのuintkeyを格納）
CREATE OPERATOR CLASS gin_kmersearch_ops
FOR TYPE kmersearch_dna2 USING gin
AS
    OPERATOR 1 =%,
    FUNCTION 1 btint8cmp(int8, int8),        -- 標準のint8比較関数を使用
    FUNCTION 2 kmersearch_extract_value_dna2,
    FUNCTION 3 kmersearch_extract_query,
    FUNCTION 4 kmersearch_consistent,
    FUNCTION 5 kmersearch_compare_partial,
    STORAGE int8;  -- uint64として格納
```

### 9. マッチング関数の修正

#### kmersearch_dna2_match() / kmersearch_dna4_match()
```c
// 現在の実装（VarBit使用）
query_keys = get_cached_query_kmer(pattern_string, kmersearch_kmer_size, &query_nkeys);

// 新しい実装（uintkey使用）
void *query_uintkey, *seq_uintkey;
int query_count, seq_count;

// キャッシュを使用してクエリのuintkeyを取得
query_uintkey = get_cached_query_uintkey(pattern_string, kmersearch_kmer_size, &query_count);

// シーケンスからuintkey抽出
kmersearch_extract_uintkey_from_dna2(sequence, &seq_uintkey, &seq_count);

// uintkey配列で直接マッチング計算
int shared_count = count_matching_uintkey(seq_uintkey, seq_count, query_uintkey, query_count);
```

## 実装手順

### Phase 1: キャッシュ関連の改訂
1. `QueryPatternCacheEntry`構造体の改訂（uintkey対応）
2. `get_cached_query_uintkey()`の実装（新規）
3. `init_query_pattern_cache_manager()`の改訂（uintkey対応）
4. `copy_uintkey_array()`の実装
5. `get_cached_actual_min_score()`の改訂（uintkey対応）

### Phase 2: 補助関数の実装
1. `filter_highfreq_uintkey()`の実装
2. `kmersearch_is_uintkey_highfreq()`の実装
3. `convert_uintkey_to_datum()`の実装
4. `filter_uintkey_and_set_actual_min_score()`の実装

### Phase 3: GIN関数の修正
1. `kmersearch_extract_value_dna2()`の簡略化
2. `kmersearch_extract_value_dna4()`の簡略化
3. `kmersearch_extract_query()`の修正
4. `kmersearch_consistent()`の修正
5. `kmersearch_compare_partial()`の修正

### Phase 4: マッチング・スコア関数の修正
1. `kmersearch_dna2_match()`の修正
2. `kmersearch_dna4_match()`の修正
3. `kmersearch_correctedscore_dna2()`の修正
4. `kmersearch_correctedscore_dna4()`の修正

### Phase 5: クリーンアップ
1. VarBit専用関数の削除
2. 旧構造体の削除
3. SQLオペレータクラスの更新

### Phase 6: テスト
1. キャッシュ機能のテスト
2. GINインデックスのテスト
3. パフォーマンステスト

## メリット

1. **キャッシュ機能の維持**: クエリパターンキャッシュを維持し、パフォーマンスを保持
2. **既存資産の活用**: 実装済みのuintkey抽出関数を最大限活用
3. **高速化**: VarBit変換オーバーヘッドの削除
4. **メモリ効率**: VarBit構造のオーバーヘッドがなくなる
5. **保守性向上**: データ形式が統一され、理解しやすい

## リスクと対策

### 1. GINインデックスの再構築
**リスク**: 既存のインデックスが使用不可
**対策**: 新規インストールのみ対応（CLAUDE.mdに記載済み）

### 2. キャッシュメカニズムの移行
**リスク**: キャッシュデータの互換性
**対策**: 新しいキャッシュ構造体で再実装、既存のキャッシュは自動的にクリア

### 3. 高頻度判定の精度
**リスク**: uintkeyベースの判定により、オカレンスカウント別の除外が発生
**対策**: これは仕様として受け入れる（より精密な制御が可能）

## テスト計画

### 1. 単体テスト
```sql
-- キャッシュ機能のテスト
SELECT kmersearch_test_query_cache_uintkey();

-- uintkey抽出テスト
SELECT pg_typeof(unnest) FROM (
    SELECT unnest(kmersearch_test_extract_value('ACGTACGT'::kmersearch_dna2))
) t;

-- 高頻度フィルタリングテスト
SELECT COUNT(*) FROM test_table WHERE seq_column =% 'ACGTACGT';
```

### 2. パフォーマンステスト
- キャッシュヒット率の測定
- インデックス構築時間: VarBit版 vs uintkey版
- 検索速度の比較
- メモリ使用量の比較

## 期待される効果

1. **キャッシュ効率の向上**: uintkey形式により、キャッシュサイズが削減
2. **処理速度**: VarBit変換の削除により20-30%の高速化
3. **メモリ使用量**: VarBit構造のオーバーヘッド（20-30%）の削減
4. **保守性**: シンプルで理解しやすいコード

## 結論

クエリパターンキャッシュ機能を維持しながら、uintkey形式への完全移行を実現します。既存のキャッシュメカニズムをuintkey対応に改訂することで、パフォーマンスを維持しつつ、コードの簡略化とメモリ効率の向上を達成できます。