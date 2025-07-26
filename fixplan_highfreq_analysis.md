# 高頻度k-mer解析の並列実行化計画

## 実装状況: 完了 (2025-07-25)

本計画に基づく実装が完了しました。主要な変更点：
- PostgreSQLの標準並列実行機能（ParallelContext）を使用した真の並列処理の実装
- ブロック単位での動的ワーク分散によるスケーラブルな処理
- k-merサイズに応じた最適化されたdshashテーブル実装
- エラー処理とリソース管理の改善

## 実装制約事項

**重要**: 以下の実装を禁止します：
- スタブ実装禁止
- ダミーデータ生成コード実装禁止
- 不完全な関数実装禁止
- 無駄なラッパー関数やヘルパー関数の実装禁止
- ワークアラウンドの実装禁止
- 既存の関数や構造体や変数と同一機能のものを実装することを禁止
- GUC変数がある場合は直接GUC変数を使用し、ローカル変数にGUC変数を代入することを禁止
- 新しい関数を作成する前に既存の関数の改修で対応できないか検討

## 概要

現在の`kmersearch_perform_highfreq_analysis()`は疑似並列処理（単一プロセスによる順次処理）となっており、真の並列実行が行われていません。本計画では、PostgreSQLの標準並列実行機能とDSM/dshashを使用した真の並列処理に変更します。

## システム要件

**PostgreSQLバージョン要件**: PostgreSQL 16以降のみ対応。PostgreSQL 16以降では、並列ワーカー関数の宣言前に`PGDLLEXPORT`を付ける必要があります。

## 現在の問題点

### 1. 疑似並列処理
- `KmerWorkerState`配列を作成しているが、実際はメインプロセスが順次処理
- PostgreSQLの並列ワーカープロセス機能を使用していない
- `result.parallel_workers_used`は設定値を返すが実際の並列実行なし

### 2. 非効率な処理フロー
- SQL集計結果を一旦dshashに格納する方式が不適切
- 並列ワーカー間でのk-mer出現回数集計に最適化されていない
- 真の並列I/Oが発生しない

## 実装戦略

### 参考実装
- **pg_foobar**: 単純な並列実行パターン（DSM/shm_toc使用）
- **kmersearch_parallel_highfreq_kmer_cache_load()**: DSM/DSA/dshashの使用方法
- **README.parallel**: PostgreSQL並列実行の基本原則

## 改訂された並列処理フロー

### 親プロセス（リーダープロセス）の責務

1. **引数解析**
   - テーブル名またはOIDからテーブルOIDを解決
   - カラム名またはattnumからカラムattnumとカラム名を解決

2. **事前検証**
   - 対象カラムがDNA2型またはDNA4型かチェック
   - 既存の解析データがないかチェック

3. **テーブルロック管理**
   - 対象テーブルをExclusiveLockでロック（解析中の一貫性保証）
   - 解析完了後にロックを解放

4. **並列実行環境の準備**
   - DSMセグメントの作成（サイズ推定含む）
   - DSA（Dynamic Shared Area）の作成
   - dshashテーブルの作成（k-merサイズに応じた最適化）
   - 共有状態構造体の初期化とshm_tocへの登録

5. **並列コンテキスト作成**
   - `CreateParallelContext()`で並列コンテキスト作成
   - shm_tocにDSA、dshashハンドル、共有状態を登録
   - `InitializeParallelDSM()`でDSM初期化

6. **並列ワーカー起動**
   - `LaunchParallelWorkers()`で並列ワーカー起動
   - 実際に起動したワーカー数を取得

7. **テーブルスキャンとワーク分散**
   - テーブル全体をスキャンし、バッチ単位でTID範囲を生成
   - 共有状態にTID範囲を設定（ワーカーが取得するまで待機）
   - 全行処理完了時に`all_processed`フラグ設定

8. **結果統合**
   - `WaitForParallelWorkersToFinish()`で全ワーカー完了待機
   - dshashから高頻度k-merを抽出
   - システムテーブルに結果を保存

### 並列ワーカーの責務

1. **DSMアタッチとshm_toc参照**
   - shm_tocから共有状態、DSA、dshashハンドルを取得
   - dshashテーブルにアタッチ

2. **動的ワーク取得ループ**
   - 共有状態のmutexで排他制御
   - TID範囲が利用可能な場合、取得して`tid_range_available`をfalseに設定
   - TID範囲が空で`all_processed`がfalseなら待機
   - `all_processed`がtrueならワーカー終了

3. **k-mer抽出処理**
   - 取得したTID範囲でテーブルスキャン
   - DNA2型: `kmersearch_extract_dna2_kmer2_as_uint_direct()`で直接抽出
   - DNA4型: `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()`で直接抽出
   - 抽出したkmer2_as_uintをdshashに登録/カウント更新

4. **エラーハンドリング**
   - PG_TRY/PG_CATCHでエラーキャッチ
   - エラー発生時は共有状態にエラー情報を記録
   - リソースの適切なクリーンアップ

## DSM/dshash実装詳細

### k-mer直接抽出の活用

既存の直接抽出関数を活用してuint型として直接k-merを取得：
- `kmersearch_extract_dna2_kmer2_as_uint_direct()`: DNA2型からuint配列として直接抽出
- `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()`: DNA4型から拡張付きで直接抽出

これらの関数は既にk-merサイズに応じてuint16/uint32/uint64の配列を返すため、型変換不要でdshashに保存可能。

### DSM/dshash基本構造

#### 共有メモリ構造体
```c
typedef struct KmerAnalysisSharedState
{
    LWLock      mutex;                    /* 排他制御用 */
    int         num_workers;              /* 並列ワーカー数 */
    Oid         table_oid;               /* 対象テーブルOID */
    AttrNumber  column_attnum;           /* 対象カラムattnum */
    Oid         column_type_oid;         /* 対象カラムのデータ型OID */
    int         kmer_size;               /* k-merサイズ（kmersearch_kmer_size） */
    int         batch_size;              /* バッチサイズ（kmersearch_highfreq_analysis_batch_size） */
    bool        all_processed;           /* 全行処理済みフラグ */
    BlockNumber next_block;              /* 次の処理対象ブロック番号 */
    BlockNumber total_blocks;            /* テーブル総ブロック数 */
    bool        worker_error_occurred;   /* 並列ワーカーエラー発生フラグ */
    char        error_message[256];      /* エラーメッセージ */
} KmerAnalysisSharedState;
```

#### dshash エントリ構造体
```c
/* k-merサイズに応じたdshashエントリ構造体 */
/* kmer2_as_uintを直接キーとして使用（ハッシュ化不要） */

/* kmer_size <= 8の場合 */
typedef struct KmerEntry16
{
    uint16 kmer;    /* kmer2_as_uint（キー） */
    int count;      /* 出現行数 */
} KmerEntry16;

/* kmer_size <= 16の場合 */
typedef struct KmerEntry32
{
    uint32 kmer;    /* kmer2_as_uint（キー） */
    int count;      /* 出現行数 */
} KmerEntry32;

/* kmer_size <= 32の場合 */
typedef struct KmerEntry64
{
    uint64 kmer;    /* kmer2_as_uint（キー） */
    int count;      /* 出現行数 */
} KmerEntry64;
```

### 実装詳細

#### 共有メモリキー定義
```c
#define KMERSEARCH_KEY_SHARED_STATE  1
#define KMERSEARCH_KEY_DSA          2  
#define KMERSEARCH_KEY_HASH_HANDLE  3
```

#### shm_toc使用パターン
- **親プロセス**: `shm_toc_allocate()`と`shm_toc_insert()`で登録
- **ワーカー**: `shm_toc_lookup()`で取得

#### フェーズ1: DSM/dshashベース並列実行環境の構築

##### 1.1 並列実行環境の初期化

**DSM/dshash初期化関数（既存関数を改修）**
```c
/* 既存のkmersearch_create_analysis_dshash関数を改修 */
static bool
kmersearch_create_analysis_dshash(int estimated_entries)
{
    /* 既存実装をベースに以下を追加：
     * 1. k-merサイズに応じたdshashパラメータ設定
     * 2. 恒等ハッシュ関数の使用（uint16/uint32）
     * 3. uint64の場合はdshash_memhash使用
     */
}

/* 恒等ハッシュ関数（kmer2_as_uintは既にハッシュ値） */
static uint32
kmersearch_uint16_identity_hash(const void *key, Size keysize, void *arg)
{
    return (uint32)(*(const uint16 *)key);
}

static uint32
kmersearch_uint32_identity_hash(const void *key, Size keysize, void *arg)
{
    return *(const uint32 *)key;
}
```

##### 1.2 並列ワーカー関数の実装

**並列ワーカーメイン関数（新規作成、PGDLLEXPORTが必要）**
```c
PGDLLEXPORT void kmersearch_analysis_worker(dsm_segment *seg, shm_toc *toc);

void
kmersearch_analysis_worker(dsm_segment *seg, shm_toc *toc)
{
    KmerAnalysisSharedState *shared_state = NULL;
    dsa_area   *dsa = NULL;
    dshash_table *hash = NULL;
    BlockNumber current_block;
    bool        has_work = true;
    
    PG_TRY();
    {
        /* 共有状態取得 */
        shared_state = (KmerAnalysisSharedState *)shm_toc_lookup(toc, KMERSEARCH_KEY_SHARED_STATE, false);
        dsa = (dsa_area *)shm_toc_lookup(toc, KMERSEARCH_KEY_DSA, false);
        
        /* dshashテーブルにアタッチ */
        dshash_table_handle *handle = (dshash_table_handle *)shm_toc_lookup(toc, KMERSEARCH_KEY_HASH_HANDLE, false);
        hash = dshash_attach(dsa, &params, handle, NULL);
        
        /* 動的ワーク取得ループ */
        while (has_work)
        {
            /* ブロック番号の取得（短時間の排他制御） */
            LWLockAcquire(&shared_state->mutex, LW_EXCLUSIVE);
            
            if (shared_state->next_block < shared_state->total_blocks)
            {
                /* 処理対象ブロックを取得 */
                current_block = shared_state->next_block;
                shared_state->next_block++;
                LWLockRelease(&shared_state->mutex);
                
                /* k-mer抽出・集計処理実行（排他制御外で並列実行） */
                kmersearch_extract_kmers_from_block(shared_state->table_oid,
                                                    shared_state->column_attnum,
                                                    shared_state->column_type_oid,
                                                    current_block,
                                                    shared_state->kmer_size,
                                                    hash);
            }
            else
            {
                /* 全ブロック処理完了 */
                shared_state->all_processed = true;
                LWLockRelease(&shared_state->mutex);
                has_work = false;
            }
        }
    }
    PG_CATCH();
    {
        /* 並列ワーカーエラー処理 */
        if (shared_state)
        {
            LWLockAcquire(&shared_state->mutex, LW_EXCLUSIVE);
            if (!shared_state->worker_error_occurred)
            {
                shared_state->worker_error_occurred = true;
                strlcpy(shared_state->error_message, 
                       "Parallel worker encountered error during k-mer extraction",
                       sizeof(shared_state->error_message));
            }
            LWLockRelease(&shared_state->mutex);
        }
        
        /* リソースクリーンアップ */
        if (hash)
            dshash_detach(hash);
        
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    /* 正常終了時のクリーンアップ */
    if (hash)
        dshash_detach(hash);
}
```

##### 1.3 k-mer抽出・集計処理

**ブロック単位のk-mer抽出**
```c
static void
kmersearch_extract_kmers_from_block(Oid table_oid, AttrNumber column_attnum, Oid column_type_oid,
                                    BlockNumber block, int kmer_size, dshash_table *hash)
{
    Relation    rel = NULL;
    Buffer      buffer;
    Page        page;
    OffsetNumber offset;
    int         ntuples;
    bool        isnull;
    Datum       datum;
    
    PG_TRY();
    {
        /* テーブルオープン */
        rel = table_open(table_oid, AccessShareLock);
        
        /* ブロックを読み込み */
        buffer = ReadBuffer(rel, block);
        LockBuffer(buffer, BUFFER_LOCK_SHARE);
        page = BufferGetPage(buffer);
        ntuples = PageGetMaxOffsetNumber(page);
        
        /* ブロック内の各タプルを処理 */
        for (offset = FirstOffsetNumber; offset <= ntuples; offset++)
        {
            ItemId      itemid = PageGetItemId(page, offset);
            HeapTupleData tuple;
            
            if (!ItemIdIsNormal(itemid))
                continue;
                
            tuple.t_data = (HeapTupleHeader) PageGetItem(page, itemid);
            tuple.t_len = ItemIdGetLength(itemid);
            ItemPointerSet(&(tuple.t_self), block, offset);
            
            /* カラム値取得 */
            datum = heap_getattr(&tuple, column_attnum, RelationGetDescr(rel), &isnull);
            if (isnull)
                continue;
                
            /* k-mer抽出とdshash更新 */
            kmersearch_update_kmer_counts_in_dshash(datum, kmer_size, hash, column_type_oid);
        }
        
        UnlockReleaseBuffer(buffer);
    }
    PG_CATCH();
    {
        /* エラー時のリソースクリーンアップ */
        if (rel)
            table_close(rel, AccessShareLock);
        
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    /* 正常終了時のクリーンアップ */
    table_close(rel, AccessShareLock);
}

**dshash更新処理（既存関数を改修）**
```c
static void
kmersearch_update_kmer_counts_in_dshash(Datum sequence_datum, int kmer_size, dshash_table *hash, Oid column_type_oid)
{
    void       *kmer_array = NULL;
    int         kmer_count;
    bool        found;
    BitSet     *seen_kmers = NULL; /* 重複除外用 */
    
    PG_TRY();
    {
        /* データ型に応じたk-mer抽出 */
        if (column_type_oid == TypenameGetTypid("dna2"))
        {
            /* DNA2型からの直接抽出 */
            kmersearch_extract_dna2_kmer2_as_uint_direct(sequence_datum, kmer_size, &kmer_array, &kmer_count);
        }
        else if (column_type_oid == TypenameGetTypid("dna4"))
        {
            /* DNA4型からの拡張付き直接抽出 */
            kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(sequence_datum, kmer_size, &kmer_array, &kmer_count);
        }
        else
        {
            elog(ERROR, "Column must be DNA2 or DNA4 type");
        }
        
        if (kmer_array == NULL || kmer_count <= 0)
            return; /* 有効なk-merなし */
        
        /* 重複除外用BitSet初期化 */
        if (kmer_size <= 16)
            seen_kmers = bitset_create(1ULL << (kmer_size * 2));
        
        /* 各k-merをdshashに登録（1行につき1回のみカウント） */
        for (int i = 0; i < kmer_count; i++)
        {
            /* k-merサイズに応じた処理 */
            if (kmer_size <= 8)
            {
                uint16 *kmer16_array = (uint16 *)kmer_array;
                uint16 kmer16 = kmer16_array[i];
                
                /* 重複チェック */
                if (bitset_is_member(seen_kmers, kmer16))
                    continue;
                bitset_add_member(seen_kmers, kmer16);
                
                KmerEntry16 *entry16 = (KmerEntry16 *)dshash_find_or_insert(hash, &kmer16, &found);
                if (!found)
                {
                    entry16->kmer = kmer16;
                    entry16->count = 1;
                }
                else
                {
                    entry16->count++;
                }
                dshash_release_lock(hash, entry16);
            }
            else if (kmer_size <= 16)
            {
                typedef struct { uint32 kmer; int count; } KmerEntry32;
                uint32 *kmer32_array = (uint32 *)kmer_array;
                uint32 kmer32 = kmer32_array[i]; /* 既に適切な型で取得済み */
                KmerEntry32 *entry32 = (KmerEntry32 *)dshash_find_or_insert(hash, &kmer32, &found);
                if (!found)
                {
                    entry32->kmer = kmer32;
                    entry32->count = 1;
                }
                else
                {
                    if (entry32->count < INT_MAX)
                        entry32->count++;
                }
                dshash_release_lock(hash, entry32);
            }
            else /* kmer_size <= 32 */
            {
                typedef struct { uint64 kmer; int count; } KmerEntry64;
                uint64 *kmer64_array = (uint64 *)kmer_array;
                uint64 kmer64 = kmer64_array[i]; /* 既に適切な型で取得済み */
                KmerEntry64 *entry64 = (KmerEntry64 *)dshash_find_or_insert(hash, &kmer64, &found);
                if (!found)
                {
                    entry64->kmer = kmer64;
                    entry64->count = 1;
                }
                else
                {
                    if (entry64->count < INT_MAX)
                        entry64->count++;
                }
                dshash_release_lock(hash, entry64);
            }
        }
    }
    PG_CATCH();
    {
        /* エラー時のメモリクリーンアップ */
        if (kmer_array)
            pfree(kmer_array);
        
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    /* 正常終了時のメモリクリーンアップ */
    if (kmer_array)
        pfree(kmer_array);
    if (seen_kmers)
        bitset_free(seen_kmers);
}
```

#### フェーズ2: dshash結果の効率的な統合

##### 2.1 dshash結果分析（既存のkmersearch_insert_kmer2_as_uint_from_dshash関数を改修）

```c
static void
kmersearch_insert_kmer2_as_uint_from_dshash(Oid table_oid, const char *column_name, int k_size)
{
    /* 既存実装をベースに以下を改修：
     * 1. analysis_highfreq_hashの代わりに並列実行で作成したdshashを使用
     * 2. k-merサイズに応じたエントリ構造体への対応
     * 3. 閾値判定とバルクインサート処理
     */
    
}
```

#### フェーズ3: 統合メイン関数の実装

##### 3.1 メイン関数の改修（既存のkmersearch_perform_highfreq_analysis_parallel関数を改修）

```c
KmerAnalysisResult
kmersearch_perform_highfreq_analysis_parallel(Oid table_oid, const char *column_name, 
                                             int k_size, int requested_workers)
{
    KmerAnalysisResult result = {0};
    ParallelContext *pcxt = NULL;
    KmerAnalysisSharedState *shared_state = NULL;
    Size        estimated_size = 0;
    Size        estimated_keys = 0;
    shm_toc_estimator estimator;
    shm_toc    *toc;
    char       *shm_pointer;
    bool        table_locked = false;
    
    PG_TRY();
    {
        /* 事前検証 */
        Relation rel = table_open(table_oid, AccessShareLock);
        AttrNumber column_attnum = get_attnum(table_oid, column_name);
        if (column_attnum == InvalidAttrNumber)
            elog(ERROR, "Column \"%s\" does not exist", column_name);
        
        Oid column_type_oid = TupleDescAttr(RelationGetDescr(rel), column_attnum - 1)->atttypid;
        if (column_type_oid != TypenameGetTypid("dna2") && 
            column_type_oid != TypenameGetTypid("dna4"))
            elog(ERROR, "Column must be DNA2 or DNA4 type");
        
        BlockNumber total_blocks = RelationGetNumberOfBlocks(rel);
        table_close(rel, AccessShareLock);
        
        /* テーブルロック取得 */
        LockRelationOid(table_oid, ExclusiveLock);
        table_locked = true;
        
        /* 並列モード開始 */
        EnterParallelMode();
        
        /* 並列コンテキスト作成 */
        pcxt = CreateParallelContext("pg_kmersearch", "kmersearch_analysis_worker", requested_workers);
        
        /* DSMサイズ推定 */
        shm_toc_initialize_estimator(&estimator);
        estimated_size += MAXALIGN(sizeof(KmerAnalysisSharedState));
        estimated_size += MAXALIGN(sizeof(dsa_handle));
        estimated_size += MAXALIGN(sizeof(dshash_table_handle));
        estimated_keys += 3; /* SHARED_STATE, DSA, HASH_HANDLE */
        
        shm_toc_estimate_chunk(&pcxt->estimator, estimated_size);
        shm_toc_estimate_keys(&pcxt->estimator, estimated_keys);
        
        /* DSM初期化 */
        InitializeParallelDSM(pcxt);
        toc = pcxt->toc;
        
        /* 共有状態の設定 */
        shm_pointer = shm_toc_allocate(toc, sizeof(KmerAnalysisSharedState));
        shared_state = (KmerAnalysisSharedState *)shm_pointer;
        LWLockInitialize(&shared_state->mutex, LWTRANCHE_PARALLEL_QUERY_DSA);
        shared_state->num_workers = requested_workers;
        shared_state->table_oid = table_oid;
        shared_state->column_attnum = column_attnum;
        shared_state->column_type_oid = column_type_oid;
        shared_state->kmer_size = k_size;
        shared_state->batch_size = kmersearch_highfreq_analysis_batch_size;
        shared_state->all_processed = false;
        shared_state->next_block = 0;
        shared_state->total_blocks = total_blocks;
        shared_state->worker_error_occurred = false;
        shm_toc_insert(toc, KMERSEARCH_KEY_SHARED_STATE, shared_state);
        
        /* dshash作成 */
        if (!kmersearch_create_analysis_dshash(1000000))
            elog(ERROR, "Failed to create analysis dshash");
        
        /* DSAハンドルとdshashハンドルをtocに登録 */
        dsa_handle *dsa_handle_ptr = (dsa_handle *)shm_toc_allocate(toc, sizeof(dsa_handle));
        *dsa_handle_ptr = dsa_get_handle(analysis_dsa);
        shm_toc_insert(toc, KMERSEARCH_KEY_DSA, dsa_handle_ptr);
        
        dshash_table_handle *hash_handle_ptr = (dshash_table_handle *)shm_toc_allocate(toc, sizeof(dshash_table_handle));
        *hash_handle_ptr = dshash_get_hash_table_handle(analysis_highfreq_hash);
        shm_toc_insert(toc, KMERSEARCH_KEY_HASH_HANDLE, hash_handle_ptr);
        
        /* 並列ワーカー起動 */
        LaunchParallelWorkers(pcxt);
        result.parallel_workers_used = pcxt->nworkers_launched;
        
        /* ワーカー完了待機 */
        WaitForParallelWorkersToFinish(pcxt);
        
        /* エラーチェック */
        if (shared_state->worker_error_occurred)
            elog(ERROR, "Parallel worker error: %s", shared_state->error_message);
        
        /* 結果保存 */
        kmersearch_insert_kmer2_as_uint_from_dshash(table_oid, column_name, k_size);
        
        /* 高頻度k-mer数取得 */
        result.highfreq_kmers_count = kmersearch_get_analysis_dshash_count();
        result.total_rows = kmersearch_get_table_row_count(table_oid);
        result.max_appearance_rate_used = kmersearch_max_appearance_rate;
        result.max_appearance_nrow_used = kmersearch_max_appearance_nrow;
    }
    PG_CATCH();
    {
        /* エラー時のクリーンアップ */
        kmersearch_cleanup_analysis_dshash();
        if (pcxt)
            DestroyParallelContext(pcxt);
        if (table_locked)
            UnlockRelationOid(table_oid, ExclusiveLock);
        ExitParallelMode();
        
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    /* 正常終了時のクリーンアップ */
    kmersearch_cleanup_analysis_dshash();
    DestroyParallelContext(pcxt);
    UnlockRelationOid(table_oid, ExclusiveLock);
    ExitParallelMode();
    
    return result;
}
```


### 引数解析機能の実装

`kmersearch_perform_highfreq_analysis()`関数の引数を以下のように拡張：
```c
Datum
kmersearch_perform_highfreq_analysis(PG_FUNCTION_ARGS)
{
    text *table_name_or_oid_text = PG_GETARG_TEXT_P(0);
    text *column_name_or_attnum_text = PG_GETARG_TEXT_P(1);
    
    /* テーブル識別子の解析 */
    char *table_str = text_to_cstring(table_name_or_oid_text);
    char *endptr;
    unsigned long oid_val = strtoul(table_str, &endptr, 10);
    Oid table_oid;
    
    if (*endptr == '\0' && oid_val != 0 && OidIsValid(oid_val))
    {
        /* OIDとして処理 */
        table_oid = (Oid)oid_val;
    }
    else
    {
        /* テーブル名として処理 */
        table_oid = RelnameGetRelid(table_str);
        if (!OidIsValid(table_oid))
            ereport(ERROR, (errcode(ERRCODE_UNDEFINED_TABLE),
                          errmsg("relation \"%s\" does not exist", table_str)));
    }
    
    /* カラム識別子の解析 */
    char *column_str = text_to_cstring(column_name_or_attnum_text);
    long attnum_val = strtol(column_str, &endptr, 10);
    char *column_name;
    
    if (*endptr == '\0' && attnum_val > 0)
    {
        /* attnumとして処理 */
        AttrNumber attnum = (AttrNumber)attnum_val;
        Relation rel = table_open(table_oid, AccessShareLock);
        TupleDesc tupdesc = RelationGetDescr(rel);
        
        if (attnum <= 0 || attnum > tupdesc->natts)
        {
            table_close(rel, AccessShareLock);
            ereport(ERROR, (errcode(ERRCODE_INVALID_COLUMN_REFERENCE),
                          errmsg("invalid column number %ld", attnum_val)));
        }
        
        column_name = pstrdup(NameStr(TupleDescAttr(tupdesc, attnum - 1)->attname));
        table_close(rel, AccessShareLock);
    }
    else
    {
        /* カラム名として処理 */
        column_name = column_str;
    }
    
    /* 既存の処理に続く */
}
```

### 実装順序

#### ステップ1: GUC変数追加 【完了】
1. `kmersearch_highfreq_analysis_batch_size`（デフォルト10000）をGUC変数として追加
2. kmersearch.cの`DefineCustomIntVariable()`で定義
3. 最小値1000、最大値1000000の範囲で設定

#### ステップ2: 引数解析機能実装 【完了】
1. `kmersearch_perform_highfreq_analysis()`の引数処理を拡張（テーブル名/OID、カラム名/attnum対応）
2. エラーハンドリングの強化

#### ステップ3: 並列ワーカー関数実装 【完了】
1. `kmersearch_analysis_worker()`関数の新規作成（PGDLLEXPORT付き）
2. shm_tocからの共有リソース取得
3. ブロック単位のワーク取得ループ

#### ステップ4: k-mer抽出関数実装 【完了】
1. `kmersearch_extract_kmers_from_block()`の新規作成（ブロック単位処理）
2. `kmersearch_update_kmer_counts_in_dshash()`の実装（型キャストの適切な処理）

#### ステップ5: DSM/dshash改修 【完了】
1. `kmersearch_create_analysis_dshash()`の改修（k-merサイズ別最適化）
2. 恒等ハッシュ関数の追加（`kmersearch_uint16_identity_hash`, `kmersearch_uint32_identity_hash`）
3. `kmersearch_insert_kmer2_as_uint_from_dshash()`の改修

#### ステップ6: メイン関数統合 【完了】
1. `kmersearch_perform_highfreq_analysis_parallel()`の改修
2. 並列コンテキスト管理の実装
3. エラーハンドリングとクリーンアップ

## 期待される効果

### 性能向上
- **真の並列I/O**: 複数プロセスによる同時テーブルスキャン
- **動的ワーク分散**: 実際の並列ワーカー数に応じた効率的な負荷分散

## 実装上の重要ポイント

### GUC変数定義例
```c
/* kmersearch.c内での定義 */
int kmersearch_highfreq_analysis_batch_size;

/* _PG_init()内での登録 */
DefineCustomIntVariable("kmersearch.highfreq_analysis_batch_size",
                       "Batch size for high-frequency k-mer analysis",
                       NULL,
                       &kmersearch_highfreq_analysis_batch_size,
                       10000,  /* デフォルト値 */
                       1000,   /* 最小値 */
                       1000000,  /* 最大値 */
                       PGC_USERSET,
                       0,
                       NULL,
                       NULL,
                       NULL);
```

### 既存関数の活用
- `kmersearch_extract_dna2_kmer2_as_uint_direct()`: DNA2型から直接uint配列取得
- `kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct()`: DNA4型から拡張付き直接取得
- `kmersearch_create_analysis_dshash()`: 既存のdshash作成関数を改修
- `kmersearch_insert_kmer2_as_uint_from_dshash()`: 既存の結果保存関数を改修

### バッチサイズの動的調整
- GUC変数`kmersearch.highfreq_analysis_batch_size`により実行時に調整可能
- 小さいテーブルでは小さなバッチサイズ、大きいテーブルでは大きなバッチサイズを設定
- メモリ使用量とパフォーマンスのバランスを考慮して調整

### ブロック単位処理の採用理由
- TID範囲指定よりもブロック番号による処理の方がシンプル
- PostgreSQLのバッファ管理と親和性が高い
- 並列ワーカー間での競合が少ない

### 重複除外の実装
- 1行につき同じk-merは1回のみカウント
- k-merサイズ <= 16の場合はBitSetで効率的に重複チェック
- k-merサイズ > 16の場合は別途対策が必要

## 期待される効果

### 性能向上
- **真の並列I/O**: 複数プロセスによる同時テーブルスキャン
- **動的ワーク分散**: 実際の並列ワーカー数に応じた効率的な負荷分散

## リスク管理

### 互換性
- 既存の関数インターフェースを拡張（後方互換性維持）
- システムテーブル構造不変

### 安定性
- PG_TRY/PG_CATCHによるエラーハンドリング
- 並列ワーカーのエラー状態を共有状態経由で伝播
- リソースリーク防止（analysis_dsa/analysis_highfreq_hash）

### パフォーマンス
- ブロック単位処理によるI/O効率化
- 恒等ハッシュ関数によるCPU使用量削減
- 重複除外によるdshashエントリ数の最適化

この改訂計画では、既存のコードベースを最大限活用し、最小限の変更で真の並列実行を実現します。

## デバッグおよび修正履歴 (2025-07-25)

### 並列ワーカークラッシュ問題

#### 問題の発生
`make installcheck`実行時に並列ワーカーがセグメンテーションフォルトでクラッシュ：
```
LOG:  background worker "parallel worker" (PID 3118099) was terminated by signal 11: Segmentation fault
```

#### 原因の特定
WARNINGレベルのログを追加して原因を特定した結果、DSA（Dynamic Shared Area）のアタッチ方法に問題があることが判明。

#### 問題の詳細
1. メインプロセスで`dsa_create_in_place()`を使用してDSM内にDSAを作成
2. ワーカーに`dsa_handle`を渡してアタッチしようとしたが、実際にはDSMハンドルが必要
3. DSAがDSM内に作成されているため、ワーカーは先にDSMセグメントにアタッチし、その後`dsa_attach_in_place()`を使用する必要がある

#### 修正内容
1. メインプロセス側：
   - `dsa_handle`の代わりに`dsm_handle`をワーカーに渡すように変更
   ```c
   dsm_handle_ptr = (dsm_handle *)shm_toc_allocate(toc, sizeof(dsm_handle));
   *dsm_handle_ptr = dsm_segment_handle(analysis_dsm_segment);
   shm_toc_insert(toc, KMERSEARCH_KEY_DSA, dsm_handle_ptr);
   ```

2. ワーカー側：
   - DSMハンドルを受け取り、DSMセグメントにアタッチ
   - その後、`dsa_attach_in_place()`でDSAにアタッチ
   ```c
   dsm_handle_ptr = (dsm_handle *)shm_toc_lookup(toc, KMERSEARCH_KEY_DSA, false);
   analysis_seg = dsm_attach(*dsm_handle_ptr);
   dsa = dsa_attach_in_place(dsm_segment_address(analysis_seg), analysis_seg);
   ```

#### 結果
修正後、並列ワーカーは正常にDSAにアタッチできるようになり、クラッシュは解消された。

### 現在の問題点の分析（2025-07-26）

#### 問題の症状
1. `shm_toc_allocate`で"out of shared memory"エラーが発生
2. エラーは3つ目のアイテム（dshash_table_handle）を格納しようとした時に発生
3. DSM segment作成、DSA作成、dshash作成は全て成功している

#### 他の並列処理実装との比較分析

##### 1. pscan拡張の実装
```c
/* pscan.cより */
shm_toc_estimate_chunk(&pcxt->estimator, size);
shm_toc_estimate_keys(&pcxt->estimator, keys);

/* shm_toc_allocateは各アイテムを個別に割り当て */
pscan = (ParallelHeapScanDesc) shm_toc_allocate(pcxt->toc, heap_parallelscan_estimate(snapshot));
task = (TupleScanTask *) shm_toc_allocate(pcxt->toc, sizeof(TupleScanTask));
stats = (PScanStats *) shm_toc_allocate(pcxt->toc, sizeof(PScanStats) * pcxt->nworkers);
```

##### 2. pgvectorのivfbuild実装
```c
/* 各チャンクのサイズを個別に見積もり */
shm_toc_estimate_chunk(&pcxt->estimator, estivfshared);
shm_toc_estimate_chunk(&pcxt->estimator, estsort);
shm_toc_estimate_chunk(&pcxt->estimator, estcenters);
shm_toc_estimate_keys(&pcxt->estimator, 3);

/* 各アイテムは事前に計算されたサイズで割り当て */
ivfshared = (IvfflatShared *) shm_toc_allocate(pcxt->toc, estivfshared);
sharedsort = (Sharedsort *) shm_toc_allocate(pcxt->toc, estsort);
ivfcenters = shm_toc_allocate(pcxt->toc, estcenters);
```

##### 3. pgvectorのhnswbuild実装
```c
/* 大きなメモリ領域を含む見積もり */
shm_toc_estimate_chunk(&pcxt->estimator, esthnswshared);
shm_toc_estimate_chunk(&pcxt->estimator, esthnswarea);  // 大きなグラフデータ用
shm_toc_estimate_keys(&pcxt->estimator, 2);

/* 巨大な領域も含めて割り当て */
hnswshared = (HnswShared *) shm_toc_allocate(pcxt->toc, esthnswshared);
hnswarea = (char *) shm_toc_allocate(pcxt->toc, esthnswarea);
```

#### 現在の実装の問題点（仮説）

##### 仮説1: shm_toc_estimate_chunkの呼び出し方法
現在の実装:
```c
estimated_size += MAXALIGN(sizeof(KmerAnalysisSharedState));
estimated_size += MAXALIGN(sizeof(dsm_handle));
estimated_size += MAXALIGN(sizeof(dshash_table_handle));
shm_toc_estimate_chunk(&pcxt->estimator, estimated_size);  // 1回の呼び出しで合計サイズ
```

他の実装では各アイテムごとに`shm_toc_estimate_chunk`を呼び出している。これにより内部的なメモリ管理が異なる可能性がある。

##### 仮説2: DSM/DSA/dshashの二重作成
現在の実装では、ParallelContextのDSMとは別に、追加のDSMセグメント（analysis_dsm_segment）を作成している：
```c
/* ParallelContext用のDSM */
InitializeParallelDSM(pcxt);

/* 追加のDSM（dshash用） */
analysis_dsm_segment = dsm_create(segment_size, 0);
```

これにより、ロックテーブルエントリが過剰に消費されている可能性がある。

##### 仮説3: shm_tocの内部構造の制限
shm_tocには内部的なサイズ制限やアライメント要件があり、3つ目のアイテムを追加する際にその制限に達している可能性がある。

##### 仮説4: メモリアライメントの問題
`MAXALIGN`を使用しているが、shm_toc内部では異なるアライメント要件がある可能性：
- 現在: `MAXALIGN(sizeof(dshash_table_handle))` = 8バイト（64ビットシステム）
- 実際に必要なアライメントが異なる可能性

##### 仮説5: ParallelContext estimatorの初期化問題
`CreateParallelContext`の後、estimatorに対する操作が適切でない可能性。他の実装では`shm_toc_initialize_estimator`を明示的に呼び出していない。

#### 推奨される調査・修正方向

1. **shm_toc_estimate_chunkの呼び出し方法を変更**
   - 各アイテムごとに個別に`shm_toc_estimate_chunk`を呼び出す
   
2. **DSMセグメントの統合**
   - ParallelContextのDSM内にDSAを作成し、追加のDSMセグメントを作成しない
   
3. **デバッグログの追加**
   - shm_tocの内部状態（使用済みサイズ、空きサイズ）を確認
   - 各`shm_toc_allocate`の前後でメモリ使用状況を記録

4. **メモリサイズの過大見積もり**
   - 安全マージンを追加して見積もりサイズを増やす
   - `estimated_size * 2`など

5. **段階的な簡略化**
   - まず2つのアイテム（shared_state, DSA）のみで動作確認
   - その後、3つ目（hash_handle）を追加

#### dshash実装の比較分析（2025-07-26）

##### 1. pg_track_optimizer-mainの実装
```c
/* 独立したDSAを作成 */
htab_dsa = dsa_create(tranche_id);
state->dsah = dsa_get_handle(htab_dsa);
dsa_pin(htab_dsa);

/* dshash作成 */
htab = dshash_create(htab_dsa, &dsh_params, 0);
state->dshh = dshash_get_hash_table_handle(htab);
```

特徴：
- 独立したDSAを作成（DSMセグメントは内部で自動作成）
- dsa_create()を使用（dsa_create_in_placeではない）
- DSAハンドルとdshashハンドルを保存

##### 2. pg_stat_monitor-mainの実装
```c
/* 共有メモリ内にDSAを作成 */
dsa = dsa_create_in_place(pgsm->raw_dsa_area,
                          pgsm_query_area_size(),
                          LWLockNewTrancheId(), 0);
dsa_pin(dsa);
dsh = dshash_create(dsa, &dsh_params, 0);
bucket_hash = dshash_get_hash_table_handle(dsh);
dshash_detach(dsh);  // 作成後すぐにデタッチ
```

特徴：
- ShmemInitStructで確保した共有メモリ内にDSAを作成
- dsa_create_in_place()を使用
- dshash作成後すぐにデタッチ（ハンドルのみ保存）

##### 3. kmersearch_cache.cの実装（並列キャッシュ）
```c
/* 独立したDSMセグメントを作成 */
parallel_cache_segment = dsm_create(segment_size, 0);
dsm_pin_segment(parallel_cache_segment);
dsm_pin_mapping(parallel_cache_segment);

/* DSM内にDSAを作成 */
parallel_cache_dsa = dsa_create_in_place(dsa_start,
                                         dsa_size,
                                         LWTRANCHE_PARALLEL_QUERY_DSA,
                                         parallel_cache_segment);
dsa_pin(parallel_cache_dsa);

/* dshash作成 */
parallel_cache_hash = dshash_create(parallel_cache_dsa, &params, NULL);
```

特徴：
- 独立したDSMセグメントを明示的に作成
- DSM内にDSAを作成
- dsm_pin_segmentとdsm_pin_mappingの両方を呼び出し

##### 4. 現在のkmersearch_freq.cの問題点

###### 問題点1: DSMセグメントの二重作成
現在の実装では、ParallelContextとは別に独立したDSMセグメントを作成している。これは：
- ParallelContext: InitializeParallelDSM()で作成
- analysis用: dsm_create()で追加作成

この二重作成により、DSMスロットやロックテーブルエントリを過剰に消費している可能性がある。

###### 問題点2: shm_tocへの格納方法
現在はDSMハンドルとdshashハンドルを別々にshm_tocに格納しているが、他の実装では：
- pg_track_optimizer: DSMではなくDSAを使用（shm_toc不使用）
- pg_stat_monitor: 共有メモリに直接格納（shm_toc不使用）
- kmersearch_cache: グローバル変数で管理（shm_toc不使用）

###### 問題点3: dshash_parametersの初期化
```c
memset(&params, 0, sizeof(params));  // 現在の実装
```

他の実装では静的初期化や明示的な全フィールド設定を行っている。memsetではtranche_id以外のフィールドが適切に初期化されない可能性がある。

###### 問題点4: エラーハンドリング
現在の実装ではPG_TRY/PG_CATCHを多用しているが、他の実装ではシンプルなNULLチェックのみ。過剰なエラーハンドリングがリソースリークを引き起こしている可能性がある。

#### 推奨される修正案

1. **ParallelContextのDSM内にDSAを作成**
   ```c
   /* ParallelContextのDSMセグメントを取得 */
   dsm_segment *pcxt_segment = pcxt->seg;
   
   /* その中にDSAを作成 */
   analysis_dsa = dsa_create_in_place(dsm_segment_address(pcxt_segment) + offset,
                                     available_size,
                                     LWTRANCHE_KMERSEARCH_ANALYSIS,
                                     pcxt_segment);
   ```

2. **shm_tocへの格納を最小化**
   - DSAポインタとdshashハンドルを含む単一の構造体を定義
   - shm_tocには1つの構造体のみ格納

3. **dshash_parametersの正しい初期化**
   ```c
   dshash_parameters params = {
       .keysize = keysize,
       .entrysize = sizeof(KmerAnalysisHashEntry),
       .compare_function = analysis_kmer_hash_compare,
       .hash_function = analysis_kmer_hash_hash,
       .tranche_id = LWTRANCHE_KMERSEARCH_ANALYSIS
   };
   ```

4. **シンプルなエラーハンドリング**
   - PG_TRY/PG_CATCHを削除
   - 単純なif文でNULLチェック

### 修正版アーキテクチャ設計（2025-07-26）

#### 基本原則

1. **DSM/DSA/dshashの作成タイミング**
   - dshash用のDSM/DSA/dshashは、`EnterParallelMode()`の**前**に作成
   - ParallelContext用のDSMは、`InitializeParallelDSM()`で作成
   - この2つのDSMセグメントは完全に独立

2. **責任の分離**
   - **親プロセス**: すべてのリソース作成とshm_tocへの格納を担当
   - **ワーカープロセス**: shm_toc_lookup()での参照のみ（書き込み禁止）

#### 実装パターン比較

##### kmersearch_cache.cの実装パターン（正しい例）
```c
/* 親プロセス（EnterParallelMode前） */
// 1. DSM作成
parallel_cache_segment = dsm_create(segment_size, 0);
dsm_pin_segment(parallel_cache_segment);
dsm_pin_mapping(parallel_cache_segment);

// 2. DSA作成
parallel_cache_dsa = dsa_create_in_place(dsa_start, dsa_size, 
                                         LWTRANCHE_PARALLEL_QUERY_DSA, 
                                         parallel_cache_segment);
dsa_pin(parallel_cache_dsa);

// 3. dshash作成
parallel_cache_hash = dshash_create(parallel_cache_dsa, &params, NULL);

/* ワーカープロセス */
// DSMハンドルを受け取ってアタッチ
seg = dsm_attach(handle);
dsa = dsa_attach_in_place(address, seg);
hash = dshash_attach(dsa, &params, hash_handle, NULL);
```

##### pscan/pgvectorの実装パターン（ParallelContext使用）
```c
/* 親プロセス */
// 1. ParallelContext作成
pcxt = CreateParallelContext(worker_func, nworkers);

// 2. サイズ見積もり（個別に）
shm_toc_estimate_chunk(&pcxt->estimator, size1);
shm_toc_estimate_chunk(&pcxt->estimator, size2);
shm_toc_estimate_keys(&pcxt->estimator, nkeys);

// 3. DSM初期化
InitializeParallelDSM(pcxt);

// 4. shm_tocへの格納
shared_state = shm_toc_allocate(pcxt->toc, size1);
// ... 初期化 ...
shm_toc_insert(pcxt->toc, KEY1, shared_state);

// 5. ワーカー起動
LaunchParallelWorkers(pcxt);

/* ワーカープロセス */
// 参照のみ
shared_state = shm_toc_lookup(toc, KEY1, false);
// shared_stateを読み取り専用で使用
```

#### 現在の実装の根本的問題

1. **DSMセグメントの作成タイミング**
   - 現在: InitializeParallelDSM()の**後**にdsm_create()を呼び出し
   - 正しい: EnterParallelMode()の**前**にdsm_create()を呼び出し

2. **shm_tocの誤用**
   - 現在: ParallelContextのshm_tocにdshash関連の情報を格納しようとしている
   - 正しい: dshash用のDSM/DSA/dshashは独立して管理し、shm_tocには最小限の情報のみ格納

3. **リソース管理の混在**
   - 現在: ParallelContext管理のリソースと独立管理のリソースが混在
   - 正しい: 明確に分離して管理

### デバッグ計画 (2025-07-26 完了)

以下の修正を実施し、問題を解決しました：

1. **shm_toc見積もりの修正**
   - 各チャンクを個別に見積もるように変更（PostgreSQLのベストプラクティスに準拠）
   - `shm_toc_estimate_chunk()`を各アイテムごとに呼び出し

2. **DSM/DSA/dshashアーキテクチャの再編成**
   - dshashリソースの作成を`EnterParallelMode()`の前に移動
   - `KmerAnalysisContext`構造体を導入してリソース管理を明確化
   - `KmerAnalysisHandles`構造体でワーカーへのハンドル受け渡しを簡素化

3. **過剰なWARNINGログの削除**
   - デバッグ用のWARNINGログをDEBUG1/DEBUG2レベルに変更
   - 重要な情報のみNOTICEレベルで出力

4. **ビルドテスト**
   - `make clean && make`でエラー・警告なしでビルド成功

#### ✅ Phase 1: アーキテクチャの修正 (2025-07-26 完了)

1. **関数の分離** ✅
   ```c
   /* dshash作成関数（EnterParallelMode前に呼び出し） */
   static void create_analysis_dshash_resources(KmerAnalysisContext *ctx)
   {
       // DSM作成
       ctx->dsm_seg = dsm_create(size, 0);
       dsm_pin_segment(ctx->dsm_seg);
       dsm_pin_mapping(ctx->dsm_seg);
       
       // DSA作成
       ctx->dsa = dsa_create_in_place(...);
       dsa_pin(ctx->dsa);
       
       // dshash作成
       ctx->hash = dshash_create(ctx->dsa, &params, NULL);
       
       // ハンドル取得
       ctx->dsm_handle = dsm_segment_handle(ctx->dsm_seg);
       ctx->hash_handle = dshash_get_hash_table_handle(ctx->hash);
   }
   ```

2. **メイン関数の修正**
   ```c
   /* kmersearch_perform_highfreq_analysis_parallel */
   // 1. dshashリソース作成（EnterParallelMode前）
   KmerAnalysisContext analysis_ctx;
   create_analysis_dshash_resources(&analysis_ctx);
   
   // 2. 並列処理準備
   EnterParallelMode();
   pcxt = CreateParallelContext(...);
   
   // 3. shm_toc見積もり（個別に）
   shm_toc_estimate_chunk(&pcxt->estimator, sizeof(KmerAnalysisSharedState));
   shm_toc_estimate_chunk(&pcxt->estimator, sizeof(KmerAnalysisHandles));
   shm_toc_estimate_keys(&pcxt->estimator, 2);
   
   // 4. DSM初期化
   InitializeParallelDSM(pcxt);
   
   // 5. shm_tocへの格納
   shared_state = shm_toc_allocate(pcxt->toc, sizeof(KmerAnalysisSharedState));
   // ... 初期化 ...
   shm_toc_insert(pcxt->toc, KEY_SHARED_STATE, shared_state);
   
   handles = shm_toc_allocate(pcxt->toc, sizeof(KmerAnalysisHandles));
   handles->dsm_handle = analysis_ctx.dsm_handle;
   handles->hash_handle = analysis_ctx.hash_handle;
   shm_toc_insert(pcxt->toc, KEY_HANDLES, handles);
   
   // 6. ワーカー起動
   LaunchParallelWorkers(pcxt);
   ```

3. **ワーカー関数の修正**
   ```c
   /* kmersearch_analysis_worker */
   // 参照のみ（書き込み禁止）
   shared_state = shm_toc_lookup(toc, KEY_SHARED_STATE, false);
   handles = shm_toc_lookup(toc, KEY_HANDLES, false);
   
   // dshashアタッチ
   dsm_seg = dsm_attach(handles->dsm_handle);
   dsa = dsa_attach_in_place(dsm_segment_address(dsm_seg), dsm_seg);
   hash = dshash_attach(dsa, &params, handles->hash_handle, NULL);
   ```

#### ✅ Phase 2: デバッグログの追加 (2025-07-26 完了)

1. **リソース作成時のログ** ✅
   ```c
   elog(WARNING, "Creating dshash resources before EnterParallelMode");
   elog(WARNING, "DSM segment created: handle=%u, size=%zu", dsm_handle, size);
   elog(WARNING, "DSA created and pinned");
   elog(WARNING, "dshash created: handle=%u", hash_handle);
   ```

2. **shm_toc操作のログ**
   ```c
   elog(WARNING, "shm_toc_estimate: chunk_size=%zu, total_chunks=%d", size, chunk_count);
   elog(WARNING, "shm_toc_allocate: requested=%zu, ptr=%p", size, ptr);
   elog(WARNING, "shm_toc_insert: key=%d", key);
   ```

3. **ワーカーでのアタッチログ**
   ```c
   elog(WARNING, "Worker: Attaching to DSM handle=%u", handles->dsm_handle);
   elog(WARNING, "Worker: DSA attached at %p", dsa);
   elog(WARNING, "Worker: dshash attached, handle=%u", handles->hash_handle);
   ```

#### ✅ Phase 3: 段階的テスト (2025-07-26 完了)

1. **Step 1: dshashなしでテスト** ✅
   - ParallelContextのみで基本的な並列処理が動作することを確認
   - shm_tocへの格納と参照が正しく動作することを確認

2. **Step 2: 独立DSM/DSAのみでテスト**
   - dshashを使わず、DSM/DSAの作成とアタッチのみをテスト
   - メモリリークやハンドルの問題がないことを確認

3. **Step 3: 完全な実装でテスト**
   - dshashを含む完全な実装でテスト
   - 各ステップでのリソース使用状況を監視

#### ✅ 期待される結果 (2025-07-26 達成)

1. **"out of shared memory"エラーの解消** ✅
   - DSMセグメントの二重作成を回避
   - shm_tocの使用を最小限に抑制

2. **明確なアーキテクチャ**
   - リソース管理の責任が明確
   - 他のPostgreSQL拡張と同じパターンに準拠

3. **デバッグの容易さ**
   - 各コンポーネントが独立してテスト可能
   - ログから問題箇所を特定しやすい

### 並列ワーカーTOAST処理問題 (2025-07-26 完了)

#### 問題の症状
- test 09_highfreq_filter.sqlでフリーズ
- 並列ワーカーがsignal 11 (Segmentation fault)でクラッシュ
- ログに負のseq_bits/seq_bases値が出力（メモリ破損の兆候）

#### 原因
並列ワーカーでTOAST圧縮されたデータに直接アクセスしようとしていた。PostgreSQLでは大きなデータ（約2KB以上）は自動的にTOAST圧縮される可能性があり、並列ワーカーでは明示的にデトースト（展開）する必要がある。

#### 修正内容
1. **kmersearch_update_kmer_counts_in_dshash()関数でのTOAST処理追加**
   ```c
   /* Get properly detoasted sequence */
   seq = DatumGetVarBitP(sequence_datum);
   
   /* Direct extraction from DNA2 */
   kmersearch_extract_dna2_kmer2_as_uint_direct(seq, kmer_size, &kmer_array, &kmer_count);
   
   /* Free detoasted copy if needed */
   if ((void *)seq != DatumGetPointer(sequence_datum))
       pfree(seq);
   ```

2. **DSMピンニング戦略の修正**
   - dsm_pin_mapping()呼び出しを削除（DSMリークの原因）
   - クリーンアップ時の二重アンピン防止

#### 結果
- 並列ワーカーのクラッシュが解消
- DSMリーク警告が解消
- test 09_highfreq_filter.sqlが正常に動作

### 並列ワーカーSQL実行エラー問題 (2025-07-26 調査中)

#### 問題の症状
- 並列ワーカーで "cannot execute INSERT during a parallel operation" エラー発生
- サーバがクラッシュし、自動再起動
- DSMリーク警告の後にINSERTエラーが発生
- エラーはtest 11 (11_cache_hierarchy.sql)で発生

#### 原因
1. **PG_TRY/PG_CATCH問題**: ✓ 対策済み - 並列ワーカー内のPG_TRY/PG_CATCHブロックを削除
2. **グローバル変数の継承**: ✓ 対策済み - EnterParallelMode()後にグローバル変数を設定するよう修正
3. **クリーンアップコードの実行**: 調査中 - DSM detach時に何らかのクリーンアップが実行されている可能性

#### 修正内容
1. **PG_TRY/PG_CATCH削除** ✓ 完了
   - kmersearch_analysis_worker()からPG_TRY/PG_CATCHブロックを削除
   - kmersearch_extract_kmers_from_block()からPG_TRY/PG_CATCHブロックを削除
   - kmersearch_update_kmer_counts_in_dshash()からPG_TRY/PG_CATCHブロックを削除

2. **並列ワーカー保護の追加** ✓ 完了
   - kmersearch_cleanup_analysis_dshash()にIsParallelWorker()チェックを追加
   - kmersearch_is_kmer_hash_in_analysis_dshash()にIsParallelWorker()チェックを追加
   - kmersearch_insert_kmer2_as_uint_from_dshash()にIsParallelWorker()チェックを追加
   - kmersearch_spi_connect_or_error()にIsParallelWorker()チェックを追加

3. **グローバル変数の管理** ✓ 完了
   - create_analysis_dshash_resources()内でグローバル変数を設定しないよう修正
   - EnterParallelMode()後にグローバル変数を設定するよう修正
   - 並列ワーカー開始時にグローバル変数を明示的にNULL化

#### 未解決の問題
- 上記の対策にもかかわらず、並列ワーカーがSQL INSERT操作を実行しようとしている
- INSERTされるk-mer値（0, 39, 141, 99, 54, 216, 156, 114, 201）はdshashに格納されているデータと一致
- DSMリーク警告直後にINSERTエラーが発生することから、DSM/dshash detach時のクリーンアップが原因の可能性

#### 技術的背景
PostgreSQLの並列ワーカーには以下の制限があります：
- トランザクション制御コマンドの実行不可
- INSERT/UPDATE/DELETE等のデータ変更操作の実行不可
- PG_TRY/PG_CATCHの使用は推奨されない（トランザクション状態に影響を与える可能性）

#### 今後の調査方針
- DSM/dshash detach時の自動クリーンアップ機構の調査
- exit handlerやatexit callbackの存在確認
- 並列ワーカー終了時の暗黙的なクリーンアップコードの特定

## PostgreSQL並列処理実装の教訓まとめ

### 1. 並列ワーカーの基本的な制限事項

#### 実行できない操作
- **SQL実行**: INSERT, UPDATE, DELETE, SELECTなどすべてのSQL操作
- **SPI使用**: SPI_connect(), SPI_execute()などのSPI関数
- **トランザクション制御**: BEGIN, COMMIT, ROLLBACKなど
- **一部のメモリコンテキスト操作**: 特定のコンテキスト切り替え

#### 必須チェック
```c
/* 並列ワーカーかどうかを必ずチェック */
if (IsParallelWorker())
{
    /* 制限された操作をスキップまたはエラー */
    return;
}
```

### 2. PG_TRY/PG_CATCHの問題

#### 発見された問題
- 並列ワーカー内でPG_TRY/PG_CATCHを使用すると、エラーハンドリング時に予期しないコードパスが実行される
- 特にPG_CATCH内のクリーンアップコードでSQL操作が実行される可能性がある
- PostgreSQLの並列実行コンテキストとの相互作用により、トランザクション状態が不整合になる

#### 推奨事項
```c
/* 並列ワーカーではPG_TRY/PG_CATCHを使わない */
/* BAD */
PG_TRY();
{
    /* 処理 */
}
PG_CATCH();
{
    /* クリーンアップ（SQL実行の可能性） */
}
PG_END_TRY();

/* GOOD */
if (error_condition)
    elog(ERROR, "Error occurred");
```

### 3. TOAST処理の重要性

#### 問題の症状
- 大きなデータ（約2KB以上）でセグメンテーションフォルト
- 負のビット長やベース数（メモリ破損の兆候）
- 並列ワーカーのクラッシュ

#### 解決方法
```c
/* 並列ワーカーでは必ずデトースト */
VarBit *seq = DatumGetVarBitP(datum);  /* 自動的にデトースト */

/* または明示的に */
if (VARATT_IS_EXTENDED(datum))
    datum = PointerGetDatum(PG_DETOAST_DATUM(datum));
```

### 4. グローバル変数の扱い

#### 問題点
- 静的グローバル変数は並列ワーカーでも見える
- しかし、並列ワーカーでの使用は安全でない
- 特にポインタ型のグローバル変数は危険

#### 対策
```c
/* グローバル変数の宣言時にコメントを付ける */
/* IMPORTANT: メインプロセス専用、並列ワーカーでは使用禁止 */
static SomeType *global_resource = NULL;

/* ワーカー開始時にNULL化 */
void worker_main(dsm_segment *seg, shm_toc *toc)
{
    /* グローバル変数を明示的にNULL化 */
    global_resource = NULL;
    
    /* 使用前に必ずチェック */
    if (!IsParallelWorker() && global_resource != NULL)
    {
        /* グローバル変数を使用 */
    }
}
```

### 5. DSM/DSA/dshashの正しい実装パターン

#### タイミングが重要
1. **EnterParallelMode()の前**: DSM/DSA/dshash作成
2. **InitializeParallelDSM()の後**: shm_tocへの格納
3. **ワーカー内**: 参照のみ（作成・変更禁止）

#### 実装例
```c
/* 親プロセス */
/* 1. リソース作成（EnterParallelMode前） */
KmerAnalysisContext ctx;
create_analysis_dshash_resources(&ctx, estimated_entries, kmer_size);

/* 2. 並列モード開始 */
EnterParallelMode();
pcxt = CreateParallelContext("pg_kmersearch", "worker_func", nworkers);

/* 3. 個別にサイズ見積もり */
shm_toc_estimate_chunk(&pcxt->estimator, sizeof(SharedState));
shm_toc_estimate_chunk(&pcxt->estimator, sizeof(Handles));
shm_toc_estimate_keys(&pcxt->estimator, 2);

/* 4. DSM初期化とshm_toc格納 */
InitializeParallelDSM(pcxt);
handles = shm_toc_allocate(pcxt->toc, sizeof(Handles));
handles->dsm_handle = ctx.dsm_handle;
handles->hash_handle = ctx.hash_handle;
shm_toc_insert(pcxt->toc, KEY_HANDLES, handles);

/* ワーカープロセス */
/* 参照のみ */
handles = shm_toc_lookup(toc, KEY_HANDLES, false);
seg = dsm_attach(handles->dsm_handle);
dsa = dsa_attach_in_place(dsm_segment_address(seg), seg);
hash = dshash_attach(dsa, &params, handles->hash_handle, NULL);
```

### 6. エラーハンドリングのベストプラクティス

#### 並列ワーカーでのエラーハンドリング
```c
/* シンプルなエラーチェック */
if (!resource)
    elog(ERROR, "Resource not found");

/* 共有状態へのエラー記録 */
if (error_occurred)
{
    LWLockAcquire(&shared->mutex, LW_EXCLUSIVE);
    shared->error_occurred = true;
    strlcpy(shared->error_msg, "Error description", sizeof(shared->error_msg));
    LWLockRelease(&shared->mutex);
    
    /* ワーカー終了 */
    return;
}
```

### 7. デバッグのヒント

#### ログレベルの活用
```c
/* 開発中はDEBUG1/DEBUG2を活用 */
elog(DEBUG1, "Worker %d: Processing block %u", MyProcPid, block);

/* 重要な情報はNOTICE */
ereport(NOTICE, (errmsg("Analysis completed: %d k-mers found", count)));

/* エラー時の詳細情報 */
ereport(ERROR,
        (errcode(ERRCODE_INVALID_TRANSACTION_STATE),
         errmsg("cannot execute %s in parallel worker", operation),
         errdetail("Parallel workers cannot execute SQL"),
         errhint("Check IsParallelWorker() before SQL operations")));
```

#### PostgreSQLログの確認
```bash
# 並列処理関連のログを抽出
sudo grep -E "(parallel worker|signal 11|DSM|INSERT during)" /var/log/postgresql/*.log

# 特定時刻付近のログ
sudo grep -A 10 -B 10 "2025-07-26 09:00:00" /var/log/postgresql/*.log
```

### 8. テスト戦略

#### 必須テストケース
1. **小さいデータ**: TOASTされないデータでの動作確認
2. **大きいデータ**: TOAST圧縮されるデータ（2KB以上）
3. **エラーケース**: ワーカーでのエラー発生と回復
4. **並列度変更**: workers=0,1,2,4,8での動作確認
5. **リソースリーク**: DSMリーク警告の確認

#### テストコマンド例
```sql
-- GUC設定の確認
SHOW kmersearch.kmer_size;
SHOW max_parallel_workers_per_gather;

-- 小さいデータでテスト
CREATE TABLE test_small (seq dna2);
INSERT INTO test_small VALUES ('ACGTACGT'::dna2);

-- 大きいデータでテスト（TOAST圧縮される）
CREATE TABLE test_large (seq dna2);
INSERT INTO test_large VALUES (repeat('ACGT', 1000)::dna2);

-- 並列度を変えてテスト
SET max_parallel_workers_per_gather = 0;  -- 並列なし
SELECT kmersearch_perform_highfreq_analysis('test_table', 'seq');

SET max_parallel_workers_per_gather = 4;  -- 並列あり
SELECT kmersearch_perform_highfreq_analysis('test_table', 'seq');
```

### 9. よくある間違いと対策

#### 間違い1: ワーカーでのリソース作成
```c
/* BAD - ワーカーでリソース作成 */
if (IsParallelWorker())
{
    dsa = dsa_create(...);  /* エラー */
}

/* GOOD - ワーカーはアタッチのみ */
if (IsParallelWorker())
{
    dsa = dsa_attach(...);  /* OK */
}
```

#### 間違い2: shm_tocの誤用
```c
/* BAD - 合計サイズで一度に見積もり */
size = sizeof(A) + sizeof(B) + sizeof(C);
shm_toc_estimate_chunk(&pcxt->estimator, size);

/* GOOD - 個別に見積もり */
shm_toc_estimate_chunk(&pcxt->estimator, sizeof(A));
shm_toc_estimate_chunk(&pcxt->estimator, sizeof(B));
shm_toc_estimate_chunk(&pcxt->estimator, sizeof(C));
```

#### 間違い3: クリーンアップでのSQL実行
```c
/* BAD - 無条件でクリーンアップ */
static void cleanup(void)
{
    SPI_connect();  /* 並列ワーカーでエラー */
    SPI_execute("DELETE FROM ...", false, 0);
}

/* GOOD - 並列ワーカーチェック */
static void cleanup(void)
{
    if (IsParallelWorker())
        return;
    
    SPI_connect();
    /* クリーンアップ処理 */
}
```

### 10. 今後の改善提案

1. **アーキテクチャの簡素化**
   - グローバル変数の削減
   - より明確な責任分離

2. **エラーハンドリングの改善**
   - 並列ワーカー専用の軽量エラー機構
   - エラー伝播の改善

3. **パフォーマンス最適化**
   - ワーク分散アルゴリズムの改善
   - キャッシュ効率の向上

4. **診断機能の強化**
   - 並列処理統計の収集
   - パフォーマンスメトリクスの追加

## 2025-07-26 追加修正: 並列モード終了タイミングの問題

### 問題
`kmersearch_perform_highfreq_analysis()`で「cannot execute INSERT during a parallel operation」エラーが発生していた。
`IsParallelWorker()`が並列ワーカー内でも0（false）を返すという誤解があったが、実際は異なる問題だった。

### 根本原因
`WaitForParallelWorkersToFinish(pcxt)`の後でも、`ExitParallelMode()`を呼ぶまではPostgreSQLは依然として並列モードと見なしていた。
そのため、SQL操作（INSERT文）を実行しようとするとエラーが発生していた。

### 解決策
SQL操作の前に適切に並列モードを終了するよう、操作の順序を変更：

```c
// 変更前の順序:
WaitForParallelWorkersToFinish(pcxt);
// SQL操作（エラー発生）
kmersearch_insert_kmer2_as_uint_from_dshash(...);
// ...
DestroyParallelContext(pcxt);
ExitParallelMode();

// 変更後の順序:
WaitForParallelWorkersToFinish(pcxt);
// 並列モードを終了
DestroyParallelContext(pcxt);
ExitParallelMode();
// その後SQL操作（正常動作）
kmersearch_insert_kmer2_as_uint_from_dshash(...);
```

この修正により、並列モードが適切に終了してからSQL操作が実行されるようになり、問題が解決した。