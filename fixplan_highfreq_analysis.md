# 高頻度k-mer解析の並列実行化計画

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

#### ステップ1: GUC変数追加
1. `kmersearch_highfreq_analysis_batch_size`（デフォルト10000）をGUC変数として追加
2. kmersearch.cの`DefineCustomIntVariable()`で定義
3. 最小値1000、最大値1000000の範囲で設定

#### ステップ2: 引数解析機能実装
1. `kmersearch_perform_highfreq_analysis()`の引数処理を拡張（テーブル名/OID、カラム名/attnum対応）
2. エラーハンドリングの強化

#### ステップ3: 並列ワーカー関数実装
1. `kmersearch_analysis_worker()`関数の新規作成（PGDLLEXPORT付き）
2. shm_tocからの共有リソース取得
3. ブロック単位のワーク取得ループ

#### ステップ4: k-mer抽出関数実装
1. `kmersearch_extract_kmers_from_block()`の新規作成（ブロック単位処理）
2. `kmersearch_update_kmer_counts_in_dshash()`の改修（重複除外機能追加）

#### ステップ5: DSM/dshash改修
1. `kmersearch_create_analysis_dshash()`の改修（k-merサイズ別最適化）
2. 恒等ハッシュ関数の追加
3. `kmersearch_insert_kmer2_as_uint_from_dshash()`の改修

#### ステップ6: メイン関数統合
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