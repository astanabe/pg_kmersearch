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

## 概要

現在の`kmersearch_perform_highfreq_analysis()`は疑似並列処理（単一プロセスによる順次処理）となっており、真の並列実行が行われていません。本計画では、PostgreSQLの標準並列実行機能とDSM/dshashを使用した真の並列処理に変更し、真の並列ワーカー間の処理フローを整理します。

## システム要件

**PostgreSQLバージョン要件**: この並列実行化計画は**PostgreSQL 16以降**にのみ対応します。PostgreSQL 15以前のバージョンに対する互換性コードは実装しません。

## 現在の問題点

### 1. 疑似並列処理
- `KmerWorkerState`配列を作成しているが、実際はメインプロセスが順次処理
- PostgreSQLの並列ワーカープロセス機能を使用していない
- `result.parallel_workers_used`は設定値を返すが実際の並列実行なし

### 2. 非効率な処理フロー
- SQL集計結果を一旦dshashに格納する方式が不適切
- 並列ワーカー間でのk-mer出現回数集計に最適化されていない
- 真の並列I/Oが発生しない

## 改訂された並列処理フロー

### 親プロセス（リーダープロセス）の責務

1. **テーブルロック管理**
   - 対象テーブルを`EXCLUSIVE MODE`でロック
   - 解析完了後にロックを解放

2. **事前検証**
   - 対象カラムがDNA2型またはDNA4型かチェック
   - 違う型の場合はエラーを出力して停止

3. **並列実行環境の準備**
   - 並列ワーカーにパラメータを伝えるためのDSM（Dynamic Shared Memory）を作成
   - k-mer出現回数を集計するためのdshashテーブルを作成
   - dshashテーブルには`kmer2_as_uint`と、その`kmer2_as_uint`の出現した行数を保存

4. **並列ワーカー管理**
   - PostgreSQLによって決定される並列ワーカープロセス数を取得（計画数であり実際の起動数と異なる場合あり）
   - 並列ワーカーの担当行数決定：最大で`kmersearch_highfreq_analysis_batch_size`行（デフォルト10000行）

5. **動的ワーク分散システムの構築**
   - DSM上に以下の動的分散制御パラメータを配置：
     - 対象テーブルOID（引数から解決済み）
     - 対象カラムのattnum（引数から解決済み）
     - 次の処理対象となる最初の行のTID（動的更新）
     - 次の処理対象となる最後の行のTID（動的更新）
     - 全行処理済みフラグ（boolean、初期値false）

6. **PostgreSQL標準並列実行機能での並列ワーカー起動**
   - 準備が完了したら標準の並列実行機能で並列ワーカーを起動

7. **動的ワーク分散の実行**
   - 全テーブル行を`kmersearch_highfreq_analysis_batch_size`サイズの範囲に分割
   - 最初のTID範囲をDSM上に設定
   - 並列ワーカーがTID範囲を取得したことを検知（DSM上のTID範囲が空になる）
   - 次のTID範囲が存在する場合、DSM上に新しいTID範囲を設定
   - 全範囲の分散が完了したら`全行処理済みフラグ`をtrueに設定

8. **結果統合と最終処理**
   - 全並列ワーカーの処理完了を待機
   - dshashテーブルの内容を分析して高頻度出現k-merを特定
   - `kmersearch_highfreq_kmer`テーブルに結果を書き出し
   - `kmersearch_highfreq_kmer_meta`テーブルにメタデータを書き込み

### 並列ワーカーの責務

1. **動的ワーク取得ループ**
   - DSM上の共有パラメータを確認：
     - 対象テーブルOID（関数引数から解決済み）
     - 対象カラムのattnum（関数引数から解決済み）
   - 動的ワーク取得処理を繰り返し実行

2. **TID範囲の取得と確保**
   - DSM上の最初の行のTID、最後の行のTIDを読み取り
   - TID範囲が設定されている場合、DSM上のTID範囲を空に設定（短時間の排他制御）
   - 取得したTID範囲で独立k-mer抽出処理を実行（排他制御外で並列実行）

3. **独立k-mer抽出処理**
   - 取得したTID範囲の対象テーブル・対象カラムから`kmer2_as_uint`を抽出
   - 重複除外した上で、dshashテーブルの`kmer2_as_uint`出現行数をカウントアップ

4. **次ワーク確認と待機**
   - 処理完了後、DSM上のTID範囲を再確認
   - TID範囲が空で`全行処理済みフラグ`がfalseの場合：1秒sleep後に再確認
   - TID範囲が空で`全行処理済みフラグ`がtrueの場合：ワーカー終了

5. **結果共有**
   - 全ての抽出結果をdshashテーブルに反映
   - 並列ワーカー間での同期はdshashテーブル経由で実現

## DSM/dshash実装参考

`kmersearch_parallel_highfreq_kmer_cache_load()`関数の実装を参考にしてDSM/dshash機能を実装します：

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
    int         kmer_size;               /* k-merサイズ */
    int         batch_size;              /* バッチサイズ（kmersearch_highfreq_analysis_batch_size） */
    bool        all_processed;           /* 全行処理済みフラグ */
    ItemPointerData next_start_tid;      /* 次の処理対象開始TID */
    ItemPointerData next_end_tid;        /* 次の処理対象終了TID */
    bool        tid_range_available;     /* TID範囲が利用可能かどうか */
    bool        worker_error_occurred;   /* 並列ワーカーエラー発生フラグ */
    char        error_message[256];      /* エラーメッセージ */
} KmerAnalysisSharedState;
```

#### dshash エントリ構造体
```c
/* 既存のKmerData構造体を活用（k-merサイズに応じた主キー） */
/* CompactKmerFreq構造体を活用（k-mer + 出現回数） */

/* k-merサイズに応じたdshashエントリ型の決定 */
typedef struct AnalysisKmerHashEntry
{
    KmerData    kmer_data;               /* k-mer（主キー、既存構造体） */
    int         occurrence_count;        /* 出現行数 */
} AnalysisKmerHashEntry;

/* 
 * 注意：実際の実装では以下のようにk-merサイズに応じて構造体を選択
 * - kmer_size <= 8:  uint16をキーとして使用
 * - kmer_size <= 16: uint32をキーとして使用  
 * - kmer_size <= 32: uint64をキーとして使用
 * KmerData構造体内のkmer2_as_uintフィールドを適切なサイズにキャストして使用
 */
```

### 変更計画

#### フェーズ1: DSM/dshashベース並列実行環境の構築

##### 1.1 並列実行環境の初期化

**DSM/dshash初期化関数**
```c
static dshash_table *
kmersearch_create_analysis_dshash(dsm_segment **dsm_seg, dsa_area **dsa, int kmer_size)
{
    Size        segsize;
    dsm_segment *seg;
    dsa_area   *area;
    dshash_table *hash;
    dshash_parameters params;
    Size        key_size;
    
    /* k-merサイズに応じたキーサイズ決定 */
    if (kmer_size <= 8)
        key_size = sizeof(uint16);
    else if (kmer_size <= 16)
        key_size = sizeof(uint32);
    else if (kmer_size <= 32)
        key_size = sizeof(uint64);
    else
        elog(ERROR, "Invalid kmer_size: %d (must be <= 32)", kmer_size);
    
    /* DSMセグメントサイズ計算 */
    segsize = MAXALIGN(sizeof(dsa_area_control)) + (1024 * 1024); /* 1MB初期サイズ */
    
    /* DSMセグメント作成 */
    seg = dsm_create(segsize, 0);
    area = dsa_create_in_place(dsm_segment_address(seg), segsize,
                              LWTRANCHE_KMERSEARCH_DSA, dsm_segment_handle(seg));
    
    /* dshashパラメータ設定 */
    memset(&params, 0, sizeof(params));
    params.key_size = key_size;
    
    /* k-merサイズに応じてエントリサイズと関数を設定 */
    if (kmer_size <= 8)
    {
        params.entry_size = sizeof(uint16) + sizeof(int); /* kmer16 + occurrence_count */
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint16_identity_hash; /* 恒等関数 */
    }
    else if (kmer_size <= 16)
    {
        params.entry_size = sizeof(uint32) + sizeof(int); /* kmer32 + occurrence_count */
        params.compare_function = dshash_memcmp;
        params.hash_function = kmersearch_uint32_identity_hash; /* 恒等関数 */
    }
    else /* kmer_size <= 32 */
    {
        params.entry_size = sizeof(uint64) + sizeof(int); /* kmer64 + occurrence_count */
        params.compare_function = dshash_memcmp;
        params.hash_function = dshash_memhash; /* uint64の場合はPostgreSQL標準ハッシュ関数を使用 */
    }
    params.tranche_id = LWTRANCHE_KMERSEARCH_HASH;
    
    /* dshashテーブル作成 */
    hash = dshash_create(area, &params, NULL);
    
    *dsm_seg = seg;
    *dsa = area;
    return hash;
}

/* kmer2_as_uint恒等ハッシュ関数群（既にハッシュ化済みなのでそのまま返す） */
static inline uint32
kmersearch_uint16_identity_hash(const void *key, Size keysize, void *arg)
{
    return (uint32)(*(const uint16 *)key);
}

static inline uint32
kmersearch_uint32_identity_hash(const void *key, Size keysize, void *arg)
{
    return *(const uint32 *)key;
}

/* 
 * 注意：kmer_size > 16の場合は、dshash_memhashを使用するため、
 * kmersearch_kmer64_identity_hash関数は不要
 */
```

##### 1.2 並列ワーカー関数の実装

**並列ワーカーメイン関数**
```c
void kmersearch_parallel_worker_main(dsm_segment *seg, shm_toc *toc)
{
    KmerAnalysisSharedState *shared_state = NULL;
    dsa_area   *dsa = NULL;
    dshash_table *hash = NULL;
    ItemPointerData current_start_tid, current_end_tid;
    bool        has_work = true;
    
    PG_TRY();
    {
        /* 共有状態取得 */
        shared_state = shm_toc_lookup(toc, KMERSEARCH_KEY_SHARED_STATE, false);
        dsa = shm_toc_lookup(toc, KMERSEARCH_KEY_DSA, false);
        
        /* dshashテーブルにアタッチ */
        hash = dshash_attach(dsa, shm_toc_lookup(toc, KMERSEARCH_KEY_HASH_HANDLE, false));
        
        /* 動的ワーク取得ループ */
        while (has_work)
        {
            /* TID範囲の取得を試行（短時間の排他制御） */
            LWLockAcquire(&shared_state->mutex, LW_EXCLUSIVE);
            
            if (shared_state->tid_range_available)
            {
                /* TID範囲を取得 */
                current_start_tid = shared_state->next_start_tid;
                current_end_tid = shared_state->next_end_tid;
                
                /* DSM上のTID範囲を空に設定 */
                shared_state->tid_range_available = false;
                LWLockRelease(&shared_state->mutex);
                
                /* k-mer抽出・集計処理実行（排他制御外で並列実行） */
                kmersearch_extract_kmers_in_range(shared_state->table_oid,
                                                 shared_state->column_attnum,
                                                 shared_state->column_type_oid,
                                                 &current_start_tid,
                                                 &current_end_tid,
                                                 shared_state->kmer_size,
                                                 hash);
            }
            else if (shared_state->all_processed)
            {
                /* 全処理完了、ワーカー終了 */
                LWLockRelease(&shared_state->mutex);
                has_work = false;
            }
            else
            {
                /* TID範囲未設定だが処理継続中、1秒待機 */
                LWLockRelease(&shared_state->mutex);
                pg_usleep(1000000L); /* 1秒 */
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

**TID範囲ベースのk-mer抽出**
```c
static void
kmersearch_extract_kmers_in_range(Oid table_oid, AttrNumber column_attnum, Oid column_type_oid,
                                 ItemPointer start_tid, ItemPointer end_tid,
                                 int kmer_size, dshash_table *hash)
{
    Relation    rel = NULL;
    TableScanDesc scan = NULL;
    TupleTableSlot *slot = NULL;
    TupleDesc   tupdesc;
    bool        isnull;
    Datum       datum;
    
    PG_TRY();
    {
        /* テーブルオープン */
        rel = table_open(table_oid, AccessShareLock);
        tupdesc = RelationGetDescr(rel);
        
        /* スキャンスロット作成 */
        slot = table_slot_create(rel, NULL);
        
        /* TID範囲スキャン開始 */
        scan = table_beginscan_tid(rel, GetActiveSnapshot());
        
        /* TID範囲の妥当性チェック */
        if (!ItemPointerIsValid(start_tid) || !ItemPointerIsValid(end_tid))
            elog(ERROR, "Invalid TID range for k-mer extraction");
        
        while (table_scan_getnextslot_tid(scan, ForwardScanDirection, slot))
        {
            /* TID範囲チェック */
            if (ItemPointerCompare(&slot->tts_tid, start_tid) < 0 ||
                ItemPointerCompare(&slot->tts_tid, end_tid) > 0)
                continue;
            
            /* カラム値取得 */
            datum = slot_getattr(slot, column_attnum, &isnull);
            if (isnull)
                continue;
                
            /* k-mer抽出とdshash更新 */
            kmersearch_update_kmer_counts_in_dshash(datum, kmer_size, hash, column_type_oid);
        }
    }
    PG_CATCH();
    {
        /* エラー時のリソースクリーンアップ */
        if (scan)
            table_endscan(scan);
        if (slot)
            ExecDropSingleTupleTableSlot(slot);
        if (rel)
            table_close(rel, AccessShareLock);
        
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    /* 正常終了時のクリーンアップ */
    table_endscan(scan);
    ExecDropSingleTupleTableSlot(slot);
    table_close(rel, AccessShareLock);
}

**dshash更新処理（エラーハンドリング付き）**
```c
static void
kmersearch_update_kmer_counts_in_dshash(Datum sequence_datum, int kmer_size, dshash_table *hash, Oid column_type_oid)
{
    void       *kmer_array = NULL; /* k-merサイズに応じた型の配列 */
    int         kmer_count;
    bool        found;
    
    PG_TRY();
    {
        /* データ型に応じたk-mer抽出（k-merサイズに最適化された型で取得） */
        if (column_type_oid == DNA2_TYPE_OID)
        {
            /* DNA2型からの直接抽出 */
            kmersearch_extract_dna2_kmer2_as_uint_direct(sequence_datum, kmer_size, &kmer_array, &kmer_count);
        }
        else if (column_type_oid == DNA4_TYPE_OID)
        {
            /* DNA4型からの拡張付き直接抽出 */
            kmersearch_extract_dna4_kmer2_as_uint_with_expansion_direct(sequence_datum, kmer_size, &kmer_array, &kmer_count);
        }
        else
        {
            elog(ERROR, "Unsupported column type for k-mer extraction");
        }
        
        if (kmer_array == NULL || kmer_count <= 0)
            return; /* 有効なk-merなし */
        
        /* 各k-merをdshashに直接登録/更新（k-merサイズに応じた型で処理） */
        for (int i = 0; i < kmer_count; i++)
        {
            /* k-merサイズに応じた直接保存（型変換不要） */
            if (kmer_size <= 8)
            {
                typedef struct { uint16 kmer; int count; } KmerEntry16;
                uint16 *kmer16_array = (uint16 *)kmer_array;
                uint16 kmer16 = kmer16_array[i]; /* 既に適切な型で取得済み */
                KmerEntry16 *entry16 = (KmerEntry16 *)dshash_find_or_insert(hash, &kmer16, &found);
                if (!found)
                {
                    entry16->kmer = kmer16;
                    entry16->count = 1;
                }
                else
                {
                    if (entry16->count < INT_MAX)
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
}
```

#### フェーズ2: dshash結果の効率的な統合

##### 2.1 dshash結果分析

**高頻度k-mer特定とテーブル挿入**
```c
static void
kmersearch_finalize_analysis_results(dshash_table *hash, 
                                    Oid table_oid, const char *column_name, int kmer_size,
                                    float max_appearance_rate, int max_appearance_nrow)
{
    dshash_seq_status status;
    void *entry; /* k-merサイズに応じて異なるエントリ構造体 */
    StringInfoData query;
    int threshold_rows;
    uint64 kmer_as_uint64;
    int occurrence_count;
    
    /* 閾値計算 */
    threshold_rows = calculate_threshold_rows(table_oid, max_appearance_rate, max_appearance_nrow);
    
    /* バルクインサート用クエリ準備 */
    initStringInfo(&query);
    appendStringInfo(&query,
        "INSERT INTO kmersearch_highfreq_kmer "
        "(table_oid, column_name, kmer2_as_uint, detection_reason) VALUES ");
    
    /* dshash全エントリスキャン */
    dshash_seq_init(&status, hash, false);
    bool first_entry = true;
    
    while ((entry = dshash_seq_next(&status)) != NULL)
    {
        /* k-merサイズに応じてエントリから値を抽出 */
        if (kmer_size <= 8)
        {
            typedef struct { uint16 kmer; int count; } KmerEntry16;
            KmerEntry16 *entry16 = (KmerEntry16 *)entry;
            kmer_as_uint64 = (uint64)entry16->kmer;
            occurrence_count = entry16->count;
        }
        else if (kmer_size <= 16)
        {
            typedef struct { uint32 kmer; int count; } KmerEntry32;
            KmerEntry32 *entry32 = (KmerEntry32 *)entry;
            kmer_as_uint64 = (uint64)entry32->kmer;
            occurrence_count = entry32->count;
        }
        else /* kmer_size <= 32 */
        {
            typedef struct { uint64 kmer; int count; } KmerEntry64;
            KmerEntry64 *entry64 = (KmerEntry64 *)entry;
            kmer_as_uint64 = entry64->kmer;
            occurrence_count = entry64->count;
        }
        
        if (occurrence_count > threshold_rows)
        {
            if (!first_entry)
                appendStringInfoString(&query, ", ");
            
            appendStringInfo(&query, "(%u, %s, " UINT64_FORMAT ", 'parallel_frequency_analysis')",
                           table_oid, quote_literal_cstr(column_name), kmer_as_uint64);
            first_entry = false;
        }
    }
    
    dshash_seq_term(&status);
    
    /* 高頻度k-merが存在する場合のみ実行 */
    if (!first_entry)
    {
        SPI_connect();
        SPI_exec(query.data, 0);
        SPI_finish();
    }
    
    pfree(query.data);
}
```

#### フェーズ3: 統合メイン関数の実装

##### 3.1 新しいメイン関数

**改訂されたメイン関数**
```c
KmerAnalysisResult
kmersearch_perform_highfreq_analysis(text *table_name_or_oid, text *column_name_or_attnum)
{
    KmerAnalysisResult result = {0};
    ParallelContext *pcxt = NULL;
    dsm_segment *dsm_seg = NULL;
    dsa_area   *dsa = NULL;
    dshash_table *hash = NULL;
    KmerAnalysisSharedState *shared_state = NULL;
    Oid         table_oid;
    AttrNumber  column_attnum;
    char       *resolved_column_name = NULL;
    int         num_workers;
    bool        table_locked = false;
    
    PG_TRY();
    {
        /* 引数解析とOID/attnum解決 */
        table_oid = kmersearch_resolve_table_identifier(table_name_or_oid);
        column_attnum = kmersearch_resolve_column_identifier(table_oid, column_name_or_attnum, &resolved_column_name);
        
        /* テーブルロック取得 */
        LockRelationOid(table_oid, ExclusiveLock);
        table_locked = true;
        
        /* 事前検証とカラム型情報取得 */
        Oid column_type_oid = kmersearch_get_column_type_oid(table_oid, column_attnum);
        if (!kmersearch_validate_column_type(column_type_oid))
            elog(ERROR, "Column must be DNA2 or DNA4 type");
        
        /* 既存解析データのチェック */
        if (kmersearch_analysis_exists(table_oid, resolved_column_name))
            elog(ERROR, "High-frequency k-mer analysis already exists for table %u column %s", 
                 table_oid, resolved_column_name);
        
        /* 並列ワーカー数決定 */
        num_workers = kmersearch_determine_optimal_workers(table_oid);
        
        /* DSM/dshash環境構築 */
        hash = kmersearch_create_analysis_dshash(&dsm_seg, &dsa, kmersearch_kmer_size);
        shared_state = kmersearch_initialize_shared_state(dsm_seg, table_oid, column_attnum, column_type_oid);
        
        /* 並列コンテキスト作成 */
        pcxt = kmersearch_setup_parallel_context(num_workers, dsm_seg, hash);
        
        /* 並列ワーカー起動 */
        LaunchParallelWorkers(pcxt);
        result.parallel_workers_used = pcxt->nworkers_launched;
        
        /* 動的ワーク分散実行 */
        kmersearch_distribute_work_dynamically(shared_state, table_oid);
        
        /* ワーカー完了待機（エラー検知含む） */
        WaitForParallelWorkersToFinish(pcxt);
        
        /* 並列ワーカーのエラー状態確認 */
        if (shared_state->worker_error_occurred)
            elog(ERROR, "Parallel worker encountered error during analysis");
        
        /* dshash結果の分析・保存 */
        kmersearch_finalize_analysis_results(hash, table_oid, resolved_column_name, kmersearch_kmer_size,
                                            kmersearch_max_appearance_rate,
                                            kmersearch_max_appearance_nrow);
        
        /* メタデータ保存 */
        kmersearch_save_analysis_metadata(table_oid, resolved_column_name, &result);
    }
    PG_CATCH();
    {
        /* エラー時のリソースクリーンアップ */
        if (hash)
            dshash_destroy(hash);
        if (dsa)
            dsa_detach(dsa);
        if (dsm_seg)
            dsm_detach(dsm_seg);
        if (pcxt)
            DestroyParallelContext(pcxt);
        if (table_locked)
            UnlockRelationOid(table_oid, ExclusiveLock);
        if (resolved_column_name)
            pfree(resolved_column_name);
        
        PG_RE_THROW();
    }
    PG_END_TRY();
    
    /* 正常終了時のクリーンアップ */
    dshash_destroy(hash);
    dsa_detach(dsa);
    dsm_detach(dsm_seg);
    DestroyParallelContext(pcxt);
    UnlockRelationOid(table_oid, ExclusiveLock);
    pfree(resolved_column_name);
    
    return result;
}
```

**動的ワーク分散関数**
```c
static void
kmersearch_distribute_work_dynamically(KmerAnalysisSharedState *shared_state, Oid table_oid)
{
    Relation    rel;
    TableScanDesc scan;
    ItemPointerData current_tid, batch_start_tid;
    int         batch_count = 0;
    bool        has_more_rows = true;
    
    rel = table_open(table_oid, AccessShareLock);
    scan = table_beginscan(rel, GetActiveSnapshot(), 0, NULL);
    
    /* 最初のTID取得 */
    if (table_scan_getnextslot(scan, ForwardScanDirection, slot))
    {
        batch_start_tid = slot->tts_tid;
        batch_count = 1;
    }
    else
    {
        /* 空テーブルの場合 */
        LWLockAcquire(&shared_state->mutex, LW_EXCLUSIVE);
        shared_state->all_processed = true;
        LWLockRelease(&shared_state->mutex);
        table_endscan(scan);
        table_close(rel, AccessShareLock);
        return;
    }
    
    while (has_more_rows)
    {
        /* バッチサイズまで行を読み進める */
        while (batch_count < shared_state->batch_size && 
               table_scan_getnextslot(scan, ForwardScanDirection, slot))
        {
            current_tid = slot->tts_tid;
            batch_count++;
        }
        
        /* バッチ範囲をワーカーに分散 */
        LWLockAcquire(&shared_state->mutex, LW_EXCLUSIVE);
        
        /* 前のTID範囲が取得されるまで待機 */
        while (shared_state->tid_range_available)
        {
            LWLockRelease(&shared_state->mutex);
            pg_usleep(10000L); /* 10ms */
            LWLockAcquire(&shared_state->mutex, LW_EXCLUSIVE);
        }
        
        /* 新しいTID範囲を設定 */
        shared_state->next_start_tid = batch_start_tid;
        shared_state->next_end_tid = current_tid;
        shared_state->tid_range_available = true;
        
        LWLockRelease(&shared_state->mutex);
        
        /* 次のバッチ準備 */
        if (table_scan_getnextslot(scan, ForwardScanDirection, slot))
        {
            batch_start_tid = slot->tts_tid;
            batch_count = 1;
        }
        else
        {
            has_more_rows = false;
        }
    }
    
    /* 全行処理完了をマーク */
    LWLockAcquire(&shared_state->mutex, LW_EXCLUSIVE);
    shared_state->all_processed = true;
    LWLockRelease(&shared_state->mutex);
    
    table_endscan(scan);
    table_close(rel, AccessShareLock);
}
```

**引数解析ユーティリティ関数**
```c
/* テーブル名またはOIDからOIDを解決 */
static Oid
kmersearch_resolve_table_identifier(text *table_name_or_oid)
{
    char *str = text_to_cstring(table_name_or_oid);
    char *endptr;
    unsigned long oid_val;
    
    /* 数値（OID）として解析を試行 */
    oid_val = strtoul(str, &endptr, 10);
    if (*endptr == '\0' && oid_val != 0)
    {
        /* 数値として解析成功、OIDとして扱う */
        if (!OidIsValid(oid_val))
            elog(ERROR, "Invalid table OID: %lu", oid_val);
        return (Oid)oid_val;
    }
    else
    {
        /* テーブル名として解析 */
        RangeVar *rv = makeRangeVarFromNameList(textToQualifiedNameList(table_name_or_oid));
        return RangeVarGetRelid(rv, NoLock, false);
    }
}

/* カラム名またはattnumからattnumを解決 */
static AttrNumber
kmersearch_resolve_column_identifier(Oid table_oid, text *column_name_or_attnum, char **resolved_name)
{
    char *str = text_to_cstring(column_name_or_attnum);
    char *endptr;
    long attnum_val;
    
    /* 数値（attnum）として解析を試行 */
    attnum_val = strtol(str, &endptr, 10);
    if (*endptr == '\0' && attnum_val > 0)
    {
        /* 数値として解析成功、attnumとして扱う */
        AttrNumber attnum = (AttrNumber)attnum_val;
        Relation rel = table_open(table_oid, AccessShareLock);
        TupleDesc tupdesc = RelationGetDescr(rel);
        
        if (attnum <= 0 || attnum > tupdesc->natts)
        {
            table_close(rel, AccessShareLock);
            elog(ERROR, "Invalid column attnum: %ld", attnum_val);
        }
        
        /* カラム名を解決 */
        *resolved_name = pstrdup(NameStr(TupleDescAttr(tupdesc, attnum - 1)->attname));
        table_close(rel, AccessShareLock);
        
        return attnum;
    }
    else
    {
        /* カラム名として解析 */
        Relation rel = table_open(table_oid, AccessShareLock);
        AttrNumber attnum = get_attnum(table_oid, str);
        
        if (attnum == InvalidAttrNumber)
        {
            table_close(rel, AccessShareLock);
            elog(ERROR, "Column \"%s\" does not exist", str);
        }
        
        *resolved_name = pstrdup(str);
        table_close(rel, AccessShareLock);
        
        return attnum;
    }
}
```

### 実装順序とポイント

#### ステップ1: GUC変数追加 ❌ 未完
1. `kmersearch_highfreq_analysis_batch_size`（デフォルト10000）をGUC変数として追加 ❌
2. 関連するdefine定数とvalidation関数の実装 ❌

#### ステップ2: 引数解析機能実装 ❌ 未完
1. `kmersearch_resolve_table_identifier()`実装（テーブル名/OID解決） ❌
2. `kmersearch_resolve_column_identifier()`実装（カラム名/attnum解決） ❌
3. 関数シグネチャの変更（text型引数への対応） ❌

#### ステップ3: DSM/dshash基盤実装 ❌ 未完
1. `kmersearch_create_analysis_dshash()`の改良（DSMセグメント・dshashテーブル作成、k-merサイズ別対応） ❌
2. k-merサイズ別最適化エントリ構造（uint16/uint32/uint64対応） ❌
3. `dshash_memcmp`比較関数と恒等ハッシュ関数の実装 ❌
4. 動的分散対応`KmerAnalysisSharedState`構造体実装 ❌

#### ステップ4: 動的ワーク分散機能実装 ❌ 未完
1. `kmersearch_distribute_work_dynamically()`実装（親プロセスのワーク分散） ❌
2. `kmersearch_parallel_worker_main()`実装（動的ワーク取得ループ、エラーハンドリング付き） ❌
3. `kmersearch_extract_kmers_in_range()`実装（TID範囲ベースk-mer抽出、エラーハンドリング付き） ❌
4. `kmersearch_update_kmer_counts_in_dshash()`実装（dshash更新、エラーハンドリング付き） ❌

#### ステップ5: エラーハンドリング・検証機能実装 ❌ 未完
1. `kmersearch_analysis_exists()`実装（重複解析チェック） ❌
2. `kmersearch_validate_column_type()`実装（カラム型検証） ❌
3. TID範囲妥当性検証機能 ❌
4. 共有状態エラー伝播機能 ❌

#### ステップ6: 結果統合・保存機能実装 ❌ 未完
1. `kmersearch_finalize_analysis_results()`実装（dshash結果分析・テーブル挿入） ❌
2. `kmersearch_setup_parallel_context()`実装（並列コンテキスト作成） ❌
3. メイン関数統合（`kmersearch_perform_highfreq_analysis_parallel()`）（PG_TRY/PG_CATCH付き） ❌

## 期待される効果

### 性能向上
- **真の並列I/O**: 複数プロセスによる同時テーブルスキャン
- **動的ワーク分散**: 実際の並列ワーカー数に応じた効率的な負荷分散

## ❌ 実装ステータス - 並列実行未完成

**開始日**: 2025-07-25
**現在の状況**: 計画段階、実装開始前の状態にリセット

### ❌ 未実装機能
- ❌ DSM (Dynamic Shared Memory) / dshash 基盤実装
- ❌ **真の並列実行**: PostgreSQL標準並列実行機能を使用した並列ワーカー起動
- ❌ **並列ワーカー起動**: `LaunchParallelWorkers()`による実際の並列ワーカー起動
- ❌ **動的ワーク分散**: TID範囲ベースの動的ワーク分散機能
- ❌ **並列I/O**: 複数プロセスによる同時テーブルスキャンとk-mer抽出
- ❌ 並列ワーカーメイン関数 (`kmersearch_parallel_worker_main`) の実装
- ❌ TID範囲ベースk-mer抽出機能の実装
- ❌ エラーハンドリング・検証機能実装

### ❌ 未実装関数
- `kmersearch_perform_highfreq_analysis_parallel()` - **現在の実装は疑似並列**
  - 現在は順次処理による疑似並列のみ
  - 真のPostgreSQL並列実行機能が未実装
  - DSM/dshashによる並列ワーカー間のデータ共有が未実装
- `kmersearch_parallel_worker_main()` - **未実装**
- `kmersearch_setup_parallel_context()` - **未実装**
- `kmersearch_distribute_work_dynamically()` - **未実装**
- `kmersearch_extract_kmers_in_range()` - **未実装**
- `kmersearch_finalize_analysis_results()` - **改良が必要**

### ❌ 残存技術的問題
- **疑似並列実行**: 現在の実装は単一プロセスによる順次処理
- **並列基盤未活用**: DSM/dshash等の並列インフラが未実装
- **性能問題**: 大規模テーブルでの性能向上が未実現
- **カウント重複**: 現在の実装で高頻度k-mer数が正しくカウントされない問題

**重要**: 本計画の核心である**真の並列実行は未実装です**。現在の実装は疑似並列処理に留まります。

## 期待される効果

### 性能向上
- **真の並列I/O**: 複数プロセスによる同時テーブルスキャン
- **動的ワーク分散**: 実際の並列ワーカー数に応じた効率的な負荷分散

## リスク管理

### 互換性
- 既存の関数インターフェース維持
- GUC変数の互換性保持
- システムテーブル構造不変

### 安定性
- DSM/dshash管理の厳密化
- 並列実行時のエラーハンドリング強化
- リソースリーク防止（DSMセグメント・dshashテーブル）
- 並列ワーカー異常終了時の適切な処理

### パフォーマンス
- 小テーブルでの並列化オーバーヘッド回避
- dshashメモリ使用量監視
- TID範囲分散によるI/O競合最小化

この改訂された計画により、`kmersearch_perform_highfreq_analysis()`は`kmersearch_parallel_highfreq_kmer_cache_load()`関数と同様のDSM/dshash技術を活用した真の並列実行機能を獲得し、大幅な性能向上が期待されます。