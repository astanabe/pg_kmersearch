# SQLite3を使用した低メモリ消費型高頻度k-mer分析の実装計画

## 1. 背景と問題点

### 現在の実装の問題
- **メモリ消費問題**: 現在のdshashベースの実装では、全てのuintkey（k-mer+オカレンスカウント）をメモリ上で管理するため、巨大データベースではメモリ不足が発生
- **PostgreSQLの制限**: 並列ワーカーからテーブル・一時テーブルの作成が不可能なため、メモリベースのdshashを使用せざるを得ない状況

### 解決方針
- SQLite3のWALモードを使用した一時ファイルベースの実装に変更
- 並列ワーカーごとに独立したSQLite3データベースファイルを使用
- メモリ使用量を大幅に削減しつつ、並列処理を維持

## 2. 新アーキテクチャの概要

### 2.1 全体構成
```
親プロセス
├── SQLite3一時DB (pg_kmersearch_親PID_XXXXXX)
│   └── kmersearch_highfreq_kmer テーブル (全体集計用)
│
並列ワーカー1
├── SQLite3一時DB (pg_kmersearch_ワーカー1PID_XXXXXX)
│   └── kmersearch_highfreq_kmer テーブル (ワーカー1集計用)
│
並列ワーカー2
├── SQLite3一時DB (pg_kmersearch_ワーカー2PID_XXXXXX)
│   └── kmersearch_highfreq_kmer テーブル (ワーカー2集計用)
...
```

### 2.2 テーブル構造
```sql
CREATE TABLE kmersearch_highfreq_kmer (
    uintkey INTEGER PRIMARY KEY,  -- k-mer + オカレンスカウント（SQLite3が自動的に1-8バイトを選択）
    nline INTEGER                  -- 出現行数（SQLite3が値に応じて1-8バイトを自動選択）
);
```

## 3. 実装詳細

### 3.1 親プロセスでの準備フェーズ

#### 3.1.1 一時ディレクトリパスの決定と共有
```c
// 親プロセスでの処理
typedef struct KmerAnalysisSharedState {
    /* 既存のフィールド */
    BlockNumber next_block;
    BlockNumber total_blocks;
    /* ... */
    
    /* 追加: 一時ディレクトリパス */
    char temp_dir_path[MAXPGPATH];
} KmerAnalysisSharedState;

// 親プロセスで一時ディレクトリパスを決定
void prepare_temp_directory(KmerAnalysisSharedState *shared_state)
{
    Oid tablespace_oid;
    char temp_path[MAXPGPATH];
    
    /* temp_tablespacesから次のテーブル空間を取得 */
    tablespace_oid = GetNextTempTableSpace();
    
    /* 一時ディレクトリパスを構築 */
    TempTablespacePath(temp_path, tablespace_oid);
    
    /* 共有メモリにコピー */
    strlcpy(shared_state->temp_dir_path, temp_path, MAXPGPATH);
    
    elog(DEBUG1, "Selected temp directory: %s (tablespace OID: %u)",
         temp_path, tablespace_oid);
}
```

### 3.2 並列ワーカーの処理フロー

#### 3.2.1 初期化フェーズ
```c
// ワーカー起動時の処理
1. 親プロセスから一時ディレクトリパスを受け取り、mkstemp()で一時ファイルを作成
   - 共有メモリから一時ディレクトリパスを取得
   - パステンプレート: <shared_state->temp_dir_path>/pg_kmersearch_<PID>_XXXXXX
   - mkstemp()で一時ファイルを作成（自動削除されないため明示的管理が必要）

2. SQLite3データベースを作成・オープン
   - 各ワーカーが独立したファイルを使用するためWALモード不要
   - 同期モード設定: PRAGMA synchronous=NORMAL;（オプション）
   - キャッシュサイズ設定: PRAGMA cache_size=10000;（オプション）

3. kmersearch_highfreq_kmerテーブルを作成
   - CREATE TABLE文の実行

4. プリペアドステートメントを準備
   - INSERT OR UPDATE用のSQL文を事前にコンパイル
```

#### 3.2.2 詳細な実装例（ワーカー側）
```c
// ワーカー側での一時ファイル作成の実装例
void worker_create_temp_sqlite(KmerAnalysisSharedState *shared_state)
{
    char sqlite_db_path[MAXPGPATH];
    int fd;
    
    /* 共有メモリから一時ディレクトリパスを取得 */
    /* shared_state->temp_dir_pathには親プロセスが設定したパスが格納されている */
    
    /* 一時ファイルのテンプレートを作成 */
    snprintf(sqlite_db_path, MAXPGPATH, "%s/pg_kmersearch_%d_XXXXXX", 
             shared_state->temp_dir_path, getpid());
    
    /* mkstemp()で一時ファイルを作成 */
    fd = mkstemp(sqlite_db_path);
    if (fd < 0) {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not create temporary file \"%s\": %m",
                        sqlite_db_path)));
    }
    close(fd);  // SQLite3が独自にファイルを開くため一旦閉じる
    
    /* SQLite3データベースとして使用 */
    sqlite3 *db;
    int rc = sqlite3_open(sqlite_db_path, &db);
    
    /* ワーカー終了時にファイルパスを共有メモリに記録 */
    /* 親プロセスが後で読み込めるようにする */
}
```

#### 3.2.3 ワーカー関数の実装構造
```c
PGDLLEXPORT void
kmersearch_analysis_worker(dsm_segment *seg, shm_toc *toc)
{
    KmerAnalysisSharedState *shared_state;
    sqlite3 *db;
    sqlite3_stmt *stmt;
    char sqlite_db_path[MAXPGPATH];
    HTAB *batch_hash = NULL;
    int batch_count = 0;
    
    /* 共有メモリから情報取得 */
    shared_state = (KmerAnalysisSharedState *)shm_toc_lookup(toc, 
                                KMERSEARCH_KEY_SHARED_STATE, false);
    
    /* SQLite3データベース作成 */
    snprintf(sqlite_db_path, MAXPGPATH, "%s/pg_kmersearch_%d_XXXXXX",
             shared_state->temp_dir_path, getpid());
    int fd = mkstemp(sqlite_db_path);
    close(fd);
    sqlite3_open(sqlite_db_path, &db);
    
    /* SQLite3初期設定（各ワーカー独立ファイルなのでWAL不要） */
    sqlite3_exec(db, "CREATE TABLE kmersearch_highfreq_kmer "
                     "(uintkey INTEGER PRIMARY KEY, nline INTEGER)", 
                 NULL, NULL, NULL);
    
    /* プリペアドステートメント準備 */
    sqlite3_prepare_v2(db,
        "INSERT INTO kmersearch_highfreq_kmer (uintkey, nline) "
        "VALUES (?, ?) "
        "ON CONFLICT(uintkey) DO UPDATE SET nline = nline + excluded.nline",
        -1, &stmt, NULL);
    
    /* メインループ：ブロック処理 */
    while (has_work) {
        BlockNumber current_block = get_next_block(shared_state);
        if (current_block == InvalidBlockNumber) break;
        
        /* ブロック内の全行を処理 */
        process_block_with_batch(current_block, shared_state, 
                                &batch_hash, &batch_count, db, stmt);
    }
    
    /* 残りのバッチをフラッシュ */
    if (batch_hash && batch_count > 0) {
        flush_batch_to_sqlite(batch_hash, db, stmt);
    }
    
    /* クリーンアップとファイルパス登録 */
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    register_worker_temp_file(shared_state, sqlite_db_path, worker_id);
}
```

#### 3.2.4 ブロック処理とバッチ管理
```c
/* 一時的なk-mer頻度エントリ構造体（total_bitsに応じて使い分け） */
typedef struct TempKmerFreqEntry16 {
    uint16 uintkey;
    uint64 nline;  /* 出現行数（将来の拡張性のためuint64） */
} TempKmerFreqEntry16;

typedef struct TempKmerFreqEntry32 {
    uint32 uintkey;
    uint64 nline;  /* 出現行数（将来の拡張性のためuint64） */
} TempKmerFreqEntry32;

typedef struct TempKmerFreqEntry64 {
    uint64 uintkey;
    uint64 nline;  /* 出現行数（将来の拡張性のためuint64） */
} TempKmerFreqEntry64;

static void
process_block_with_batch(BlockNumber block, 
                        KmerAnalysisSharedState *shared_state,
                        HTAB **batch_hash, int *batch_count,
                        sqlite3 *db, sqlite3_stmt *stmt)
{
    Relation rel;
    Buffer buffer;
    Page page;
    OffsetNumber maxoff;
    int total_bits = shared_state->kmer_size * 2 + shared_state->occur_bitlen;
    
    /* バッチハッシュテーブルの初期化（必要に応じて） */
    if (*batch_hash == NULL) {
        HASHCTL hashctl;
        memset(&hashctl, 0, sizeof(hashctl));
        
        /* total_bitsに基づいてキーサイズとエントリサイズを決定 */
        if (total_bits <= 16) {
            hashctl.keysize = sizeof(uint16);
            hashctl.entrysize = sizeof(TempKmerFreqEntry16);
        } else if (total_bits <= 32) {
            hashctl.keysize = sizeof(uint32);
            hashctl.entrysize = sizeof(TempKmerFreqEntry32);
        } else {
            hashctl.keysize = sizeof(uint64);
            hashctl.entrysize = sizeof(TempKmerFreqEntry64);
        }
        
        *batch_hash = hash_create("KmerBatchHash",
                                 kmersearch_highfreq_analysis_batch_size,
                                 &hashctl, HASH_ELEM | HASH_BLOBS);
    }
    
    /* ブロックを読み込み */
    rel = table_open(shared_state->table_oid, AccessShareLock);
    buffer = ReadBuffer(rel, block);
    LockBuffer(buffer, BUFFER_LOCK_SHARE);
    page = BufferGetPage(buffer);
    maxoff = PageGetMaxOffsetNumber(page);
    
    /* ページ内の全タプルを処理 */
    for (OffsetNumber offnum = FirstOffsetNumber; 
         offnum <= maxoff; 
         offnum = OffsetNumberNext(offnum)) {
        
        ItemId itemid = PageGetItemId(page, offnum);
        if (!ItemIdIsNormal(itemid)) continue;
        
        HeapTupleData tuple;
        tuple.t_data = (HeapTupleHeader) PageGetItem(page, itemid);
        
        /* シーケンスデータ取得 */
        bool isnull;
        Datum datum = heap_getattr(&tuple, shared_state->column_attnum,
                                  RelationGetDescr(rel), &isnull);
        if (isnull) continue;
        
        /* k-mer抽出と集計 */
        void *kmer_array = NULL;
        int kmer_count;
        
        if (shared_state->column_type_oid == dna2_oid) {
            VarBit *seq = DatumGetVarBitP(datum);
            kmersearch_extract_uintkey_from_dna2(seq, &kmer_array, &kmer_count);
        } else {
            VarBit *seq = DatumGetVarBitP(datum);
            kmersearch_extract_uintkey_from_dna4(seq, &kmer_array, &kmer_count);
        }
        
        /* バッチハッシュテーブルに追加 */
        for (int i = 0; i < kmer_count; i++) {
            bool found;
            void *entry;
            
            if (total_bits <= 16) {
                uint16 uintkey = ((uint16 *)kmer_array)[i];
                TempKmerFreqEntry16 *freq_entry;
                
                entry = hash_search(*batch_hash, &uintkey, HASH_ENTER, &found);
                freq_entry = (TempKmerFreqEntry16 *)entry;
                if (!found) {
                    freq_entry->uintkey = uintkey;
                    freq_entry->nline = 0;
                }
                freq_entry->nline++;
                
            } else if (total_bits <= 32) {
                uint32 uintkey = ((uint32 *)kmer_array)[i];
                TempKmerFreqEntry32 *freq_entry;
                
                entry = hash_search(*batch_hash, &uintkey, HASH_ENTER, &found);
                freq_entry = (TempKmerFreqEntry32 *)entry;
                if (!found) {
                    freq_entry->uintkey = uintkey;
                    freq_entry->nline = 0;
                }
                freq_entry->nline++;
                
            } else {
                uint64 uintkey = ((uint64 *)kmer_array)[i];
                TempKmerFreqEntry64 *freq_entry;
                
                entry = hash_search(*batch_hash, &uintkey, HASH_ENTER, &found);
                freq_entry = (TempKmerFreqEntry64 *)entry;
                if (!found) {
                    freq_entry->uintkey = uintkey;
                    freq_entry->nline = 0;
                }
                freq_entry->nline++;
            }
        }
        
        if (kmer_array) pfree(kmer_array);
        (*batch_count)++;
        
        /* バッチサイズに達したらSQLite3に書き込み */
        if (*batch_count >= kmersearch_highfreq_analysis_batch_size) {
            flush_batch_to_sqlite(*batch_hash, db, stmt);
            hash_destroy(*batch_hash);
            *batch_hash = NULL;
            *batch_count = 0;
        }
    }
    
    UnlockReleaseBuffer(buffer);
    table_close(rel, AccessShareLock);
}
```

#### 3.2.5 バッチフラッシュ関数
```c
static void
flush_batch_to_sqlite(HTAB *batch_hash, sqlite3 *db, sqlite3_stmt *stmt,
                     int total_bits)
{
    HASH_SEQ_STATUS status;
    void *entry;
    
    /* トランザクション開始 */
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, NULL);
    
    /* ハッシュテーブルをスキャンしてSQLite3に書き込み */
    hash_seq_init(&status, batch_hash);
    
    while ((entry = hash_seq_search(&status)) != NULL) {
        /* total_bitsに応じて適切な型とbind関数を使用してSQLite3に書き込み */
        if (total_bits <= 16) {
            TempKmerFreqEntry16 *e = (TempKmerFreqEntry16 *)entry;
            /* uint16をint16として扱い、SQLite3に2バイトで保存 */
            sqlite3_bind_int(stmt, 1, (int16)e->uintkey);
            sqlite3_bind_int64(stmt, 2, (int64)e->nline);
        } else if (total_bits <= 32) {
            TempKmerFreqEntry32 *e = (TempKmerFreqEntry32 *)entry;
            /* uint32をint32として扱い、SQLite3に4バイトで保存 */
            sqlite3_bind_int(stmt, 1, (int32)e->uintkey);
            sqlite3_bind_int64(stmt, 2, (int64)e->nline);
        } else {
            TempKmerFreqEntry64 *e = (TempKmerFreqEntry64 *)entry;
            /* uint64をint64として扱い、SQLite3に8バイトで保存 */
            sqlite3_bind_int64(stmt, 1, (int64)e->uintkey);
            sqlite3_bind_int64(stmt, 2, (int64)e->nline);
        }
        
        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }
    
    /* トランザクションコミット */
    sqlite3_exec(db, "COMMIT", NULL, NULL, NULL);
}
```

#### 3.2.6 データ処理フェーズ（簡略版）
```c
// バッチ処理ループの要点
while (担当ブロックが残っている) {
    // 1. メモリ上のハッシュテーブルを初期化
    HTAB *batch_hash = hash_create(
        "KmerBatchHash",
        kmersearch.highfreq_analysis_batch_size,
        ...
    );
    
    // 2. batch_size分の行を処理
    for (int i = 0; i < batch_size && has_more_rows; i++) {
        Datum sequence_datum;
        void *kmer_array = NULL;
        int kmer_count;
        
        // シーケンスデータを取得
        sequence_datum = heap_getattr(...);
        
        // 既存の関数でk-mer（uintkey）を抽出
        if (column_type_oid == dna2_oid) {
            VarBit *seq = DatumGetVarBitP(sequence_datum);
            kmersearch_extract_uintkey_from_dna2(seq, &kmer_array, &kmer_count);
        } else if (column_type_oid == dna4_oid) {
            VarBit *seq = DatumGetVarBitP(sequence_datum);
            kmersearch_extract_uintkey_from_dna4(seq, &kmer_array, &kmer_count);
        }
        
        // メモリ上のハッシュテーブルで集計
        for (int j = 0; j < kmer_count; j++) {
            uint64 uintkey;
            if (total_bits <= 16) {
                uintkey = ((uint16 *)kmer_array)[j];
            } else if (total_bits <= 32) {
                uintkey = ((uint32 *)kmer_array)[j];
            } else {
                uintkey = ((uint64 *)kmer_array)[j];
            }
            
            // ハッシュテーブルで出現行数をカウント
            update_batch_hash(batch_hash, uintkey);
        }
        
        if (kmer_array) pfree(kmer_array);
    }
    
    // 3. SQLite3トランザクション開始
    sqlite3_exec(db, "BEGIN TRANSACTION", ...);
    
    // 4. バッチの内容をSQLite3に反映
    for (each entry in batch_hash) {
        sqlite3_bind_int(stmt, 1, entry->uintkey);
        sqlite3_bind_int(stmt, 2, entry->nline);
        sqlite3_step(stmt);  // INSERT OR UPDATE実行
        sqlite3_reset(stmt);
    }
    
    // 5. トランザクションコミット
    sqlite3_exec(db, "COMMIT", ...);
    
    // 6. メモリ上のハッシュテーブルを破棄
    hash_destroy(batch_hash);
}
```

#### 3.2.4 SQL文の詳細
```sql
-- uintkeyの出現行数を更新するSQL
INSERT INTO kmersearch_highfreq_kmer (uintkey, nline) 
VALUES (?, ?)
ON CONFLICT(uintkey) DO UPDATE 
SET nline = nline + excluded.nline;
```

### 3.3 親プロセスの処理フロー

#### 3.3.1 ワーカーのファイルパス収集
```c
// 共有メモリにワーカーのファイルパス格納領域を追加
typedef struct KmerAnalysisSharedState {
    /* 既存のフィールド */
    char temp_dir_path[MAXPGPATH];
    
    /* 追加: ワーカーの一時ファイルパス配列 */
    int num_workers;
    char worker_temp_files[MAX_PARALLEL_WORKERS][MAXPGPATH];
    SpinLock worker_file_lock;  // ファイルパス登録時の排他制御
} KmerAnalysisSharedState;

// ワーカー側: 一時ファイルパスを登録
void register_worker_temp_file(KmerAnalysisSharedState *shared_state, 
                               const char *file_path, int worker_id)
{
    SpinLockAcquire(&shared_state->worker_file_lock);
    strlcpy(shared_state->worker_temp_files[worker_id], 
            file_path, MAXPGPATH);
    SpinLockRelease(&shared_state->worker_file_lock);
}
```

#### 3.3.2 並列ワーカー完了後の集計
```c
// 1. 自身の一時SQLite3データベースを作成
char parent_db_path[MAXPGPATH];
snprintf(parent_db_path, MAXPGPATH, "%s/pg_kmersearch_%d_XXXXXX",
         shared_state->temp_dir_path, getpid());
int fd = mkstemp(parent_db_path);
close(fd);
sqlite3 *parent_db;
sqlite3_open(parent_db_path, &parent_db);

// 2. 各ワーカーの結果を集計
for (int i = 0; i < shared_state->num_workers; i++) {
    char *worker_db_path = shared_state->worker_temp_files[i];
    
    if (strlen(worker_db_path) == 0)
        continue;  // ワーカーがファイルを作成しなかった場合
    // ワーカーのDBをアタッチ
    sqlite3_exec(parent_db, 
        "ATTACH DATABASE ? AS worker_db", 
        worker_db_path, ...);
    
    // ワーカーのデータを親DBに集計
    sqlite3_exec(parent_db,
        "INSERT INTO main.kmersearch_highfreq_kmer (uintkey, nline) "
        "SELECT uintkey, nline FROM worker_db.kmersearch_highfreq_kmer "
        "ON CONFLICT(uintkey) DO UPDATE "
        "SET nline = main.kmersearch_highfreq_kmer.nline + excluded.nline",
        ...);
    
    // ワーカーDBをデタッチ
    sqlite3_exec(parent_db, "DETACH DATABASE worker_db", ...);
    
    // ワーカーの一時ファイルを削除
    unlink(worker_db_path);
}
```

#### 3.3.3 高頻度k-merの特定と保存
```c
// 1. 閾値の計算
uint64 threshold_rows = calculate_threshold(
    total_rows, 
    kmersearch_max_appearance_rate,
    kmersearch_max_appearance_nrow
);

// 2. 高頻度k-merの抽出
sqlite3_prepare_v2(parent_db,
    "SELECT uintkey FROM kmersearch_highfreq_kmer "
    "WHERE nline > ?",
    ...);
sqlite3_bind_int64(stmt, 1, (int64)threshold_rows);

// 3. PostgreSQLのkmersearch_highfreq_kmerテーブルに保存
SPI_connect();
while (sqlite3_step(stmt) == SQLITE_ROW) {
    uint64 uintkey;
    
    /* total_bitsに応じて適切な型で読み取り、uint64に変換 */
    if (total_bits <= 32) {
        /* sqlite3_column_int()はint型を返すので、uint64にキャスト */
        uintkey = (uint64)sqlite3_column_int(stmt, 0);
    } else {
        uintkey = (uint64)sqlite3_column_int64(stmt, 0);
    }
    
    // PostgreSQLテーブルに挿入
    SPI_exec(
        "INSERT INTO kmersearch_highfreq_kmer "
        "(table_oid, column_name, uintkey, detection_reason) "
        "VALUES (...)", 
        ...);
}
SPI_finish();

// 4. 親の一時ファイルを削除
sqlite3_close(parent_db);
unlink(parent_db_path);
```

## 4. メモリ管理戦略

### 4.1 メモリ使用量の制御
- **バッチサイズ**: `kmersearch.highfreq_analysis_batch_size`で制御（デフォルト: 10000行）
- **メモリ上限**: バッチ処理ごとにハッシュテーブルを作成・破棄することで、メモリ使用量を一定以下に維持
- **予想メモリ使用量（total_bitsに応じて最適化）**: 
  - total_bits ≤ 16: バッチサイズ × (2 + 8 + オーバーヘッド) ≈ 150KB〜1MB
  - total_bits ≤ 32: バッチサイズ × (4 + 8 + オーバーヘッド) ≈ 180KB〜1.2MB
  - total_bits ≤ 64: バッチサイズ × (8 + 8 + オーバーヘッド) ≈ 240KB〜1.6MB
  - 全体: ワーカー数 × ワーカーあたりのメモリ使用量

### 4.2 ストレージ使用量
- **一時ファイルサイズ**: データ量とユニークなuintkeyの数に依存
- **SQLite3のINTEGER型の特性**:
  - 値に応じて自動的に1〜8バイトを使用（可変長ストレージ）
  - 小さい値は少ないバイト数で効率的に保存
- **最悪ケース**: 全行が異なるk-merを持つ場合
  - ファイルサイズ ≈ ユニークuintkey数 × (実際の値サイズ + B-treeオーバーヘッド）
- **削除タイミング**: 
  - ワーカーファイル: 親プロセスでの集計完了後即座に削除
  - 親ファイル: 高頻度k-mer特定・保存完了後即座に削除

## 5. エラーハンドリング

### 5.1 SQLite3エラー処理
```c
// SQLite3操作のエラーチェック
int rc = sqlite3_exec(db, sql, ...);
if (rc != SQLITE_OK) {
    // エラーメッセージを取得
    const char *errmsg = sqlite3_errmsg(db);
    
    // PostgreSQLエラーとして報告
    ereport(ERROR,
        (errcode(ERRCODE_INTERNAL_ERROR),
         errmsg("SQLite3 error: %s", errmsg)));
}
```

### 5.2 クリーンアップ処理
```c
// PG_TRY/PG_CATCHブロックでの確実なクリーンアップ
PG_TRY();
{
    // 通常処理
}
PG_CATCH();
{
    // SQLite3データベースのクローズ
    if (db) sqlite3_close(db);
    
    // 一時ファイルの削除
    if (temp_file_path) unlink(temp_file_path);
    
    PG_RE_THROW();
}
PG_END_TRY();
```

## 6. パフォーマンス最適化

### 6.1 SQLite3の最適化設定
```sql
-- ジャーナルモード（各ワーカー独立ファイルのためデフォルトで可）
-- PRAGMA journal_mode = DELETE;  /* デフォルト */

-- 同期モード（パフォーマンス重視）
PRAGMA synchronous = NORMAL;

-- ページキャッシュサイズ
PRAGMA cache_size = 10000;

-- 一時ストレージ
PRAGMA temp_store = MEMORY;

-- メモリマップI/O
PRAGMA mmap_size = 268435456;  -- 256MB
```

### 6.2 バッチ処理の最適化
- トランザクションごとに複数のINSERT/UPDATEをまとめて実行
- プリペアドステートメントの再利用
- インデックスの自動作成（PRIMARY KEY）

## 7. 実装上の注意点

### 7.0 dshash削除時の注意事項
- **削除対象は高頻度分析専用のdshashのみ**:
  - `kmersearch_perform_highfreq_analysis_parallel()`内で使用されるdshash
  - `kmersearch_analysis_worker()`内で使用されるdshash
- **以下のdshash使用箇所は削除しない**:
  - GINインデックス構築（`gin_extract_value_dna2/4`など）
  - 並列キャッシュ管理（`kmersearch_parallel_highfreq_kmer_cache_*`）
  - その他のキャッシュ管理機能
- **削除前に必ずgrep等で使用箇所を確認**

### 7.1 型変換の注意事項
- **SQLite3のINTEGER型への書き込み（容量最適化）**:
  - uint16: `(int16)`にキャストして`sqlite3_bind_int()`で書き込み（2バイト保存）
  - uint32: `(int32)`にキャストして`sqlite3_bind_int()`で書き込み（4バイト保存）
  - uint64: `(int64)`にキャストして`sqlite3_bind_int64()`で書き込み（8バイト保存）
  - SQLite3のINTEGER型は可変長なので、適切なbind関数を使用することで容量を最適化
- **SQLite3からの読み込み**:
  - 読み込み時は`sqlite3_column_int()`または`sqlite3_column_int64()`を使用
  - 必要に応じて元のuint型に再解釈キャスト
- **注意**: uint型をint型にキャストすると最上位ビットが1の場合は負数になるが、ビットパターンは保持される

### 7.2 既存関数の活用
- **k-mer抽出には既存関数を使用**:
  - `kmersearch_extract_uintkey_from_dna2()`: DNA2シーケンスからuintkey抽出
  - `kmersearch_extract_uintkey_from_dna4()`: DNA4シーケンスからuintkey抽出
  - これらの関数は既にオカレンスカウント込みのuintkeyを返すため、そのまま使用可能
- **新規関数の作成は最小限に**:
  - SQLite3インターフェース関数のみ新規作成
  - 既存のユーティリティ関数を最大限活用

### 7.3 ファイルパスの管理
- PostgreSQLの一時ファイル配置ルールに従う
  - `temp_tablespaces`が設定されている場合: 指定されたテーブル空間の`pgsql_tmp`ディレクトリ
  - 未設定の場合: デフォルトテーブル空間（通常は`$PGDATA/base/pgsql_tmp`）
- プロセスIDを含むユニークなファイル名
- セキュアな一時ファイル作成（`mkstemp()`使用）
- PostgreSQL内部APIを使用して適切なパスを取得
  - `GetNextTempTableSpace()`でテーブル空間OIDを取得
  - `TempTablespacePath()`で一時ディレクトリパスを構築
  - 重要: PostgreSQLの自動削除機能を使わず、独自に一時ファイルを管理
  - 並列ワーカー終了時に自動削除されないよう、明示的なファイル管理が必要

### 7.4 並行性の考慮
- 各ワーカーが独立したファイルを使用するため、ロック競合なし
- ファイル間の同期や排他制御が不要

### 7.5 PostgreSQL統合
- SPIインターフェースとSQLite3 APIの適切な使い分け
- エラーレポートの統一（PostgreSQLのereport使用）
- メモリコンテキストの適切な管理

## 8. 一時ファイルクリーンアップ機能

### 8.1 kmersearch_delete_tempfiles()関数の設計

#### 8.1.1 関数シグネチャ
```sql
-- SQL関数定義
CREATE OR REPLACE FUNCTION kmersearch_delete_tempfiles()
RETURNS TABLE(
    deleted_count integer,
    deleted_size bigint,
    error_count integer
) AS 'MODULE_PATHNAME', 'kmersearch_delete_tempfiles'
LANGUAGE C STRICT;

-- 使用例
SELECT * FROM kmersearch_delete_tempfiles();
```

#### 8.1.2 実装詳細
```c
Datum
kmersearch_delete_tempfiles(PG_FUNCTION_ARGS)
{
    int deleted_count = 0;
    int64 deleted_size = 0;
    int error_count = 0;
    List *tablespace_oids = NIL;
    ListCell *lc;
    
    /* 1. temp_tablespacesから全てのテーブル空間OIDを取得 */
    int num_tablespaces;
    Oid *tablespace_array = palloc(sizeof(Oid) * 256);  // 最大256個
    
    num_tablespaces = GetTempTablespaces(tablespace_array, 256);
    
    /* 2. デフォルトテーブル空間も追加（InvalidOid） */
    tablespace_oids = lappend_oid(tablespace_oids, InvalidOid);
    
    /* 3. temp_tablespacesのテーブル空間を追加 */
    for (int i = 0; i < num_tablespaces; i++) {
        tablespace_oids = lappend_oid(tablespace_oids, tablespace_array[i]);
    }
    
    /* 4. 各テーブル空間のpgsql_tmpディレクトリをスキャン */
    foreach(lc, tablespace_oids) {
        Oid tablespace_oid = lfirst_oid(lc);
        char temp_path[MAXPGPATH];
        DIR *dir;
        struct dirent *de;
        
        /* テーブル空間の一時ディレクトリパスを取得 */
        TempTablespacePath(temp_path, tablespace_oid);
        
        /* ディレクトリをオープン */
        dir = opendir(temp_path);
        if (dir == NULL) {
            /* ディレクトリが存在しない場合は警告のみ */
            ereport(WARNING,
                    (errmsg("could not open temp directory \"%s\": %m", temp_path)));
            continue;
        }
        
        /* ディレクトリ内のファイルをスキャン */
        while ((de = readdir(dir)) != NULL) {
            char full_path[MAXPGPATH];
            struct stat st;
            
            /* pg_kmersearch_で始まるファイル名をチェック */
            if (strncmp(de->d_name, "pg_kmersearch_", 14) != 0)
                continue;
            
            /* フルパスを構築 */
            snprintf(full_path, MAXPGPATH, "%s/%s", temp_path, de->d_name);
            
            /* ファイル情報を取得 */
            if (stat(full_path, &st) < 0) {
                ereport(WARNING,
                        (errmsg("could not stat file \"%s\": %m", full_path)));
                error_count++;
                continue;
            }
            
            /* 通常ファイルかチェック */
            if (!S_ISREG(st.st_mode))
                continue;
            
            /* ファイルの作成時間をチェック（オプション） */
            /* 現在時刻から一定時間以上古いファイルのみ削除する場合 */
            time_t current_time = time(NULL);
            if (current_time - st.st_mtime < 60) {  // 60秒以内のファイルはスキップ
                ereport(INFO,
                        (errmsg("skipping recently created file \"%s\"", full_path)));
                continue;
            }
            
            /* ファイルサイズを記録 */
            int64 file_size = st.st_size;
            
            /* ファイルを削除 */
            if (unlink(full_path) < 0) {
                ereport(WARNING,
                        (errmsg("could not delete file \"%s\": %m", full_path)));
                error_count++;
            } else {
                deleted_count++;
                deleted_size += file_size;
                ereport(INFO,
                        (errmsg("deleted temporary file \"%s\" (size: %ld bytes)",
                                full_path, file_size)));
            }
        }
        
        closedir(dir);
    }
    
    /* 5. 結果を返す */
    TupleDesc tupdesc;
    Datum values[3];
    bool nulls[3] = {false};
    HeapTuple tuple;
    
    if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE) {
        ereport(ERROR,
                (errmsg("function returning record called in context "
                        "that cannot accept a record")));
    }
    
    values[0] = Int32GetDatum(deleted_count);
    values[1] = Int64GetDatum(deleted_size);
    values[2] = Int32GetDatum(error_count);
    
    tuple = heap_form_tuple(tupdesc, values, nulls);
    
    PG_RETURN_DATUM(HeapTupleGetDatum(tuple));
}
```

#### 8.1.3 安全性の考慮事項
1. **ファイル名パターンのチェック**
   - `pg_kmersearch_`で始まるファイルのみを対象とする
   - 他のPostgreSQLプロセスのファイルを誤って削除しない

2. **タイムスタンプチェック**
   - 作成から一定時間（例: 60秒）経過したファイルのみ削除
   - 現在実行中の分析のファイルを削除しない

3. **エラーハンドリング**
   - ファイル削除に失敗してもプロセス全体は継続
   - 各エラーをWARNINGレベルで記録

4. **権限チェック**
   - ファイル削除権限がない場合はエラーとして記録
   - スーパーユーザー権限を要求するかは要検討

### 8.2 自動クリーンアップの仕組み

#### 8.2.1 分析開始時の自動クリーンアップ
```c
// kmersearch_perform_highfreq_analysis()の冒頭に追加
static void
kmersearch_cleanup_old_tempfiles(void)
{
    /* 前回の異常終了で残った可能性のある古い一時ファイルを削除 */
    DirectFunctionCall0(kmersearch_delete_tempfiles);
}
```

#### 8.2.2 PostgreSQL起動時のクリーンアップ
```c
// _PG_init()に追加
void
_PG_init(void)
{
    /* ... 既存の初期化処理 ... */
    
    /* 起動時に古い一時ファイルをクリーンアップ */
    if (!IsBootstrapProcessingMode() && !IsUnderPostmaster) {
        kmersearch_cleanup_old_tempfiles();
    }
}
```

### 8.3 手動クリーンアップコマンド
```sql
-- 手動で一時ファイルをクリーンアップ
SELECT * FROM kmersearch_delete_tempfiles();

-- 結果例:
-- deleted_count | deleted_size | error_count
-- --------------+--------------+-------------
--             5 |     10485760 |           0
```

## 9. 移行計画

### 9.1 既存コードの変更箇所
1. `kmersearch_freq.c`
   - `kmersearch_perform_highfreq_analysis_parallel()`の大幅改修
   - `kmersearch_analysis_worker()`内にSQLite3処理を直接実装
     - `kmersearch_extract_kmers_from_block()`の処理をワーカー内に統合
     - メモリバッチ処理とSQLite3更新を直接実装
   - **高頻度分析専用のdshash関連コードのみ削除**:
     - `analysis_dsm_segment`
     - `analysis_dsa`
     - `analysis_highfreq_hash`
     - `kmersearch_update_kmer_counts_in_dshash()` - 高頻度分析専用
     - `kmersearch_extract_kmers_from_block()` - ワーカー内に統合
   - **既存のk-mer抽出関数を再利用**:
     - `kmersearch_extract_uintkey_from_dna2()`
     - `kmersearch_extract_uintkey_from_dna4()`

2. `kmersearch.h`
   - SQLite3関連の構造体・関数宣言追加
   - **重要: 以下のdshash関連は削除しない（他機能で使用）**:
     - GINインデックス構築用のdshash関数
     - `kmersearch_parallel_highfreq_kmer_cache_load()`関連
     - `kmersearch_parallel_highfreq_kmer_cache_free()`関連
     - その他のキャッシュ管理用dshash構造体
   - 新規関数の追加は最小限に抑える

3. `Makefile`
   - SQLite3ライブラリのリンク追加
   - 必要に応じてコンパイルフラグ追加

### 9.2 テスト計画
1. 単体テスト
   - SQLite3操作の正常性確認
   - バッチ処理の境界値テスト
   - エラーハンドリングのテスト

2. 統合テスト
   - 小規模データでの動作確認
   - 大規模データでのメモリ使用量測定
   - 並列度を変えた性能測定

3. 回帰テスト
   - 既存のmake installcheckが引き続き成功することを確認

## 10. 期待される効果

### 10.1 メモリ使用量の削減
- **現在**: データサイズに比例してメモリ使用量が増加
- **改善後**: バッチサイズで制限された一定量のメモリのみ使用
- **削減率**: 大規模データで90%以上のメモリ削減が期待

### 10.2 スケーラビリティの向上
- 事実上無制限のデータサイズに対応可能
- ストレージ容量が許す限り処理可能

### 10.3 トレードオフ
- ディスクI/Oの増加による処理時間の増加
- 一時的なストレージ使用量の増加
- SQLite3依存性の追加

## 11. 今後の拡張可能性

### 11.1 圧縮オプション
- SQLite3の圧縮拡張機能の利用検討
- 一時ファイルサイズの削減

### 11.2 分散処理対応
- ネットワーク経由での結果集約
- 複数ノードでの並列処理

### 11.3 インクリメンタル分析
- 差分データのみの分析
- 既存の分析結果との統合