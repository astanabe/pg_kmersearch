# DNA2/DNA4型GINインデックス並列作成機能実装計画 [PARTIALLY COMPLETED - Sequential execution only]

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

既存のパーティションテーブル（DNA2型またはDNA4型カラムを含む）の各パーティションに対してCREATE INDEXをPostgreSQLの標準並列実行機能で並列ワーカーを起動して実行する`kmersearch_parallel_create_index()`関数の実装計画。

## 関数仕様

### 関数シグネチャ
```sql
kmersearch_parallel_create_index(
    table_name text,
    column_name text,
    tablespace_name text DEFAULT NULL
) RETURNS TABLE(
    partition_name text,
    index_name text,
    rows_processed bigint,
    execution_time_ms bigint,
    worker_pid int,
    success boolean,
    error_message text
)
```

### パラメータ
- `table_name`: 対象パーティションテーブル名
- `column_name`: DNA2/DNA4型の対象カラム名
- `tablespace_name`: インデックス作成先のテーブルスペース名（省略可能、デフォルト: 対象テーブルと同一のテーブルスペース）

### 戻り値
各パーティションのインデックス作成結果を返すテーブル形式

## 技術的要件  

### 前提条件
- **対象テーブル**: PostgreSQLのパーティションテーブル（`relkind = 'p'`）であること
- **対象カラム**: DNA2型またはDNA4型のカラムであること

### 高頻度k-mer除外設定（オプション）
高頻度k-mer除外を行う場合は以下の設定が必要：
- **高頻度k-mer除外**: `kmersearch.preclude_highfreq_kmer = true`
- **並列キャッシュ強制使用**: `kmersearch.force_use_parallel_highfreq_kmer_cache = true`
- **事前キャッシュロード**: `kmersearch_parallel_highfreq_kmer_cache_load(table_name, column_name)`実行済み

高頻度k-mer除外を行わない場合：
- **高頻度k-mer除外**: `kmersearch.preclude_highfreq_kmer = false`

### PostgreSQL並列実行設定への依存
- `max_parallel_workers_per_gather`: ギャザーノードあたりの最大並列ワーカー数
- `max_parallel_workers`: システム全体の最大並列ワーカー数
- `max_worker_processes`: 最大ワーカープロセス数
- `max_parallel_maintenance_workers`: メンテナンス作業の最大並列ワーカー数

## アーキテクチャ設計

### 1. パーティション検出とメタデータ取得

#### パーティションテーブル判定
対象テーブルがパーティションテーブルかどうかを判定：
```sql
SELECT relkind = 'p' FROM pg_class WHERE oid = table_oid;
```

#### パーティション一覧取得
既存パーティションの一覧を取得：
```sql
SELECT inhrelid::regclass AS partition_name,
       inhrelid AS partition_oid
FROM pg_inherits 
WHERE inhparent = table_oid;
```

#### カラム存在確認と型チェック
各パーティションでの対象カラムの存在とDNA2/DNA4型であることを確認：
```sql
SELECT attname, atttypid
FROM pg_attribute a
JOIN pg_class c ON a.attrelid = c.oid
WHERE c.oid = partition_oid
  AND attname = column_name
  AND NOT attisdropped;
```

カラム型の検証：
```c
// 既存の関数を使用してDNA2/DNA4型を判定
if (attr->atttypid == get_dna4_type_oid()) {
    is_dna4_type = true;
} else if (attr->atttypid == get_dna2_type_oid()) {
    is_dna4_type = false;
} else {
    ereport(ERROR, (errmsg("Column '%s' must be DNA2 or DNA4 type", column_name)));
}
```

### 2. 並列実行制御

#### ワーカー数決定ロジック
検出されたパーティション数とシステム制限に基づいてワーカー数を決定：
```c
int partition_count = detected_partitions_count;
int max_workers = Min(max_parallel_workers, max_parallel_maintenance_workers);
int effective_workers = Min(partition_count, max_workers);
```

#### バックグラウンドワーカー起動
PostgreSQLの標準バックグラウンドワーカーAPIを使用：
- `RegisterBackgroundWorker()`: ワーカープロセス登録
- `WaitForBackgroundWorkerStartup()`: 起動完了待機
- `TerminateBackgroundWorker()`: ワーカー終了処理

#### パーティション-ワーカー割り当て
各ワーカーに対して処理すべきパーティションを割り当て：
```c
typedef struct ParallelIndexWorkerArgs {
    Oid partition_oid;
    char partition_name[NAMEDATALEN];
    char column_name[NAMEDATALEN];
    int worker_id;
} ParallelIndexWorkerArgs;
```

### 3. 高頻度k-mer除外機能連携

#### GUC設定とキャッシュ状態の検証
関数開始時に高頻度k-mer除外の設定組み合わせを検証：

```c
// 高頻度k-mer除外設定の検証
if (kmersearch_preclude_highfreq_kmer) {
    // preclude_highfreq_kmer=trueの場合、force_use_parallel_highfreq_kmer_cache=trueが必須
    if (!kmersearch_force_use_parallel_highfreq_kmer_cache) {
        ereport(ERROR, (errmsg("force_use_parallel_highfreq_kmer_cache must be true when preclude_highfreq_kmer is true")));
    }
    
    // 並列キャッシュのロード状態検証
    if (!kmersearch_is_parallel_highfreq_cache_loaded() ||
        !kmersearch_parallel_highfreq_kmer_cache_is_valid(table_oid, column_name, k_value)) {
        ereport(ERROR, (errmsg("Parallel high-frequency k-mer cache not loaded or invalid"),
                       errhint("Run kmersearch_parallel_highfreq_kmer_cache_load('%s', '%s') first", 
                              table_name, column_name)));
    }
} else {
    // preclude_highfreq_kmer=falseの場合、高頻度k-mer除外は行わない
    ereport(NOTICE, (errmsg("High-frequency k-mer exclusion disabled (preclude_highfreq_kmer=false)")));
}
```

#### 各ワーカーでのキャッシュ利用
- **preclude_highfreq_kmer=false**: 高頻度k-mer除外なし、通常のGINインデックス作成
- **preclude_highfreq_kmer=true**: 既存のGIN構築処理で高頻度k-mer除外が自動実行
  - `force_use_parallel_highfreq_kmer_cache=true`により並列キャッシュを使用
  - 各パーティションで同一のキャッシュ内容を参照してk-mer除外を実施

### 4. エラーハンドリングと監視

#### エラー処理戦略
- **パーティション単位の失敗処理**: 一部パーティションが失敗しても他は継続
- **非パーティションテーブルエラー**: 対象テーブルがパーティションテーブルでない場合は即座にエラー出力
- **非DNA2/DNA4型カラムエラー**: 対象カラムがDNA2/DNA4型でない場合はエラー出力
- **GUC設定不整合エラー**: 
  - `preclude_highfreq_kmer=true` + `force_use_parallel_highfreq_kmer_cache=false` → エラー停止
  - `preclude_highfreq_kmer=true` + キャッシュ未ロード → エラー停止
- **ロールバック**: 失敗時は作成済みインデックスをクリーンアップ

#### 進捗監視
- 各ワーカーの進捗を共有メモリで管理
- パーティション単位での実行状況追跡
- 戻り値テーブルでの詳細結果提供

## 実装詳細

### 1. 主要関数構成

#### メイン関数
```c
Datum kmersearch_parallel_create_index(PG_FUNCTION_ARGS)
{
    // 1. パラメータ検証とテーブル存在確認
    // 2. パーティションテーブル判定（非パーティションはエラー）
    // 3. 高頻度k-mer除外設定検証（GUC組み合わせとキャッシュ状態）
    // 4. パーティション一覧取得とカラム存在確認
    // 5. カラム型検証（DNA2/DNA4型以外はエラー）
    // 6. 全条件満たす場合のみワーカー数決定と並列実行開始
    // 7. 結果収集と戻り値構築
}
```

#### パーティション検出関数
```c
static List* get_table_partitions(Oid table_oid, const char* column_name)
{
    // パーティション一覧取得
    // 各パーティションでのカラム存在確認
    // カラム型検証（DNA2/DNA4型判定）
    // PartitionInfo構造体のリストを返却
}
```

#### パーティションテーブル判定関数
```c
static bool is_partitioned_table(Oid table_oid)
{
    // pg_classでrelkind='p'を確認
}
```

#### ワーカー実行関数
```c
static void parallel_index_worker_main(Datum main_arg)
{
    // 1. パーティションでのGINインデックス作成（CONCURRENTLYなし）
    //    - CREATE INDEXが自動的にACCESS EXCLUSIVE LOCKを取得
    //    - tablespace_nameが指定されている場合はTABLESPACE句を追加
    // 2. 高頻度k-mer除外処理（自動実行、既存キャッシュ利用）
    // 3. 実行結果の共有メモリへの書き込み
}
```

### 2. 共有メモリ構造

#### パーティション情報構造体
```c
typedef struct PartitionInfo {
    Oid partition_oid;
    char partition_name[NAMEDATALEN];
    bool has_target_column;
    Oid column_type_oid;
    bool is_dna4_type;  // DNA2 or DNA4 type identification
} PartitionInfo;
```

#### ワーカー状態管理
```c
typedef struct ParallelIndexState {
    int total_partitions;
    int completed_partitions;
    bool error_occurred;
    ParallelIndexPartitionResult results[FLEXIBLE_ARRAY_MEMBER];
} ParallelIndexState;

typedef struct ParallelIndexPartitionResult {
    char partition_name[NAMEDATALEN];
    char index_name[NAMEDATALEN];
    int64 rows_processed;
    int64 execution_time_ms;
    int worker_pid;
    bool success;
    char error_message[1024];
} ParallelIndexPartitionResult;
```

### 3. インデックス命名規則

各パーティションのインデックス名：
```
{partition_name}_{column_name}_gin_idx
```
例: `sales_2023_q1_dna_sequence_gin_idx`

### 4. テーブルスペース処理

インデックス作成時のテーブルスペース決定ロジック：
```c
/* テーブルスペース名の決定 */
if (tablespace_name != NULL && strlen(tablespace_name) > 0) {
    /* 明示的に指定されたテーブルスペースを使用 */
    snprintf(create_index_sql, sizeof(create_index_sql),
        "CREATE INDEX %s ON %s USING gin (%s) TABLESPACE %s",
        index_name, partition_name, column_name, tablespace_name);
} else {
    /* テーブルスペース句なし（パーティションと同じ場所に作成） */
    snprintf(create_index_sql, sizeof(create_index_sql),
        "CREATE INDEX %s ON %s USING gin (%s)",
        index_name, partition_name, column_name);
}
```

### 5. 実行シーケンス

1. **前処理**
   - パラメータ検証
   - テーブル存在確認
   - パーティションテーブル判定（非パーティションはエラー）
   - 高頻度k-mer除外設定検証（GUC組み合わせとキャッシュ状態）

2. **パーティション検出と検証**
   - 既存パーティション一覧取得
   - 各パーティションでのカラム存在確認
   - カラム型検証（DNA2/DNA4型以外はエラー）
   - 処理対象パーティションリスト作成

3. **並列実行**
   - 共有メモリセグメント初期化
   - ワーカー数決定（パーティション数とシステム制限考慮）
   - バックグラウンドワーカー起動
   - 各ワーカーで以下を実行：
     - GINインデックス作成（CONCURRENTLYなし、自動的にACCESS EXCLUSIVE LOCK取得）
     - 高頻度k-mer除外処理（既存キャッシュ利用）

4. **結果収集**
   - ワーカー完了待機
   - パーティション単位の実行結果収集
   - 成功/失敗状況の集約

5. **後処理**
   - エラー時の部分ロールバック
   - 統計情報更新
   - 詳細結果のテーブル形式返却

## 制限事項と注意点

### 技術的制限
- **パーティションテーブル必須**: 対象テーブルはPostgreSQLのパーティションテーブルである必要がある
- **同時実行制限**: PostgreSQLの並列ワーカー制限に依存
- **GUC設定制限**: 高頻度k-mer除外機能の有効・無効による設定要件の違い

### 運用上の注意点
- **事前パーティション設定**: テーブルが適切にパーティション化されている必要がある
- **カラム型制限**: 対象カラムはDNA2型またはDNA4型である必要がある
- **GUC設定の整合性**: 
  - `preclude_highfreq_kmer=true`の場合、`force_use_parallel_highfreq_kmer_cache=true`が必須
  - 高頻度k-mer除外を行う場合は事前にキャッシュロードが必要
- **カラム一貫性**: 全パーティションで対象カラムが同一型で存在する必要がある
- **排他的アクセス**: CREATE INDEX（CONCURRENTLYなし）により各パーティションは作成中完全に排他アクセスとなる

### パフォーマンス考慮事項
- **パーティション数と並列度**: パーティション数がワーカー数より多い場合、順次実行される
- **テーブルサイズ**: 小さなパーティションでは並列化のオーバーヘッドが性能劣化要因
- **排他ロックの影響**: CREATE INDEXのACCESS EXCLUSIVE LOCKによりパーティション単位で完全ブロック
- **I/Oボトルネック**: 大量の並列I/Oによるディスク性能への影響
- **CONCURRENTLYなし**: 高速なインデックス作成だが、テーブルアクセス完全ブロック

## テスト戦略

### 単体テスト
- パーティション検出ロジックの正確性
- 非パーティションテーブルでのエラーハンドリング
- 非DNA2/DNA4型カラムでのエラーハンドリング
- GUC設定不整合時のエラーハンドリング：
  - `preclude_highfreq_kmer=true` + `force_use_parallel_highfreq_kmer_cache=false`
  - `preclude_highfreq_kmer=true` + キャッシュ未ロード
- カラム存在確認の適切性
- CREATE INDEXによる自動ロック取得の確認
- メモリリーク検出

### 統合テスト  
- 高頻度k-mer除外機能の有効・無効での動作確認
- `preclude_highfreq_kmer=false`での通常インデックス作成
- `preclude_highfreq_kmer=true`での高頻度k-mer除外インデックス作成
- 様々なパーティション数での動作確認
- DNA2型とDNA4型の混在パーティションでの動作確認
- 並列度を変えた性能測定
- 部分失敗時の動作確認
- CREATE INDEX実行中の同時アクセス制限テスト

### 回帰テスト
- 既存のGINインデックス機能への影響確認
- 並列キャッシュシステムとの互換性確認
- 従来の単一テーブルインデックス作成への影響

## 実装優先度

### Phase 1: 基本機能実装
- パーティション検出ロジック
- パーティションテーブル判定機能
- 基本的な並列実行制御
- エラーハンドリング

### Phase 2: 高度な機能
- 動的ワーカー数調整
- 詳細な進捗監視
- 性能最適化

### Phase 3: 運用機能
- ログ出力強化
- 統計情報収集
- 監視機能拡張

## 期待効果

### 性能向上
- 大規模パーティションテーブルでのインデックス作成時間短縮
- 複数パーティションの並列処理による高速化
- CPU使用率とI/O並列化の最適化

### 運用性向上
- パーティション化されたDNA2/DNA4テーブルのメンテナンス効率化
- バッチ処理での統一的なインデックス作成
- システムリソースの有効活用

### スケーラビリティ
- パーティション数の増加に対するリニアな性能向上
- 大規模データセットでの実用的なインデックス構築時間実現

この実装により、パーティション化されたDNA2/DNA4型の大規模テーブルに対して効率的なGINインデックス作成が可能となり、k-mer検索システムでのパーティショニング戦略を活用した高性能データ処理基盤の構築が実現される。

## 実装状況 - 2025年7月現在

### 実装済み機能（約60-70%完了）

1. **kmersearch_partition_table関数**
   - 第3引数として`tablespace_name`を追加（省略可能）
   - 省略時は元のテーブルと同じテーブルスペースを使用
   - 指定時は指定されたテーブルスペースにパーティションを作成

2. **kmersearch_parallel_create_index関数**
   - 第3引数として`tablespace_name`を追加（省略可能）
   - 省略時はパーティションと同じテーブルスペースを使用
   - 指定時は指定されたテーブルスペースにインデックスを作成

### 未実装機能（主要機能）

1. **並列実行機能** - 核心機能が未実装
   - 現在は順次実行のみ（`create_partition_indexes()`内のコメント参照）
   - バックグラウンドワーカーAPIの使用は未実装
   - パーティション-ワーカー割り当てロジックなし
   - 並列実行制御機構が存在しない

2. **詳細な実行結果追跡**
   - 個別パーティションの実行時間、処理行数などは未追跡
   - worker_pidは常に0（プレースホルダー）
   - 実際の並列実行統計情報なし

3. **高度な機能**
   - 動的ワーカー数調整
   - 共有メモリを使用した進捗監視
   - 詳細なパフォーマンス統計収集

### 実装の相違点

1. **データ移行方式（kmersearch_partition_table）**
   - 計画：maintenance_work_memベースの動的バッチサイズ計算
   - 実装：シンプルなINSERT SELECT + TRUNCATE（一括処理）

2. **並列インデックス作成（kmersearch_parallel_create_index）**
   - 計画：PostgreSQLの標準並列実行機能で複数ワーカー起動
   - 実装：単一プロセスでの順次実行

3. **エラーハンドリング**
   - 計画：パーティション単位の詳細な失敗処理
   - 実装：基本的なエラーハンドリングのみ

### 使用例

```sql
-- デフォルトテーブルスペース（元のテーブルと同じ）
SELECT kmersearch_partition_table('my_table', 4);

-- 明示的なテーブルスペース指定
SELECT kmersearch_partition_table('my_table', 4, 'fast_ssd');

-- インデックス作成（デフォルト）
SELECT kmersearch_parallel_create_index('my_table', 'sequence');

-- インデックス作成（テーブルスペース指定）
SELECT kmersearch_parallel_create_index('my_table', 'sequence', 'fast_nvme');
```

### 注意事項

PostgreSQLの制約により、パーティションテーブル作成時に`pg_default`テーブルスペースを明示的に指定することはできません。`pg_default`を使用したい場合は、テーブルスペースパラメータを省略するか、NULLを指定してください。

```sql
-- 正しい使い方（pg_defaultを使用）
SELECT kmersearch_partition_table('my_table', 4);        -- 省略
SELECT kmersearch_partition_table('my_table', 4, NULL);  -- NULL指定

-- エラーになる使い方
SELECT kmersearch_partition_table('my_table', 4, 'pg_default');  -- ERROR
```

## PostgreSQL並列処理実装ガイドライン

`kmersearch_perform_highfreq_analysis()`の実装経験から得られた、PostgreSQL拡張で並列処理を実装する際に従うべきガイドライン：

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

### 5. DSM/DSA/dshashアーキテクチャ

#### 正しいパターン
1. DSM/DSA/dshash資源をEnterParallelMode()の**前**に作成
2. shm_tocを通じてハンドルを渡す
3. ワーカーはアタッチと読み取りのみ、作成や構造変更は禁止

```c
/* 親プロセス - EnterParallelMode()前 */
create_dshash_resources(&ctx);

/* InitializeParallelDSM()後 */
handles->dsm_handle = ctx.dsm_handle;
shm_toc_insert(toc, KEY, handles);

/* ワーカー - 読み取り専用 */
handles = shm_toc_lookup(toc, KEY, false);
dsm_attach(handles->dsm_handle);
```

### 6. リソースクリーンアップ

#### 問題
クリーンアップ関数がSQLやその他の制限された操作を実行する可能性がある。

#### 解決策
すべてのクリーンアップをIsParallelWorker()チェックでガード：
```c
static void cleanup_resources(void)
{
    /* 並列ワーカーでは決してクリーンアップを実行しない */
    if (IsParallelWorker())
        return;
    
    /* クリーンアップコード */
}
```

### 7. エラーハンドリングのベストプラクティス
- ワーカーでは単純なNULLチェックとelog(ERROR)を使用
- SQLをトリガーする可能性のある複雑なエラーハンドリングを避ける
- 共有状態を通じてPostgreSQLにワーカーエラーを処理させる

### 8. 並列モード終了タイミング

#### 重大な問題
PostgreSQLは、WaitForParallelWorkersToFinish()の後でもExitParallelMode()が呼ばれるまで並列モードのままである。

#### 問題
ワーカーが終了してもExitParallelMode()前にSQL操作を試みると「cannot execute INSERT during a parallel operation」エラーが発生。

#### 解決策
SQL操作を実行する前に必ず並列モードを終了：
```c
/* BAD - 並列モードのままSQL操作 */
WaitForParallelWorkersToFinish(pcxt);
SPI_connect();
SPI_execute("INSERT ...", false, 0);  /* エラー！ */
DestroyParallelContext(pcxt);
ExitParallelMode();

/* GOOD - 先に並列モードを終了 */
WaitForParallelWorkersToFinish(pcxt);
DestroyParallelContext(pcxt);
ExitParallelMode();
/* これでSQL操作が安全に実行可能 */
SPI_connect();
SPI_execute("INSERT ...", false, 0);
```

これはよくある落とし穴で、メインプロセスでIsParallelWorker()がfalseを返しても、ExitParallelMode()が呼ばれるまでPostgreSQLは依然として並列モードと見なすためである。

### 9. 一般的な並列処理パターン

#### 正しい並列実行フロー
```c
/* 1. 共有リソースを作成（並列モード前） */
create_shared_resources();

/* 2. 並列モードに入る */
EnterParallelMode();

/* 3. 並列コンテキストを作成・設定 */
pcxt = CreateParallelContext("worker_function", nworkers);
InitializeParallelDSM(pcxt);

/* 4. ワーカーを起動 */
LaunchParallelWorkers(pcxt);

/* 5. 完了を待つ */
WaitForParallelWorkersToFinish(pcxt);

/* 6. SQL操作の前に並列モードを終了 */
DestroyParallelContext(pcxt);
ExitParallelMode();

/* 7. これでSQL操作が安全に実行可能 */
SPI_connect();
/* SQL操作 */
SPI_finish();

/* 8. 共有リソースをクリーンアップ */
cleanup_shared_resources();
```

### 10. 並列コードのテスト

以下の条件で必ずテスト：
- 小さいデータ（非TOAST）
- 大きいデータ（TOAST圧縮）
- 複数ワーカー
- ワーカー障害
- PostgreSQLログでDSMリークやSQL実行試行の警告をチェック
- 「cannot execute ... during a parallel operation」エラーが発生しないことを確認

### 11. よくある間違いと対策

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

これらのガイドラインに従うことで、PostgreSQL拡張における並列処理実装の一般的な落とし穴を回避し、安定した高性能な実装を実現できる。

## 補助関数: kmersearch_partition_table [COMPLETED]

### 概要

既存の非パーティションテーブルをハッシュパーティションテーブルに変換する関数。DNA2/DNA4型カラムを基準にハッシュ分割を行い、データを保持したままテーブル構造を変換する。

### 関数仕様

#### 関数シグネチャ
```sql
kmersearch_partition_table(
    table_name text,             -- 対象テーブル名（OIDも可）
    partition_count int,         -- 分割数
    tablespace_name text DEFAULT NULL  -- テーブルスペース名（省略可能）
) RETURNS void
```

#### パラメータ
- `table_name`: 変換対象の非パーティションテーブル名（regclass型でOID指定も可能）
- `partition_count`: ハッシュパーティションの分割数（1以上の整数）
- `tablespace_name`: パーティションテーブル作成先のテーブルスペース名（省略可能、デフォルト: 元のテーブルと同一のテーブルスペース）

### 実行条件

以下のすべての条件を満たす必要がある：
1. 対象テーブルがパーティションテーブルではない（`relkind != 'p'`）
2. 対象テーブルにDNA2型またはDNA4型カラムがちょうど1つ存在する
3. 対象テーブルが存在し、アクセス権限がある

条件を満たさない場合はエラーで停止する。

### 処理フロー

#### 1. 事前検証
```sql
-- パーティションテーブルかどうかを確認
SELECT relkind FROM pg_class WHERE oid = table_oid;
-- エラー: "Table '%s' is already a partitioned table"

-- DNA2/DNA4型カラムの検出
SELECT attname, atttypid 
FROM pg_attribute 
WHERE attrelid = table_oid 
  AND NOT attisdropped 
  AND atttypid IN (dna2_type_oid, dna4_type_oid);
-- エラー: "Table must have exactly one DNA2 or DNA4 column"
```

#### 2. 一時パーティションテーブル作成
```sql
-- タイムスタンプを含む一時テーブル名生成
temp_table_name = sprintf("%s_temp%ld", table_name, time(NULL));

-- テーブルスペースの決定
if (tablespace_name != NULL && strlen(tablespace_name) > 0) {
    /* 指定されたテーブルスペースを使用 */
    target_tablespace = tablespace_name;
} else {
    /* 元のテーブルと同じテーブルスペースを使用 */
    target_tablespace = get_table_tablespace(table_oid);
}

-- 親テーブル作成
CREATE TABLE {temp_table_name} (
    LIKE {table_name} INCLUDING ALL
) PARTITION BY HASH ({dna_column_name})
TABLESPACE {target_tablespace};
```

#### 3. パーティション（子テーブル）作成
```sql
-- 各パーティションを作成（N = 0 から partition_count-1）
CREATE TABLE {table_name}_{N} 
PARTITION OF {temp_table_name} 
FOR VALUES WITH (modulus {partition_count}, remainder {N})
TABLESPACE {target_tablespace};
```

#### 4. データ移行（バッチ処理）
```c
/* maintenance_work_memに基づいてバッチサイズを動的に決定 */
static int calculate_partition_batch_size(Oid table_oid, AttrNumber attnum)
{
    int64 maintenance_work_mem_bytes = maintenance_work_mem * 1024L;
    int64 avg_row_size;
    int batch_size;
    
    /* 平均行サイズを推定（統計情報から） */
    avg_row_size = estimate_avg_row_size(table_oid);
    if (avg_row_size <= 0)
        avg_row_size = 1024;  /* デフォルト1KB */
    
    /* maintenance_work_memの1/4をバッチに使用 */
    batch_size = (maintenance_work_mem_bytes / 4) / avg_row_size;
    
    /* 最小1,000行、最大100,000行に制限 */
    if (batch_size < 1000)
        batch_size = 1000;
    else if (batch_size > 100000)
        batch_size = 100000;
        
    elog(NOTICE, "Using batch size %d based on maintenance_work_mem=%dMB",
         batch_size, maintenance_work_mem);
         
    return batch_size;
}

/* バッチ単位でのデータ移行 */
int batch_size = calculate_partition_batch_size(table_oid, attnum);

while (true) {
    /* トランザクション開始 */
    SPI_execute("BEGIN", false, 0);
    
    /* バッチデータ取得 */
    snprintf(query, sizeof(query),
        "SELECT * FROM %s LIMIT %d FOR UPDATE",
        table_name, batch_size);
    SPI_execute(query, false, 0);
    
    if (SPI_processed == 0) {
        /* 全データ移行完了 */
        SPI_execute("COMMIT", false, 0);
        break;
    }
    
    /* パーティションテーブルへの挿入 */
    for (i = 0; i < SPI_processed; i++) {
        /* 各行をINSERT */
        insert_into_partition_table(temp_table_name, SPI_tuptable->vals[i]);
    }
    
    /* 元テーブルから削除 */
    snprintf(query, sizeof(query),
        "DELETE FROM %s WHERE ctid IN (SELECT ctid FROM %s LIMIT %d)",
        table_name, table_name, BATCH_SIZE);
    SPI_execute(query, false, 0);
    
    /* トランザクションコミット */
    SPI_execute("COMMIT", false, 0);
    
    /* メモリ解放とVACUUM考慮 */
    if (total_rows_processed % (BATCH_SIZE * 10) == 0) {
        /* 定期的なVACUUMの検討（オプション） */
        check_vacuum_threshold(table_oid);
    }
}
```

#### 5. テーブル置換
```sql
-- 元テーブルが空であることを確認
SELECT COUNT(*) FROM {table_name};
-- COUNT = 0でなければエラー

-- 元テーブルを削除
DROP TABLE {table_name};

-- パーティションテーブルを元の名前に変更
ALTER TABLE {temp_table_name} RENAME TO {table_name};

-- 各パーティションの名前も調整（オプション）
-- {temp_table_name}_N → {table_name}_N
```

### エラーハンドリング

#### トランザクション管理
- バッチ処理の各イテレーションは独立したトランザクションで実行
- エラー発生時は現在のバッチのみロールバック
- 全体的な整合性は一時テーブルの存在により保証

#### リカバリ処理
```c
/* エラー時のクリーンアップ */
PG_CATCH();
{
    /* 一時テーブルとパーティションの削除 */
    cleanup_temp_partitions(temp_table_name, partition_count);
    
    /* エラーを再スロー */
    PG_RE_THROW();
}
PG_END_TRY();
```

### パフォーマンス考慮事項

#### ストレージ使用量
- 移行中の最大ストレージ使用量: 元のテーブルサイズ + 現在のバッチサイズ
- バッチ処理により、完全な2倍のストレージを必要としない

#### ロック戦略
- 元テーブル: バッチ単位でのFOR UPDATE行ロック
- パーティションテーブル: 通常のINSERTロック
- 長時間の排他ロックを回避

#### 最適化オプション
```c
/* PostgreSQL標準のmaintenance_work_memを使用 */
/* バッチサイズは動的に計算され、追加のGUC変数は不要 */

/* オプション: VACUUM制御用GUC変数 */
bool kmersearch_partition_vacuum_enabled = true;  /* 自動VACUUM有効/無効 */
int kmersearch_partition_vacuum_threshold = 100000;  /* VACUUM実行閾値 */
```

### 実装例

```c
Datum
kmersearch_partition_table(PG_FUNCTION_ARGS)
{
    text *table_name_text = PG_GETARG_TEXT_PP(0);
    int32 partition_count = PG_GETARG_INT32(1);
    text *tablespace_name_text = PG_ARGISNULL(2) ? NULL : PG_GETARG_TEXT_PP(2);
    char *table_name;
    char *tablespace_name = NULL;
    Oid table_oid;
    Oid dna_column_type;
    char *dna_column_name;
    char temp_table_name[NAMEDATALEN];
    
    /* パラメータ検証 */
    if (partition_count < 1)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
                 errmsg("partition_count must be at least 1")));
    
    /* テーブル名/OID解決 */
    table_name = text_to_cstring(table_name_text);
    table_oid = get_table_oid(table_name, false);
    
    /* テーブルスペース名取得 */
    if (tablespace_name_text != NULL)
        tablespace_name = text_to_cstring(tablespace_name_text);
    
    /* 事前検証 */
    validate_table_for_partitioning(table_oid, &dna_column_name, &dna_column_type);
    
    /* 一時テーブル名生成 */
    snprintf(temp_table_name, sizeof(temp_table_name), 
             "%s_temp%ld", table_name, (long)time(NULL));
    
    /* パーティションテーブル作成 */
    create_partition_table(temp_table_name, table_name, dna_column_name, 
                          partition_count, tablespace_name);
    
    /* データ移行（バッチサイズは自動計算） */
    migrate_data_in_batches(table_name, temp_table_name, table_oid);
    
    /* テーブル置換 */
    replace_table_with_partition(table_name, temp_table_name);
    
    PG_RETURN_VOID();
}
```

### 使用例

```sql
-- 単純な使用例（デフォルトテーブルスペース）
SELECT kmersearch_partition_table('sequences', 4);

-- テーブルスペースを指定した例
SELECT kmersearch_partition_table('sequences', 4, 'fast_ssd');

-- OIDを使用した例
SELECT kmersearch_partition_table(16384::regclass, 8, 'large_storage');

-- 大規模テーブルの場合（maintenance_work_memを調整）
SET maintenance_work_mem = '256MB';  -- より大きなバッチサイズが自動計算される
SELECT kmersearch_partition_table('large_sequences', 16, 'fast_nvme');
```

### 制限事項と注意点

1. **実行時間**: 大規模テーブルの場合、データ量に比例した時間が必要
2. **ディスク容量**: 移行中は追加のディスク容量が必要（最大でバッチサイズ分）
3. **外部キー制約**: 外部キー制約は自動的に再作成されない
4. **トリガー**: LIKE INCLUDING ALLにより基本的なトリガーは複製される
5. **同時アクセス**: 移行中のテーブルへの同時アクセスは制限される
6. **インデックス**: 既存のインデックスは削除され、パーティション化後に再作成が必要

### 後続作業

パーティション化完了後、以下の作業が推奨される：

1. **インデックス再作成**: 
   ```sql
   SELECT kmersearch_parallel_create_index('sequences', 'dna_sequence');
   ```

2. **統計情報更新**:
   ```sql
   ANALYZE sequences;
   ```

3. **高頻度k-mer分析**（必要に応じて）:
   ```sql
   SELECT kmersearch_perform_highfreq_analysis('sequences', 'dna_sequence');
   ```

この補助関数により、既存の非パーティションテーブルを効率的にパーティション化し、`kmersearch_parallel_create_index()`関数と組み合わせることで、大規模データセットでの高速なk-mer検索インフラストラクチャを構築できる。