# SQLite3からファイルベースハッシュテーブルへの移行実装計画

## 概要

`kmersearch_perform_highfreq_analysis()`における並列ワーカーの一時ファイル保存方式を、SQLite3からシンプルなファイルベースハッシュテーブルに変更する。

### 現状の課題

現在の実装では、並列ワーカーがSQLite3を使用して「uintkey - appearance_nrow」をファイルに保存しているが：
- 複数プロセスからの並列アクセスは発生しない（各ワーカーが独自のファイルを使用）
- 保存するのは単純なkey-value型情報のみ
- SQLite3のオーバーヘッドが大きい（SQL解析、トランザクション管理、インデックス管理等）

### 設計方針

1. **ハッシュ生成不要**: uintkey(uint16/uint32/uint64)をそのままキーとして使用可能
2. **削除処理不要**: 分析完了後にファイル全体を削除するのみ
3. **キャッシュ不要**: OSのディスクキャッシュに任せる（mmap不使用）
4. **3パターン実装**: uint16/uint32/uint64それぞれ専用の実装で速度を最優先

---

## データ構造設計

### 1. uint16 key - uint64 value（配列実装）

uint16のキー空間は65,536エントリのみであるため、単純な配列で実装する。

```
ファイル構造:
+------------------+
| Header (32 bytes)|
+------------------+
| Array[65536]     | <- uint64 * 65536 = 512KB
| (appearance_nrow)|
+------------------+

Header構造:
- magic: uint32 (0x4B4D5231 = "KMR1")
- version: uint32 (1)
- key_type: uint32 (16 = uint16)
- entry_count: uint64 (有効エントリ数)
- reserved: uint64 (予約)
- checksum: uint32 (簡易チェックサム)
```

**実装詳細**:
- ファイルサイズ固定: 32 + 65536 * 8 = 524,320 bytes
- 配列インデックス = uintkey
- appearance_nrow = 0 は未使用エントリを意味
- 読み書きはpread()/pwrite()で直接アクセス（キャッシュはOSに任せる）

**add()の動作**:
- 別の行での出現: 既存値に加算（preadで読み取り、加算、pwriteで書き戻し）
- 同一行内での重複出現: 呼び出し元で予め重複除外されているため考慮不要
- 上書きが必要になる事態は発生しない

### 2. uint32 key - uint64 value（ハッシュテーブル実装）

uint32のキー空間は最大約43億エントリであり、実際に使用されるエントリ数は限定的なためハッシュテーブルで実装。

```
ファイル構造:
+------------------+
| Header (64 bytes)|
+------------------+
| Bucket Directory | <- uint64 * bucket_count (各バケットの開始オフセット)
+------------------+
| Entry Area       | <- 可変長、エントリを追加していく
+------------------+

Header構造:
- magic: uint32 (0x4B4D5232 = "KMR2")
- version: uint32 (1)
- key_type: uint32 (32 = uint32)
- bucket_count: uint32 (バケット数、2のべき乗)
- entry_count: uint64 (有効エントリ数)
- next_entry_offset: uint64 (次のエントリ書き込み位置)
- reserved: uint64 * 4
- checksum: uint32

Entry構造 (20 bytes):
- uintkey: uint32 (4 bytes)
- appearance_nrow: uint64 (8 bytes)
- next_offset: uint64 (同一バケット内の次エントリへのオフセット、0=終端)
```

**ハッシュ関数**: MurmurHash3の32ビット版

```c
static inline uint32
kmersearch_murmurhash32(uint32 key)
{
    uint32 h = key;
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}
```

**バケット数の決定**:
- 初期バケット数: kmersearch_highfreq_analysis_hashtable_size / 4（ロードファクター0.25目標）
- 最小: 4096
- 最大: 16,777,216

**add()の動作**:
- 別の行での出現: 既存エントリを検索し、見つかれば加算、なければ新規エントリ追加
- 同一行内での重複出現: 呼び出し元で予め重複除外されているため考慮不要
- 上書きが必要になる事態は発生しない

### 3. uint64 key - uint64 value（ハッシュテーブル実装）

uint32と同様の構造だが、エントリサイズが異なる。

```
ファイル構造:
+------------------+
| Header (64 bytes)|
+------------------+
| Bucket Directory | <- uint64 * bucket_count (各バケットの開始オフセット)
+------------------+
| Entry Area       | <- 可変長、エントリを追加していく
+------------------+

Header構造:
- magic: uint32 (0x4B4D5233 = "KMR3")
- version: uint32 (1)
- key_type: uint32 (64 = uint64)
- bucket_count: uint32 (バケット数、2のべき乗)
- entry_count: uint64 (有効エントリ数)
- next_entry_offset: uint64 (次のエントリ書き込み位置)
- reserved: uint64 * 4
- checksum: uint32

Entry構造 (24 bytes):
- uintkey: uint64 (8 bytes)
- appearance_nrow: uint64 (8 bytes)
- next_offset: uint64 (8 bytes)
```

**ハッシュ関数**: MurmurHash3の64ビット版

```c
static inline uint64
kmersearch_murmurhash64(uint64 key)
{
    uint64 h = key;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}
```

**add()の動作**: uint32と同様

---

## マージ処理の詳細

並列ワーカー終了後、複数の一時ファイルを集約する際のマージ処理：

1. **マージ対象**: ファイルA（マージ先）とファイルB（マージ元）
2. **処理内容**: ファイルBの全エントリを走査し、各uintkeyのappearance_nrowをファイルAに加算
3. **終了処理**: マージ完了後、ファイルBを削除（unlink）
4. **マージ順序**: 現在のSQLite3実装と同様、ペアワイズで並列マージを繰り返す

---

## API設計

### 新規関数

```c
/* ファイルハッシュテーブルのコンテキスト */
typedef struct FileHashTable16Context
{
    int         fd;
    char        path[MAXPGPATH];
    uint64      entry_count;
} FileHashTable16Context;

typedef struct FileHashTable32Context
{
    int         fd;
    char        path[MAXPGPATH];
    uint32      bucket_count;
    uint64      entry_count;
    uint64      next_entry_offset;
} FileHashTable32Context;

typedef struct FileHashTable64Context
{
    int         fd;
    char        path[MAXPGPATH];
    uint32      bucket_count;
    uint64      entry_count;
    uint64      next_entry_offset;
} FileHashTable64Context;

/* uint16用API */
FileHashTable16Context *kmersearch_fht16_create(const char *path);
FileHashTable16Context *kmersearch_fht16_open(const char *path);
void kmersearch_fht16_close(FileHashTable16Context *ctx);
void kmersearch_fht16_add(FileHashTable16Context *ctx, uint16 uintkey, uint64 appearance_nrow);
uint64 kmersearch_fht16_get(FileHashTable16Context *ctx, uint16 uintkey);
void kmersearch_fht16_flush(FileHashTable16Context *ctx);

/* uint32用API */
FileHashTable32Context *kmersearch_fht32_create(const char *path, uint32 bucket_count);
FileHashTable32Context *kmersearch_fht32_open(const char *path);
void kmersearch_fht32_close(FileHashTable32Context *ctx);
void kmersearch_fht32_add(FileHashTable32Context *ctx, uint32 uintkey, uint64 appearance_nrow);
uint64 kmersearch_fht32_get(FileHashTable32Context *ctx, uint32 uintkey);
void kmersearch_fht32_flush(FileHashTable32Context *ctx);

/* uint64用API */
FileHashTable64Context *kmersearch_fht64_create(const char *path, uint32 bucket_count);
FileHashTable64Context *kmersearch_fht64_open(const char *path);
void kmersearch_fht64_close(FileHashTable64Context *ctx);
void kmersearch_fht64_add(FileHashTable64Context *ctx, uint64 uintkey, uint64 appearance_nrow);
uint64 kmersearch_fht64_get(FileHashTable64Context *ctx, uint64 uintkey);
void kmersearch_fht64_flush(FileHashTable64Context *ctx);

/* マージ用API */
void kmersearch_fht16_merge(const char *source_path, const char *target_path);
void kmersearch_fht32_merge(const char *source_path, const char *target_path);
void kmersearch_fht64_merge(const char *source_path, const char *target_path);

/* イテレータAPI（結果取得用） */
typedef struct FileHashTableIterator16
{
    FileHashTable16Context *ctx;
    uint32      current_index;
} FileHashTableIterator16;

typedef struct FileHashTableIterator32
{
    FileHashTable32Context *ctx;
    uint32      current_bucket;
    uint64      current_offset;
} FileHashTableIterator32;

typedef struct FileHashTableIterator64
{
    FileHashTable64Context *ctx;
    uint32      current_bucket;
    uint64      current_offset;
} FileHashTableIterator64;

/* イテレータ初期化API */
void kmersearch_fht16_iterator_init(FileHashTableIterator16 *iter, FileHashTable16Context *ctx);
void kmersearch_fht32_iterator_init(FileHashTableIterator32 *iter, FileHashTable32Context *ctx);
void kmersearch_fht64_iterator_init(FileHashTableIterator64 *iter, FileHashTable64Context *ctx);

/* イテレータ走査API（次のエントリを取得、エントリがなければfalseを返す） */
bool kmersearch_fht16_iterate(FileHashTableIterator16 *iter, uint16 *uintkey, uint64 *appearance_nrow);
bool kmersearch_fht32_iterate(FileHashTableIterator32 *iter, uint32 *uintkey, uint64 *appearance_nrow);
bool kmersearch_fht64_iterate(FileHashTableIterator64 *iter, uint64 *uintkey, uint64 *appearance_nrow);
```

---

## 型選択ロジック

既存のkmersearch_freq.cと同様のロジックで、使用するファイルハッシュテーブルの型を決定する。

```c
int total_bits = kmersearch_kmer_size * 2 + kmersearch_occur_bitlen;

if (total_bits <= 16) {
    /* uint16配列を使用 */
} else if (total_bits <= 32) {
    /* uint32ハッシュテーブルを使用 */
} else {
    /* uint64ハッシュテーブルを使用 */
}
```

**参考: total_bitsの具体例**:
| kmer_size | occur_bitlen | total_bits | 使用する型 |
|-----------|--------------|------------|-----------|
| 4 | 0 | 8 | uint16 |
| 8 | 0 | 16 | uint16 |
| 8 | 8 | 24 | uint32 |
| 16 | 0 | 32 | uint32 |
| 16 | 8 | 40 | uint64 |
| 32 | 0 | 64 | uint64 |

---

## 実装手順

### フェーズ1: 新規ファイル作成とデータ構造実装

1. **新規ファイル作成**: `kmersearch_tmpfile.c`（一時ファイルベースハッシュテーブル）
2. **ヘッダーファイル更新**: `kmersearch.h`に新しい構造体とAPI宣言を追加
3. **Makefile更新**: `kmersearch_tmpfile.c`を追加

### フェーズ2: uint16配列実装

1. `kmersearch_fht16_create()` - ファイル作成と初期化
2. `kmersearch_fht16_add()` - appearance_nrowの加算
3. `kmersearch_fht16_get()` - appearance_nrowの取得
4. `kmersearch_fht16_flush()` - バッファのフラッシュ
5. `kmersearch_fht16_close()` - ファイルクローズ
6. `kmersearch_fht16_merge()` - 2つのファイルをマージ
7. `kmersearch_fht16_iterate()` - 全エントリの走査

### フェーズ3: uint32ハッシュテーブル実装

1. MurmurHash32実装
2. `kmersearch_fht32_create()` - ファイル作成、バケットディレクトリ初期化
3. `kmersearch_fht32_add()` - エントリ追加（チェイン法）
4. `kmersearch_fht32_get()` - エントリ検索
5. `kmersearch_fht32_flush()` - ヘッダー更新
6. `kmersearch_fht32_close()` - ファイルクローズ
7. `kmersearch_fht32_merge()` - マージ処理
8. `kmersearch_fht32_iterate()` - 全エントリの走査

### フェーズ4: uint64ハッシュテーブル実装

1. MurmurHash64実装
2. uint32と同様のAPI実装

### フェーズ5: kmersearch_freq.cの修正

1. **SQLiteWorkerContext構造体の置き換え**:
   ```c
   typedef struct FileHashWorkerContext
   {
       char        file_path[MAXPGPATH];
       HTAB        *batch_hash;
       int         batch_count;
       int         total_bits;
       Oid         dna2_oid;
       Oid         dna4_oid;
       Oid         column_type_oid;
       MemoryContext batch_memory_context;
       BufferAccessStrategy strategy;
       /* SQLite3関連を削除し、FileHashTable*Contextに置き換え */
       void        *fht_ctx;  /* FileHashTable16/32/64Contextのいずれか */
   } FileHashWorkerContext;
   ```

2. **kmersearch_flush_batch_to_sqlite()をkmersearch_flush_batch_to_fht()に置き換え**

3. **kmersearch_analysis_worker()の修正**:
   - SQLite3初期化コードを削除
   - FileHashTable API呼び出しに置き換え

4. **kmersearch_parallel_merge_worker()の修正**:
   - SQLite3マージ処理をFileHashTableマージ処理に置き換え

5. **結果集約処理の修正**:
   - SQLite3からの読み取りをFileHashTableイテレータに置き換え

### フェーズ6: kmersearch_highfreq_kmerテーブルの変更

1. **pg_kmersearch--1.0.sqlの修正**:
   ```sql
   CREATE TABLE kmersearch_highfreq_kmer (
       table_oid oid NOT NULL,
       column_name name NOT NULL,
       uintkey bigint NOT NULL,
       appearance_nrow bigint NOT NULL,  -- 新規追加
       detection_reason text,
       created_at timestamp with time zone DEFAULT now(),
       PRIMARY KEY (table_oid, column_name, uintkey)
   );
   ```

2. **kmersearch_freq.cの結果保存処理修正**:
   - PostgreSQLテーブルへのINSERT時に`appearance_nrow`を含める

### フェーズ7: SQLite3依存の削除

1. `#include <sqlite3.h>`を削除
2. Makefileから`-lsqlite3`を削除
3. SQLite3関連の設定コード削除（`sqlite3_config`、`sqlite3_soft_heap_limit64`等）

---

## テスト計画

### 単体テスト

1. **uint16配列テスト**:
   - 空の配列に値を追加
   - 同一キーへの複数回加算
   - 全キー空間（0-65535）への書き込み
   - マージ処理の検証

2. **uint32ハッシュテーブルテスト**:
   - 衝突のないケース
   - 同一バケットへの複数エントリ追加（チェイン）
   - マージ処理の検証
   - 大量エントリ（100万件）のパフォーマンス

3. **uint64ハッシュテーブルテスト**:
   - uint32と同様のテスト

### 統合テスト

1. 既存の回帰テスト（09_highfreq_filter.sql等）が全てパス
2. 並列ワーカー処理の検証
3. パーティションテーブルでの動作確認
4. メモリ使用量の確認（SQLite3比較）
5. 処理時間の確認（SQLite3比較）

---

## 期待される改善

1. **パフォーマンス向上**:
   - SQL解析オーバーヘッドの削除
   - トランザクション管理オーバーヘッドの削除
   - 直接ファイルアクセスによる高速化

2. **メモリ使用量削減**:
   - SQLite3の内部バッファ不要
   - より単純なデータ構造

3. **依存関係削減**:
   - SQLite3ライブラリへの依存を削除
   - ビルド環境の簡素化

4. **コードの簡素化**:
   - SQLite3のエラーハンドリング削除
   - より直接的なデータアクセス

---

## リスクと対策

### リスク1: ディスクI/Oパフォーマンス

**対策**: OSのディスクキャッシュを活用。バッチ書き込みでI/O回数を最小化。

### リスク2: ファイル破損

**対策**: ヘッダーにマジックナンバーとバージョン、チェックサムを含める。分析途中でエラーが発生した場合は一時ファイルを削除。

### リスク3: 大規模データでのメモリ不足

**対策**: uint32/uint64のバケット数を適切に設定。メモリ上のバッチハッシュテーブルサイズはGUC変数で制御可能。

---

## 補足: nrowからappearance_nrowへの名称変更

現在のコードベースで`nrow`として使用されている「並列ワーカーの担当範囲におけるそのuintkeyの出現行数」を、より明確な`appearance_nrow`に変更する。

### 変更箇所

1. **kmersearch_freq.c**:
   - `TempKmerFreqEntry16/32/64`構造体の`nrow`メンバーを`appearance_nrow`に変更
   - SQLite3テーブル作成SQLの`nrow`カラムを`appearance_nrow`に変更
   - 関連する全ての変数名・コメントを更新

2. **pg_kmersearch--1.0.sql**:
   - `kmersearch_highfreq_kmer`テーブルに`appearance_nrow`カラムを追加

3. **回帰テスト**:
   - `appearance_nrow`カラムを参照するテストケースを追加
