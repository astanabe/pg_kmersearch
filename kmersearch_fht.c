/*
 * kmersearch_fht.c - File-based hash table for temporary k-mer storage
 *
 * This module provides efficient file-based storage for k-mer frequency data
 * during high-frequency k-mer analysis. It replaces SQLite3 with direct file
 * I/O for improved performance.
 *
 * Three implementations are provided for different key sizes:
 * - uint16: Direct array (65536 entries, 512KB fixed size)
 * - uint32: Chain-based hash table with MurmurHash3
 * - uint64: Chain-based hash table with MurmurHash3 64-bit
 */

#include "kmersearch.h"
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

/* Magic numbers for file format validation */
#define FHT16_MAGIC  0x4B4D5231  /* "KMR1" */
#define FHT32_MAGIC  0x4B4D5232  /* "KMR2" */
#define FHT64_MAGIC  0x4B4D5233  /* "KMR3" */
#define FHT_VERSION  1

/* File structure sizes */
#define FHT16_HEADER_SIZE  32
#define FHT32_HEADER_SIZE  64
#define FHT64_HEADER_SIZE  64

/* Array size for uint16 (fixed) */
#define FHT16_ARRAY_SIZE   65536

/* Entry sizes */
#define FHT32_ENTRY_SIZE   20  /* uint32 key + uint64 value + uint64 next_offset */
#define FHT64_ENTRY_SIZE   24  /* uint64 key + uint64 value + uint64 next_offset */

/* Default bucket counts */
#define FHT_MIN_BUCKET_COUNT       4096
#define FHT_MAX_BUCKET_COUNT       16777216
#define FHT_DEFAULT_LOAD_FACTOR    4

/* Memory overhead factor for HTAB (approximately 2x for safety) */
#define FHT_HTAB_OVERHEAD_FACTOR   2

/*
 * File header structures
 */
typedef struct FileHashTable16Header
{
    uint32      magic;              /* Magic number: FHT16_MAGIC */
    uint32      version;            /* Version: 1 */
    uint32      key_type;           /* Key type: 16 */
    uint32      reserved1;          /* Reserved */
    uint64      entry_count;        /* Number of non-zero entries */
    uint64      reserved2;          /* Reserved */
    uint32      checksum;           /* Simple checksum */
    uint32      reserved3;          /* Reserved */
} FileHashTable16Header;

typedef struct FileHashTable32Header
{
    uint32      magic;              /* Magic number: FHT32_MAGIC */
    uint32      version;            /* Version: 1 */
    uint32      key_type;           /* Key type: 32 */
    uint32      bucket_count;       /* Number of buckets (power of 2) */
    uint64      entry_count;        /* Number of entries */
    uint64      next_entry_offset;  /* Next entry write position */
    uint64      reserved[4];        /* Reserved for future use */
    uint32      checksum;           /* Simple checksum */
    uint32      reserved2;          /* Reserved */
} FileHashTable32Header;

typedef struct FileHashTable64Header
{
    uint32      magic;              /* Magic number: FHT64_MAGIC */
    uint32      version;            /* Version: 1 */
    uint32      key_type;           /* Key type: 64 */
    uint32      bucket_count;       /* Number of buckets (power of 2) */
    uint64      entry_count;        /* Number of entries */
    uint64      next_entry_offset;  /* Next entry write position */
    uint64      reserved[4];        /* Reserved for future use */
    uint32      checksum;           /* Simple checksum */
    uint32      reserved2;          /* Reserved */
} FileHashTable64Header;

/*
 * On-disk entry structures for hash tables
 */
typedef struct FileHashTable32Entry
{
    uint32      uintkey;            /* K-mer key */
    uint64      appearance_nrow;    /* Number of rows containing this k-mer */
    uint64      next_offset;        /* Offset to next entry in chain (0 = end) */
} FileHashTable32Entry;

typedef struct FileHashTable64Entry
{
    uint64      uintkey;            /* K-mer key */
    uint64      appearance_nrow;    /* Number of rows containing this k-mer */
    uint64      next_offset;        /* Offset to next entry in chain (0 = end) */
} FileHashTable64Entry;

/*
 * MurmurHash3 finalization mix for 32-bit keys
 */
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

/*
 * MurmurHash3 finalization mix for 64-bit keys
 */
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

/*
 * Calculate optimal bucket count based on expected entries
 */
static inline uint32
kmersearch_calculate_bucket_count(uint64 expected_entries)
{
    uint32 bucket_count;

    bucket_count = (uint32)(expected_entries / FHT_DEFAULT_LOAD_FACTOR);

    if (bucket_count < FHT_MIN_BUCKET_COUNT)
        bucket_count = FHT_MIN_BUCKET_COUNT;
    if (bucket_count > FHT_MAX_BUCKET_COUNT)
        bucket_count = FHT_MAX_BUCKET_COUNT;

    /* Round up to next power of 2 */
    bucket_count--;
    bucket_count |= bucket_count >> 1;
    bucket_count |= bucket_count >> 2;
    bucket_count |= bucket_count >> 4;
    bucket_count |= bucket_count >> 8;
    bucket_count |= bucket_count >> 16;
    bucket_count++;

    return bucket_count;
}

/*
 * ============================================================================
 * uint16 Array Implementation (FHT16)
 * ============================================================================
 */

/*
 * Create a new uint16 file hash table
 */
FileHashTable16Context *
kmersearch_fht16_create(const char *path)
{
    FileHashTable16Context *ctx;
    FileHashTable16Header header;
    uint64 *array;
    ssize_t written;

    ctx = palloc0(sizeof(FileHashTable16Context));
    strlcpy(ctx->path, path, MAXPGPATH);

    /* Create and open file */
    ctx->fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0600);
    if (ctx->fd < 0)
    {
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not create file hash table \"%s\": %m", path)));
    }

    /* Initialize header */
    memset(&header, 0, sizeof(header));
    header.magic = FHT16_MAGIC;
    header.version = FHT_VERSION;
    header.key_type = 16;
    header.entry_count = 0;

    /* Write header */
    written = write(ctx->fd, &header, sizeof(header));
    if (written != sizeof(header))
    {
        close(ctx->fd);
        unlink(path);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write file hash table header: %m")));
    }

    /* Initialize array with zeros */
    array = palloc0(FHT16_ARRAY_SIZE * sizeof(uint64));
    written = write(ctx->fd, array, FHT16_ARRAY_SIZE * sizeof(uint64));
    pfree(array);

    if (written != (ssize_t)(FHT16_ARRAY_SIZE * sizeof(uint64)))
    {
        close(ctx->fd);
        unlink(path);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not initialize file hash table array: %m")));
    }

    ctx->entry_count = 0;

    return ctx;
}

/*
 * Open an existing uint16 file hash table
 */
FileHashTable16Context *
kmersearch_fht16_open(const char *path)
{
    FileHashTable16Context *ctx;
    FileHashTable16Header header;
    ssize_t bytes_read;

    ctx = palloc0(sizeof(FileHashTable16Context));
    strlcpy(ctx->path, path, MAXPGPATH);

    ctx->fd = open(path, O_RDWR);
    if (ctx->fd < 0)
    {
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not open file hash table \"%s\": %m", path)));
    }

    /* Read and validate header */
    bytes_read = pread(ctx->fd, &header, sizeof(header), 0);
    if (bytes_read != sizeof(header))
    {
        close(ctx->fd);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read file hash table header: %m")));
    }

    if (header.magic != FHT16_MAGIC || header.version != FHT_VERSION)
    {
        close(ctx->fd);
        pfree(ctx);
        ereport(ERROR,
                (errcode(ERRCODE_DATA_CORRUPTED),
                 errmsg("invalid file hash table format")));
    }

    ctx->entry_count = header.entry_count;

    return ctx;
}

/*
 * Close a uint16 file hash table
 */
void
kmersearch_fht16_close(FileHashTable16Context *ctx)
{
    if (ctx == NULL)
        return;

    if (ctx->fd >= 0)
    {
        /* Update entry count in header before closing */
        FileHashTable16Header header;

        if (pread(ctx->fd, &header, sizeof(header), 0) == sizeof(header))
        {
            header.entry_count = ctx->entry_count;
            if (pwrite(ctx->fd, &header, sizeof(header), 0) != sizeof(header))
                ereport(WARNING,
                        (errcode_for_file_access(),
                         errmsg("could not update file hash table header: %m")));
        }

        fsync(ctx->fd);
        close(ctx->fd);
    }

    pfree(ctx);
}

/*
 * Add appearance_nrow to a uint16 key
 */
void
kmersearch_fht16_add(FileHashTable16Context *ctx, uint16 uintkey, uint64 appearance_nrow)
{
    off_t offset;
    uint64 current_value;
    ssize_t bytes_read;
    ssize_t bytes_written;

    offset = FHT16_HEADER_SIZE + (off_t)uintkey * sizeof(uint64);

    /* Read current value */
    bytes_read = pread(ctx->fd, &current_value, sizeof(uint64), offset);
    if (bytes_read != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read from file hash table: %m")));
    }

    /* Track new entries */
    if (current_value == 0 && appearance_nrow > 0)
        ctx->entry_count++;

    /* Add to current value */
    current_value += appearance_nrow;

    /* Write back */
    bytes_written = pwrite(ctx->fd, &current_value, sizeof(uint64), offset);
    if (bytes_written != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write to file hash table: %m")));
    }
}

/*
 * Get appearance_nrow for a uint16 key
 */
uint64
kmersearch_fht16_get(FileHashTable16Context *ctx, uint16 uintkey)
{
    off_t offset;
    uint64 value;
    ssize_t bytes_read;

    offset = FHT16_HEADER_SIZE + (off_t)uintkey * sizeof(uint64);

    bytes_read = pread(ctx->fd, &value, sizeof(uint64), offset);
    if (bytes_read != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read from file hash table: %m")));
    }

    return value;
}

/*
 * Flush any pending writes
 */
void
kmersearch_fht16_flush(FileHashTable16Context *ctx)
{
    if (ctx && ctx->fd >= 0)
    {
        /* Update entry count in header */
        FileHashTable16Header header;

        if (pread(ctx->fd, &header, sizeof(header), 0) == sizeof(header))
        {
            header.entry_count = ctx->entry_count;
            if (pwrite(ctx->fd, &header, sizeof(header), 0) != sizeof(header))
                ereport(WARNING,
                        (errcode_for_file_access(),
                         errmsg("could not update file hash table header: %m")));
        }

        fsync(ctx->fd);
    }
}

/*
 * Bulk add from memory array to file hash table
 * Reads current file values, adds memory array values, writes back in one operation
 */
void
kmersearch_fht16_bulk_add(FileHashTable16Context *ctx, uint64 *memory_array)
{
    uint64 file_array[FHT16_ARRAY_SIZE];
    ssize_t bytes;
    uint64 entry_count = 0;
    FileHashTable16Header header;

    bytes = pread(ctx->fd, file_array, sizeof(file_array), FHT16_HEADER_SIZE);
    if (bytes != sizeof(file_array))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read file hash table array: %m")));
    }

    for (int i = 0; i < FHT16_ARRAY_SIZE; i++)
    {
        file_array[i] += memory_array[i];
        if (file_array[i] > 0)
            entry_count++;
    }

    bytes = pwrite(ctx->fd, file_array, sizeof(file_array), FHT16_HEADER_SIZE);
    if (bytes != sizeof(file_array))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write file hash table array: %m")));
    }

    bytes = pread(ctx->fd, &header, sizeof(header), 0);
    if (bytes == sizeof(header))
    {
        header.entry_count = entry_count;
        if (pwrite(ctx->fd, &header, sizeof(header), 0) != sizeof(header))
            ereport(WARNING,
                    (errcode_for_file_access(),
                     errmsg("could not update file hash table header: %m")));
    }

    ctx->entry_count = entry_count;
    fsync(ctx->fd);
}

/*
 * Merge source file into target file (target += source)
 */
void
kmersearch_fht16_merge(const char *source_path, const char *target_path)
{
    int source_fd;
    int target_fd;
    uint64 source_array[FHT16_ARRAY_SIZE];
    uint64 target_array[FHT16_ARRAY_SIZE];
    FileHashTable16Header header;
    ssize_t bytes;

    source_fd = open(source_path, O_RDONLY);
    if (source_fd < 0)
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not open source file \"%s\": %m", source_path)));
    }

    target_fd = open(target_path, O_RDWR);
    if (target_fd < 0)
    {
        close(source_fd);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not open target file \"%s\": %m", target_path)));
    }

    /* Read source array */
    bytes = pread(source_fd, source_array, sizeof(source_array), FHT16_HEADER_SIZE);
    if (bytes != sizeof(source_array))
    {
        close(source_fd);
        close(target_fd);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read source array: %m")));
    }

    /* Read target array */
    bytes = pread(target_fd, target_array, sizeof(target_array), FHT16_HEADER_SIZE);
    if (bytes != sizeof(target_array))
    {
        close(source_fd);
        close(target_fd);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read target array: %m")));
    }

    /* Merge: add source values to target */
    for (int i = 0; i < FHT16_ARRAY_SIZE; i++)
    {
        target_array[i] += source_array[i];
    }

    /* Write merged array back to target */
    bytes = pwrite(target_fd, target_array, sizeof(target_array), FHT16_HEADER_SIZE);
    if (bytes != sizeof(target_array))
    {
        close(source_fd);
        close(target_fd);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write merged array: %m")));
    }

    /* Update entry count in target header */
    bytes = pread(target_fd, &header, sizeof(header), 0);
    if (bytes == sizeof(header))
    {
        uint64 count = 0;
        for (int i = 0; i < FHT16_ARRAY_SIZE; i++)
        {
            if (target_array[i] > 0)
                count++;
        }
        header.entry_count = count;
        if (pwrite(target_fd, &header, sizeof(header), 0) != sizeof(header))
            ereport(WARNING,
                    (errcode_for_file_access(),
                     errmsg("could not update file hash table header: %m")));
    }

    fsync(target_fd);
    close(source_fd);
    close(target_fd);

    /* Delete source file after successful merge */
    unlink(source_path);
}

/*
 * Initialize iterator for uint16 file hash table
 */
void
kmersearch_fht16_iterator_init(FileHashTableIterator16 *iter, FileHashTable16Context *ctx)
{
    iter->ctx = ctx;
    iter->current_index = 0;
}

/*
 * Get next entry from uint16 file hash table iterator
 */
bool
kmersearch_fht16_iterate(FileHashTableIterator16 *iter, uint16 *uintkey, uint64 *appearance_nrow)
{
    off_t offset;
    uint64 value;
    ssize_t bytes_read;

    while (iter->current_index < FHT16_ARRAY_SIZE)
    {
        offset = FHT16_HEADER_SIZE + (off_t)iter->current_index * sizeof(uint64);

        bytes_read = pread(iter->ctx->fd, &value, sizeof(uint64), offset);
        if (bytes_read != sizeof(uint64))
        {
            return false;
        }

        if (value > 0)
        {
            *uintkey = (uint16)iter->current_index;
            *appearance_nrow = value;
            iter->current_index++;
            return true;
        }

        iter->current_index++;
    }

    return false;
}

/*
 * ============================================================================
 * uint32 Hash Table Implementation (FHT32)
 * ============================================================================
 */

/*
 * Create a new uint32 file hash table
 */
FileHashTable32Context *
kmersearch_fht32_create(const char *path, uint32 bucket_count)
{
    FileHashTable32Context *ctx;
    FileHashTable32Header header;
    uint64 *bucket_array;
    off_t data_start;
    ssize_t written;

    if (bucket_count == 0)
        bucket_count = kmersearch_calculate_bucket_count(kmersearch_highfreq_analysis_hashtable_size);

    ctx = palloc0(sizeof(FileHashTable32Context));
    strlcpy(ctx->path, path, MAXPGPATH);
    ctx->bucket_count = bucket_count;

    /* Create and open file */
    ctx->fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0600);
    if (ctx->fd < 0)
    {
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not create file hash table \"%s\": %m", path)));
    }

    /* Initialize header */
    memset(&header, 0, sizeof(header));
    header.magic = FHT32_MAGIC;
    header.version = FHT_VERSION;
    header.key_type = 32;
    header.bucket_count = bucket_count;
    header.entry_count = 0;

    /* Calculate data start position (after header and bucket directory) */
    data_start = FHT32_HEADER_SIZE + (off_t)bucket_count * sizeof(uint64);
    header.next_entry_offset = data_start;
    ctx->next_entry_offset = data_start;

    /* Write header */
    written = write(ctx->fd, &header, sizeof(header));
    if (written != sizeof(header))
    {
        close(ctx->fd);
        unlink(path);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write file hash table header: %m")));
    }

    /* Initialize bucket directory with zeros (0 = empty bucket) */
    bucket_array = palloc0(bucket_count * sizeof(uint64));
    written = write(ctx->fd, bucket_array, bucket_count * sizeof(uint64));
    pfree(bucket_array);

    if (written != (ssize_t)(bucket_count * sizeof(uint64)))
    {
        close(ctx->fd);
        unlink(path);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not initialize bucket directory: %m")));
    }

    ctx->entry_count = 0;

    return ctx;
}

/*
 * Open an existing uint32 file hash table
 */
FileHashTable32Context *
kmersearch_fht32_open(const char *path)
{
    FileHashTable32Context *ctx;
    FileHashTable32Header header;
    ssize_t bytes_read;

    ctx = palloc0(sizeof(FileHashTable32Context));
    strlcpy(ctx->path, path, MAXPGPATH);

    ctx->fd = open(path, O_RDWR);
    if (ctx->fd < 0)
    {
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not open file hash table \"%s\": %m", path)));
    }

    /* Read and validate header */
    bytes_read = pread(ctx->fd, &header, sizeof(header), 0);
    if (bytes_read != sizeof(header))
    {
        close(ctx->fd);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read file hash table header: %m")));
    }

    if (header.magic != FHT32_MAGIC || header.version != FHT_VERSION)
    {
        close(ctx->fd);
        pfree(ctx);
        ereport(ERROR,
                (errcode(ERRCODE_DATA_CORRUPTED),
                 errmsg("invalid file hash table format")));
    }

    ctx->bucket_count = header.bucket_count;
    ctx->entry_count = header.entry_count;
    ctx->next_entry_offset = header.next_entry_offset;

    return ctx;
}

/*
 * Close a uint32 file hash table
 */
void
kmersearch_fht32_close(FileHashTable32Context *ctx)
{
    if (ctx == NULL)
        return;

    if (ctx->fd >= 0)
    {
        kmersearch_fht32_flush(ctx);
        close(ctx->fd);
    }

    pfree(ctx);
}

/*
 * Add appearance_nrow to a uint32 key
 */
void
kmersearch_fht32_add(FileHashTable32Context *ctx, uint32 uintkey, uint64 appearance_nrow)
{
    uint32 bucket_idx;
    off_t bucket_offset;
    uint64 entry_offset;
    FileHashTable32Entry entry;
    ssize_t bytes;

    bucket_idx = kmersearch_murmurhash32(uintkey) & (ctx->bucket_count - 1);
    bucket_offset = FHT32_HEADER_SIZE + (off_t)bucket_idx * sizeof(uint64);

    /* Read bucket head offset */
    bytes = pread(ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read bucket: %m")));
    }

    /* Search for existing entry in chain */
    while (entry_offset != 0)
    {
        bytes = pread(ctx->fd, &entry, sizeof(entry), entry_offset);
        if (bytes != sizeof(entry))
        {
            ereport(ERROR,
                    (errcode_for_file_access(),
                     errmsg("could not read entry: %m")));
        }

        if (entry.uintkey == uintkey)
        {
            /* Found existing entry, add to its value */
            entry.appearance_nrow += appearance_nrow;
            bytes = pwrite(ctx->fd, &entry, sizeof(entry), entry_offset);
            if (bytes != sizeof(entry))
            {
                ereport(ERROR,
                        (errcode_for_file_access(),
                         errmsg("could not update entry: %m")));
            }
            return;
        }

        entry_offset = entry.next_offset;
    }

    /* Entry not found, create new one */
    entry.uintkey = uintkey;
    entry.appearance_nrow = appearance_nrow;

    /* Read current bucket head to set as next */
    bytes = pread(ctx->fd, &entry.next_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read bucket head: %m")));
    }

    /* Write new entry at next available position */
    entry_offset = ctx->next_entry_offset;
    bytes = pwrite(ctx->fd, &entry, sizeof(entry), entry_offset);
    if (bytes != sizeof(entry))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write new entry: %m")));
    }

    /* Update bucket head to point to new entry */
    bytes = pwrite(ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not update bucket head: %m")));
    }

    ctx->next_entry_offset += sizeof(entry);
    ctx->entry_count++;
}

/*
 * Get appearance_nrow for a uint32 key
 */
uint64
kmersearch_fht32_get(FileHashTable32Context *ctx, uint32 uintkey)
{
    uint32 bucket_idx;
    off_t bucket_offset;
    uint64 entry_offset;
    FileHashTable32Entry entry;
    ssize_t bytes;

    bucket_idx = kmersearch_murmurhash32(uintkey) & (ctx->bucket_count - 1);
    bucket_offset = FHT32_HEADER_SIZE + (off_t)bucket_idx * sizeof(uint64);

    /* Read bucket head offset */
    bytes = pread(ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        return 0;
    }

    /* Search chain for entry */
    while (entry_offset != 0)
    {
        bytes = pread(ctx->fd, &entry, sizeof(entry), entry_offset);
        if (bytes != sizeof(entry))
        {
            return 0;
        }

        if (entry.uintkey == uintkey)
        {
            return entry.appearance_nrow;
        }

        entry_offset = entry.next_offset;
    }

    return 0;
}

/*
 * Flush uint32 file hash table
 */
void
kmersearch_fht32_flush(FileHashTable32Context *ctx)
{
    FileHashTable32Header header;
    ssize_t bytes;

    if (ctx == NULL || ctx->fd < 0)
        return;

    /* Read current header */
    bytes = pread(ctx->fd, &header, sizeof(header), 0);
    if (bytes == sizeof(header))
    {
        header.entry_count = ctx->entry_count;
        header.next_entry_offset = ctx->next_entry_offset;
        if (pwrite(ctx->fd, &header, sizeof(header), 0) != sizeof(header))
            ereport(WARNING,
                    (errcode_for_file_access(),
                     errmsg("could not update file hash table header: %m")));
    }

    fsync(ctx->fd);
}

/*
 * Bulk add from batch hash table to file hash table
 * Reads FHT file into memory, merges with batch_hash, writes back
 */
void
kmersearch_fht32_bulk_add(FileHashTable32Context *ctx, HTAB *batch_hash)
{
    HTAB *merge_htab;
    HASHCTL hash_ctl;
    MemoryContext merge_context;
    MemoryContext old_context;
    FileHashTableIterator32 iter;
    uint32 uintkey;
    uint64 appearance_nrow;
    KmerFreqEntry32 *entry;
    bool found;
    HASH_SEQ_STATUS hash_seq;
    HASH_SEQ_STATUS batch_seq;
    void *batch_entry;
    uint32 saved_bucket_count;
    char saved_path[MAXPGPATH];
    long max_entries;
    FileHashTable32Context *new_ctx;

    merge_context = AllocSetContextCreate(CurrentMemoryContext,
                                          "FHT32 BulkAdd Context",
                                          ALLOCSET_DEFAULT_SIZES);
    old_context = MemoryContextSwitchTo(merge_context);

    max_entries = ctx->entry_count + hash_get_num_entries(batch_hash);
    if (max_entries < 1024)
        max_entries = 1024;

    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint32);
    hash_ctl.entrysize = sizeof(KmerFreqEntry32);
    hash_ctl.hcxt = merge_context;

    merge_htab = hash_create("FHT32 BulkAdd Hash",
                             max_entries,
                             &hash_ctl,
                             HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);

    kmersearch_fht32_iterator_init(&iter, ctx);
    while (kmersearch_fht32_iterate(&iter, &uintkey, &appearance_nrow))
    {
        entry = (KmerFreqEntry32 *) hash_search(merge_htab, &uintkey, HASH_ENTER, &found);
        if (found)
            entry->appearance_nrow += appearance_nrow;
        else
            entry->appearance_nrow = appearance_nrow;
    }

    hash_seq_init(&batch_seq, batch_hash);
    while ((batch_entry = hash_seq_search(&batch_seq)) != NULL)
    {
        /*
         * Cast to KmerFreqEntry32 which has the same memory layout as
         * TempKmerFreqEntry32 for the first two fields (uint32 + uint64).
         * Direct offset calculation with sizeof(uint32) is incorrect due to
         * struct alignment (4-byte padding between uint32 and uint64).
         */
        KmerFreqEntry32 *freq_entry = (KmerFreqEntry32 *)batch_entry;

        entry = (KmerFreqEntry32 *) hash_search(merge_htab, &freq_entry->uintkey, HASH_ENTER, &found);
        if (found)
            entry->appearance_nrow += freq_entry->appearance_nrow;
        else
            entry->appearance_nrow = freq_entry->appearance_nrow;
    }

    saved_bucket_count = ctx->bucket_count;
    strlcpy(saved_path, ctx->path, MAXPGPATH);

    if (ftruncate(ctx->fd, 0) < 0)
    {
        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(merge_context);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not truncate file hash table: %m")));
    }
    lseek(ctx->fd, 0, SEEK_SET);
    kmersearch_fht32_close(ctx);

    new_ctx = kmersearch_fht32_create(saved_path, saved_bucket_count);

    hash_seq_init(&hash_seq, merge_htab);
    while ((entry = (KmerFreqEntry32 *) hash_seq_search(&hash_seq)) != NULL)
    {
        kmersearch_fht32_add(new_ctx, entry->uintkey, entry->appearance_nrow);
    }

    kmersearch_fht32_close(new_ctx);

    MemoryContextSwitchTo(old_context);
    MemoryContextDelete(merge_context);
}

/*
 * Merge source file into target file (target += source)
 * Uses in-memory merge when sufficient memory is available,
 * falls back to entry-by-entry processing otherwise.
 */
void
kmersearch_fht32_merge(const char *source_path, const char *target_path)
{
    FileHashTable32Context *source_ctx;
    FileHashTable32Context *target_ctx;
    FileHashTableIterator32 iter;
    uint32 uintkey;
    uint64 appearance_nrow;
    Size memory_required;
    Size memory_available;
    bool use_inmemory_merge;
    HTAB *merge_htab = NULL;
    HASHCTL hash_ctl;
    MemoryContext merge_context = NULL;
    MemoryContext old_context;

    source_ctx = kmersearch_fht32_open(source_path);
    target_ctx = kmersearch_fht32_open(target_path);

    memory_available = (Size)maintenance_work_mem * 1024L;
    memory_required = (source_ctx->entry_count + target_ctx->entry_count) *
                      sizeof(KmerFreqEntry32) * FHT_HTAB_OVERHEAD_FACTOR;

    use_inmemory_merge = (memory_required < memory_available / 2);

    if (use_inmemory_merge)
    {
        long max_entries;
        KmerFreqEntry32 *entry;
        bool found;
        HASH_SEQ_STATUS hash_seq;

        merge_context = AllocSetContextCreate(CurrentMemoryContext,
                                              "FHT32 Merge Context",
                                              ALLOCSET_DEFAULT_SIZES);
        old_context = MemoryContextSwitchTo(merge_context);

        max_entries = source_ctx->entry_count + target_ctx->entry_count;
        if (max_entries < 1024)
            max_entries = 1024;

        memset(&hash_ctl, 0, sizeof(hash_ctl));
        hash_ctl.keysize = sizeof(uint32);
        hash_ctl.entrysize = sizeof(KmerFreqEntry32);
        hash_ctl.hcxt = merge_context;

        merge_htab = hash_create("FHT32 Merge Hash",
                                 max_entries,
                                 &hash_ctl,
                                 HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);

        kmersearch_fht32_iterator_init(&iter, target_ctx);
        while (kmersearch_fht32_iterate(&iter, &uintkey, &appearance_nrow))
        {
            entry = (KmerFreqEntry32 *) hash_search(merge_htab, &uintkey,
                                                    HASH_ENTER, &found);
            if (found)
                entry->appearance_nrow += appearance_nrow;
            else
                entry->appearance_nrow = appearance_nrow;
        }

        kmersearch_fht32_iterator_init(&iter, source_ctx);
        while (kmersearch_fht32_iterate(&iter, &uintkey, &appearance_nrow))
        {
            entry = (KmerFreqEntry32 *) hash_search(merge_htab, &uintkey,
                                                    HASH_ENTER, &found);
            if (found)
                entry->appearance_nrow += appearance_nrow;
            else
                entry->appearance_nrow = appearance_nrow;
        }

        {
            uint32 saved_bucket_count = target_ctx->bucket_count;
            int truncate_fd;

            kmersearch_fht32_close(source_ctx);
            kmersearch_fht32_close(target_ctx);

            truncate_fd = open(target_path, O_RDWR | O_TRUNC);
            if (truncate_fd < 0)
            {
                MemoryContextSwitchTo(old_context);
                MemoryContextDelete(merge_context);
                ereport(ERROR,
                        (errcode_for_file_access(),
                         errmsg("could not truncate target file: %m")));
            }
            close(truncate_fd);

            target_ctx = kmersearch_fht32_create(target_path, saved_bucket_count);
        }

        hash_seq_init(&hash_seq, merge_htab);
        while ((entry = (KmerFreqEntry32 *) hash_seq_search(&hash_seq)) != NULL)
        {
            kmersearch_fht32_add(target_ctx, entry->uintkey, entry->appearance_nrow);
        }

        kmersearch_fht32_close(target_ctx);

        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(merge_context);
    }
    else
    {
        kmersearch_fht32_iterator_init(&iter, source_ctx);
        while (kmersearch_fht32_iterate(&iter, &uintkey, &appearance_nrow))
        {
            kmersearch_fht32_add(target_ctx, uintkey, appearance_nrow);
        }

        kmersearch_fht32_close(source_ctx);
        kmersearch_fht32_close(target_ctx);
    }

    unlink(source_path);
}

/*
 * Initialize iterator for uint32 file hash table
 */
void
kmersearch_fht32_iterator_init(FileHashTableIterator32 *iter, FileHashTable32Context *ctx)
{
    iter->ctx = ctx;
    iter->current_bucket = 0;
    iter->current_offset = 0;
}

/*
 * Get next entry from uint32 file hash table iterator
 */
bool
kmersearch_fht32_iterate(FileHashTableIterator32 *iter, uint32 *uintkey, uint64 *appearance_nrow)
{
    FileHashTable32Entry entry;
    ssize_t bytes;

    /* If we have a current chain to follow */
    if (iter->current_offset != 0)
    {
        bytes = pread(iter->ctx->fd, &entry, sizeof(entry), iter->current_offset);
        if (bytes == sizeof(entry))
        {
            *uintkey = entry.uintkey;
            *appearance_nrow = entry.appearance_nrow;
            iter->current_offset = entry.next_offset;
            return true;
        }
        iter->current_offset = 0;
    }

    /* Find next non-empty bucket */
    while (iter->current_bucket < iter->ctx->bucket_count)
    {
        off_t bucket_offset;
        uint64 entry_offset;

        bucket_offset = FHT32_HEADER_SIZE + (off_t)iter->current_bucket * sizeof(uint64);
        bytes = pread(iter->ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);

        iter->current_bucket++;

        if (bytes == sizeof(uint64) && entry_offset != 0)
        {
            bytes = pread(iter->ctx->fd, &entry, sizeof(entry), entry_offset);
            if (bytes == sizeof(entry))
            {
                *uintkey = entry.uintkey;
                *appearance_nrow = entry.appearance_nrow;
                iter->current_offset = entry.next_offset;
                return true;
            }
        }
    }

    return false;
}

/*
 * ============================================================================
 * uint64 Hash Table Implementation (FHT64)
 * ============================================================================
 */

/*
 * Create a new uint64 file hash table
 */
FileHashTable64Context *
kmersearch_fht64_create(const char *path, uint32 bucket_count)
{
    FileHashTable64Context *ctx;
    FileHashTable64Header header;
    uint64 *bucket_array;
    off_t data_start;
    ssize_t written;

    if (bucket_count == 0)
        bucket_count = kmersearch_calculate_bucket_count(kmersearch_highfreq_analysis_hashtable_size);

    ctx = palloc0(sizeof(FileHashTable64Context));
    strlcpy(ctx->path, path, MAXPGPATH);
    ctx->bucket_count = bucket_count;

    /* Create and open file */
    ctx->fd = open(path, O_RDWR | O_CREAT | O_TRUNC, 0600);
    if (ctx->fd < 0)
    {
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not create file hash table \"%s\": %m", path)));
    }

    /* Initialize header */
    memset(&header, 0, sizeof(header));
    header.magic = FHT64_MAGIC;
    header.version = FHT_VERSION;
    header.key_type = 64;
    header.bucket_count = bucket_count;
    header.entry_count = 0;

    /* Calculate data start position (after header and bucket directory) */
    data_start = FHT64_HEADER_SIZE + (off_t)bucket_count * sizeof(uint64);
    header.next_entry_offset = data_start;
    ctx->next_entry_offset = data_start;

    /* Write header */
    written = write(ctx->fd, &header, sizeof(header));
    if (written != sizeof(header))
    {
        close(ctx->fd);
        unlink(path);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write file hash table header: %m")));
    }

    /* Initialize bucket directory with zeros (0 = empty bucket) */
    bucket_array = palloc0(bucket_count * sizeof(uint64));
    written = write(ctx->fd, bucket_array, bucket_count * sizeof(uint64));
    pfree(bucket_array);

    if (written != (ssize_t)(bucket_count * sizeof(uint64)))
    {
        close(ctx->fd);
        unlink(path);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not initialize bucket directory: %m")));
    }

    ctx->entry_count = 0;

    return ctx;
}

/*
 * Open an existing uint64 file hash table
 */
FileHashTable64Context *
kmersearch_fht64_open(const char *path)
{
    FileHashTable64Context *ctx;
    FileHashTable64Header header;
    ssize_t bytes_read;

    ctx = palloc0(sizeof(FileHashTable64Context));
    strlcpy(ctx->path, path, MAXPGPATH);

    ctx->fd = open(path, O_RDWR);
    if (ctx->fd < 0)
    {
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not open file hash table \"%s\": %m", path)));
    }

    /* Read and validate header */
    bytes_read = pread(ctx->fd, &header, sizeof(header), 0);
    if (bytes_read != sizeof(header))
    {
        close(ctx->fd);
        pfree(ctx);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read file hash table header: %m")));
    }

    if (header.magic != FHT64_MAGIC || header.version != FHT_VERSION)
    {
        close(ctx->fd);
        pfree(ctx);
        ereport(ERROR,
                (errcode(ERRCODE_DATA_CORRUPTED),
                 errmsg("invalid file hash table format")));
    }

    ctx->bucket_count = header.bucket_count;
    ctx->entry_count = header.entry_count;
    ctx->next_entry_offset = header.next_entry_offset;

    return ctx;
}

/*
 * Close a uint64 file hash table
 */
void
kmersearch_fht64_close(FileHashTable64Context *ctx)
{
    if (ctx == NULL)
        return;

    if (ctx->fd >= 0)
    {
        kmersearch_fht64_flush(ctx);
        close(ctx->fd);
    }

    pfree(ctx);
}

/*
 * Add appearance_nrow to a uint64 key
 */
void
kmersearch_fht64_add(FileHashTable64Context *ctx, uint64 uintkey, uint64 appearance_nrow)
{
    uint32 bucket_idx;
    off_t bucket_offset;
    uint64 entry_offset;
    FileHashTable64Entry entry;
    ssize_t bytes;

    bucket_idx = (uint32)(kmersearch_murmurhash64(uintkey) & (ctx->bucket_count - 1));
    bucket_offset = FHT64_HEADER_SIZE + (off_t)bucket_idx * sizeof(uint64);

    /* Read bucket head offset */
    bytes = pread(ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read bucket: %m")));
    }

    /* Search for existing entry in chain */
    while (entry_offset != 0)
    {
        bytes = pread(ctx->fd, &entry, sizeof(entry), entry_offset);
        if (bytes != sizeof(entry))
        {
            ereport(ERROR,
                    (errcode_for_file_access(),
                     errmsg("could not read entry: %m")));
        }

        if (entry.uintkey == uintkey)
        {
            /* Found existing entry, add to its value */
            entry.appearance_nrow += appearance_nrow;
            bytes = pwrite(ctx->fd, &entry, sizeof(entry), entry_offset);
            if (bytes != sizeof(entry))
            {
                ereport(ERROR,
                        (errcode_for_file_access(),
                         errmsg("could not update entry: %m")));
            }
            return;
        }

        entry_offset = entry.next_offset;
    }

    /* Entry not found, create new one */
    entry.uintkey = uintkey;
    entry.appearance_nrow = appearance_nrow;

    /* Read current bucket head to set as next */
    bytes = pread(ctx->fd, &entry.next_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not read bucket head: %m")));
    }

    /* Write new entry at next available position */
    entry_offset = ctx->next_entry_offset;
    bytes = pwrite(ctx->fd, &entry, sizeof(entry), entry_offset);
    if (bytes != sizeof(entry))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not write new entry: %m")));
    }

    /* Update bucket head to point to new entry */
    bytes = pwrite(ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not update bucket head: %m")));
    }

    ctx->next_entry_offset += sizeof(entry);
    ctx->entry_count++;
}

/*
 * Get appearance_nrow for a uint64 key
 */
uint64
kmersearch_fht64_get(FileHashTable64Context *ctx, uint64 uintkey)
{
    uint32 bucket_idx;
    off_t bucket_offset;
    uint64 entry_offset;
    FileHashTable64Entry entry;
    ssize_t bytes;

    bucket_idx = (uint32)(kmersearch_murmurhash64(uintkey) & (ctx->bucket_count - 1));
    bucket_offset = FHT64_HEADER_SIZE + (off_t)bucket_idx * sizeof(uint64);

    /* Read bucket head offset */
    bytes = pread(ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);
    if (bytes != sizeof(uint64))
    {
        return 0;
    }

    /* Search chain for entry */
    while (entry_offset != 0)
    {
        bytes = pread(ctx->fd, &entry, sizeof(entry), entry_offset);
        if (bytes != sizeof(entry))
        {
            return 0;
        }

        if (entry.uintkey == uintkey)
        {
            return entry.appearance_nrow;
        }

        entry_offset = entry.next_offset;
    }

    return 0;
}

/*
 * Flush uint64 file hash table
 */
void
kmersearch_fht64_flush(FileHashTable64Context *ctx)
{
    FileHashTable64Header header;
    ssize_t bytes;

    if (ctx == NULL || ctx->fd < 0)
        return;

    /* Read current header */
    bytes = pread(ctx->fd, &header, sizeof(header), 0);
    if (bytes == sizeof(header))
    {
        header.entry_count = ctx->entry_count;
        header.next_entry_offset = ctx->next_entry_offset;
        if (pwrite(ctx->fd, &header, sizeof(header), 0) != sizeof(header))
            ereport(WARNING,
                    (errcode_for_file_access(),
                     errmsg("could not update file hash table header: %m")));
    }

    fsync(ctx->fd);
}

/*
 * Bulk add from batch hash table to file hash table
 * Reads FHT file into memory, merges with batch_hash, writes back
 */
void
kmersearch_fht64_bulk_add(FileHashTable64Context *ctx, HTAB *batch_hash)
{
    HTAB *merge_htab;
    HASHCTL hash_ctl;
    MemoryContext merge_context;
    MemoryContext old_context;
    FileHashTableIterator64 iter;
    uint64 uintkey;
    uint64 appearance_nrow;
    KmerFreqEntry64 *entry;
    bool found;
    HASH_SEQ_STATUS hash_seq;
    HASH_SEQ_STATUS batch_seq;
    void *batch_entry;
    uint32 saved_bucket_count;
    char saved_path[MAXPGPATH];
    long max_entries;
    FileHashTable64Context *new_ctx;

    merge_context = AllocSetContextCreate(CurrentMemoryContext,
                                          "FHT64 BulkAdd Context",
                                          ALLOCSET_DEFAULT_SIZES);
    old_context = MemoryContextSwitchTo(merge_context);

    max_entries = ctx->entry_count + hash_get_num_entries(batch_hash);
    if (max_entries < 1024)
        max_entries = 1024;

    memset(&hash_ctl, 0, sizeof(hash_ctl));
    hash_ctl.keysize = sizeof(uint64);
    hash_ctl.entrysize = sizeof(KmerFreqEntry64);
    hash_ctl.hcxt = merge_context;

    merge_htab = hash_create("FHT64 BulkAdd Hash",
                             max_entries,
                             &hash_ctl,
                             HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);

    kmersearch_fht64_iterator_init(&iter, ctx);
    while (kmersearch_fht64_iterate(&iter, &uintkey, &appearance_nrow))
    {
        entry = (KmerFreqEntry64 *) hash_search(merge_htab, &uintkey, HASH_ENTER, &found);
        if (found)
            entry->appearance_nrow += appearance_nrow;
        else
            entry->appearance_nrow = appearance_nrow;
    }

    hash_seq_init(&batch_seq, batch_hash);
    while ((batch_entry = hash_seq_search(&batch_seq)) != NULL)
    {
        /*
         * Cast to KmerFreqEntry64 which has the same memory layout as
         * TempKmerFreqEntry64 for the first two fields (uint64 + uint64).
         * While TempKmerFreqEntry64 happens to have no padding, using struct
         * member access is more consistent and maintainable.
         */
        KmerFreqEntry64 *freq_entry = (KmerFreqEntry64 *)batch_entry;

        entry = (KmerFreqEntry64 *) hash_search(merge_htab, &freq_entry->uintkey, HASH_ENTER, &found);
        if (found)
            entry->appearance_nrow += freq_entry->appearance_nrow;
        else
            entry->appearance_nrow = freq_entry->appearance_nrow;
    }

    saved_bucket_count = ctx->bucket_count;
    strlcpy(saved_path, ctx->path, MAXPGPATH);

    if (ftruncate(ctx->fd, 0) < 0)
    {
        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(merge_context);
        ereport(ERROR,
                (errcode_for_file_access(),
                 errmsg("could not truncate file hash table: %m")));
    }
    lseek(ctx->fd, 0, SEEK_SET);
    kmersearch_fht64_close(ctx);

    new_ctx = kmersearch_fht64_create(saved_path, saved_bucket_count);

    hash_seq_init(&hash_seq, merge_htab);
    while ((entry = (KmerFreqEntry64 *) hash_seq_search(&hash_seq)) != NULL)
    {
        kmersearch_fht64_add(new_ctx, entry->uintkey, entry->appearance_nrow);
    }

    kmersearch_fht64_close(new_ctx);

    MemoryContextSwitchTo(old_context);
    MemoryContextDelete(merge_context);
}

/*
 * Merge source file into target file (target += source)
 * Uses in-memory merge when sufficient memory is available,
 * falls back to entry-by-entry processing otherwise.
 */
void
kmersearch_fht64_merge(const char *source_path, const char *target_path)
{
    FileHashTable64Context *source_ctx;
    FileHashTable64Context *target_ctx;
    FileHashTableIterator64 iter;
    uint64 uintkey;
    uint64 appearance_nrow;
    Size memory_required;
    Size memory_available;
    bool use_inmemory_merge;
    HTAB *merge_htab = NULL;
    HASHCTL hash_ctl;
    MemoryContext merge_context = NULL;
    MemoryContext old_context;

    source_ctx = kmersearch_fht64_open(source_path);
    target_ctx = kmersearch_fht64_open(target_path);

    memory_available = (Size)maintenance_work_mem * 1024L;
    memory_required = (source_ctx->entry_count + target_ctx->entry_count) *
                      sizeof(KmerFreqEntry64) * FHT_HTAB_OVERHEAD_FACTOR;

    use_inmemory_merge = (memory_required < memory_available / 2);

    if (use_inmemory_merge)
    {
        long max_entries;
        KmerFreqEntry64 *entry;
        bool found;
        HASH_SEQ_STATUS hash_seq;

        merge_context = AllocSetContextCreate(CurrentMemoryContext,
                                              "FHT64 Merge Context",
                                              ALLOCSET_DEFAULT_SIZES);
        old_context = MemoryContextSwitchTo(merge_context);

        max_entries = source_ctx->entry_count + target_ctx->entry_count;
        if (max_entries < 1024)
            max_entries = 1024;

        memset(&hash_ctl, 0, sizeof(hash_ctl));
        hash_ctl.keysize = sizeof(uint64);
        hash_ctl.entrysize = sizeof(KmerFreqEntry64);
        hash_ctl.hcxt = merge_context;

        merge_htab = hash_create("FHT64 Merge Hash",
                                 max_entries,
                                 &hash_ctl,
                                 HASH_ELEM | HASH_BLOBS | HASH_CONTEXT);

        kmersearch_fht64_iterator_init(&iter, target_ctx);
        while (kmersearch_fht64_iterate(&iter, &uintkey, &appearance_nrow))
        {
            entry = (KmerFreqEntry64 *) hash_search(merge_htab, &uintkey,
                                                    HASH_ENTER, &found);
            if (found)
                entry->appearance_nrow += appearance_nrow;
            else
                entry->appearance_nrow = appearance_nrow;
        }

        kmersearch_fht64_iterator_init(&iter, source_ctx);
        while (kmersearch_fht64_iterate(&iter, &uintkey, &appearance_nrow))
        {
            entry = (KmerFreqEntry64 *) hash_search(merge_htab, &uintkey,
                                                    HASH_ENTER, &found);
            if (found)
                entry->appearance_nrow += appearance_nrow;
            else
                entry->appearance_nrow = appearance_nrow;
        }

        {
            uint32 saved_bucket_count = target_ctx->bucket_count;
            int truncate_fd;

            kmersearch_fht64_close(source_ctx);
            kmersearch_fht64_close(target_ctx);

            truncate_fd = open(target_path, O_RDWR | O_TRUNC);
            if (truncate_fd < 0)
            {
                MemoryContextSwitchTo(old_context);
                MemoryContextDelete(merge_context);
                ereport(ERROR,
                        (errcode_for_file_access(),
                         errmsg("could not truncate target file: %m")));
            }
            close(truncate_fd);

            target_ctx = kmersearch_fht64_create(target_path, saved_bucket_count);
        }

        hash_seq_init(&hash_seq, merge_htab);
        while ((entry = (KmerFreqEntry64 *) hash_seq_search(&hash_seq)) != NULL)
        {
            kmersearch_fht64_add(target_ctx, entry->uintkey, entry->appearance_nrow);
        }

        kmersearch_fht64_close(target_ctx);

        MemoryContextSwitchTo(old_context);
        MemoryContextDelete(merge_context);
    }
    else
    {
        kmersearch_fht64_iterator_init(&iter, source_ctx);
        while (kmersearch_fht64_iterate(&iter, &uintkey, &appearance_nrow))
        {
            kmersearch_fht64_add(target_ctx, uintkey, appearance_nrow);
        }

        kmersearch_fht64_close(source_ctx);
        kmersearch_fht64_close(target_ctx);
    }

    unlink(source_path);
}

/*
 * Initialize iterator for uint64 file hash table
 */
void
kmersearch_fht64_iterator_init(FileHashTableIterator64 *iter, FileHashTable64Context *ctx)
{
    iter->ctx = ctx;
    iter->current_bucket = 0;
    iter->current_offset = 0;
}

/*
 * Get next entry from uint64 file hash table iterator
 */
bool
kmersearch_fht64_iterate(FileHashTableIterator64 *iter, uint64 *uintkey, uint64 *appearance_nrow)
{
    FileHashTable64Entry entry;
    ssize_t bytes;

    /* If we have a current chain to follow */
    if (iter->current_offset != 0)
    {
        bytes = pread(iter->ctx->fd, &entry, sizeof(entry), iter->current_offset);
        if (bytes == sizeof(entry))
        {
            *uintkey = entry.uintkey;
            *appearance_nrow = entry.appearance_nrow;
            iter->current_offset = entry.next_offset;
            return true;
        }
        iter->current_offset = 0;
    }

    /* Find next non-empty bucket */
    while (iter->current_bucket < iter->ctx->bucket_count)
    {
        off_t bucket_offset;
        uint64 entry_offset;

        bucket_offset = FHT64_HEADER_SIZE + (off_t)iter->current_bucket * sizeof(uint64);
        bytes = pread(iter->ctx->fd, &entry_offset, sizeof(uint64), bucket_offset);

        iter->current_bucket++;

        if (bytes == sizeof(uint64) && entry_offset != 0)
        {
            bytes = pread(iter->ctx->fd, &entry, sizeof(entry), entry_offset);
            if (bytes == sizeof(entry))
            {
                *uintkey = entry.uintkey;
                *appearance_nrow = entry.appearance_nrow;
                iter->current_offset = entry.next_offset;
                return true;
            }
        }
    }

    return false;
}
