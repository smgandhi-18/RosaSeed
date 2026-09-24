/*************************************************************************************
                           The MIT License

   RosaSeed (RosaSeed: Faster and Accurate Short Read Alignment Using a Configurable Seeding Strategy),
   Copyright (C) 2026  University of Alberta, Gandhi Shyama.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Contacts: Shyama Gandhi <smgandhi@ualberta.ca>

*****************************************************************************************/

#include "load_data.h"
#include "file_dec.h"
#include "helper_functions.h"
#include <ctype.h>

#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>

uint32_t no_of_jumps, _jtable_string_length;
int64_t sentinel_index = -1;

/* =========================================================================
   mmap-based parallel loader, 1-step version
   Same principle as 2-step: mmap(MAP_POPULATE) 
   ========================================================================= */

/* mmap region registry for cleanup */
typedef struct { void *addr; size_t len; } MmapRegion;
#define MAX_MMAP_REGIONS 8
static MmapRegion g_mmap_regions[MAX_MMAP_REGIONS];
static int        g_n_mmap_regions = 0;

static void register_mmap(void *addr, size_t len) {
    if (g_n_mmap_regions < MAX_MMAP_REGIONS) {
        g_mmap_regions[g_n_mmap_regions].addr = addr;
        g_mmap_regions[g_n_mmap_regions].len  = len;
        g_n_mmap_regions++;
    }
}

void unmap_all_index_regions(void) {
    for (int i = 0; i < g_n_mmap_regions; ++i) {
        if (g_mmap_regions[i].addr && g_mmap_regions[i].addr != MAP_FAILED)
            munmap(g_mmap_regions[i].addr, g_mmap_regions[i].len);
        g_mmap_regions[i].addr = NULL;
        g_mmap_regions[i].len  = 0;
    }
    g_n_mmap_regions = 0;
}

/* -------------------------------------------------------------------------
   mmap_file_populate(): map a file region directly into address space.
   ------------------------------------------------------------------------- */
static void *mmap_file_populate(const char *path,
                                 size_t      skip_bytes,
                                 size_t      data_bytes,
                                 const char *label,
                                 double     *elapsed_out)
{
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int fd = open(path, O_RDONLY);
    if (fd < 0) { perror(path); exit(EXIT_FAILURE); }

    long   page_sz  = sysconf(_SC_PAGESIZE);
    size_t align_off = skip_bytes & ~((size_t)(page_sz - 1));
    size_t extra     = skip_bytes - align_off;
    size_t map_len   = data_bytes + extra;

    void *map = mmap(NULL, map_len,
                     PROT_READ,
                     MAP_SHARED | MAP_POPULATE,
                     fd, (off_t)align_off);
    close(fd);

    if (map == MAP_FAILED) {
        perror("mmap"); exit(EXIT_FAILURE);
    }

    madvise(map, map_len, MADV_SEQUENTIAL | MADV_WILLNEED);
    register_mmap(map, map_len);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) +
                     (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    if (elapsed_out) *elapsed_out = elapsed;

    fprintf(stderr, "[LOAD] %-28s  %6.2f GiB  %5.2f s  (%4.1f GiB/s)\n",
            label,
            (double)data_bytes / (1024.0*1024.0*1024.0),
            elapsed,
            elapsed > 0 ? (double)data_bytes/(1024.0*1024.0*1024.0)/elapsed : 0.0);

    return (uint8_t *)map + extra;
}

/* Parallel mmap worker */
typedef struct {
    const char *path;
    size_t      skip_bytes;
    size_t      data_bytes;
    const char *label;
    void      **dest_ptr;
    int         done;
    double      elapsed_sec;
} MmapTask;

static void *mmap_worker(void *arg)
{
    MmapTask *t = (MmapTask *)arg;
    double elapsed = 0.0;
    void *addr = mmap_file_populate(t->path, t->skip_bytes, t->data_bytes,
                                    t->label, &elapsed);
    *(t->dest_ptr) = addr;
    t->elapsed_sec = elapsed;
    t->done = 1;
    return NULL;
}

#define RS1_REF4_MAGIC 0x5253315245463401ULL
typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t bases_per_byte;
    uint64_t n_bases;
    uint64_t packed_bytes;
} rs1_ref4_header_t;

#define RS1_OCC_MAGIC 0x5253314F43433332ULL
typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t occ_interval;
    uint32_t alphabet_size;
    uint32_t reserved;
    uint64_t bwt_len_no_term;
    uint64_t bwt_len_total;
    uint64_t num_blocks;
    int64_t  sentinel_index;
    uint64_t c_array[4];
} rs1_occ_header_t;

/* -------------------------------------------------------------------------
   rosaseed_compact_read_index_metadata()
   Reads the cp_occ_compact.bin header before any allocation and sets the
   genome-dependent sizes, so the aligner works for any genome the index
   builder accepts instead of a compiled-in one.
   ------------------------------------------------------------------------- */
void rosaseed_compact_read_index_metadata(void)
{
    FILE *f = fopen(OCC_BIN_FILE, "rb");
    if (!f) { perror(OCC_BIN_FILE); exit(EXIT_FAILURE); }
    rs1_occ_header_t hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1) {
        fprintf(stderr, "Failed to read cp_occ header: %s\n", OCC_BIN_FILE);
        fclose(f); exit(EXIT_FAILURE);
    }
    if (fseeko(f, 0, SEEK_END) != 0) { perror(OCC_BIN_FILE); fclose(f); exit(EXIT_FAILURE); }
    off_t file_size = ftello(f);
    fclose(f);

    if (hdr.magic != RS1_OCC_MAGIC) {
        fprintf(stderr, "Error: %s is not a RosaSeed-Compact cp_occ file (bad magic).\n"
                        "       Build the index with index-builder/build_compact_pipeline.sh.\n",
                OCC_BIN_FILE);
        exit(EXIT_FAILURE);
    }
    if (hdr.occ_interval != OCC_INTERVAL || hdr.alphabet_size != ALPHABET_SIZE) {
        fprintf(stderr, "Error: %s has occ_interval=%u alphabet=%u; this binary expects %d / %d.\n",
                OCC_BIN_FILE, hdr.occ_interval, hdr.alphabet_size, OCC_INTERVAL, ALPHABET_SIZE);
        exit(EXIT_FAILURE);
    }
    const uint64_t N = hdr.bwt_len_total;
    if (N < 3 || (N & 1ULL) == 0 || hdr.bwt_len_no_term + 1 != N) {
        fprintf(stderr, "Error: %s: BWT length %llu is not forward + reverse complement + terminator.\n",
                OCC_BIN_FILE, (unsigned long long)N);
        exit(EXIT_FAILURE);
    }
    if (N > RS_MAX_BWT_LEN) {
        fprintf(stderr, "\n[FATAL] BWT length %llu exceeds the 2^33 limit of the index format\n"
                        "        (reference genomes up to ~4.29 Gbp). See README, 'Genome size limits'.\n",
                (unsigned long long)N);
        exit(EXIT_FAILURE);
    }
    const uint64_t expected_blocks = N / OCC_INTERVAL + 1;
    if (hdr.num_blocks != expected_blocks ||
        (uint64_t)file_size != (uint64_t)sizeof(hdr) + hdr.num_blocks * (uint64_t)sizeof(cp_occ32_t)) {
        fprintf(stderr, "Error: %s has %llu blocks / %lld bytes; expected %llu blocks / %llu bytes.\n"
                        "       The index may be truncated. Regenerate it with build_compact_pipeline.sh.\n",
                OCC_BIN_FILE, (unsigned long long)hdr.num_blocks, (long long)file_size,
                (unsigned long long)expected_blocks,
                (unsigned long long)(sizeof(hdr) + expected_blocks * sizeof(cp_occ32_t)));
        exit(EXIT_FAILURE);
    }

    BWT_SIZE_REFERENCE_SIZE = N;
    rosaseed_L = (N - 1) / 2;
    fprintf(stderr, "[INFO] RosaSeed-Compact index: N=%llu  genome L=%llu bp\n",
            (unsigned long long)N, (unsigned long long)rosaseed_L);
}

static void load_cp_occ_1step_mmap(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(EXIT_FAILURE); }

    rs1_occ_header_t hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1) {
        fprintf(stderr, "Failed to read 1-step cp_occ header: %s\n", path);
        exit(EXIT_FAILURE);
    }
    fclose(f);

    if (hdr.magic != RS1_OCC_MAGIC) {
        fprintf(stderr, "Bad 1-step cp_occ magic in %s\n", path);
        exit(EXIT_FAILURE);
    }
    if (hdr.occ_interval != OCC_INTERVAL) {
        fprintf(stderr, "cp_occ OCC_INTERVAL mismatch: file=%u code=%d\n",
                hdr.occ_interval, OCC_INTERVAL);
        exit(EXIT_FAILURE);
    }
    if (hdr.alphabet_size != ALPHABET_SIZE) {
        fprintf(stderr, "cp_occ alphabet mismatch: file=%u code=%d\n",
                hdr.alphabet_size, ALPHABET_SIZE);
        exit(EXIT_FAILURE);
    }
    if (hdr.bwt_len_total != BWT_SIZE_REFERENCE_SIZE) {
        fprintf(stderr, "WARNING: BWT length mismatch: file=%llu code=%llu\n",
                (unsigned long long)hdr.bwt_len_total,
                (unsigned long long)BWT_SIZE_REFERENCE_SIZE);
    }
    uint64_t expected_blocks = BWT_SIZE_REFERENCE_SIZE / OCC_INTERVAL + 1;
    if (hdr.num_blocks != expected_blocks) {
        fprintf(stderr, "cp_occ block count mismatch: file=%llu expected=%llu\n",
                (unsigned long long)hdr.num_blocks,
                (unsigned long long)expected_blocks);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < ALPHABET_SIZE; ++i)
        c_vec[i] = hdr.c_array[i];
    c_vec[ALPHABET_SIZE] = hdr.bwt_len_total;
    sentinel_index = hdr.sentinel_index;

    /* mmap the body, skipping the header */
    size_t body_bytes = hdr.num_blocks * sizeof(cp_occ32_t);
    void *mapped = mmap_file_populate(path,
                                      sizeof(hdr),
                                      body_bytes,
                                      "cp_occ (mmap+MAP_POPULATE)",
                                      NULL);

    /* Free the _mm_malloc'd buffer, redirect pointer to mmap region */
    if (cp_occ) { _mm_free(cp_occ); cp_occ = NULL; }
    cp_occ = (cp_occ32_t *)mapped;

    fprintf(stderr, "1-step cp_occ loaded: %llu blocks\n",
            (unsigned long long)hdr.num_blocks);
    fprintf(stderr, "BWT length total    : %llu\n",
            (unsigned long long)hdr.bwt_len_total);
    fprintf(stderr, "BWT length no-term  : %llu\n",
            (unsigned long long)hdr.bwt_len_no_term);
    fprintf(stderr, "Sentinel index      : %ld\n", sentinel_index);
    fprintf(stderr, "C[A]=%llu C[C]=%llu C[G]=%llu C[T]=%llu\n",
            (unsigned long long)c_vec[0], (unsigned long long)c_vec[1],
            (unsigned long long)c_vec[2], (unsigned long long)c_vec[3]);
}

/* -------------------------------------------------------------------------
   load_bwt_data_structures(), 1-step mmap version

   Thread 1: mmap sa_ls_word
   Thread 2: mmap sa_ms_byte
   Thread 3: mmap jump_pointers
   Thread 4: mmap ref4_packed / ref_4 (skip 32-byte header)
   Main:     mmap cp_occ body (header parsed first → fills c_vec+sentinel)

   All five mmap(MAP_POPULATE) calls run simultaneously 
   ------------------------------------------------------------------------- */
void load_bwt_data_structures(void)
{
    struct timespec wall_t0, wall_t1;
    clock_gettime(CLOCK_MONOTONIC, &wall_t0);

    fprintf(stderr, "\n[LOAD] Starting parallel 1-step mmap index load (MAP_POPULATE)...\n");

    /* Read ref4 header to get exact packed byte count */
    size_t ref_data_bytes = 0;
    {
        FILE *fp = fopen(REFERENCE_BIN_FILE, "rb");
        if (!fp) { perror(REFERENCE_BIN_FILE); exit(EXIT_FAILURE); }

        rs1_ref4_header_t hdr;
        if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
            fprintf(stderr, "ERROR: failed to read ref4 header\n");
            exit(EXIT_FAILURE);
        }
        fclose(fp);

        if (hdr.magic != RS1_REF4_MAGIC) {
            fprintf(stderr, "ERROR: bad ref4 binary magic\n"); exit(EXIT_FAILURE);
        }
        if (hdr.version != 1 || hdr.bases_per_byte != 4) {
            fprintf(stderr, "ERROR: unsupported ref4 format\n"); exit(EXIT_FAILURE);
        }
        uint64_t expected_bases = 2ULL * rosaseed_L;
        uint64_t expected_bytes = (expected_bases + 3ULL) >> 2;
        if (hdr.n_bases != expected_bases || hdr.packed_bytes != expected_bytes) {
            fprintf(stderr,
                "ERROR: ref4 size mismatch: n_bases=%llu/%llu packed=%llu/%llu\n",
                (unsigned long long)hdr.n_bases, (unsigned long long)expected_bases,
                (unsigned long long)hdr.packed_bytes, (unsigned long long)expected_bytes);
            exit(EXIT_FAILURE);
        }
        ref_data_bytes = (size_t)hdr.packed_bytes;
        fprintf(stderr, "[LOAD] ref4: n_bases=%llu packed_bytes=%zu (%.3f GiB)\n",
                (unsigned long long)hdr.n_bases, ref_data_bytes,
                (double)ref_data_bytes / (1024.0*1024.0*1024.0));
    }

    const size_t sa_n = (size_t)SA_SAMPLED_SIZE;

    /* Set up mmap tasks */
    MmapTask task_sa_ls = {
        .path       = SA_LS_WORD_BIN_FILE,
        .skip_bytes = 0,
        .data_bytes = sa_n * sizeof(uint32_t),
        .label      = "sa_ls_word (mmap)",
        .dest_ptr   = (void **)&sa_ls_word,
        .done       = 0,
    };
    MmapTask task_sa_ms = {
        .path       = SA_MSB_BIT_BIN_FILE,
        .skip_bytes = 0,
        .data_bytes = sa_n * sizeof(uint8_t),
        .label      = "sa_ms_byte (mmap)",
        .dest_ptr   = (void **)&sa_ms_byte,
        .done       = 0,
    };
    MmapTask task_jt = {
        .path       = JUMP_TABLE_FILE,
        .skip_bytes = 0,
        .data_bytes = (size_t)JUMP_TABLE_ENTRIES * sizeof(uint64_t),
        .label      = "jump_pointers (mmap)",
        .dest_ptr   = (void **)&jump_pointers,
        .done       = 0,
    };

#if PACK_REF4
    void **ref_dest_ptr = (void **)&ref4_packed;
    const char *ref_label = "ref4_packed (mmap)";
#else
    void **ref_dest_ptr = (void **)&ref_4;
    const char *ref_label = "ref_4 (mmap)";
#endif

    MmapTask task_ref = {
        .path       = REFERENCE_BIN_FILE,
        .skip_bytes = sizeof(rs1_ref4_header_t),   /* skip 32-byte header */
        .data_bytes = ref_data_bytes,
        .label      = ref_label,
        .dest_ptr   = ref_dest_ptr,
        .done       = 0,
    };

    /* Free posix_memalign/mm_malloc buffers */
    free(sa_ls_word);       sa_ls_word    = NULL;
    free(sa_ms_byte);       sa_ms_byte    = NULL;
    free(jump_pointers);    jump_pointers = NULL;
#if PACK_REF4
    if (ref4_packed) { free(ref4_packed); ref4_packed = NULL; }
#else
    if (ref_4) { _mm_free(ref_4); ref_4 = NULL; }
#endif

    /* Launch 4 mmap threads */
    pthread_t thr_sa_ls, thr_sa_ms, thr_jt, thr_ref;
    pthread_create(&thr_sa_ls, NULL, mmap_worker, &task_sa_ls);
    pthread_create(&thr_sa_ms, NULL, mmap_worker, &task_sa_ms);
    pthread_create(&thr_jt,    NULL, mmap_worker, &task_jt);
    pthread_create(&thr_ref,   NULL, mmap_worker, &task_ref);

    /* Main thread: mmap cp_occ (largest) in parallel with threads.
       Also fills c_vec and sentinel_index from the header. */
    load_cp_occ_1step_mmap(OCC_BIN_FILE);

    /* Join threads */
    pthread_join(thr_sa_ls, NULL);
    pthread_join(thr_sa_ms, NULL);
    pthread_join(thr_jt,    NULL);
    pthread_join(thr_ref,   NULL);

    if (!task_sa_ls.done || !task_sa_ms.done ||
        !task_jt.done    || !task_ref.done) {
        fprintf(stderr, "[ERROR] One or more 1-step mmap threads failed\n");
        exit(EXIT_FAILURE);
    }

    /* SA sentinel scan */
    for (size_t v = 0; v < sa_n; ++v) {
        if (sa_ls_word[v] == 0 && sa_ms_byte[v] == 0) {
            int64_t sentinel_from_sa = (int)((uint64_t)v << SA_COMPRESSION_FACTOR_POWER);
            if (sentinel_from_sa != sentinel_index) {
                fprintf(stderr,
                    "NOTE: sentinel from SA scan (%ld) differs from "
                    "cp_occ header (%ld), using cp_occ value\n",
                    sentinel_from_sa, sentinel_index);
            }
            fprintf(stderr, "SA loaded: %zu sampled entries (CF=%d)\n",
                    sa_n, 1 << SA_COMPRESSION_FACTOR_POWER);
            fprintf(stderr, "Sentinel index from SA scan: %ld\n", sentinel_from_sa);
            break;
        }
    }

    {
        uint64_t unique = 0, empty = 0, multi = 0;
        for (uint64_t i = 0; i < JUMP_TABLE_ENTRIES; i++) {
            uint64_t jp = jump_pointers[i];
            if (jp == 0) empty++;
            else if ((jp >> 63) & 1ULL) unique++;
            else multi++;
        }
        fprintf(stderr,
            "Jump table loaded: unique=%llu multi=%llu empty=%llu\n",
            (unsigned long long)unique,
            (unsigned long long)multi,
            (unsigned long long)empty);
    }

    clock_gettime(CLOCK_MONOTONIC, &wall_t1);
    double wall_elapsed = (wall_t1.tv_sec  - wall_t0.tv_sec) +
                          (wall_t1.tv_nsec - wall_t0.tv_nsec) * 1e-9;

    fprintf(stderr,
        "\n[LOAD] All 1-step structures mmap'd in %.2f sec wall time\n"
        "       (second run will be faster, pages remain in page cache)\n\n",
        wall_elapsed);
}

static void load_sa_binary(const char *path_ls, const char *path_msb)
{
    FILE *f32 = fopen(path_ls, "rb");
    FILE *f8  = fopen(path_msb, "rb");
    if (!f32) { perror(path_ls);  exit(EXIT_FAILURE); }
    if (!f8)  { perror(path_msb); exit(EXIT_FAILURE); }
    const size_t sa_n = (size_t)SA_SAMPLED_SIZE;
    size_t r32 = fread(sa_ls_word, sizeof(uint32_t), sa_n, f32);
    size_t r8  = fread(sa_ms_byte, sizeof(uint8_t),  sa_n, f8);
    if (r32 != sa_n || r8 != sa_n) {
        fprintf(stderr,
            "SA binary read mismatch: expected %zu, got ls=%zu ms=%zu\n",
            sa_n, r32, r8);
        exit(EXIT_FAILURE);
    }
    fclose(f32); fclose(f8);
    for (size_t v = 0; v < sa_n; v++) {
        if (sa_ls_word[v] == 0 && sa_ms_byte[v] == 0) {
            sentinel_index = (int64_t)(v << SA_COMPRESSION_FACTOR_POWER);
            break;
        }
    }
    fprintf(stderr, "SA loaded: %zu sampled entries (CF=%d)\n",
            sa_n, 1 << SA_COMPRESSION_FACTOR_POWER);
    fprintf(stderr, "Sentinel index: %ld\n", sentinel_index);
}

static void load_jump_table_binary(const char *path)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) { perror(path); exit(EXIT_FAILURE); }
    size_t got = fread(jump_pointers, sizeof(uint64_t), JUMP_TABLE_ENTRIES, fp);
    if (got != JUMP_TABLE_ENTRIES) {
        fprintf(stderr, "ERROR: JT binary read mismatch\n"); exit(EXIT_FAILURE);
    }
    fclose(fp);
}
