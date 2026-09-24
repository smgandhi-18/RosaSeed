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
#include "read_init.h"

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

void rosaseed_read_index_metadata(const char *index_dir) {
    char path[4096];
    snprintf(path, sizeof(path), "%s/cp_occ_full.bin", index_dir);

    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(EXIT_FAILURE); }

    rs_occ_full_header_t hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1) {
        fprintf(stderr, "Failed to read cp_occ header\n");
        fclose(f); exit(EXIT_FAILURE);
    }

    if (fseeko(f, 0, SEEK_END) != 0) { perror(path); fclose(f); exit(EXIT_FAILURE); }
    off_t file_size = ftello(f);
    fclose(f);
    if (file_size < 0) { perror(path); exit(EXIT_FAILURE); }

    /* ---- validation ---- */
    if (hdr.magic != RS_OCC_MAGIC) {
        fprintf(stderr, "Error: %s is not a RosaSeed cp_occ file (bad magic).\n", path);
        exit(EXIT_FAILURE);
    }
    if (hdr.occ_interval != OCC_INTERVAL || hdr.alphabet_size != ALPHABET_SIZE) {
        fprintf(stderr,
            "Error: index built with occ_interval=%u alphabet=%u, "
            "this binary expects %d / %d.\n",
            hdr.occ_interval, hdr.alphabet_size, OCC_INTERVAL, ALPHABET_SIZE);
        exit(EXIT_FAILURE);
    }
    {
        uint64_t expected_size =
            (uint64_t)sizeof(hdr) + hdr.num_blocks * (uint64_t)sizeof(cp_occ32_t);
        if ((uint64_t)file_size != expected_size) {
            fprintf(stderr,
                "Error: %s is %llu bytes, expected %llu.\n"
                "       The index may be truncated, or built by a different\n"
                "       version of make_cp_occ_2step.  Regenerate with:\n"
                "         make_cp_occ_2step <bwt> cp_occ_full.bin c_vector.txt\n",
                path,
                (unsigned long long)file_size,
                (unsigned long long)expected_size);
            exit(EXIT_FAILURE);
        }
    }
    
    if (hdr.bwt_len_total_with_dollar > RS_MAX_BWT_LEN) {
        fprintf(stderr,
            "\n[FATAL] BWT length %llu exceeds the 2^33 limit of the current index\n"
            "        format (reference genomes up to ~4.29 Gbp).\n"
            "        See README, 'Genome size limits'.\n",
            (unsigned long long)hdr.bwt_len_total_with_dollar);
        exit(EXIT_FAILURE);
    }    

    BWT_SIZE_REFERENCE_SIZE = hdr.bwt_len_total_with_dollar;
    rosaseed_L              = (BWT_SIZE_REFERENCE_SIZE + 1) / 2;

    fprintf(stderr, "[RosaSeed] Index metadata: BWT_SIZE=%llu  L=%llu\n",
            (unsigned long long)BWT_SIZE_REFERENCE_SIZE,
            (unsigned long long)rosaseed_L);
}

/* =========================================================================
   mmap-based parallel loader
   =========================================================================*/

typedef struct {
    void   *addr;
    size_t  len;
} MmapRegion;

#define MAX_MMAP_REGIONS 8
static MmapRegion g_mmap_regions[MAX_MMAP_REGIONS];
static int        g_n_mmap_regions = 0;

static pthread_mutex_t g_mmap_lock = PTHREAD_MUTEX_INITIALIZER;

static void register_mmap(void *addr, size_t len) {
    pthread_mutex_lock(&g_mmap_lock);
    if (g_n_mmap_regions < MAX_MMAP_REGIONS) {
        g_mmap_regions[g_n_mmap_regions].addr = addr;
        g_mmap_regions[g_n_mmap_regions].len  = len;
        g_n_mmap_regions++;
    }
    pthread_mutex_unlock(&g_mmap_lock);
}

/* this instead of free() for mmap'd pointers */
void unmap_all_index_regions(void) {
    for (int i = 0; i < g_n_mmap_regions; ++i) {
        if (g_mmap_regions[i].addr && g_mmap_regions[i].addr != MAP_FAILED)
            munmap(g_mmap_regions[i].addr, g_mmap_regions[i].len);
        g_mmap_regions[i].addr = NULL;
        g_mmap_regions[i].len  = 0;
    }
    g_n_mmap_regions = 0;
}

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

    long page_sz      = sysconf(_SC_PAGESIZE);
    size_t align_off  = skip_bytes & ~((size_t)(page_sz - 1));  /* page-aligned offset */
    size_t extra      = skip_bytes - align_off;                  /* bytes before data */
    size_t map_len    = data_bytes + extra;                      /* total map size */

    void *map = mmap(NULL, map_len,
                     PROT_READ,
                     MAP_SHARED | MAP_POPULATE,   /* MAP_POPULATE = prefault all pages */
                     fd, (off_t)align_off);

    close(fd);   

    if (map == MAP_FAILED) {
        perror("mmap");
        fprintf(stderr, "[ERROR] mmap failed for %s (%zu bytes)\n",
                path, map_len);
        exit(EXIT_FAILURE);
    }

    madvise(map, map_len, MADV_RANDOM | MADV_WILLNEED);

    register_mmap(map, map_len);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) +
                     (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    if (elapsed_out) *elapsed_out = elapsed;

    fprintf(stderr, "[LOAD] %-28s  %6.2f GiB  %5.2f s  (%4.1f GiB/s)\n",
            label,
            (double)data_bytes / (1024.0*1024.0*1024.0),
            elapsed,
            elapsed > 0 ? (double)data_bytes / (1024.0*1024.0*1024.0) / elapsed : 0.0);

    return (uint8_t *)map + extra;   
}

/* -------------------------------------------------------------------------
   Parallel mmap worker 
   ------------------------------------------------------------------------- */
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

/* -------------------------------------------------------------------------
   read_C_vector()
   ------------------------------------------------------------------------- */
static uint64_t read_C_vector(FILE *fp, uint64_t *Cout)
{
    for (int i = 0; i < ALPHABET_SIZE; ++i) {
        if (fscanf(fp, "%lu", &c_vec[i]) != 1) {
            fprintf(stderr, "C-vector read error at col %d\n", i);
            exit(1);
        }
    }
    c_vec[ALPHABET_SIZE] = BWT_SIZE_REFERENCE_SIZE;
    return 0;
}

/* -------------------------------------------------------------------------
   load_occ_full_binary_mmap()
   ------------------------------------------------------------------------- */
static uint64_t load_occ_full_binary_mmap(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(EXIT_FAILURE); }

    rs_occ_full_header_t hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1) {
        fprintf(stderr, "Error: failed to read OCC full header\n");
        fclose(f);                                   
        exit(EXIT_FAILURE);
    }

    if (fseeko(f, 0, SEEK_END) != 0) { perror(path); fclose(f); exit(EXIT_FAILURE); }
    off_t file_size = ftello(f);
    fclose(f);
    if (file_size < 0) { perror(path); exit(EXIT_FAILURE); }

    if (hdr.magic != RS_OCC_MAGIC) {
        fprintf(stderr,
            "Error: %s is not a RosaSeed cp_occ file (bad magic).\n", path);
        exit(EXIT_FAILURE);
    }

    if (hdr.occ_interval != OCC_INTERVAL || hdr.alphabet_size != ALPHABET_SIZE) {
        fprintf(stderr, "Error: OCC full binary metadata mismatch\n");
        exit(EXIT_FAILURE);
    }

    /* floor(N/32) + 1 blocks, so GET_OCC32 can read block (pos >> 5) for
       every pos in 0..N. Must match make_cp_occ_2step.c. */
    uint64_t expected_blocks =
        hdr.bwt_len_total_with_dollar / OCC_INTERVAL + 1;
    if (hdr.num_blocks != expected_blocks) {
        fprintf(stderr,
            "Error: OCC block mismatch: header=%llu expected=%llu\n",
            (unsigned long long)hdr.num_blocks,
            (unsigned long long)expected_blocks);
        exit(EXIT_FAILURE);
    }

    {
        uint64_t expected_size =
            (uint64_t)sizeof(hdr) + hdr.num_blocks * (uint64_t)sizeof(cp_occ32_t);

        if ((uint64_t)file_size != expected_size) {
            fprintf(stderr,
                "Error: %s is %llu bytes, expected %llu.\n"
                "       The index may be truncated, or built by a different\n"
                "       version of make_cp_occ_2step.  Regenerate with:\n"
                "         make_cp_occ_2step <bwt> cp_occ_full.bin c_vector.txt\n",
                path,
                (unsigned long long)file_size,
                (unsigned long long)expected_size);
            exit(EXIT_FAILURE);
        }
    }

    size_t body_bytes = hdr.num_blocks * sizeof(cp_occ32_t);

    void *mapped = mmap_file_populate(path,
                                      sizeof(hdr),   
                                      body_bytes,
                                      "cp_occ (mmap+MAP_POPULATE)",
                                      NULL);

    if (cp_occ) { _mm_free(cp_occ); cp_occ = NULL; }
    cp_occ = (cp_occ32_t *)mapped;

    sentinel_index = hdr.sentinel_index;

    fprintf(stderr,
        "Full OCC+BWT bitmaps loaded: bwt_non_dollar=%llu total=%llu "
        "blocks=%llu sentinel=%lld\n",
        (unsigned long long)hdr.bwt_len_non_dollar,
        (unsigned long long)hdr.bwt_len_total_with_dollar,
        (unsigned long long)hdr.num_blocks,
        (long long)hdr.sentinel_index);

    return hdr.bwt_len_non_dollar;
}

/* -------------------------------------------------------------------------
   load_bwt_data_structures() 

   Strategy:
     Thread 1: mmap sa_ls_word  with MAP_POPULATE (~11 GiB at CF=2)
     Thread 2: mmap sa_ms_byte  with MAP_POPULATE (~3 GiB)
     Thread 3: mmap jump_pointers with MAP_POPULATE (~2 GiB)
     Thread 4: mmap ref16_packed with MAP_POPULATE (~3 GiB)
     Main:     mmap cp_occ body  with MAP_POPULATE (~23 GiB)
               + c_vec (tiny text file)
               + SA sentinel scan
   ------------------------------------------------------------------------- */
void load_bwt_data_structures(void)
{
    struct timespec wall_t0, wall_t1;
    clock_gettime(CLOCK_MONOTONIC, &wall_t0);

    fprintf(stderr, "\n[LOAD] Starting parallel mmap index load (MAP_POPULATE)...\n");

    /* c_vec: tiny text file, load instantly on main */
    {
        FILE *fc = fopen(C_VEC_FILE, "r");
        if (!fc) { perror("open C-vector"); exit(1); }
        read_C_vector(fc, c_vec);
        fclose(fc);
    }

    /* Read ref16 header to get exact packed byte count */
    uint64_t ref16_symbols_from_file = 0;
    {
        FILE *fp_ref_hdr = fopen(REFERENCE_GENOME_FILE, "rb");
        if (!fp_ref_hdr) { perror(REFERENCE_GENOME_FILE); exit(EXIT_FAILURE); }
        if (fread(&ref16_symbols_from_file, sizeof(uint64_t), 1, fp_ref_hdr) != 1) {
            fprintf(stderr, "Error: failed to read ref16 symbol count\n");
            exit(EXIT_FAILURE);
        }
        fclose(fp_ref_hdr);
    }
    const size_t ref_packed_bytes = (ref16_symbols_from_file + 1) >> 1;
    fprintf(stderr,
        "[LOAD] ref16: symbols=%llu packed_bytes=%zu (%.3f GiB)\n",
        (unsigned long long)ref16_symbols_from_file,
        ref_packed_bytes,
        (double)ref_packed_bytes / (1024.0*1024.0*1024.0));

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
    MmapTask task_ref = {
        .path       = REFERENCE_GENOME_FILE,
        .skip_bytes = sizeof(uint64_t),   /* skip 8-byte symbol count header */
        .data_bytes = ref_packed_bytes,
        .label      = "ref16_packed (mmap)",
        .dest_ptr   = (void **)&ref16_packed,
        .done       = 0,
    };

    free(sa_ls_word);   sa_ls_word   = NULL;
    free(sa_ms_byte);   sa_ms_byte   = NULL;
    free(jump_pointers); jump_pointers = NULL;
    if (ref16_packed) { free(ref16_packed); ref16_packed = NULL; }

    /* Launch 4 mmap threads simultaneously */
    pthread_t thr_sa_ls, thr_sa_ms, thr_jt, thr_ref;
    pthread_create(&thr_sa_ls, NULL, mmap_worker, &task_sa_ls);
    pthread_create(&thr_sa_ms, NULL, mmap_worker, &task_sa_ms);
    pthread_create(&thr_jt,    NULL, mmap_worker, &task_jt);
    pthread_create(&thr_ref,   NULL, mmap_worker, &task_ref);

    /* Main thread: mmap cp_occ (largest, 23 GiB) in parallel with threads */
    fprintf(stderr,
        "[RS_LOAD_DATA_DEBUG] SA_COMPRESSION_FACTOR_POWER=%d CF=%d\n",
        SA_COMPRESSION_FACTOR_POWER,
        1 << SA_COMPRESSION_FACTOR_POWER);

    uint64_t bwt_len = load_occ_full_binary_mmap(OCC_BIN_FILE);

    fprintf(stderr,
        "BWT(non-$) length = %llu (total = %llu), OCC interval = %d\n",
        (unsigned long long)bwt_len,
        (unsigned long long)(bwt_len + 1),
        OCC_INTERVAL);

    pthread_join(thr_sa_ls, NULL);
    pthread_join(thr_sa_ms, NULL);
    pthread_join(thr_jt,    NULL);
    pthread_join(thr_ref,   NULL);

    if (!task_sa_ls.done || !task_sa_ms.done ||
        !task_jt.done    || !task_ref.done) {
        fprintf(stderr, "[ERROR] One or more mmap threads failed\n");
        exit(EXIT_FAILURE);
    }

    /* Sentinel scan */
    if(sentinel_index <0){
        for (size_t v = 0; v < sa_n; ++v) {
            if (sa_ls_word[v] == 0 && sa_ms_byte[v] == 0) {
                sentinel_index = (int64_t)((uint64_t)v << SA_COMPRESSION_FACTOR_POWER);
                break;
            }
        }
    }
    
    fprintf(stderr, "SA loaded: %zu sampled entries (CF=%d), sentinel=%ld\n",
            sa_n, 1 << SA_COMPRESSION_FACTOR_POWER, sentinel_index);

    clock_gettime(CLOCK_MONOTONIC, &wall_t1);
    double wall_elapsed = (wall_t1.tv_sec - wall_t0.tv_sec) +
                          (wall_t1.tv_nsec - wall_t0.tv_nsec) * 1e-9;

    fprintf(stderr,
        "\n[LOAD] All structures mmap'd in %.2f sec wall time\n"
        "       (second run will be faster: pages remain in page cache)\n\n",
        wall_elapsed);

    fprintf(stderr, "Jump table loaded: %zu entries\n",
            (size_t)JUMP_TABLE_ENTRIES);
}

static void load_sa_binary(const char *path_ls, const char *path_msb)
{
    FILE *f32 = fopen(path_ls,  "rb");
    fprintf(stderr,
        "[RS_LOAD_DATA_DEBUG] SA_COMPRESSION_FACTOR_POWER=%d CF=%d\n",
        SA_COMPRESSION_FACTOR_POWER, 1 << SA_COMPRESSION_FACTOR_POWER);
    FILE *f8  = fopen(path_msb, "rb");
    if (!f32) { perror(path_ls);  exit(1); }
    if (!f8)  { perror(path_msb); exit(1); }
    const size_t sa_n = (size_t)SA_SAMPLED_SIZE;
    size_t r32 = fread(sa_ls_word, sizeof(uint32_t), sa_n, f32);
    size_t r8  = fread(sa_ms_byte, sizeof(uint8_t),  sa_n, f8);
    if (r32 != sa_n || r8 != sa_n) {
        fprintf(stderr,
            "SA binary read mismatch: expected %zu, got ls=%zu ms=%zu\n",
            sa_n, r32, r8);
        exit(1);
    }
    fclose(f32); fclose(f8);
    for (size_t v = 0; v < sa_n; v++) {
        if (sa_ls_word[v] == 0 && sa_ms_byte[v] == 0) {
            sentinel_index = (int64_t)(v << SA_COMPRESSION_FACTOR_POWER);
            break;
        }
    }
    fprintf(stderr, "SA loaded from binary: %zu sampled entries (CF=%d)\n",
            sa_n, 1 << SA_COMPRESSION_FACTOR_POWER);
    fprintf(stderr, "Sentinel index: %ld\n", sentinel_index);
}

static uint64_t load_occ_full_binary(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(EXIT_FAILURE); }
    rs_occ_full_header_t hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1) {
        fprintf(stderr, "Error: failed to read OCC full header\n");
        exit(EXIT_FAILURE);
    }
    if (hdr.occ_interval != OCC_INTERVAL || hdr.alphabet_size != ALPHABET_SIZE) {
        fprintf(stderr, "Error: OCC full binary metadata mismatch\n");
        exit(EXIT_FAILURE);
    }
    /* floor(N/32) + 1 blocks, so GET_OCC32 can read block (pos >> 5) for
       every pos in 0..N. Must match make_cp_occ_2step.c. */
    uint64_t expected_blocks =
        hdr.bwt_len_total_with_dollar / OCC_INTERVAL + 1;
    if (hdr.num_blocks != expected_blocks) {
        fprintf(stderr,
            "Error: OCC block mismatch: header blocks=%llu expected=%llu\n",
            (unsigned long long)hdr.num_blocks,
            (unsigned long long)expected_blocks);
        exit(EXIT_FAILURE);
    }
    size_t got = fread(cp_occ, sizeof(cp_occ32_t), hdr.num_blocks, f);
    if (got != hdr.num_blocks) {
        fprintf(stderr,
            "Error: cp_occ full read mismatch: expected %llu got %zu\n",
            (unsigned long long)hdr.num_blocks, got);
        exit(EXIT_FAILURE);
    }
    fclose(f);
    sentinel_index = hdr.sentinel_index;
    fprintf(stderr,
        "Full OCC+BWT bitmaps loaded: bwt_non_dollar=%llu total=%llu "
        "blocks=%llu sentinel=%lld\n",
        (unsigned long long)hdr.bwt_len_non_dollar,
        (unsigned long long)hdr.bwt_len_total_with_dollar,
        (unsigned long long)hdr.num_blocks,
        (long long)hdr.sentinel_index);
    return hdr.bwt_len_non_dollar;
}
