/*
 * convert_sa_to_bin_CF.c
 *
 * Reads a binary suffix array (int64_t per entry, as produced by gsufsort-64
 * after trimming the first entry) and writes compressed SA files for
 * compression factors CF = 1, 2, 4, 8.
 *
 * For each CF:
 *   - Samples every CF-th SA entry  (entries where row % CF == 0)
 *   - Splits each 64-bit value into:
 *       sa_ls_word_cfX.bin  —  lower 32 bits  (uint32_t, 4 bytes each)
 *       sa_ms_byte_cfX.bin  —  bits 32..39    (uint8_t,  1 byte  each)
 *   - Seeding kernel reconstructs: val = ((uint64_t)ms << 32) | ls
 *
 * Usage:
 *   ./convert_sa_to_bin_CF <sa_binary_file> <out_dir>
 *
 * Example:
 *   ./convert_sa_to_bin_CF genome.txt.8.sa /data/index/
 *
 * The SA_N (number of entries) is computed automatically from the file size.
 * No hardcoding needed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#define READ_BUF_ENTRIES  (1 << 22)     /* 4M entries = 32 MB read buffer  */
#define WRITE_BUF_ENTRIES (1 << 20)     /* 1M entries write buffer per file */
#define REPORT_EVERY      100000000ULL  /* progress every 100M rows         */

/* ── get file size in bytes ─────────────────────────────────────────────── */
static uint64_t file_size_bytes(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(EXIT_FAILURE); }
    if (fseeko(f, 0, SEEK_END) != 0) { perror("fseeko"); exit(EXIT_FAILURE); }
    off_t sz = ftello(f);
    fclose(f);
    return (uint64_t)sz;
}

/* ── process one compression factor ────────────────────────────────────── */
static void convert_one_cf(const char *sa_path,
                            const char *out_dir,
                            uint64_t    sa_n,
                            int         cf)
{
    if (cf != 1 && cf != 2 && cf != 4 && cf != 8) {
        fprintf(stderr, "Invalid CF=%d (must be 1,2,4,8)\n", cf);
        exit(EXIT_FAILURE);
    }

    int      power       = __builtin_ctz((unsigned)cf);   /* log2(cf) */
    uint64_t sample_mask = (uint64_t)cf - 1;
    uint64_t sa_sampled  = (sa_n + sample_mask) >> power; /* ceil(sa_n / cf) */

    /* ── output file paths ── */
    char path_ls[1024], path_ms[1024];
    snprintf(path_ls, sizeof(path_ls), "%s/sa_ls_word_cf%d.bin", out_dir, cf);
    snprintf(path_ms, sizeof(path_ms), "%s/sa_ms_byte_cf%d.bin", out_dir, cf);

    /* ── open files ── */
    FILE *fin = fopen(sa_path, "rb");     /* binary read */
    FILE *f32 = fopen(path_ls, "wb");
    FILE *f8  = fopen(path_ms, "wb");

    if (!fin) { perror(sa_path); exit(EXIT_FAILURE); }
    if (!f32) { perror(path_ls); exit(EXIT_FAILURE); }
    if (!f8)  { perror(path_ms); exit(EXIT_FAILURE); }

    /* ── allocate buffers ── */
    uint64_t *rbuf  = (uint64_t *)malloc(READ_BUF_ENTRIES  * sizeof(uint64_t));
    uint32_t *wls   = (uint32_t *)malloc(WRITE_BUF_ENTRIES * sizeof(uint32_t));
    uint8_t  *wms   = (uint8_t  *)malloc(WRITE_BUF_ENTRIES * sizeof(uint8_t));

    if (!rbuf || !wls || !wms) { perror("malloc"); exit(EXIT_FAILURE); }

    uint64_t row           = 0;   /* current SA row */
    uint64_t total_written = 0;
    uint64_t w_buf         = 0;   /* write buffer fill level */
    uint64_t sentinel_row  = UINT64_MAX;

    clock_t t0 = clock();

    fprintf(stderr, "[CF=%d] sa_n=%llu  sampled=%llu  paths:\n"
                    "        LS: %s\n"
                    "        MS: %s\n",
            cf, (unsigned long long)sa_n,
            (unsigned long long)sa_sampled, path_ls, path_ms);

    /* ── main loop: read large chunks, scatter into write buffers ── */
    size_t nread;
    while ((nread = fread(rbuf, sizeof(uint64_t),
                          READ_BUF_ENTRIES, fin)) > 0) {

        for (size_t i = 0; i < nread; i++, row++) {

            uint64_t x = rbuf[i];

            /* track sentinel (SA entry == 0 marks the $ position) */
            if (x == 0) sentinel_row = row;

            /* sample every CF-th entry */
            if ((row & sample_mask) == 0) {
                wls[w_buf] = (uint32_t)( x         & 0xFFFFFFFFULL);
                wms[w_buf] = (uint8_t) ((x >> 32)  & 0xFFULL);
                w_buf++;

                /* flush write buffers when full */
                if (w_buf == WRITE_BUF_ENTRIES) {
                    if (fwrite(wls, sizeof(uint32_t), w_buf, f32) != w_buf ||
                        fwrite(wms, sizeof(uint8_t),  w_buf, f8)  != w_buf) {
                        fprintf(stderr, "Write error CF=%d\n", cf);
                        exit(EXIT_FAILURE);
                    }
                    total_written += w_buf;
                    w_buf = 0;

                    double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
                    fprintf(stderr, "\r[CF=%d]  %6.2f%%  (%llu / %llu sampled)  %.0fs",
                            cf,
                            100.0 * (double)total_written / (double)sa_sampled,
                            (unsigned long long)total_written,
                            (unsigned long long)sa_sampled,
                            el);
                }
            }
        }

        /* progress after each read chunk even if no flush yet */
        if (row % REPORT_EVERY == 0) {
            double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
            fprintf(stderr, "\r[CF=%d]  row %llu / %llu  %.0fs   ",
                    cf, (unsigned long long)row,
                    (unsigned long long)sa_n, el);
        }
    }

    /* flush remaining write buffer */
    if (w_buf > 0) {
        if (fwrite(wls, sizeof(uint32_t), w_buf, f32) != w_buf ||
            fwrite(wms, sizeof(uint8_t),  w_buf, f8)  != w_buf) {
            fprintf(stderr, "Final write error CF=%d\n", cf);
            exit(EXIT_FAILURE);
        }
        total_written += w_buf;
    }

    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;

    fclose(fin); fclose(f32); fclose(f8);
    free(rbuf);  free(wls);   free(wms);

    /* ── summary ── */
    fprintf(stderr,
            "\r[CF=%d]  done in %.1fs\n"
            "  SA rows read             : %llu\n"
            "  sampled entries expected : %llu\n"
            "  sampled entries written  : %llu\n"
            "  sentinel (SA==0) row     : %llu\n"
            "  LS file size             : %llu bytes  (%.2f GB)\n"
            "  MS file size             : %llu bytes  (%.2f GB)\n",
            cf, elapsed,
            (unsigned long long)row,
            (unsigned long long)sa_sampled,
            (unsigned long long)total_written,
            (unsigned long long)sentinel_row,
            (unsigned long long)(total_written * 4),
            (double)(total_written * 4) / 1073741824.0 ,
            (unsigned long long)(total_written * 1),
            (double)(total_written * 1) / 1073741824.0 );

    if (total_written != sa_sampled) {
        fprintf(stderr,
                "[ERROR] CF=%d: wrote %llu but expected %llu\n",
                cf,
                (unsigned long long)total_written,
                (unsigned long long)sa_sampled);
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "  [OK] CF=%d count matches.\n\n", cf);
}

/* ── main ───────────────────────────────────────────────────────────────── */
int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr,
            "\nUsage: %s <sa_binary_file> <out_dir>\n\n"
            "  sa_binary_file  trimmed binary SA from gsufsort-64 (.8.sa)\n"
            "                  each entry is a little-endian int64_t\n"
            "  out_dir         directory to write CF output files into\n\n"
            "Outputs for each CF in {1,2,4,8}:\n"
            "  sa_ls_word_cfX.bin  lower 32 bits  (uint32_t per sampled entry)\n"
            "  sa_ms_byte_cfX.bin  bits 32..39    (uint8_t  per sampled entry)\n\n"
            "Kernel reconstructs: sa_val = ((uint64_t)ms << 32) | ls\n\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    const char *sa_path = argv[1];
    const char *out_dir = argv[2];

    /* auto-detect SA_N from file size */
    uint64_t sa_bytes = file_size_bytes(sa_path);
    if (sa_bytes % sizeof(uint64_t) != 0) {
        fprintf(stderr,
                "[ERROR] SA file size %llu is not a multiple of 8.\n"
                "        Is this the correct binary SA file?\n",
                (unsigned long long)sa_bytes);
        return EXIT_FAILURE;
    }
    uint64_t sa_n = sa_bytes / sizeof(uint64_t);

    fprintf(stderr,
            "──────────────────────────────────────────────\n"
            " convert_sa_to_bin_CF\n"
            "──────────────────────────────────────────────\n"
            " SA file  : %s\n"
            " SA size  : %llu bytes\n"
            " SA_N     : %llu entries\n"
            " Out dir  : %s\n"
            "──────────────────────────────────────────────\n\n",
            sa_path,
            (unsigned long long)sa_bytes,
            (unsigned long long)sa_n,
            out_dir);

    int cfs[] = {1, 2, 4, 8};
    for (int i = 0; i < 4; i++)
        convert_one_cf(sa_path, out_dir, sa_n, cfs[i]);

    fprintf(stderr,
            "──────────────────────────────────────────────\n"
            " All CF files generated.\n"
            "──────────────────────────────────────────────\n");
    return EXIT_SUCCESS;
}
