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

/*
 * make_jumptable_compact.c
 *
 * Generates binary jump tables for 14-mer, 15-mer, and/or 16-mer k-mer seeding.
 *
 * Uses cp_occ.bin for all Occ() queries, no BWT file needed.
 * Uses SA CF=1 files for unique entry SA lookup (direct lookup, no walk).
 *
 * ── Why CF=1 (not CF=4) ───────────────────────────────────────────────────────
 *
 *   Our SA CF files sample by ROW index (rows 0,4,8,...).
 *   For a walk-based lookup to be bounded, sampling must be by SA VALUE
 *   (rows where SA[row] % CF == 0), so that LF-mapping always reaches a
 *   sampled row within CF-1 steps.
 *   With row-based sampling, the walk can take thousands of steps → use CF=1.
 *   CF=1 stores every SA entry → direct O(1) lookup, no walk at all.
 *
 * ── Jump pointer encoding ─────────────────────────────────────────────────────
 *
 *   diff = hi - lo
 *   diff == 0  → jp = 0                              (no mapping)
 *   diff == 1  → jp = (1ULL<<63) | (SA[lo] & 0x1FFFFFFFFULL)   (unique)
 *   diff  > 1  → jp = lo | (diff << 34)              (non-unique)
 *
 * ── rank_excl, the key formula ───────────────────────────────────────────────
 *
 *   rank_excl(c, i) = count of c in BWT[0..i-1]  (exclusive upper bound)
 *
 *   This is what the FM-index backward search needs:
 *     lo = C[c] + rank_excl(c, lo)
 *     hi = C[c] + rank_excl(c, hi)
 *
 *   Implemented from cp_occ.bin:
 *     if i == 0: return 0
 *     blk = (i-1) >> 5
 *     off = (i-1) & 31
 *     return cp_count[blk][c] + popcount(one_hot[blk][c] >> (31-off))
 *
 * ── Memory requirements ───────────────────────────────────────────────────────
 *
 *   cp_occ.bin blocks    :  ~6.2 GB
 *   sa_ls_word_cf1.bin   : ~24.9 GB
 *   sa_ms_byte_cf1.bin   :  ~6.2 GB
 *   ─────────────────────────────────
 *   Total loaded          : ~37 GB
 *
 * ── Output file sizes ─────────────────────────────────────────────────────────
 *
 *   jumptable_14nt.bin :  4^14 × 8 =  ~2.1 GB
 *   jumptable_15nt.bin :  4^15 × 8 =  ~8.6 GB
 *   jumptable_16nt.bin :  4^16 × 8 = ~34.4 GB
 *
 * Usage:
 *   ./make_jumptable_compact <cp_occ.bin> <sa_ls_cf1.bin> <sa_ms_cf1.bin> <out_dir> [k]
 *
 *   k = 14, 15, 16, or 0 (default) for all three
 *
 * Examples:
 *   ./make_jumptable_compact genome.cp_occ.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin /data/index
 *   ./make_jumptable_compact genome.cp_occ.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin /data/index 15
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <time.h>

/* ── must match make_cp_occ.c exactly ──────────────────────────────────────── */
#define RS1_OCC_MAGIC  0x5253314F43433332ULL

typedef struct {
    uint32_t cp_count[4];   /* cumulative A/C/G/T counts before this block   */
    uint32_t one_hot[4];    /* bitvectors: bit(31-off) = BWT[blk*32+off]=base */
} cp_occ_block_t;

typedef struct {
    uint64_t magic;
    uint32_t version, occ_interval, alphabet_size, reserved;
    uint64_t bwt_len_no_term, bwt_len_total, num_blocks;
    int64_t  sentinel_index;
    uint64_t c_array[4];
} rs1_occ_header_t;

/* ── global loaded data ─────────────────────────────────────────────────────── */
static cp_occ_block_t *blocks = NULL;
static uint64_t        C[4];
static uint64_t        N;
static uint64_t        Nblocks;

static uint32_t       *sa_ls = NULL;   /* lower 32 bits of SA, CF=1 (all rows)  */
static uint8_t        *sa_ms = NULL;   /* bits 32-39   of SA, CF=1 (all rows)  */

/* ── rank_excl ───────────────────────────────────────────────────────────────
 *
 *  rank_excl(c, i) = count of c in BWT[0..i-1]
 *
 *  Proof that this (not inclusive) is correct for FM-index backward search:
 *    LF(r) = C[c] + #{j < r : BWT[j] = c}  = C[c] + rank_excl(c, r)
 *    SA[LF(r)] = SA[r] - 1  ← verified with BANANA$ example in comments above
 *
 *  With inclusive occ: LF(r) = C[c] + #{j<=r : BWT[j]=c} gives WRONG result.
 */
static inline uint64_t rank_excl(int c, uint64_t i)
{
    if (i == 0) return 0ULL;
    uint64_t pos = i - 1;
    uint64_t blk = pos >> 5;
    uint32_t off = (uint32_t)(pos & 31u);
    return (uint64_t)blocks[blk].cp_count[c]
         + (uint64_t)__builtin_popcount(blocks[blk].one_hot[c] >> (31u - off));
}

/* ── SA lookup, CF=1, direct lookup, no walk needed ────────────────────────
 *
 *  With CF=1, every row is sampled → sa_ls[row] / sa_ms[row] = SA[row].
 *  No LF-mapping walk required.
 */
static inline uint64_t get_sa(uint64_t row)
{
    return ((uint64_t)sa_ms[row] << 32) | (uint64_t)sa_ls[row];
}

/* ── jump-pointer field widths ───────────────────────────────────────────────
 * Packed layout:
 *   bit  63     : unique flag
 *   bits 62..34 : diff  (29 bits)
 *   bits 33..0  : lo, or SA position when unique (34 bits; unique SA values
 *                 only use bits 32..0, i.e. < 2^33)
 *
 * JT_DIFF_MAX is overridable at compile time so a test can exercise the
 * overflow path with a small limit (-DJT_DIFF_MAX=<n>).
 */
#ifndef JT_DIFF_MAX
#define JT_DIFF_MAX 536870911ULL   /* 2^29 - 1 */
#endif
#define JT_LO_BITS  34
#define JT_LO_MAX   ((1ULL << JT_LO_BITS) - 1ULL)        /* 17,179,869,183 */
#define JT_SA_BITS  33
#define JT_SA_MAX   ((1ULL << JT_SA_BITS) - 1ULL)        /*  8,589,934,591 */

/* ── jump pointer encoding ───────────────────────────────────────────────────
 *  Matches 2-step RosaSeed convert_jt_txt_2_bin.c exactly.
 */
static inline uint64_t encode_jp(uint64_t lo, uint64_t hi, const char *table_name,
                                  uint64_t *max_diff_seen)
{
    uint64_t diff = hi - lo;
    if (diff == 0) return 0ULL;

    if (diff > *max_diff_seen) *max_diff_seen = diff;

    if (diff == 1) {
        uint64_t sa = get_sa(lo);
        if (sa > JT_SA_MAX) {
            fprintf(stderr,
                "\n[ERROR] %s: unique SA value %llu exceeds the %d-bit field (limit %llu).\n"
                "        BWT length is too large for this jump-table encoding.\n",
                table_name, (unsigned long long)sa, JT_SA_BITS,
                (unsigned long long)JT_SA_MAX);
            exit(EXIT_FAILURE);
        }
        return (1ULL << 63) | (sa & JT_SA_MAX);
    }

    if (lo > JT_LO_MAX) {
        fprintf(stderr,
            "\n[ERROR] %s: interval start %llu exceeds the %d-bit lo field (limit %llu).\n"
            "        BWT length is too large for this jump-table encoding.\n",
            table_name, (unsigned long long)lo, JT_LO_BITS,
            (unsigned long long)JT_LO_MAX);
        exit(EXIT_FAILURE);
    }

    if (diff > JT_DIFF_MAX) {
        fprintf(stderr,
            "\n[ERROR] %s: interval width %llu exceeds the 29-bit diff field (limit %llu).\n"
            "        Encoding it would spill into bit 63 and the aligner would decode\n"
            "        this entry as a unique hit at a bogus position.\n",
            table_name, (unsigned long long)diff, (unsigned long long)JT_DIFF_MAX);
        exit(EXIT_FAILURE);
    }

    return lo | (diff << 34);
}

/* ── single k-mer table generation ──────────────────────────────────────────── */
static void make_table(int k, const char *out_path)
{
    uint64_t total = 1ULL << (2 * k);

    FILE *fout = fopen(out_path, "wb");
    if (!fout) { perror(out_path); exit(EXIT_FAILURE); }

    const size_t WBUF = 1u << 20;   /* 1M entries = 8 MB */
    uint64_t *buf = (uint64_t *)malloc(WBUF * sizeof(uint64_t));
    if (!buf) { perror("malloc"); exit(EXIT_FAILURE); }

    uint64_t n_buf = 0, written = 0;
    uint64_t cnt_unique = 0, cnt_multi = 0, cnt_none = 0;
    uint64_t max_diff_seen = 0;
    clock_t t0 = clock();

    /*
     * Enumerate all 4^k k-mers in index order:
     *   index = b[0]*4^(k-1) + ... + b[k-1]*4^0
     *   b[j] = (idx >> (2*(k-1-j))) & 3
     *
     * Backward search processes right-to-left (b[k-1] first):
     *   j = k-1 → c = (idx >> 0) & 3 = b[k-1]  (rightmost) ← first step ✓
     *   j = 0   → c = (idx >> 2*(k-1)) & 3 = b[0] (leftmost) ← last step ✓
     *
     * rank_excl(c, lo) = count of c in BWT[0..lo-1]
     *   = old occ_mat[lo][c] exactly
     */
    for (uint64_t idx = 0; idx < total; idx++) {

        uint64_t lo = 0, hi = N;

        for (int j = k - 1; j >= 0 && lo < hi; j--) {
            int c = (int)((idx >> (2u * (unsigned)(k - 1 - j))) & 3u);
            lo = C[c] + rank_excl(c, lo);
            hi = C[c] + rank_excl(c, hi);
        }

        uint64_t diff = hi - lo;
        if      (diff == 0) cnt_none++;
        else if (diff == 1) cnt_unique++;
        else                cnt_multi++;

        buf[n_buf++] = encode_jp(lo, hi, out_path, &max_diff_seen);

        if (n_buf == WBUF) {
            if (fwrite(buf, sizeof(uint64_t), n_buf, fout) != n_buf) {
                fprintf(stderr, "\n[ERROR] write failed for k=%d\n", k);
                exit(EXIT_FAILURE);
            }
            written += n_buf;
            n_buf = 0;
            if ((written & ((64u << 20) - 1)) == 0) {
                double el = (double)(clock() - t0) / CLOCKS_PER_SEC;
                fprintf(stderr,
                    "\r  k=%d  %5.1f%%  (%.0fs)  unique=%.1f%%  none=%.1f%%",
                    k, 100.0 * written / (double)total, el,
                    100.0 * cnt_unique / (double)written,
                    100.0 * cnt_none   / (double)written);
            }
        }
    }
    if (n_buf > 0) {
        if (fwrite(buf, sizeof(uint64_t), n_buf, fout) != n_buf) {
            fprintf(stderr, "\n[ERROR] final write failed k=%d\n", k);
            exit(EXIT_FAILURE);
        }
        written += n_buf;
    }

    fclose(fout);
    free(buf);

    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "\r  k=%d  100.0%%  done in %.1fs                               \n", k, elapsed);
    fprintf(stderr, "  entries written  : %llu\n",  (unsigned long long)written);
    fprintf(stderr, "  unique mapping   : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_unique, 100.0*cnt_unique/(double)total);
    fprintf(stderr, "  multi mapping    : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_multi,  100.0*cnt_multi /(double)total);
    fprintf(stderr, "  no mapping       : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_none,   100.0*cnt_none  /(double)total);
    fprintf(stderr, "  output           : %s  (%.2f GB)\n",
            out_path, (double)(written * 8) / 1073741824.0 );
    fprintf(stderr, "  max interval = %llu (limit %llu)\n",
            (unsigned long long)max_diff_seen, (unsigned long long)JT_DIFF_MAX);
}

/* ── main ────────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv)
{
    if (argc < 5 || argc > 6) {
        fprintf(stderr,
            "\nUsage: %s <cp_occ.bin> <sa_ls_cf1.bin> <sa_ms_cf1.bin> <out_dir> [k]\n\n"
            "  cp_occ.bin       combined checkpoint+bitvector from make_cp_occ\n"
            "  sa_ls_cf1.bin    lower 32 bits of full SA (sa_ls_word_cf1.bin)\n"
            "  sa_ms_cf1.bin    bits 32-39  of full SA (sa_ms_byte_cf1.bin)\n"
            "  out_dir          directory for jumptable_Xnt.bin output files\n"
            "  k                14, 15, 16, or 0/omit for all three (default)\n\n"
            "Memory required:  ~37 GB  (6.2 cp_occ + 24.9 sa_ls + 6.2 sa_ms)\n\n"
            "Examples:\n"
            "  %s genome.cp_occ.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin /data/index\n"
            "  %s genome.cp_occ.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin /data/index 15\n\n",
            argv[0], argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    int k_only = (argc == 6) ? atoi(argv[5]) : 0;
    if (k_only != 0 && k_only != 14 && k_only != 15 && k_only != 16) {
        fprintf(stderr, "[ERROR] k must be 14, 15, 16, or 0 for all.\n");
        return EXIT_FAILURE;
    }

    /* ── load cp_occ.bin ── */
    fprintf(stderr, "Loading cp_occ.bin ...\n");
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return EXIT_FAILURE; }

    rs1_occ_header_t hdr;
    if (fread(&hdr, sizeof(hdr), 1, f) != 1) {
        fprintf(stderr, "[ERROR] failed to read header\n"); return EXIT_FAILURE;
    }
    if (hdr.magic != RS1_OCC_MAGIC) {
        fprintf(stderr, "[ERROR] bad magic\n"); return EXIT_FAILURE;
    }
    if (hdr.occ_interval != 32) {
        fprintf(stderr, "[ERROR] occ_interval=%u, only 32 supported\n",
                hdr.occ_interval); return EXIT_FAILURE;
    }

    N       = hdr.bwt_len_total;
    Nblocks = hdr.num_blocks;
    C[0] = hdr.c_array[0]; C[1] = hdr.c_array[1];
    C[2] = hdr.c_array[2]; C[3] = hdr.c_array[3];

    blocks = (cp_occ_block_t *)malloc(Nblocks * sizeof(cp_occ_block_t));
    if (!blocks) { perror("malloc blocks"); return EXIT_FAILURE; }

    fprintf(stderr, "  reading %llu blocks (%.2f GB) ...\n",
            (unsigned long long)Nblocks,
            (double)(Nblocks * sizeof(cp_occ_block_t)) / 1073741824.0 );

    if (fread(blocks, sizeof(cp_occ_block_t), Nblocks, f) != Nblocks) {
        fprintf(stderr, "[ERROR] block read failed\n"); return EXIT_FAILURE;
    }
    fclose(f);

    fprintf(stderr, "  BWT length     : %llu\n",  (unsigned long long)N);
    fprintf(stderr, "  sentinel index : %lld\n",  (long long)hdr.sentinel_index);
    fprintf(stderr, "  C[A]=%-12llu C[C]=%-12llu C[G]=%-12llu C[T]=%llu\n",
            (unsigned long long)C[0], (unsigned long long)C[1],
            (unsigned long long)C[2], (unsigned long long)C[3]);

    /* ── load SA CF=1 (all entries, direct lookup) ── */
    uint64_t sa_n = N;   /* CF=1: every row stored */
    fprintf(stderr, "\nLoading SA CF=1 (%llu entries, %.2f GB) ...\n",
            (unsigned long long)sa_n, (double)(sa_n * 5) / 1073741824.0 );

    sa_ls = (uint32_t *)malloc(sa_n * sizeof(uint32_t));
    sa_ms = (uint8_t  *)malloc(sa_n * sizeof(uint8_t));
    if (!sa_ls || !sa_ms) { perror("malloc SA"); return EXIT_FAILURE; }

    FILE *fls = fopen(argv[2], "rb");
    FILE *fms = fopen(argv[3], "rb");
    if (!fls) { perror(argv[2]); return EXIT_FAILURE; }
    if (!fms) { perror(argv[3]); return EXIT_FAILURE; }

    if (fread(sa_ls, sizeof(uint32_t), sa_n, fls) != sa_n) {
        fprintf(stderr, "[ERROR] sa_ls read incomplete\n"); return EXIT_FAILURE;
    }
    if (fread(sa_ms, sizeof(uint8_t), sa_n, fms) != sa_n) {
        fprintf(stderr, "[ERROR] sa_ms read incomplete\n"); return EXIT_FAILURE;
    }
    fclose(fls); fclose(fms);
    fprintf(stderr, "  SA CF=1 loaded.\n");

    /* ── generate tables ── */
    int k_list[] = {14, 15, 16};
    int n_tables = (k_only == 0) ? 3 : 1;

    for (int ki = 0; ki < n_tables; ki++) {
        int k = (k_only == 0) ? k_list[ki] : k_only;
        char out_path[4096];
        snprintf(out_path, sizeof(out_path),
                 "%s/jumptable_%dnt.bin", argv[4], k);

        fprintf(stderr,
                "\n[ k=%d ]  4^%d = %llu entries  →  %.1f GB output\n"
                "  writing to: %s\n",
                k, k, (1ULL << (2*k)),
                (double)(8ULL << (2*k)) / 1073741824.0 , out_path);

        make_table(k, out_path);
    }

    free(blocks); free(sa_ls); free(sa_ms);
    fprintf(stderr, "\nAll done.\n");
    return EXIT_SUCCESS;
}
