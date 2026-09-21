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
 * make_cp_occ_2step.c
 *
 * Single-pass builder of cp_occ_full.bin + c_vector.txt
 * for the 2-step RosaSeed FM-index.
 *
 * Reads the raw BINARY BWT from gsufsort-64 (after trimming first entry).
 * Valid symbols: '0'-'9' (0x30-0x39), 'A'-'F' (0x41-0x46), uppercase only.
 * Terminator:    0x00 (null byte): position advances, no count incremented.
 *
 * ── Block layout (128 bytes, OCC_INTERVAL = 32) ───────────────────────────
 *
 *   cp_count[16]  uint32_t[16]  cumulative counts BEFORE this block
 *   one_hot[16]   uint32_t[16]  bitvectors: bit(31-off) = BWT[blk*32+off]==sym
 *
 *   block 0:  cp_count = [0,...,0]  ALWAYS  (nothing seen yet)
 *   block k:  cp_count = counts in BWT[0 .. k*32 - 1]
 *
 * ── Occ(sym, i) in the seeding kernel ────────────────────────────────────
 *
 *   blk = i >> 5
 *   off = i & 31
 *   return cp_count[blk][sym] + popcount(one_hot[blk][sym] >> (31 - off))
 *
 * ── Terminator (0x00) handling ────────────────────────────────────────────
 *
 *   Position `i` advances for the null byte.
 *   No bit is set in any one_hot[0..15].
 *   No cp_count is incremented.
 *   Its position is recorded as sentinel_index in the header.
 *
 * ── C vector (c_vector.txt) ───────────────────────────────────────────────
 *
 *   16 tab-separated uint64_t values read by the aligner's read_C_vector().
 *   c_vec[sym] = count of BWT chars strictly less than sym
 *   c_vec[0]   = 1   (only the null terminator is below '0')
 *   c_vec[k]   = c_vec[k-1] + count(sym k-1)
 *   c_vec[16]  = BWT_SIZE_REFERENCE_SIZE  (set by the aligner, not this file)
 *
 * Usage:
 *   ./make_cp_occ_2step <bwt_binary> <cp_occ_full.bin> <c_vector.txt>
 *
 * Example:
 *   ./make_cp_occ_2step genome_2step.bwt cp_occ_full.bin c_vector.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <time.h>

/* ── constants ── */
#define OCC_INTERVAL   32
#define ALPHABET_SIZE  16
#define READ_BUF       (64  * 1024 * 1024)   /* 64 MB read buffer  */
#define WRITE_BUF_BLKS (8192)                 /* blocks per write   */
#define REPORT_EVERY   500000000ULL           /* progress interval  */

#define RS_OCC_MAGIC   0x52534F434346554CULL  
#define RS_OCC_VERSION 1u

typedef struct {
    uint32_t cp_count[16];   /* 64 bytes: cumulative counts before block  */
    uint32_t one_hot[16];    /* 64 bytes: bitvectors                      */
} cp_occ32_t;                /* 128 bytes total                           */

typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t occ_interval;
    uint32_t alphabet_size;
    uint32_t reserved;
    uint64_t bwt_len_non_dollar;
    uint64_t bwt_len_total_with_dollar;
    uint64_t num_blocks;
    int64_t  sentinel_index;
    uint64_t reserved2;     
} rs_occ_full_header_t;

_Static_assert(sizeof(rs_occ_full_header_t) == 64,
               "header must be 64 bytes to keep the cp_occ body 64B aligned");
_Static_assert(sizeof(cp_occ32_t) == 128,
               "cp_occ block must be exactly 2 cache lines");
               
static int8_t lut[256];

static void init_lut(void) {
    memset(lut, -1, sizeof(lut));
    /* '0'-'9' → 0-9 */
    for (int i = 0; i <= 9; i++) lut['0' + i] = (int8_t)i;
    /* 'A'-'F' → 10-15 (uppercase) */
    for (int i = 0; i <= 5; i++) lut['A' + i] = (int8_t)(10 + i);
    /* lowercase a-f → same (robustness) */
    for (int i = 0; i <= 5; i++) lut['a' + i] = (int8_t)(10 + i);
    /* 0x00 = terminator → stays -1 */
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr,
            "\nUsage: %s <bwt_binary> <cp_occ_full.bin> <c_vector.txt>\n\n"
            "  bwt_binary     raw binary BWT (gsufsort output, trimmed)\n"
            "  cp_occ_full.bin  output combined checkpoint+bitvector file\n"
            "  c_vector.txt   output C[] prefix sums (16 tab-separated values)\n\n"
            "Example:\n"
            "  %s genome_2step.bwt cp_occ_full.bin c_vector.txt\n\n",
            argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    init_lut();

    FILE *fbwt = fopen(argv[1], "rb");
    FILE *fout  = fopen(argv[2], "wb");
    FILE *fcvec = fopen(argv[3], "w");

    if (!fbwt) { fprintf(stderr,"Cannot open BWT '%s': %s\n", argv[1], strerror(errno)); return 1; }
    if (!fout) { fprintf(stderr,"Cannot open output '%s': %s\n", argv[2], strerror(errno)); return 1; }
    if (!fcvec){ fprintf(stderr,"Cannot open c_vec '%s': %s\n", argv[3], strerror(errno)); return 1; }

    uint8_t      *rbuf   = (uint8_t *)malloc(READ_BUF);
    cp_occ32_t   *wbuf   = (cp_occ32_t *)calloc(WRITE_BUF_BLKS, sizeof(cp_occ32_t));
    if (!rbuf || !wbuf) { fprintf(stderr,"malloc failed\n"); return 1; }

    rs_occ_full_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    if (fwrite(&hdr, sizeof(hdr), 1, fout) != 1) {
        fprintf(stderr,"Header placeholder write failed\n"); return 1;
    }

    /* ── state ── */
    uint64_t running[16] = {0};   /* cumulative counts [0..pos-1] */
    uint64_t pos         = 0;     /* absolute BWT position         */
    uint64_t n_term      = 0;     /* terminator (0x00) count       */
    int64_t  sentinel    = -1;    /* BWT position of 0x00          */
    uint64_t blocks_written = 0;
    uint64_t wbuf_n      = 0;     /* write buffer fill level       */

    cp_occ32_t blk;               /* current block being assembled */
    memset(&blk, 0, sizeof(blk));

    clock_t t0 = clock();
    size_t nread;

    fprintf(stderr, "[make_cp_occ_2step] processing BWT...\n");
    fprintf(stderr, "  OCC_INTERVAL  = %d\n", OCC_INTERVAL);
    fprintf(stderr, "  block size    = %zu bytes\n", sizeof(cp_occ32_t));

    while ((nread = fread(rbuf, 1, READ_BUF, fbwt)) > 0) {
        for (size_t i = 0; i < nread; i++, pos++) {
            uint32_t off = (uint32_t)(pos & (OCC_INTERVAL - 1));  /* 0..31 */

            /* ── block boundary: snapshot counts, reset one_hot ── */
            if (off == 0) {
                for (int s = 0; s < ALPHABET_SIZE; s++) {
                    blk.cp_count[s] = (uint32_t)running[s];
                    blk.one_hot[s]  = 0u;
                }
                /*
                 * block 0: running[] = {0,...,0} → cp_count = {0,...,0}
                 * This is the "row of zeros" the user requires:
                 * BWT[0..31] has not been processed yet, so all counts = 0.
                 */
            }

            int8_t sym = lut[rbuf[i]];

            if (sym >= 0) {
                /* valid base-16 symbol: set bitvector bit and increment count */
                blk.one_hot[sym] |= (1u << (31u - off));
                running[sym]++;
            } else {
                /*
                 * 0x00 null terminator:
                 *   - position advances (pos++)
                 *   - NO bit set in any one_hot
                 *   - NO running count incremented
                 *   - record sentinel position (0-based)
                 */
                if (sentinel == -1) sentinel = (int64_t)pos;
                n_term++;
            }

            /* ── end of block: buffer for writing ── */
            if (off == (uint32_t)(OCC_INTERVAL - 1)) {
                wbuf[wbuf_n++] = blk;
                blocks_written++;

                if (wbuf_n == WRITE_BUF_BLKS) {
                    if (fwrite(wbuf, sizeof(cp_occ32_t), wbuf_n, fout) != wbuf_n) {
                        fprintf(stderr,"Write error at block %llu\n",
                                (unsigned long long)blocks_written);
                        return EXIT_FAILURE;
                    }
                    wbuf_n = 0;
                }
            }
        }

        if (pos % REPORT_EVERY < (uint64_t)nread) {
            double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
            fprintf(stderr, "\r[make_cp_occ_2step]  %7lluM  (%.0fs)  "
                            "0=%llu 5=%llu A=%llu F=%llu",
                    (unsigned long long)(pos/1000000ULL), el,
                    (unsigned long long)running[0],
                    (unsigned long long)running[5],
                    (unsigned long long)running[10],
                    (unsigned long long)running[15]);
        }
    }

    /* ── flush final block: always emit floor(N/32) + 1 blocks ──
     * The aligner's GET_OCC32 reads block (pos >> 5) and queries every
     * position 0..N inclusive (the last symbol range ends at C[16] = N).
     *   remainder != 0: the partially filled block; bits for positions
     *                   remainder..31 are already 0.
     *   remainder == 0: every block is full, so emit a terminal block with
     *                   the final cumulative counts and no bits set.
     */
    uint32_t remainder = (uint32_t)(pos & (OCC_INTERVAL - 1));
    if (remainder == 0) {
        for (int s = 0; s < ALPHABET_SIZE; s++) {
            blk.cp_count[s] = (uint32_t)running[s];
            blk.one_hot[s]  = 0u;
        }
    }
    wbuf[wbuf_n++] = blk;
    blocks_written++;

    /* flush write buffer */
    if (wbuf_n > 0) {
        if (fwrite(wbuf, sizeof(cp_occ32_t), wbuf_n, fout) != wbuf_n) {
            fprintf(stderr,"Final write error\n"); return EXIT_FAILURE;
        }
    }

    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
    const char *sym_names[16] = {
        "0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"
    };
    
    for (int s = 0; s < ALPHABET_SIZE; s++) {
        if (running[s] > 0xFFFFFFFFULL) {
            fprintf(stderr,
                "\n[FATAL] symbol '%s' occurs %llu times, exceeding the 32-bit\n"
                "        cp_count field (max %u).  This reference cannot be\n"
                "        represented by the current cp_occ format.\n"
                "        See README, 'Genome size limits'.\n",
                sym_names[s], (unsigned long long)running[s], 0xFFFFFFFFu);
            return EXIT_FAILURE;
        }
    }

    /* ── C vector ────────────────────────────────────────────────────────
     *
     * c_vec[sym] = count of BWT chars strictly less than sym
     *
     * Alphabet rank:  0x00($) < '0' < '1' < ... < '9' < 'A' < ... < 'F'
     *
     * c_vec[0]  = n_term          ← only null byte is below '0'
     * c_vec[1]  = n_term + count('0')
     * c_vec[k]  = c_vec[k-1] + count(sym k-1)
     *
     * The aligner's read_C_vector() reads exactly 16 values and sets
     * c_vec[16] = BWT_SIZE_REFERENCE_SIZE itself
     */
    uint64_t c_vec[ALPHABET_SIZE];
    c_vec[0] = n_term;
    for (int k = 1; k < ALPHABET_SIZE; k++)
        c_vec[k] = c_vec[k-1] + running[k-1];

    /* write tab-separated, one line */
    for (int k = 0; k < ALPHABET_SIZE; k++) {
        fprintf(fcvec, "%llu", (unsigned long long)c_vec[k]);
        if (k < ALPHABET_SIZE - 1) fputc('\t', fcvec);
    }
    fputc('\n', fcvec);

    /* ── fill in real header and seek back to write it ── */
    uint64_t bwt_non_dollar = pos - n_term;
    /* The $ row occupies a slot, so blocks cover all N rows plus the query
       position N itself. Must match the aligner's check in load_data_mmap.c. */
    uint64_t num_blocks_expected = pos / OCC_INTERVAL + 1;

    hdr.magic                    = RS_OCC_MAGIC;
    hdr.version                  = RS_OCC_VERSION;
    hdr.occ_interval             = OCC_INTERVAL;
    hdr.alphabet_size            = ALPHABET_SIZE;
    hdr.reserved                 = 0;
    hdr.bwt_len_non_dollar       = bwt_non_dollar;
    hdr.bwt_len_total_with_dollar = pos;
    hdr.num_blocks               = blocks_written;
    hdr.sentinel_index           = sentinel;

    if (fseeko(fout, 0, SEEK_SET) != 0) {
        fprintf(stderr,"fseeko failed: %s\n", strerror(errno)); return 1;
    }
    if (fwrite(&hdr, sizeof(hdr), 1, fout) != 1) {
        fprintf(stderr,"Header rewrite failed\n"); return 1;
    }

    fclose(fbwt); fclose(fout); fclose(fcvec);
    free(rbuf); free(wbuf);

    fprintf(stderr, "\n\n[make_cp_occ_2step] done in %.1fs\n", elapsed);
    fprintf(stderr, "  BWT length (total)       = %llu\n",  (unsigned long long)pos);
    fprintf(stderr, "  BWT length (non-dollar)  = %llu\n",  (unsigned long long)bwt_non_dollar);
    fprintf(stderr, "  terminator (0x00) count  = %llu  (expected 1)\n",
            (unsigned long long)n_term);
    fprintf(stderr, "  sentinel index           = %lld  (0-based)\n", (long long)sentinel);
    fprintf(stderr, "  OCC_INTERVAL             = %d\n",   OCC_INTERVAL);
    fprintf(stderr, "  blocks expected          = %llu\n",
            (unsigned long long)num_blocks_expected);
    fprintf(stderr, "  blocks written           = %llu\n",
            (unsigned long long)blocks_written);
    fprintf(stderr, "  last block remainder     = %u positions%s\n",
            remainder, remainder ? "" : "  (terminal block)");
    fprintf(stderr, "  output file size         = %llu bytes  (%.2f GB)\n",
            (unsigned long long)(sizeof(rs_occ_full_header_t) +
                                 blocks_written * sizeof(cp_occ32_t)),
            (double)(sizeof(rs_occ_full_header_t) +
                     blocks_written * sizeof(cp_occ32_t)) / 1e9);

    /* symbol counts */
    fprintf(stderr, "\n  Base-16 symbol counts:\n");

    for (int k = 0; k < ALPHABET_SIZE; k++)
        fprintf(stderr, "    %s = %llu\n", sym_names[k],
                (unsigned long long)running[k]);

    /* C vector */
    fprintf(stderr, "\n  C vector (c_vec[0..15]):\n");
    for (int k = 0; k < ALPHABET_SIZE; k++)
        fprintf(stderr, "    c_vec[%2d] (%s) = %llu\n",
                k, sym_names[k], (unsigned long long)c_vec[k]);
    fprintf(stderr, "    c_vec[16] (N)  = set by aligner to BWT_SIZE_REFERENCE_SIZE\n");

    /* sanity checks */
    fprintf(stderr, "\n  Sanity checks:\n");
    int ok = 1;
    #define CHK(label, cond) \
        fprintf(stderr, "  %s  %s\n", (cond) ? "[OK]  " : "[WARN]", label); \
        if (!(cond)) ok = 0

    CHK("terminator count == 1",         n_term == 1);
    CHK("sentinel index >= 0",           sentinel >= 0);
    CHK("blocks written == expected",    blocks_written == num_blocks_expected);
    CHK("total = non_dollar + n_term",   pos == bwt_non_dollar + n_term);

    #undef CHK

    if (!ok) {
        fprintf(stderr, "\n  [ERROR] Sanity check failed; cp_occ_full.bin is not usable.\n");
        return EXIT_FAILURE;
    }
    fprintf(stderr, "\n  All checks passed. Ready to use.\n");

    return EXIT_SUCCESS;
}
