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
 * make_cp_occ_compact.c
 *
 * Builds the combined checkpoint+bitvector OCC file for the
 * 1-step RosaSeed FM-index seeding kernel.
 *
 * Input : raw binary BWT from gsufsort-64 + trim pipeline
 *         (ACGT bytes + exactly one 0x00 null terminator)
 * Output: cp_occ.bin  =  header  +  array of cp_occ_block_t
 *
 * ── Block layout (OCC_INTERVAL = 32, 32 bytes per block) ────────────────────
 *
 *   Byte  Type         Field         Description
 *   0     uint32_t[4]  cp_count      A/C/G/T counts BEFORE this block
 *   16    uint32_t[4]  one_hot       bitvectors: bit(31−off) set if BWT[blk×32+off]==base
 *   ─────────────────────────────────────────────────────────────────────────
 *   32 bytes total  (power-of-2, naturally 32-byte aligned)
 *
 * ── Terminator (0x00) handling ──────────────────────────────────────────────
 *
 *   The 0x00 null byte output by gsufsort acts as the BWT sentinel ($).
 *   It IS a BWT position:  i increments, blk/off advance normally.
 *   It is NOT an ACGT character:
 *     - no bit is set in any one_hot[0..3]
 *     - no cp_count is incremented
 *   Its position is stored in the header as sentinel_index.
 *   To reconstruct BWT[i] from the block: find which c has bit (31-off)
 *   set in one_hot[blk][c]; if none, BWT[i] = 0x00 (terminator).
 *
 * ── Occ(c, i) in the seeding kernel ─────────────────────────────────────────
 *
 *   blk  = i >> 5                           ( = i / 32 )
 *   off  = i & 31                           ( = i % 32 )
 *   Occ  = cp_count[blk][c]
 *        + popcount( one_hot[blk][c] >> (31 - off) )
 *
 *   For off=0:  counts only BWT[blk*32+0]
 *   For off=31: counts all 32 positions in the block
 *
 * ── Remainder (last partial block) ──────────────────────────────────────────
 *
 *   When n % 32 != 0 the last block has r < 32 valid positions.
 *   Bits for positions r..31 are zero.  Occ() is correct as-is because
 *   queries never reach beyond position n-1.
 *
 * ── Memory alignment for other sample distances ──────────────────────────────
 *
 *   dist  one_hot type     cp_count type    block bytes    aligned to
 *   16    uint16_t[4]      uint32_t[4]      24 → pad 32    32 B
 *   32    uint32_t[4]      uint32_t[4]      32             32 B  ← this code
 *   64    uint64_t[4]      uint64_t[4]      64             64 B  (cache-line perfect)
 *   128   uint64_t[4][2]   uint64_t[4]      96 → pad 128   128 B
 *
 * ── This file replaces both .occ and .bwt at query time ──────────────────────
 *
 *   The one_hot bitvectors encode the BWT completely.  BWT[i] can be
 *   recovered without the original .bwt file:
 *
 *     for c in 0..3:
 *       if one_hot[blk][c] & (1u << (31-off)):  return c
 *     return 0  // terminator
 *
 * Usage:
 *   ./make_cp_occ_compact <bwt_binary_file> <cp_occ_output.bin>
 *
 * Example:
 *   ./make_cp_occ_compact genome.txt.bwt genome.cp_occ.bin
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <time.h>

/* The RosaSeed-Compact aligner requires OCC_INTERVAL == 32 (one 32-bit
   bitvector per base per block); the loader rejects any other value. */
#define OCC_INTERVAL   32
#define ALPHABET_SIZE   4           /* A=0  C=1  G=2  T=3                      */
#define READ_BUF       (64 << 20)   /* 64 MB read buffer                       */

/* ── magic / version ─────────────────────────────────────────────────────── */
#define RS1_OCC_MAGIC    0x5253314F43433332ULL  /* "RS1OCC32" */
#define RS1_OCC_VERSION  1

/* ── block struct (32 bytes, OCC_INTERVAL = 32) ───────────────────────────── */
typedef struct {
    uint32_t cp_count[4];   /* 16 bytes: cumulative A/C/G/T counts before blk   */
    uint32_t one_hot[4];    /* 16 bytes: bitvectors, bit(31-off) = BWT pos off   */
} cp_occ_block_t;            /* 32 bytes total                                   */

/* ── file header ──────────────────────────────────────────────────────────── */
typedef struct {
    uint64_t magic;                  /* RS1_OCC_MAGIC                            */
    uint32_t version;                /* RS1_OCC_VERSION                          */
    uint32_t occ_interval;           /* OCC_INTERVAL (e.g. 32)                   */
    uint32_t alphabet_size;          /* 4                                        */
    uint32_t reserved;
    uint64_t bwt_len_no_term;        /* total ACGT characters (excludes 0x00)    */
    uint64_t bwt_len_total;          /* total BWT length (ACGT + terminator)     */
    uint64_t num_blocks;             /* number of cp_occ_block_t entries          */
    int64_t  sentinel_index;         /* BWT position of 0x00 terminator           */
    uint64_t c_array[4];             /* C[A] C[C] C[G] C[T] for LF-mapping       */
} rs1_occ_header_t;

/* ── base lookup table ────────────────────────────────────────────────────── */
/* lut[byte] = 0(A) 1(C) 2(G) 3(T)  or  -1 (terminator / unexpected)        */
static int8_t lut[256];

static void init_lut(void) {
    memset(lut, -1, sizeof(lut));
    lut['A'] = lut['a'] = 0;
    lut['C'] = lut['c'] = 1;
    lut['G'] = lut['g'] = 2;
    lut['T'] = lut['t'] = 3;
}

/* ── main ─────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr,
            "\nUsage: %s <bwt_binary_file> <cp_occ_output.bin>\n\n"
            "  bwt_binary_file   raw binary BWT from gsufsort-64 + trim\n"
            "  cp_occ_output     combined checkpoint+bitvector binary\n\n"
            "  OCC_INTERVAL = %d  (recompile to change)\n"
            "  Block size   = %zu bytes\n\n",
            argv[0], OCC_INTERVAL, sizeof(cp_occ_block_t));
        return EXIT_FAILURE;
    }

    init_lut();

    FILE *fbwt = fopen(argv[1], "rb");
    if (!fbwt) {
        fprintf(stderr, "Cannot open BWT '%s': %s\n", argv[1], strerror(errno));
        return EXIT_FAILURE;
    }
    FILE *fout = fopen(argv[2], "wb");
    if (!fout) {
        fprintf(stderr, "Cannot open output '%s': %s\n", argv[2], strerror(errno));
        fclose(fbwt);
        return EXIT_FAILURE;
    }

    uint8_t *rbuf = (uint8_t *)malloc(READ_BUF);
    if (!rbuf) { fprintf(stderr, "malloc failed\n"); return EXIT_FAILURE; }

    /* write placeholder header, we'll seek back and overwrite it at the end */
    rs1_occ_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    if (fwrite(&hdr, sizeof(hdr), 1, fout) != 1) {
        fprintf(stderr, "Header write failed\n"); return EXIT_FAILURE;
    }

    /* ── state ─────────────────────────────────────────────────────────── */
    uint64_t running[4] = {0, 0, 0, 0};  /* cumulative A/C/G/T counts       */
    cp_occ_block_t blk;                   /* block being assembled            */
    memset(&blk, 0, sizeof(blk));

    uint64_t pos      = 0;                /* absolute BWT position            */
    uint64_t n_blocks = 0;                /* blocks written                   */
    int64_t  sentinel = -1;               /* position of 0x00 terminator      */
    uint64_t n_term   = 0;                /* count of terminator bytes seen   */

    clock_t t0 = clock();

    /*
     * ── single-pass main loop ────────────────────────────────────────────
     *
     * For each BWT byte at absolute position pos:
     *
     *   1. If pos is a block boundary (off == 0):
     *        snapshot running counts into cp_count  (BEFORE this block)
     *        reset one_hot to zero
     *        → block 0: cp_count = [0,0,0,0] always
     *
     *   2. Classify the byte:
     *        ACGT → set bit (31−off) in one_hot[base],  increment running[base]
     *        0x00 → record sentinel_index, increment n_term
     *               NO bit set, NO running count incremented
     *        other → should not appear; silently skipped
     *
     *   3. If pos is the last position in a block (off == OCC_INTERVAL−1):
     *        write the completed block to fout
     */

    size_t nread;
    while ((nread = fread(rbuf, 1, READ_BUF, fbwt)) > 0) {
        for (size_t i = 0; i < nread; i++, pos++) {

            uint32_t off = (uint32_t)(pos & (OCC_INTERVAL - 1));  /* 0 .. 31 */

            /* ── block boundary: snapshot counts, clear one_hot ── */
            if (off == 0) {
                blk.cp_count[0] = (uint32_t)running[0];
                blk.cp_count[1] = (uint32_t)running[1];
                blk.cp_count[2] = (uint32_t)running[2];
                blk.cp_count[3] = (uint32_t)running[3];
                blk.one_hot[0]  = 0;
                blk.one_hot[1]  = 0;
                blk.one_hot[2]  = 0;
                blk.one_hot[3]  = 0;
            }

            int8_t base = lut[rbuf[i]];

            if (base >= 0) {
                /* ACGT: set MSB-first bit for this position */
                blk.one_hot[base] |= (1u << (31 - off));
                running[base]++;
            } else if (rbuf[i] == 0x00) {
                /* 0x00 terminator: record position, set no bit, count no base */
                if (sentinel == -1) sentinel = (int64_t)pos;
                n_term++;
            }
            /* any other byte: should not occur, silently skip */

            /* ── end of block: write it ── */
            if (off == (uint32_t)(OCC_INTERVAL - 1)) {
                if (fwrite(&blk, sizeof(cp_occ_block_t), 1, fout) != 1) {
                    fprintf(stderr, "Write error at block %llu\n",
                            (unsigned long long)n_blocks);
                    return EXIT_FAILURE;
                }
                n_blocks++;
            }
        }

        /* progress every ~500M characters */
        if (pos % 500000000ULL < (uint64_t)nread) {
            double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
            fprintf(stderr, "\r[make_cp_occ_compact] %7lluM (%.0fs) A=%llu C=%llu G=%llu T=%llu",
                    (unsigned long long)(pos/1000000ULL), el,
                    (unsigned long long)running[0], (unsigned long long)running[1],
                    (unsigned long long)running[2], (unsigned long long)running[3]);
        }
    }

    /* ── flush final block: always emit N/32 + 1 blocks ──
     * The aligner reads block (pos >> 5) for every pos in 0..N inclusive.
     *   remainder != 0: the partially filled block (unused bits are 0).
     *   remainder == 0: every block is full, so emit a terminal block with
     *                   the final cumulative counts and no bits set.
     */
    uint32_t remainder = (uint32_t)(pos & (OCC_INTERVAL - 1));
    if (remainder == 0) {
        for (int b = 0; b < ALPHABET_SIZE; b++) {
            blk.cp_count[b] = (uint32_t)running[b];
            blk.one_hot[b]  = 0;
        }
    }
    if (fwrite(&blk, sizeof(cp_occ_block_t), 1, fout) != 1) {
        fprintf(stderr, "Final block write failed\n");
        return EXIT_FAILURE;
    }
    n_blocks++;

    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;

    /* ── C[] array: number of BWT chars strictly less than each base ── */
    /* alphabet order:  0x00($) < A < C < G < T                         */
    uint64_t C[4];
    C[0] = (uint64_t)n_term;               /* C[A]: only the terminator(s)    */
    C[1] = C[0] + running[0];              /* C[C]                            */
    C[2] = C[1] + running[1];              /* C[G]                            */
    C[3] = C[2] + running[2];              /* C[T]                            */

    /* ── fill in real header and seek back to overwrite placeholder ── */
    hdr.magic            = RS1_OCC_MAGIC;
    hdr.version          = RS1_OCC_VERSION;
    hdr.occ_interval     = OCC_INTERVAL;
    hdr.alphabet_size    = ALPHABET_SIZE;
    hdr.bwt_len_no_term  = pos - n_term;
    hdr.bwt_len_total    = pos;
    hdr.num_blocks       = n_blocks;
    hdr.sentinel_index   = sentinel;
    hdr.c_array[0]       = C[0];
    hdr.c_array[1]       = C[1];
    hdr.c_array[2]       = C[2];
    hdr.c_array[3]       = C[3];

    if (fseeko(fout, 0, SEEK_SET) != 0) {
        fprintf(stderr, "fseeko failed: %s\n", strerror(errno));
        return EXIT_FAILURE;
    }
    if (fwrite(&hdr, sizeof(hdr), 1, fout) != 1) {
        fprintf(stderr, "Header rewrite failed\n");
        return EXIT_FAILURE;
    }

    fclose(fbwt); fclose(fout); free(rbuf);

    /* ── summary ─────────────────────────────────────────────────────── */
    uint64_t out_bytes = sizeof(rs1_occ_header_t) + n_blocks * sizeof(cp_occ_block_t);

    fprintf(stderr, "\n\n[make_cp_occ_compact] done in %.1fs\n", elapsed);
    fprintf(stderr, "  OCC_INTERVAL         = %d\n",       OCC_INTERVAL);
    fprintf(stderr, "  block size           = %zu bytes\n", sizeof(cp_occ_block_t));
    fprintf(stderr, "  BWT length (total)   = %llu\n",     (unsigned long long)pos);
    fprintf(stderr, "  BWT length (no term) = %llu\n",     (unsigned long long)(pos - n_term));
    fprintf(stderr, "  terminator count     = %llu (expected 1)\n", (unsigned long long)n_term);
    fprintf(stderr, "  sentinel index       = %lld\n",     (long long)sentinel);
    fprintf(stderr, "  blocks written       = %llu\n",     (unsigned long long)n_blocks);
    fprintf(stderr, "  last block remainder = %u positions%s\n",
            remainder, remainder ? "" : "  (terminal block)");
    fprintf(stderr, "  output file size     = %llu bytes  (%.2f GB)\n",
            (unsigned long long)out_bytes, (double)out_bytes / 1073741824.0 );

    fprintf(stderr, "\n  A=%-14llu C=%-14llu G=%-14llu T=%llu\n",
            (unsigned long long)running[0], (unsigned long long)running[1],
            (unsigned long long)running[2], (unsigned long long)running[3]);
    fprintf(stderr, "  C[A]=%-12llu C[C]=%-12llu C[G]=%-12llu C[T]=%llu\n",
            (unsigned long long)C[0], (unsigned long long)C[1],
            (unsigned long long)C[2], (unsigned long long)C[3]);

    /* sanity checks, any failure means the index is unusable */
    int ok = 1;
    #define CHK(label, cond) do { \
        fprintf(stderr, "  %s  %s\n", (cond) ? "[OK]  " : "[FAIL]", label); \
        if (!(cond)) ok = 0; } while (0)
    fprintf(stderr, "\n  Sanity checks:\n");
    CHK("terminator count == 1",             n_term == 1);
    CHK("sentinel index found",              sentinel >= 0);
    CHK("A count == T count (fwd+RC)",       running[0] == running[3]);
    CHK("C count == G count (fwd+RC)",       running[1] == running[2]);
    CHK("blocks written == N/32 + 1",        n_blocks == pos / OCC_INTERVAL + 1);
    CHK("per-base counts fit in uint32_t",   running[0] < 4294967295ULL && running[1] < 4294967295ULL);
    #undef CHK
    if (!ok) {
        fprintf(stderr, "\n  [ERROR] Sanity check failed; %s is not usable.\n", argv[2]);
        return EXIT_FAILURE;
    }
    fprintf(stderr, "\n  All checks passed.\n");

    return EXIT_SUCCESS;
}
