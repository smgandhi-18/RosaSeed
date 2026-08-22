/*
 * make_jumptable_2step.c
 *
 * Generates binary jump tables for 14-mer, 15-mer, and/or 16-mer seeding
 * for the 2-step RosaSeed FM-index.
 *
 * ── Key differences from 1-step ──────────────────────────────────────────
 *
 *   Alphabet: 16 base-16 symbols (0-9, A-F encoding NT pairs)
 *   cp_occ block: 128 bytes (cp_count[16] + one_hot[16], vs 32 bytes in 1-step)
 *   C[] vector: loaded from c_vector.txt (NOT embedded in cp_occ header)
 *   N = bwt_len_total_with_dollar = 6,234,551,001
 *
 * ── Three k-mer encoding schemes ─────────────────────────────────────────
 *
 *   14-mer (7 full base-16 pairs):
 *     Index: 0 .. 16^7-1  (268,435,456 entries)
 *     idx = sym[0]*16^6 + ... + sym[6]*16^0
 *     sym[j] = (idx >> 4*j) & 0xF    (j=0 = rightmost symbol)
 *     7 backward steps, each consuming one base-16 symbol
 *
 *   15-mer (7 pairs + 1 rightmost single NT):
 *     Index: 0 .. 4^15-1  (1,073,741,824 entries, base-4 nucleotide ordering)
 *     b[i] = (idx >> 2*(14-i)) & 3   (b[14] = idx & 3 = rightmost NT)
 *
 *     Step 1 — rightmost single nucleotide b[14]:
 *       lo = C[4*b14]
 *       hi = (4*b14+4 < 16) ? C[4*b14+4] : N
 *       Rationale: b[14] selects a GROUP of 4 base-16 symbols (those starting
 *       with b[14]). This gives the SA range for all suffixes whose b[14]
 *       matches, without committing to b[13] yet.
 *
 *     Steps 2-8 — 7 pairs right to left: (b[12],b[13]),...,(b[0],b[1]):
 *       For pair step j = 0..6:
 *         b_right = (idx >> (2 + 4*j)) & 3
 *         b_left  = (idx >> (4 + 4*j)) & 3
 *         sym     = 4*b_left + b_right
 *         lo = C[sym] + rank_excl(sym, lo)
 *         hi = C[sym] + rank_excl(sym, hi)
 *
 *   16-mer (8 full base-16 pairs):
 *     Index: 0 .. 4^16-1  (4,294,967,296 entries, base-4 nucleotide ordering)
 *     b[i] = (idx >> 2*(15-i)) & 3   (b[15] = idx & 3 = rightmost NT)
 *
 *     8 pairs right to left: (b[14],b[15]), ..., (b[0],b[1]):
 *       For pair step j = 0..7:
 *         b_right = (idx >> (4*j)) & 3
 *         b_left  = (idx >> (4*j+2)) & 3
 *         sym     = 4*b_left + b_right
 *         lo = C[sym] + rank_excl(sym, lo)
 *         hi = C[sym] + rank_excl(sym, hi)
 *
 * ── rank_excl ─────────────────────────────────────────────────────────────
 *
 *   rank_excl(sym, i) = count of sym in BWT[0..i-1]  (exclusive)
 *
 *   if i==0: return 0
 *   blk = (i-1) >> 5
 *   off = (i-1) & 31
 *   return cp_count[blk][sym] + popcount(one_hot[blk][sym] >> (31-off))
 *
 *   Backward search:
 *     lo = C[sym] + rank_excl(sym, lo)
 *     hi = C[sym] + rank_excl(sym, hi)
 *
 * ── Initial pointers ──────────────────────────────────────────────────────
 *
 *   lo = 0 (inclusive)
 *   hi = N = bwt_len_total = 6,234,551,001 (exclusive)
 *   These are 0-based row indices into the BWT/SA.
 *
 * ── Jump pointer encoding (identical to 1-step) ───────────────────────────
 *
 *   diff = hi - lo
 *   diff==0 → jp = 0                                  (no mapping)
 *   diff==1 → jp = (1ULL<<63) | (SA[lo] & 0x1FFFFFFFF) (unique)
 *   diff>1  → jp = lo | (diff << 34)                 (non-unique)
 *
 * ── Memory requirements ───────────────────────────────────────────────────
 *
 *   cp_occ_full.bin  : ~24.9 GB  (194,829,719 × 128 bytes)
 *   sa_ls_word_cf1   : ~23.2 GB  (6,234,551,001 × 4 bytes)
 *   sa_ms_byte_cf1   :  ~5.8 GB  (6,234,551,001 × 1 byte)
 *   ────────────────────────────────────────────────
 *   Total            : ~54 GB
 *
 * ── Output file sizes ─────────────────────────────────────────────────────
 *
 *   jumptable_14nt.bin :  16^7 × 8 =  2.0 GiB
 *   jumptable_15nt.bin :  4^15 × 8 =  8.0 GiB
 *   jumptable_16nt.bin :  4^16 × 8 = 32.0 GiB
 *
 * Usage:
 *   ./make_jumptable_2step <cp_occ_full.bin> <sa_ls_cf1.bin> <sa_ms_cf1.bin>
 *                          <c_vector.txt> <out_dir> [k]
 *   k = 14, 15, 16, or 0 for all three (default)
 *
 * Examples:
 *   ./make_jumptable_2step cp_occ_full.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin
 *                          c_vector.txt . 15
 *   ./make_jumptable_2step cp_occ_full.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin
 *                          c_vector.txt /data/index
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <time.h>

/* ── constants ── */
#define ALPHABET_SIZE  16
#define OCC_INTERVAL   32

/* ── cp_occ structures (must match make_cp_occ_2step.c) ── */
#define RS_OCC_MAGIC   0x52534F434346554CULL   /* "RSOCCFUL" */

typedef struct {
    uint32_t cp_count[16];   /* 64 bytes */
    uint32_t one_hot[16];    /* 64 bytes */
} cp_occ32_t;                /* 128 bytes */

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
} rs_occ_full_header_t;

/* ── global data ── */
static cp_occ32_t *blocks = NULL;
static uint64_t    C[ALPHABET_SIZE + 1];   /* C[0..15] from file, C[16]=N */
static uint64_t    N;                       /* bwt_len_total_with_dollar */
static uint64_t    Nblocks;

static uint32_t   *sa_ls = NULL;
static uint8_t    *sa_ms = NULL;

/* ── rank_excl ────────────────────────────────────────────────────────────
 *
 *  rank_excl(sym, i) = count of sym in BWT[0..i-1]  (exclusive upper bound)
 *
 *  Block 0 always has cp_count = [0,...,0], so rank_excl(*,0) = 0. ✓
 *
 *  This is what the FM-index backward search needs:
 *    lo = C[sym] + rank_excl(sym, lo)
 *    hi = C[sym] + rank_excl(sym, hi)
 */
static inline uint64_t rank_excl(int sym, uint64_t i)
{
    if (i == 0) return 0ULL;
    uint64_t pos = i - 1;
    uint64_t blk = pos >> 5;
    uint32_t off = (uint32_t)(pos & 31u);
    return (uint64_t)blocks[blk].cp_count[sym]
         + (uint64_t)__builtin_popcount(blocks[blk].one_hot[sym] >> (31u - off));
}

/* ── SA lookup (CF=1 direct) ── */
static inline uint64_t get_sa(uint64_t row)
{
    return ((uint64_t)sa_ms[row] << 32) | (uint64_t)sa_ls[row];
}

/* ── jump-pointer field widths ───────────────────────────────────────────
 * Packed layout, must match extract_jump_bounds() in helper_functions.c:
 *   bit  63     : unique flag
 *   bits 62..34 : diff  (29 bits)
 *   bit  33     : unused (reserved)
 *   bits 32..0  : lo, or SA position when unique (33 bits)
 */
#define JT_LO_BITS    33
#define JT_DIFF_BITS  29
#define JT_LO_MAX     ((1ULL << JT_LO_BITS)   - 1ULL)   /* 8,589,934,591 */
#define JT_DIFF_MAX   ((1ULL << JT_DIFF_BITS) - 1ULL)   /*   536,870,911 */

static uint64_t g_max_diff_seen = 0;
static uint64_t g_max_lo_seen   = 0;

static inline uint64_t encode_jp(uint64_t lo, uint64_t hi)
{
    uint64_t diff = hi - lo;

    if (diff == 0) return 0ULL;

    if (diff > g_max_diff_seen) g_max_diff_seen = diff;
    if (lo   > g_max_lo_seen)   g_max_lo_seen   = lo;

    if (diff == 1) {
        uint64_t sa = get_sa(lo);
        if (sa > JT_LO_MAX) {
            fprintf(stderr,
                "\n[FATAL] SA value %llu exceeds the %d-bit field (max %llu).\n"
                "        BWT length is >= 2^%d; this reference is too large.\n"
                "        See README, 'Genome size limits'.\n",
                (unsigned long long)sa, JT_LO_BITS,
                (unsigned long long)JT_LO_MAX, JT_LO_BITS);
            exit(EXIT_FAILURE);
        }
        return (1ULL << 63) | (sa & 0x1FFFFFFFFULL);
    }

    if (diff > JT_DIFF_MAX) {
        fprintf(stderr,
            "\n[FATAL] interval width %llu exceeds the %d-bit diff field (max %llu).\n"
            "        Encoding it would set bit 63, and the aligner would decode\n"
            "        this entry as a UNIQUE seed at a bogus position.\n"
            "        Use a longer k-mer (-j 15 / -j 16) or widen the diff field.\n",
            (unsigned long long)diff, JT_DIFF_BITS,
            (unsigned long long)JT_DIFF_MAX);
        exit(EXIT_FAILURE);
    }

    if (lo > JT_LO_MAX) {
        fprintf(stderr,
            "\n[FATAL] interval start %llu exceeds the %d-bit lo field (max %llu).\n"
            "        BWT length is >= 2^%d; this reference is too large.\n"
            "        See README, 'Genome size limits'.\n",
            (unsigned long long)lo, JT_LO_BITS,
            (unsigned long long)JT_LO_MAX, JT_LO_BITS);
        exit(EXIT_FAILURE);
    }

    return lo | (diff << 34);
}

/* ── 14-mer table ────────────────────────────────────────────────────────
 *
 *  16^7 = 268,435,456 entries
 *  Each entry corresponds to a 7-symbol base-16 k-mer.
 *
 *  Index encoding: idx = sym[0]*16^6 + ... + sym[6]*16^0
 *  sym[j] = (idx >> (4*j)) & 0xF   (j=0 = rightmost symbol)
 *
 *  Backward search: process sym[6] (leftmost) down to sym[0] (rightmost)
 *  BUT: we enumerate idx in order 0..16^7-1, and within the loop we
 *  need to process the rightmost symbol first. Since sym[j] = (idx>>(4j))&0xF
 *  gives j=0 as rightmost, we iterate j from 0 up to 6.
 *
 *  Wait — backward search goes RIGHT to LEFT:
 *  idx=0 → sym[0..6] all = 0 (rightmost to leftmost all '0')
 *  First step of backward search: sym[6] (leftmost nucleotide in time
 *  but rightmost in BWT order? No...)
 *
 *  Actually for a 7-mer base-16: the query is sym[0]sym[1]...sym[6]
 *  Backward search processes from right to left: sym[6] first, sym[0] last.
 *  sym[6] = (idx >> 0) & 0xF  (j=0 in our formula gives rightmost = sym[6])
 *
 *  So iterating j from 0 to 6 processes sym[6]→sym[0] = right to left. ✓
 */
static void make_table_14(const char *out_path)
{
    const uint64_t total = 268435456ULL;   /* 16^7 */
    g_max_diff_seen = 0;
    g_max_lo_seen   = 0;

    FILE *fout = fopen(out_path, "wb");
    if (!fout) { perror(out_path); exit(EXIT_FAILURE); }

    const size_t WBUF = 1u << 20;
    uint64_t *buf = (uint64_t *)malloc(WBUF * sizeof(uint64_t));
    if (!buf) { perror("malloc"); exit(EXIT_FAILURE); }

    uint64_t n_buf = 0, written = 0;
    uint64_t cnt_unique = 0, cnt_multi = 0, cnt_none = 0;
    clock_t t0 = clock();

    for (uint64_t idx = 0; idx < total; idx++) {

        uint64_t lo = 0, hi = N;

        /*
         * Process 7 base-16 symbols right to left.
         * sym[j] = (idx >> (4*j)) & 0xF  gives rightmost (j=0) first.
         * j=0: sym[6] of the 7-mer  ← rightmost symbol, processed first ✓
         * j=6: sym[0] of the 7-mer  ← leftmost symbol, processed last  ✓
         */
        for (int j = 0; j < 7 && lo < hi; j++) {
            int sym = (int)((idx >> (4u * j)) & 0xFu);
            lo = C[sym] + rank_excl(sym, lo);
            hi = C[sym] + rank_excl(sym, hi);
        }

        uint64_t diff = hi - lo;
        if      (diff == 0) cnt_none++;
        else if (diff == 1) cnt_unique++;
        else                cnt_multi++;

        buf[n_buf++] = encode_jp(lo, hi);
        if (n_buf == WBUF) {
            fwrite(buf, sizeof(uint64_t), n_buf, fout);
            written += n_buf; n_buf = 0;
            if ((written & ((64u << 20) - 1)) == 0) {
                double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
                fprintf(stderr, "\r  k=14  %5.1f%%  (%.0fs)  "
                                "unique=%.1f%%  none=%.1f%%",
                        100.0*written/total, el,
                        100.0*cnt_unique/written,
                        100.0*cnt_none/written);
            }
        }
    }
    if (n_buf > 0) { fwrite(buf, sizeof(uint64_t), n_buf, fout); written += n_buf; }
    fclose(fout); free(buf);

    double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
    fprintf(stderr, "\r  k=14  100.0%%  done in %.1fs                            \n", el);
    fprintf(stderr, "  entries  : %llu\n",  (unsigned long long)written);
    fprintf(stderr, "  unique   : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_unique, 100.0*cnt_unique/(double)total);
    fprintf(stderr, "  multi    : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_multi,  100.0*cnt_multi/(double)total);
    fprintf(stderr, "  none     : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_none,   100.0*cnt_none/(double)total);
    fprintf(stderr, "  max diff : %llu  (%.2f%% of %llu limit)\n",
            (unsigned long long)g_max_diff_seen,
            100.0 * (double)g_max_diff_seen / (double)JT_DIFF_MAX,
            (unsigned long long)JT_DIFF_MAX);
    fprintf(stderr, "  max lo   : %llu  (%.2f%% of %llu limit)\n",
            (unsigned long long)g_max_lo_seen,
            100.0 * (double)g_max_lo_seen / (double)JT_LO_MAX,
            (unsigned long long)JT_LO_MAX);        
    fprintf(stderr, "  output   : %s  (%.2f GiB)\n",
            out_path, (double)(written*8)/1073741824.0);
}

/* ── 15-mer table ────────────────────────────────────────────────────────
 *
 *  4^15 = 1,073,741,824 entries
 *  Index: base-4 nucleotide ordering (A=0,C=1,G=2,T=3)
 *  idx = b[0]*4^14 + ... + b[14]*4^0
 *  b[i] = (idx >> (2*(14-i))) & 3    (b[14] = idx & 3 = rightmost NT)
 *
 *  Backward search:
 *  Step 1: rightmost single nucleotide b[14] = idx & 3
 *    sym_lo = 4*b14   (group of 4 base-16 symbols starting with NT b14)
 *    lo = C[sym_lo]
 *    hi = (sym_lo+4 < 16) ? C[sym_lo+4] : N
 *
 *    Rationale: b[14] is a single nucleotide (half a base-16 pair). The
 *    2-step reference encodes NT pair (left,right) → sym = 4*left+right.
 *    So all base-16 symbols with LEFT nucleotide = b14 are {4*b14 .. 4*b14+3}.
 *    The SA range covering all of these is [C[4*b14], C[4*b14+4]).
 *
 *  Steps 2-8: 7 pairs (b[12],b[13]) through (b[0],b[1]) right to left:
 *    For pair step j = 0..6:
 *      b_right = (idx >> (2 + 4*j)) & 3  = b[13-2j+1] = b[14-2j-1] hmm
 *
 *    Let me re-derive:
 *    pair step j=0: processes pair at positions (12,13) in the 15-mer
 *      b_right = b[13] = (idx >> (2*(14-13))) & 3 = (idx >> 2) & 3
 *      b_left  = b[12] = (idx >> (2*(14-12))) & 3 = (idx >> 4) & 3
 *    pair step j=1: positions (10,11)
 *      b_right = b[11] = (idx >> 6) & 3
 *      b_left  = b[10] = (idx >> 8) & 3
 *    pair step j: positions (12-2j, 13-2j)
 *      b_right = (idx >> (2 + 4*j)) & 3
 *      b_left  = (idx >> (4 + 4*j)) & 3
 *    pair step j=6: positions (0,1)
 *      b_right = (idx >> 26) & 3 = b[1]
 *      b_left  = (idx >> 28) & 3 = b[0]
 *
 *    sym = 4*b_left + b_right
 */
static void make_table_15(const char *out_path)
{
    const uint64_t total = 1073741824ULL;   /* 4^15 */
    g_max_diff_seen = 0;
    g_max_lo_seen   = 0;

    FILE *fout = fopen(out_path, "wb");
    if (!fout) { perror(out_path); exit(EXIT_FAILURE); }

    const size_t WBUF = 1u << 20;
    uint64_t *buf = (uint64_t *)malloc(WBUF * sizeof(uint64_t));
    if (!buf) { perror("malloc"); exit(EXIT_FAILURE); }

    uint64_t n_buf = 0, written = 0;
    uint64_t cnt_unique = 0, cnt_multi = 0, cnt_none = 0;
    clock_t t0 = clock();

    for (uint64_t idx = 0; idx < total; idx++) {

        /* Step 1: rightmost single nucleotide b[14] */
        int b14 = (int)(idx & 3u);
        int sym_lo = 4 * b14;
        uint64_t lo = C[sym_lo];
        uint64_t hi = (sym_lo + 4 < ALPHABET_SIZE) ? C[sym_lo + 4] : N;
        /*
         * b14=0(A): lo=C[0],  hi=C[4]   ← all AA,AC,AG,AT suffixes
         * b14=1(C): lo=C[4],  hi=C[8]   ← all CA,CC,CG,CT suffixes
         * b14=2(G): lo=C[8],  hi=C[12]  ← all GA,GC,GG,GT suffixes
         * b14=3(T): lo=C[12], hi=N      ← all TA,TC,TG,TT suffixes
         */

        /* Steps 2-8: 7 pairs right to left */
        for (int j = 0; j < 7 && lo < hi; j++) {
            int b_right = (int)((idx >> (2u + 4u*j)) & 3u);
            int b_left  = (int)((idx >> (4u + 4u*j)) & 3u);
            int sym = 4 * b_left + b_right;
            lo = C[sym] + rank_excl(sym, lo);
            hi = C[sym] + rank_excl(sym, hi);
        }

        uint64_t diff = hi - lo;
        if      (diff == 0) cnt_none++;
        else if (diff == 1) cnt_unique++;
        else                cnt_multi++;

        buf[n_buf++] = encode_jp(lo, hi);
        if (n_buf == WBUF) {
            fwrite(buf, sizeof(uint64_t), n_buf, fout);
            written += n_buf; n_buf = 0;
            if ((written & ((64u << 20) - 1)) == 0) {
                double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
                fprintf(stderr, "\r  k=15  %5.1f%%  (%.0fs)  "
                                "unique=%.1f%%  none=%.1f%%",
                        100.0*written/total, el,
                        100.0*cnt_unique/written,
                        100.0*cnt_none/written);
            }
        }
    }
    if (n_buf > 0) { fwrite(buf, sizeof(uint64_t), n_buf, fout); written += n_buf; }
    fclose(fout); free(buf);

    double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
    fprintf(stderr, "\r  k=15  100.0%%  done in %.1fs                            \n", el);
    fprintf(stderr, "  entries  : %llu\n",  (unsigned long long)written);
    fprintf(stderr, "  unique   : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_unique, 100.0*cnt_unique/(double)total);
    fprintf(stderr, "  multi    : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_multi,  100.0*cnt_multi/(double)total);
    fprintf(stderr, "  none     : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_none,   100.0*cnt_none/(double)total);
    fprintf(stderr, "  max diff : %llu  (%.2f%% of %llu limit)\n",
            (unsigned long long)g_max_diff_seen,
            100.0 * (double)g_max_diff_seen / (double)JT_DIFF_MAX,
            (unsigned long long)JT_DIFF_MAX);
    fprintf(stderr, "  max lo   : %llu  (%.2f%% of %llu limit)\n",
            (unsigned long long)g_max_lo_seen,
            100.0 * (double)g_max_lo_seen / (double)JT_LO_MAX,
            (unsigned long long)JT_LO_MAX);        
    fprintf(stderr, "  output   : %s  (%.2f GiB)\n",
            out_path, (double)(written*8)/1073741824.0);
}

/* ── 16-mer table ────────────────────────────────────────────────────────
 *
 *  4^16 = 4,294,967,296 entries
 *  Index: base-4 nucleotide ordering
 *  idx = b[0]*4^15 + ... + b[15]*4^0
 *  b[15] = idx & 3 = rightmost NT
 *
 *  8 pairs right to left: (b[14],b[15]), (b[12],b[13]), ..., (b[0],b[1])
 *  For pair step j = 0..7:
 *    b_right = (idx >> (4*j)) & 3   = b[15-2j] (rightmost of pair)
 *    b_left  = (idx >> (4*j+2)) & 3 = b[14-2j] (leftmost of pair)
 *    sym = 4*b_left + b_right
 *
 *  Verify j=0: b_right=(idx>>0)&3=b[15], b_left=(idx>>2)&3=b[14] ✓
 *  Verify j=7: b_right=(idx>>28)&3=b[1], b_left=(idx>>30)&3=b[0] ✓
 *
 *  lo=0, hi=N (both pairs fully determined, no single-NT init needed)
 */
static void make_table_16(const char *out_path)
{
    const uint64_t total = 4294967296ULL;   /* 4^16 */
    g_max_diff_seen = 0;
    g_max_lo_seen   = 0;

    FILE *fout = fopen(out_path, "wb");
    if (!fout) { perror(out_path); exit(EXIT_FAILURE); }

    const size_t WBUF = 1u << 20;
    uint64_t *buf = (uint64_t *)malloc(WBUF * sizeof(uint64_t));
    if (!buf) { perror("malloc"); exit(EXIT_FAILURE); }

    uint64_t n_buf = 0, written = 0;
    uint64_t cnt_unique = 0, cnt_multi = 0, cnt_none = 0;
    clock_t t0 = clock();

    for (uint64_t idx = 0; idx < total; idx++) {

        uint64_t lo = 0, hi = N;

        /* 8 full pairs right to left */
        for (int j = 0; j < 8 && lo < hi; j++) {
            int b_right = (int)((idx >> (4u * j))     & 3u);
            int b_left  = (int)((idx >> (4u * j + 2u)) & 3u);
            int sym = 4 * b_left + b_right;
            lo = C[sym] + rank_excl(sym, lo);
            hi = C[sym] + rank_excl(sym, hi);
        }

        uint64_t diff = hi - lo;
        if      (diff == 0) cnt_none++;
        else if (diff == 1) cnt_unique++;
        else                cnt_multi++;

        buf[n_buf++] = encode_jp(lo, hi);
        if (n_buf == WBUF) {
            fwrite(buf, sizeof(uint64_t), n_buf, fout);
            written += n_buf; n_buf = 0;
            if ((written & ((64u << 20) - 1)) == 0) {
                double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
                fprintf(stderr, "\r  k=16  %5.1f%%  (%.0fs)  "
                                "unique=%.1f%%  none=%.1f%%",
                        100.0*written/total, el,
                        100.0*cnt_unique/written,
                        100.0*cnt_none/written);
            }
        }
    }
    if (n_buf > 0) { fwrite(buf, sizeof(uint64_t), n_buf, fout); written += n_buf; }
    fclose(fout); free(buf);

    double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
    fprintf(stderr, "\r  k=16  100.0%%  done in %.1fs                            \n", el);
    fprintf(stderr, "  entries  : %llu\n",  (unsigned long long)written);
    fprintf(stderr, "  unique   : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_unique, 100.0*cnt_unique/(double)total);
    fprintf(stderr, "  multi    : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_multi,  100.0*cnt_multi/(double)total);
    fprintf(stderr, "  none     : %llu  (%.2f%%)\n",
            (unsigned long long)cnt_none,   100.0*cnt_none/(double)total);
    fprintf(stderr, "  max diff : %llu  (%.2f%% of %llu limit)\n",
            (unsigned long long)g_max_diff_seen,
            100.0 * (double)g_max_diff_seen / (double)JT_DIFF_MAX,
            (unsigned long long)JT_DIFF_MAX);
    fprintf(stderr, "  max lo   : %llu  (%.2f%% of %llu limit)\n",
            (unsigned long long)g_max_lo_seen,
            100.0 * (double)g_max_lo_seen / (double)JT_LO_MAX,
            (unsigned long long)JT_LO_MAX);
    fprintf(stderr, "  output   : %s  (%.2f GiB)\n",
            out_path, (double)(written*8)/1073741824.0);
}

/* ── main ── */
int main(int argc, char *argv[])
{
    if (argc < 6 || argc > 7) {
        fprintf(stderr,
            "\nUsage: %s <cp_occ_full.bin> <sa_ls_cf1.bin> <sa_ms_cf1.bin>"
            " <c_vector.txt> <out_dir> [k]\n\n"
            "  cp_occ_full.bin   2-step OCC file (24.9 GB, 16-alphabet blocks)\n"
            "  sa_ls_cf1.bin     lower 32 bits of SA, CF=1\n"
            "  sa_ms_cf1.bin     bits 32-39 of SA, CF=1\n"
            "  c_vector.txt      16 tab-separated C[] prefix sums\n"
            "  out_dir           output directory\n"
            "  k                 14, 15, 16, or 0/omit for all three\n\n"
            "Memory: ~54 GB  (cp_occ 24.9 + sa_ls 23.2 + sa_ms 5.8)\n\n"
            "Examples:\n"
            "  %s cp_occ_full.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin"
            " c_vector.txt . 15\n"
            "  %s cp_occ_full.bin sa_ls_word_cf1.bin sa_ms_byte_cf1.bin"
            " c_vector.txt /data/index\n\n",
            argv[0], argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    int k_only = (argc == 7) ? atoi(argv[6]) : 0;
    if (k_only != 0 && k_only != 14 && k_only != 15 && k_only != 16) {
        fprintf(stderr, "[ERROR] k must be 14, 15, 16, or 0\n");
        return EXIT_FAILURE;
    }

    /* ── load C vector ── */
    fprintf(stderr, "Loading c_vector.txt ...\n");
    {
        FILE *f = fopen(argv[4], "r");
        if (!f) { perror(argv[4]); return EXIT_FAILURE; }
        for (int i = 0; i < ALPHABET_SIZE; i++) {
            if (fscanf(f, "%llu", (unsigned long long *)&C[i]) != 1) {
                fprintf(stderr, "[ERROR] c_vector read error at col %d\n", i);
                return EXIT_FAILURE;
            }
        }
        fclose(f);
    }
    fprintf(stderr, "  C[0]=%llu  C[4]=%llu  C[8]=%llu  C[12]=%llu\n",
            (unsigned long long)C[0],  (unsigned long long)C[4],
            (unsigned long long)C[8],  (unsigned long long)C[12]);

    /* ── load cp_occ_full.bin ── */
    fprintf(stderr, "Loading cp_occ_full.bin ...\n");
    {
        FILE *f = fopen(argv[1], "rb");
        if (!f) { perror(argv[1]); return EXIT_FAILURE; }

        rs_occ_full_header_t hdr;
        if (fread(&hdr, sizeof(hdr), 1, f) != 1) {
            fprintf(stderr, "[ERROR] header read failed\n"); return EXIT_FAILURE;
        }
        if (hdr.magic != RS_OCC_MAGIC) {
            fprintf(stderr, "[ERROR] bad magic — wrong file?\n"); return EXIT_FAILURE;
        }
        if (hdr.occ_interval != OCC_INTERVAL || hdr.alphabet_size != ALPHABET_SIZE) {
            fprintf(stderr, "[ERROR] metadata mismatch: occ_interval=%u alphabet=%u\n",
                    hdr.occ_interval, hdr.alphabet_size); return EXIT_FAILURE;
        }

        N       = hdr.bwt_len_total_with_dollar;
        Nblocks = hdr.num_blocks;
        C[ALPHABET_SIZE] = N;   /* c_vec[16] = N, used as upper bound */

        fprintf(stderr, "  N (bwt_total)     = %llu\n", (unsigned long long)N);
        fprintf(stderr, "  num_blocks        = %llu\n", (unsigned long long)Nblocks);
        fprintf(stderr, "  sentinel_index    = %lld\n", (long long)hdr.sentinel_index);

        blocks = (cp_occ32_t *)malloc(Nblocks * sizeof(cp_occ32_t));
        if (!blocks) { perror("malloc blocks"); return EXIT_FAILURE; }

        fprintf(stderr, "  reading %llu blocks (%.2f GB) ...\n",
                (unsigned long long)Nblocks,
                (double)(Nblocks * sizeof(cp_occ32_t)) / 1e9);

        if (fread(blocks, sizeof(cp_occ32_t), Nblocks, f) != Nblocks) {
            fprintf(stderr, "[ERROR] block read failed\n"); return EXIT_FAILURE;
        }
        fclose(f);
        fprintf(stderr, "  cp_occ loaded.\n");
    }

    /* ── load SA CF=1 ── */
    uint64_t sa_n = N;
    fprintf(stderr, "\nLoading SA CF=1 (%llu entries, %.2f GB) ...\n",
            (unsigned long long)sa_n, (double)(sa_n * 5) / 1e9);

    sa_ls = (uint32_t *)malloc(sa_n * sizeof(uint32_t));
    sa_ms = (uint8_t  *)malloc(sa_n * sizeof(uint8_t));
    if (!sa_ls || !sa_ms) { perror("malloc SA"); return EXIT_FAILURE; }

    {
        FILE *f = fopen(argv[2], "rb");
        if (!f) { perror(argv[2]); return EXIT_FAILURE; }
        if (fread(sa_ls, sizeof(uint32_t), sa_n, f) != sa_n) {
            fprintf(stderr, "[ERROR] sa_ls read incomplete\n"); return EXIT_FAILURE;
        }
        fclose(f);
    }
    {
        FILE *f = fopen(argv[3], "rb");
        if (!f) { perror(argv[3]); return EXIT_FAILURE; }
        if (fread(sa_ms, sizeof(uint8_t), sa_n, f) != sa_n) {
            fprintf(stderr, "[ERROR] sa_ms read incomplete\n"); return EXIT_FAILURE;
        }
        fclose(f);
    }
    fprintf(stderr, "  SA CF=1 loaded.\n");

    /* ── generate tables ── */
    const char *out_dir = argv[5];

    if (k_only == 0 || k_only == 14) {
        char path[4096];
        snprintf(path, sizeof(path), "%s/jumptable_14nt.bin", out_dir);
        fprintf(stderr, "\n[ k=14 ]  16^7 = 268,435,456 entries  →  2.0 GiB\n"
                        "  writing: %s\n", path);
        make_table_14(path);
    }

    if (k_only == 0 || k_only == 15) {
        char path[4096];
        snprintf(path, sizeof(path), "%s/jumptable_15nt.bin", out_dir);
        fprintf(stderr, "\n[ k=15 ]  4^15 = 1,073,741,824 entries  →  8.0 GiB\n"
                        "  writing: %s\n", path);
        make_table_15(path);
    }

    if (k_only == 0 || k_only == 16) {
        char path[4096];
        snprintf(path, sizeof(path), "%s/jumptable_16nt.bin", out_dir);
        fprintf(stderr, "\n[ k=16 ]  4^16 = 4,294,967,296 entries  →  32.0 GiB\n"
                        "  writing: %s\n", path);
        make_table_16(path);
    }

    free(blocks); free(sa_ls); free(sa_ms);
    fprintf(stderr, "\nAll done.\n");
    return EXIT_SUCCESS;
}
