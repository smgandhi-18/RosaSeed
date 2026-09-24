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

#include "seeding_kernel.h"
#include "file_dec.h"
#include "load_data.h"
#include "macros.h"
#include "profiling.h"
#include "helper_functions.h"
#include "bwa.h"

/* Build k-mer jump-table address from base4 read.
 *
 * seed_end_base is the rightmost base index.
 * Rightmost base goes into lowest 2 bits.
 */
static inline __attribute__((always_inline))
uint64_t compute_jumpN_from_base4(const uint8_t *pat4,
                                  int seed_end_base,
                                  int jump_len_nt)
{
    uint64_t addr = 0;
    int b = seed_end_base;

    for (int k = 0; k < jump_len_nt; ++k, --b) {
        if (b < 0 || pat4[b] >= 4) return UINT64_MAX;   /* N or out of range */
        uint8_t nt = pat4[b];
        addr |= ((uint64_t)nt) << (2 * k);
    }

    return addr;
}

static inline __attribute__((always_inline))
int jump_addr_is_valid(uint64_t addr) { return addr != UINT64_MAX; }

/* One-base FM backward step for 1-step/base4 FM-index. */
static inline __attribute__((always_inline))
int fm_step1_b4(const uint8_t *pat4,
                int *read_idx_base,
                uint64_t *l,
                uint64_t *h)
{
    int i = *read_idx_base;

    if (i < 0)
        return 0;

    uint8_t c = pat4[i];
    if (c >= 4) return 0;   /* N-HANDLING: stop extension at an ambiguous base */

    uint64_t L2 = *l;
    uint64_t H2 = *h;

    fm_index_mapping_backward_search(&L2, &H2, c);

    if (L2 >= H2)
        return 0;

    *l = L2;
    *h = H2;
    *read_idx_base = i - 1;

    return 1;
}

/* Reference walk after a seed becomes unique.
 *
 * ref_pos_left_extend = position immediately left of the already matched seed.
 * Stops early if need_at_least additional bases are matched.
 */
static inline __attribute__((always_inline))
int ref_walk_left_b4_until(const uint8_t *pat,
                           int read_idx_base,
                           uint64_t ref_pos_left_extend,
                           int need_at_least,
                           int *reached_target,
                           int *mismatch_base,
                           uint64_t *ref_after,
                           int *new_read_idx_base,
                           uint64_t *leftmost_ref_pos)
{
    *reached_target = 0;
    *mismatch_base = -1;

    uint64_t r = ref_pos_left_extend;
    int i = read_idx_base;
    int matched = 0;

    uint64_t best_ref = r + 1;

    while (i >= 0) {
        uint8_t rb = REF_AT(r);

        if (rb != pat[i]) {
            *mismatch_base = (int)r;
            *ref_after = r;
            *new_read_idx_base = i;
            *leftmost_ref_pos = best_ref;
            return matched;
        }

        best_ref = r;
        matched++;

        if (need_at_least > 0 && matched >= need_at_least) {
            *reached_target = 1;
            *ref_after = r;
            *new_read_idx_base = i - 1;
            *leftmost_ref_pos = best_ref;
            return matched;
        }

        if (r == 0) {
            i--;
            break;
        }

        r--;
        i--;
    }

    if (need_at_least > 0 && matched >= need_at_least)
        *reached_target = 1;

    *ref_after = r;
    *new_read_idx_base = i;
    *leftmost_ref_pos = best_ref;

    return matched;
}

static inline __attribute__((always_inline))
void emit_phaseB_seed(uint8_t shortread_pattern[],
                      SMEM *matchArray,
                      uint64_t *total_smem,
                      int *rid,
                      int read_len_bases,
                      int chosen_strand,
                      int seed_end_base,
                      int seed_len_bases,
                      uint64_t l,
                      uint64_t h,
                      int unique,
                      uint64_t unique_raw_start)
{
    if (*total_smem >= MATCH_ARRAY_CAPACITY) {
        fprintf(stderr,
                "[ERROR] PhaseB matchArray overflow: rid=%d cap=%d\n",
                *rid, MATCH_ARRAY_CAPACITY);
        exit(EXIT_FAILURE);
    }

    SMEM *s = &matchArray[*total_smem];

    int start_idx = seed_end_base - seed_len_bases + 1;
    int end_idx = seed_end_base;

    if (start_idx < 0) start_idx = 0;
    if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

    s->rid = *rid;

    if (chosen_strand == 0) {
        s->start_idx = start_idx;
        s->end_idx = end_idx;
    } else {
        int f_start = (read_len_bases - 1) - end_idx;
        int f_end   = (read_len_bases - 1) - start_idx;

        s->start_idx = f_start;
        s->end_idx = f_end;
    }

    s->unique_flag = unique ? true : false;
    s->unique_ref_pointer = unique ? unique_raw_start : 0;

    s->low_ptr_read = l;
    s->low_ptr_rc_read = UINT64_MAX;

    s->smem_score = unique ? 1 : (int)(h - l);
    s->seed_strand = (uint8_t)chosen_strand;

    /* 1-step has no base16 symbol offset. */
    s->unique_leftmost_base_off = 0;

#ifdef DEBUG_SANITY
    int computed_len = s->end_idx - s->start_idx + 1;
    if (computed_len != seed_len_bases) {
        fprintf(stderr,
                "[SANITY-FAIL PhaseB-1step] rid=%d start=%d end=%d "
                "expected_len=%d got_len=%d strand=%d\n",
                s->rid,
                s->start_idx,
                s->end_idx,
                seed_len_bases,
                computed_len,
                chosen_strand);
        exit(EXIT_FAILURE);
    }
#endif

    (*total_smem)++;
}

/* -----------------------------------------------------------
   Phase B / Phase II: fixed-pivot, no-merge, 1-step/base4
   ----------------------------------------------------------- */
void phaseII_routine(
    uint8_t shortread_pattern[],
    SMEM *matchArray,
    uint64_t *total_smem,
    int *rid,
    int64_t min_intv,
    int *chosen_strand,
    int read_len)
{
    const int min_seed_bc = g_min_seed_len_BC;
    const int num_pivots  = g_num_pivots_B;

    const int read_len_bases = read_len;

    if (num_pivots <= 0) return;

    const int pivot_spacing = read_len_bases / num_pivots;

    int pivots[num_pivots];
    uint64_t addr[num_pivots];

    for (int i = 0; i < num_pivots; ++i) {
        pivots[i] = (read_len_bases - 1) - i * pivot_spacing;

        if (pivots[i] < JT_LEN_NT - 1) {
            addr[i] = UINT64_MAX;
            continue;
        }

        addr[i] = compute_jumpN_from_base4(shortread_pattern,
                                           pivots[i],
                                           JT_LEN_NT);

        if (jump_addr_is_valid(addr[i]))
            __builtin_prefetch(&jump_pointers[addr[i]], 0, 1);
    }

    /* ================================================================
       COROUTINE INTERLEAVER across the num_pivots Phase II (1-step)
       pivots. Same principle as 2-step: interleave FM
       extension one step at a time across all pivots, prefetching
       cp_occ after each step before yielding. All exit conditions and
       emit logic identical to the original serial loop.
       ================================================================ */

    typedef struct {
        int      pivot_base;
        int      valid;           /* 1 = successfully initialised in PASS 2 */
        int      active;          /* 1 = still running in PASS 3            */
        uint64_t l, h, diff;
        uint64_t jp;              /* raw jump pointer, needed for (jp>>63)
                                     unique-at-jump flag check             */
        int      seed_end_base;
        int      seed_len_bases;
        int      read_idx_base;
        uint64_t unique_raw_start;
        int      phase;           /* 0=classify  1=FM step  2=emit          */
        int      unique_at_jump;  /* 1 = (jp>>63)&1 && diff==1 at jump      */
    } P1sSlot;

    P1sSlot p1slots[num_pivots];
    int n_active = 0;

    for (int i = 0; i < num_pivots; ++i) {
        p1slots[i].valid            = 0;
        p1slots[i].active           = 0;
        p1slots[i].l                = 0;
        p1slots[i].h                = 0;
        p1slots[i].diff             = 0;
        p1slots[i].jp               = 0;
        p1slots[i].seed_len_bases   = 0;
        p1slots[i].read_idx_base    = 0;
        p1slots[i].unique_raw_start = 0;
        p1slots[i].phase            = 0;
        p1slots[i].unique_at_jump   = 0;
    }

    /* PASS 2: read warm jump entries, init slots, prefetch first cp_occ */
    for (int i = 0; i < num_pivots; ++i) {
        int pivot_base = pivots[i];
        if (pivot_base < JT_LEN_NT - 1)   continue;
        if (pivot_base + 1 < min_seed_bc) continue;
        if (addr[i] == UINT64_MAX)         continue;

        uint64_t jp = jump_pointers[addr[i]];   /* warm from PASS 1 */
        uint64_t l, h, diff;
        extract_jump_bounds(jp, &l, &h, &diff);

        /* same guard as original: diff==0 also bails */
        if (diff == 0 || l >= h) continue;   /* leave valid=0 */

        p1slots[i].valid            = 1;
        p1slots[i].active           = 1;
        p1slots[i].pivot_base       = pivot_base;
        p1slots[i].l                = l;
        p1slots[i].h                = h;
        p1slots[i].diff             = diff;
        p1slots[i].jp               = jp;
        p1slots[i].seed_end_base    = pivot_base;
        p1slots[i].seed_len_bases   = JT_LEN_NT;
        p1slots[i].read_idx_base    = pivot_base - JT_LEN_NT;
        p1slots[i].unique_raw_start = 0;
        p1slots[i].phase            = 0;
        p1slots[i].unique_at_jump   = 0;

        uint64_t blk_l = l >> 5, blk_h = h >> 5;
        __builtin_prefetch((const char *)&cp_occ[blk_l],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[blk_l] + 64, 0, 1);
        if (blk_h != blk_l) {
            __builtin_prefetch((const char *)&cp_occ[blk_h],      0, 1);
            __builtin_prefetch((const char *)&cp_occ[blk_h] + 64, 0, 1);
        }
        n_active++;
    }

    /* PASS 3: round-robin FM interleaver */
    while (n_active > 0) {
        for (int i = 0; i < num_pivots; ++i) {
            P1sSlot *sl = &p1slots[i];
            if (!sl->active) continue;

            /* ---- Phase 0: classify unique-at-jump vs non-unique ---- */
            if (sl->phase == 0) {
                if (((sl->jp >> 63) & 1ULL) && sl->diff == 1) {
                    /* unique at jump, ref_walk done in PASS 4 */
                    sl->unique_at_jump = 1;
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }
                sl->unique_at_jump = 0;
                sl->phase = 1;

                /* check read_idx_base >= 0 before entering FM loop */
                if (sl->read_idx_base < 0) {
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }
                /* cp_occ prefetched in PASS 2, yield to next pivot */
                continue;
            }

            /* ---- Phase 1: one fm_step1_b4 step (1 base) ---- */
            if (sl->phase == 1) {
                uint64_t L_before   = sl->l;
                uint64_t H_before   = sl->h;
                int      idx_before = sl->read_idx_base;

                int ok = fm_step1_b4(shortread_pattern,
                                     &sl->read_idx_base, &sl->l, &sl->h);

                if (!ok) {
                    /* rollback */
                    sl->l = L_before;
                    sl->h = H_before;
                    sl->read_idx_base = idx_before;
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                sl->seed_len_bases++;   /* 1 base per step */
                sl->diff = sl->h - sl->l;

                /* min_intv early-exit: same as original while() including
                   inline unique_raw_start reconstruction when diff==1 */
                if (sl->diff < (uint64_t)min_intv &&
                    sl->seed_len_bases >= min_seed_bc) {
                    if (sl->diff == 1)
                        sl->unique_raw_start = reconstruct_ref_base_index(sl->l);
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                /* became unique mid-extension */
                if (sl->diff == 1) {
                    uint64_t ref_start = reconstruct_ref_base_index(sl->l);
                    sl->unique_raw_start = ref_start;

                    if (sl->read_idx_base >= 0 && ref_start > 0) {
                        int need = min_seed_bc - sl->seed_len_bases;
                        if (need < 0) need = 0;

                        int reached = 0, mismatch_base = -1;
                        uint64_t ref_after = 0;
                        int new_read_idx = sl->read_idx_base;
                        uint64_t leftmost_ref_start = ref_start;

                        int matched = ref_walk_left_b4_until(
                            shortread_pattern,
                            sl->read_idx_base,
                            ref_start - 1,
                            need,
                            &reached,
                            &mismatch_base,
                            &ref_after,
                            &new_read_idx,
                            &leftmost_ref_start);

                        sl->seed_len_bases   += matched;
                        sl->read_idx_base     = new_read_idx;
                        sl->unique_raw_start  = leftmost_ref_start;
                    }

                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                /* read exhausted */
                if (sl->read_idx_base < 0) {
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                /* more FM steps, prefetch next cp_occ before yielding */
                __builtin_prefetch((const char *)&cp_occ[sl->l >> 5],      0, 1);
                __builtin_prefetch((const char *)&cp_occ[sl->l >> 5] + 64, 0, 1);
                __builtin_prefetch((const char *)&cp_occ[sl->h >> 5],      0, 1);
                __builtin_prefetch((const char *)&cp_occ[sl->h >> 5] + 64, 0, 1);
                continue;
            }

            /* phase 2 is terminal, deactivated above */
            sl->active = 0; n_active--;
        }
    }

    /* PASS 4: handle unique-at-jump ref_walk, trailing unique_raw_start
       reconstruction, and emit: identical to original serial loop tail */
    for (int i = 0; i < num_pivots; ++i) {
        P1sSlot *sl = &p1slots[i];
        if (!sl->valid) continue;

        int      seed_end_base  = sl->seed_end_base;
        int      seed_len_bases = sl->seed_len_bases;
        int      read_idx_base  = sl->read_idx_base;
        uint64_t l              = sl->l;
        uint64_t h              = sl->h;
        uint64_t unique_raw_start = sl->unique_raw_start;

        if (sl->unique_at_jump) {
            uint64_t ref_start = sl->jp & 0x1FFFFFFFFULL;
            unique_raw_start = ref_start;

            if (read_idx_base >= 0 && ref_start > 0) {
                int need = min_seed_bc - seed_len_bases;
                if (need < 0) need = 0;

                int reached = 0, mismatch_base = -1;
                uint64_t ref_after = 0;
                int new_read_idx = read_idx_base;
                uint64_t leftmost_ref_start = ref_start;

                int matched = ref_walk_left_b4_until(
                    shortread_pattern,
                    read_idx_base,
                    ref_start - 1,
                    need,
                    &reached,
                    &mismatch_base,
                    &ref_after,
                    &new_read_idx,
                    &leftmost_ref_start);

                seed_len_bases     += matched;
                read_idx_base       = new_read_idx;
                unique_raw_start    = leftmost_ref_start;
            }
        } else {
            /* non-unique trailing reconstruction */
            uint64_t diff = h - l;
            if (diff == 1 && unique_raw_start == 0)
                unique_raw_start = reconstruct_ref_base_index(l);
        }

        uint64_t diff = h - l;
        (void)read_idx_base;

        if (seed_len_bases >= min_seed_bc &&
            l < h &&
            diff < (uint64_t)min_intv) {

            emit_phaseB_seed(shortread_pattern,
                             matchArray,
                             total_smem,
                             rid,
                             read_len_bases,
                             *chosen_strand,
                             seed_end_base,
                             seed_len_bases,
                             l,
                             h,
                             (diff == 1),
                             unique_raw_start);
        }
    } /* for pivots PASS 4 */
} /* phaseII_routine */
