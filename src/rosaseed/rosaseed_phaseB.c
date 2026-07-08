
/*************************************************************************************
                           The MIT License

   RosaSeed (Fast and Configurable seeding for short-read alignment),
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

#include "seeding_kernel_mt.h"
#include "file_dec.h"
#include "load_data.h"
#include "macros.h"
#include "profiling.h"
#include "read_init.h"
#include "helper_functions.h"
#include "bwa.h"

/* Build a k-mer jump-table address directly from base-4 read.
   seed_end_base is a base index (0..read_len_bases-1). */
   static inline __attribute__((always_inline))
   uint64_t compute_jumpN_from_base4(const uint8_t *pat4,
                                     int seed_end_base,
                                     int jump_len_nt)
   {
       uint64_t addr = 0;
       int b = seed_end_base;
   
       for (int k = 0; k < jump_len_nt; ++k, --b) {
           uint8_t nt = pat4[b];
           addr |= ((uint64_t)nt) << (2 * k);
       }
       return addr;
   }
   
   
   /* Universal 2-base FM backward step for base-4 read.
      - pat4 is read in base-4 (array of 0..3)
      - read_idx_base is the index of the RIGHT base (i)
      - We read pat4[i] (right) and pat4[i-1] (left)
      - Compose 4-bit base16 symbol: (left<<2) | right
   */
   static inline __attribute__((always_inline))
   int fm_step2_b4(const uint8_t *pat4,
                   int *read_idx_base,  // right base index i
                   uint64_t *l,
                   uint64_t *h)
   {
       int i = *read_idx_base;
   
       /* Must have two bases available */
       if (i <= 0)
           return 0;
   
       uint8_t right = pat4[i];
       uint8_t left  = pat4[i - 1];
   
       base16_t sym = (base16_t)((left << 2) | right);   // LS2 bits are RIGHT base
   
       uint64_t L2 = *l;
       uint64_t H2 = *h;
   
       fm_index_mapping_backward_search(&L2, &H2, sym);
   
       if (L2 >= H2)
           return 0;   // fail
   
       /* success */
       *l = L2;
       *h = H2;
   
       /* consumed TWO bases */
       *read_idx_base = i - 2;
   
       return 1;
   }
   
   
   /* Try LEFT={A,C,G,T} except the original left base.
      Fix RIGHT=pat4[i].  Construct sym=(left<<2)|right.
      Return 1 if any succeeds; update best L/H; consume 1 read base (left).
   */
   static inline __attribute__((always_inline))
   int fm_rescue_left_b4(const uint8_t *pat4,
                         int fail_idx,        // i where FM failed
                         uint64_t L_before,
                         uint64_t H_before,
                         uint64_t *bestL,
                         uint64_t *bestH)
   {
       uint8_t right = pat4[fail_idx];   // RIGHT base stays same
       uint8_t original_left = pat4[fail_idx - 1];
   
       for (uint8_t left = 0; left < 4; left++) {
           if (left == original_left) continue;
   
           base16_t sym = (base16_t)((left << 2) | right);
   
           uint64_t L2 = L_before;
           uint64_t H2 = H_before;
   
           fm_index_mapping_backward_search(&L2, &H2, sym);
   
           if (L2 < H2) {
               *bestL = L2;
               *bestH = H2;
               return 1;   // success
           }
       }
   
       return 0;   // all failed
   }
   
static inline __attribute__((always_inline)) int ref_walk_pairwise_back_b4(
    uint8_t *pat,
    int left_base_idx,
    uint64_t ref_pos_rightbase,
    int *mismatch_base,
    uint64_t *ref_after,
    int *new_read_idx_base,
    int need_at_least,
    int *reached_target,
    uint64_t *leftmost_sym,
    uint8_t  *leftmost_off)
{
    *mismatch_base = -1;
    uint64_t r = ref_pos_rightbase;
    int matched = 0;
    *reached_target = 0;

    uint64_t best_sym = r + 1;
    uint8_t  best_off = 0;

    DEBUG_PRINTF("Entering the ref-read match 'for' loop:\n");
    DEBUG_PRINTF("[pair] need_at_least=%d\n", need_at_least);

    int i = left_base_idx;
    for (; i >= 1; i -= 2)
    {
        uint8_t sym16 = REF_AT(r);

        DEBUG_PRINTF("ref sym:%d\n", sym16);

        uint8_t rbase = (sym16 & 0x3);
        uint8_t lbase = (sym16 >> 2) & 0x3;

        uint8_t right_read = pat[i];
        uint8_t left_read  = pat[i - 1];

        // RIGHT base
        if (rbase != right_read) {
            *mismatch_base = (int)r;
            *ref_after = r;
            *new_read_idx_base = i;
            *leftmost_sym = best_sym;
            *leftmost_off = best_off;
            DEBUG_PRINTF("right base matching failed!!\n");
            DEBUG_PRINTF("total bases matched so far : %d\n", matched);
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 1;

        if (need_at_least > 0 && matched >= need_at_least) {
            *reached_target = 1;
            *ref_after = r;
            *new_read_idx_base = i - 1;
            *leftmost_sym = best_sym;
            *leftmost_off = best_off;
            DEBUG_PRINTF("[pair] early-stop hit after RIGHT @r=%lu matched=%d\n", r, matched);
            return matched;
        }

        // LEFT base
        if (lbase != left_read) {
            *mismatch_base = (int)r;
            *ref_after = r;
            *new_read_idx_base = i - 1;
            *leftmost_sym = best_sym;
            *leftmost_off = best_off;
            DEBUG_PRINTF("left base matching failed!!\n");
            DEBUG_PRINTF("total bases matched so far : %d\n", matched);
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 0;

        if (need_at_least > 0 && matched >= need_at_least) {
            *reached_target = 1;
            *ref_after = r;
            *new_read_idx_base = i - 2;
            *leftmost_sym = best_sym;
            *leftmost_off = best_off;
            DEBUG_PRINTF("[pair] early-stop hit after LEFT @r=%lu matched=%d\n", r, matched);
            return matched;
        }

        --r;
    }

    if (i == -1) {
        *new_read_idx_base = i;
        *ref_after = r;
        *leftmost_sym = best_sym;
        *leftmost_off = best_off;
        if (need_at_least > 0 && matched >= need_at_least) {
            *reached_target = 1;
        }
        return matched;
    }

    if (i == 0) {
        uint8_t sym16 = REF_AT(r);
        uint8_t rbase = sym16 & 0x3;
        uint8_t right_read = pat[0];

        if (rbase != right_read) {
            *mismatch_base = (int)r;
            *ref_after = r;
            *new_read_idx_base = 0;
            *leftmost_sym = best_sym;
            *leftmost_off = best_off;
            if (need_at_least > 0 && matched >= need_at_least) {
                *reached_target = 1;
            }
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 1;
        *new_read_idx_base = -1;
    }

    *ref_after = r;
    *leftmost_sym = best_sym;
    *leftmost_off = best_off;
    if (need_at_least > 0 && matched >= need_at_least) {
        *reached_target = 1;
    }

    return matched;
}
   
   /* Decide next pivot in *base* coordinates.
      - seed_end_base: base index of rightmost base of seed
      - seed_len_bases: total bases consumed for this seed
      - emitted: 1 if we actually emitted a seed; 0 otherwise
   
      If emitted: skip exactly the seed length → disjoint seeds.
      If not emitted: skip just one base (the "culprit" base). */
   static inline __attribute__((always_inline))
   int next_pivot_after_seed_base(int seed_end_base,
                                  int seed_len_bases,
                                  int emitted)
   {
       DEBUG_PRINTF("[next_pivot] seed_end_base=%d seed_len=%d emitted=%d\n",
                    seed_end_base, seed_len_bases, emitted);
   
       int backoff_bases =
           emitted ? (seed_len_bases > 0 ? seed_len_bases : 1)
                   : 1;
   
       int next_base = seed_end_base - backoff_bases;
       if (next_base < 0) return -1;
   
       DEBUG_PRINTF("[next_pivot] backoff=%d -> next_base=%d\n",
                    backoff_bases, next_base);
       return next_base;
   }
   
   /* Try extending a leftover single base at the left edge.
      pat4[0] is the RIGHT base.
      Test all 4 possible LEFT bases.
      If any succeeds, seed_len_bases += 1 (not 2).
   */
   static inline __attribute__((always_inline))
   int fm_left_edge_rescue_b4(const uint8_t *pat4,
                              uint64_t *l,
                              uint64_t *h,
                              int *seed_len_bases)
   {
       uint8_t right = pat4[0];   // leftover base acts as RIGHT half
       uint64_t bestL = *l;
       uint64_t bestH = *h;
   
       for (uint8_t left = 0; left < 4; left++) {
           base16_t sym = (base16_t)((left << 2) | right);
   
           uint64_t L2 = *l;
           uint64_t H2 = *h;
   
           fm_index_mapping_backward_search(&L2, &H2, sym);
   
           if (L2 < H2) {
               // success
               *l = L2;
               *h = H2;
               *seed_len_bases += 1;   // Only 1 new base matched
               return 1;
           }
       }
   
       return 0;   // all 4 failed
   }


/* -----------------------------------------------------------
   Phase II: disjoint, multi-pivot, base-4 read
   shortread_pattern[] is base-4 (0..3), length read_len
   ----------------------------------------------------------- */
   void phaseII_routine(
    uint8_t shortread_pattern[],  /* base-4 read, chosen strand */
    SMEM *matchArray,
    uint64_t *total_smem,
    int *rid,
    int64_t min_intv,
    int *chosen_strand,
    int read_len)
{
    const int     min_seed_bc  = g_min_seed_len_BC;
    const int     num_pivots   = g_num_pivots_B;   /* user-configurable via -pb */

    const int read_len_bases = read_len;   // e.g. 150

    const int pivot_spacing  = read_len_bases / num_pivots;

    int pivots[num_pivots];

    /* Place pivots in BASE coordinates:
       rightmost base index = read_len_bases - 1 */
    for (int i = 0; i < num_pivots; ++i) {
        pivots[i] = (read_len_bases - 1) - i * pivot_spacing;
    }

    /* Precompute jump-table addresses and prefetch */
    uint64_t addr[num_pivots];
    for (int i = 0; i < num_pivots; ++i) {
        int seed_end_base = pivots[i];

        /* Pivot must have at least JT_LEN_NT bases to its left (inclusive) */
        if (seed_end_base < JT_LEN_NT - 1) {
            addr[i] = UINT64_MAX;   // mark as invalid
            continue;
        }

        short_read_len = seed_end_base;
        addr[i] = compute_jumpN_from_base4(shortread_pattern,
                                           seed_end_base,
                                           JT_LEN_NT);
        __builtin_prefetch(&jump_pointers[addr[i]], 0, 1);
    }

    /* ================================================================
       COROUTINE INTERLEAVER across the num_pivots Phase II pivots.
       ================================================================ */

    /*
       PASS 1 (jump_pointers prefetch) already done above.
       PASS 2: read warm jump entries, init slots, prefetch cp_occ.
       PASS 3: round-robin FM interleaver.
       PASS 4: emit seeds from completed slots.                        */

    typedef struct {
        int      pivot_base;
        int      valid;         /* 1 = successfully initialised            */
        int      active;        /* 1 = still running in PASS 3             */
        uint64_t l, h, diff;
        int      seed_end_base;
        int      seed_len_bases;
        int      read_idx_base;
        uint64_t ref_after;
        uint64_t leftmost_sym_tmp;
        uint8_t  leftmost_off_tmp;
        int      phase;         /* 0=classify  1=FM step  2=emit          */
        int      unique_at_jump;
    } P2Slot;

    P2Slot p2slots[num_pivots];
    int n_active = 0;

    for (int i = 0; i < num_pivots; ++i) {
        p2slots[i].valid           = 0;
        p2slots[i].active          = 0;
        p2slots[i].l               = 0;
        p2slots[i].h               = 0;
        p2slots[i].diff            = 0;
        p2slots[i].seed_len_bases  = 0;
        p2slots[i].read_idx_base   = 0;
        p2slots[i].ref_after       = 0;
        p2slots[i].leftmost_sym_tmp = 0;
        p2slots[i].leftmost_off_tmp = 0;
        p2slots[i].phase           = 0;
        p2slots[i].unique_at_jump  = 0;
    }

    /* PASS 2: read warm jump entries, init slots, prefetch first cp_occ */
    for (int i = 0; i < num_pivots; ++i) {
        int pivot_base = pivots[i];
        if (pivot_base < JT_LEN_NT - 1)    continue;
        if (pivot_base + 1 < min_seed_bc)  continue;
        if (addr[i] == UINT64_MAX)          continue;

        uint64_t jp = jump_pointers[addr[i]];   /* warm from PASS 1 */
        uint64_t l, h, diff;
        extract_jump_bounds(jp, &l, &h, &diff);
        if (l >= h) continue;   /* leave valid=0 */

        p2slots[i].valid           = 1;
        p2slots[i].active          = 1;
        p2slots[i].pivot_base      = pivot_base;
        p2slots[i].l               = l;
        p2slots[i].h               = h;
        p2slots[i].diff            = diff;
        p2slots[i].seed_end_base   = pivot_base;
        p2slots[i].seed_len_bases  = JT_LEN_NT;
        p2slots[i].read_idx_base   = pivot_base - JT_LEN_NT;
        p2slots[i].ref_after       = 0;
        p2slots[i].leftmost_sym_tmp = 0;
        p2slots[i].leftmost_off_tmp = 0;
        p2slots[i].phase           = 0;
        p2slots[i].unique_at_jump  = 0;

        uint64_t blk_l = l >> 5, blk_h = h >> 5;
        __builtin_prefetch((const char *)&cp_occ[blk_l],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[blk_l] + 64, 0, 1);
        if (blk_h != blk_l) {
            __builtin_prefetch((const char *)&cp_occ[blk_h],      0, 1);
            __builtin_prefetch((const char *)&cp_occ[blk_h] + 64, 0, 1);
        }
        n_active++;
    }

    /* PASS 3: round-robin interleaver across pivots */
    while (n_active > 0) {
        for (int i = 0; i < num_pivots; ++i) {
            P2Slot *sl = &p2slots[i];
            if (!sl->active) continue;

            /* ---- Phase 0: classify unique vs non-unique ---- */
            if (sl->phase == 0) {
                if (sl->diff == 1) {
                    sl->unique_at_jump = 1;
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }
                sl->unique_at_jump = 0;
                sl->phase = 1;
                /* prefetch already fired in PASS 2 — yield to next pivot */
                continue;
            }

            /* ---- Phase 1: one fm_step2_b4 step ---- */
            if (sl->phase == 1) {
                uint64_t L_before   = sl->l;
                uint64_t H_before   = sl->h;
                int      idx_before = sl->read_idx_base;

                int step_ok = fm_step2_b4(shortread_pattern,
                                          &sl->read_idx_base, &sl->l, &sl->h);

#ifdef TRY_FOUR_CASES_ON_FAIL
                if (!step_ok && idx_before > 0) {
                    uint64_t bestL = L_before, bestH = H_before;
                    int rescue = fm_rescue_left_b4(shortread_pattern, idx_before,
                                                   L_before, H_before, &bestL, &bestH);
                    if (rescue) {
                        sl->l = bestL; sl->h = bestH;
                        sl->read_idx_base = idx_before - 1;
                        sl->seed_len_bases += 1;
                        sl->phase = 2;
                        sl->active = 0; n_active--;
                        continue;
                    }
                }
#endif
                if (!step_ok) {
                    sl->l = L_before; sl->h = H_before;
                    sl->read_idx_base = idx_before;
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                sl->seed_len_bases += 2;
                sl->diff = sl->h - sl->l;

                /* specificity threshold — same as original while(1) */
                if (sl->diff < (uint64_t)min_intv && sl->seed_len_bases >= min_seed_bc) {
                    if (sl->diff == 1 && sl->ref_after == 0) {
                        sl->ref_after = reconstruct_ref_base_index(sl->l);
                        sl->leftmost_sym_tmp = sl->ref_after + 1;
                        sl->leftmost_off_tmp = 0;
                    }
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                /* became unique mid-extension */
                if (sl->diff == 1) {
                    uint64_t ref_pos = reconstruct_ref_base_index(sl->l);
                    int mismatch_base = -1, new_read_idx = sl->read_idx_base, reached = 0;
                    uint64_t ref_after_tmp = ref_pos;
                    sl->leftmost_sym_tmp = ref_pos + 1;
                    sl->leftmost_off_tmp = 0;

                    int need = min_seed_bc - sl->seed_len_bases;
                    if (need < 0) need = 0;

                    int matched = ref_walk_pairwise_back_b4(
                        shortread_pattern, sl->read_idx_base, ref_pos,
                        &mismatch_base, &ref_after_tmp, &new_read_idx,
                        need, &reached,
                        &sl->leftmost_sym_tmp, &sl->leftmost_off_tmp);

                    sl->seed_len_bases += matched;
                    sl->read_idx_base   = new_read_idx;
                    sl->ref_after       = ref_after_tmp;
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                if (sl->read_idx_base < 0) {
                    sl->phase = 2;
                    sl->active = 0; n_active--;
                    continue;
                }

                /* more FM steps — prefetch next cp_occ before yielding */
                __builtin_prefetch((const char *)&cp_occ[sl->l >> 5],      0, 1);
                __builtin_prefetch((const char *)&cp_occ[sl->l >> 5] + 64, 0, 1);
                __builtin_prefetch((const char *)&cp_occ[sl->h >> 5],      0, 1);
                __builtin_prefetch((const char *)&cp_occ[sl->h >> 5] + 64, 0, 1);
                /* stay in phase 1 — continue round-robin */
                continue;
            }

            /* phase 2 should never be active — deactivated above */
            sl->active = 0; n_active--;
        }
    }

    /* PASS 4: handle unique-at-jump ref_walk and emit */
    for (int i = 0; i < num_pivots; ++i) {
        P2Slot *sl = &p2slots[i];
        if (!sl->valid) continue;

        int pivot_base     = sl->pivot_base;
        int seed_end_base  = sl->seed_end_base;
        int seed_len_bases = sl->seed_len_bases;
        int read_idx_base  = sl->read_idx_base;
        uint64_t l         = sl->l;
        uint64_t h         = sl->h;
        uint64_t ref_after = sl->ref_after;
        uint64_t leftmost_sym_tmp = sl->leftmost_sym_tmp;
        uint8_t  leftmost_off_tmp = sl->leftmost_off_tmp;

        if (sl->unique_at_jump) {
            uint64_t ref_pos = l - 1;
            int mismatch_base = -1, new_read_idx = read_idx_base, reached = 0;
            uint64_t ref_after_tmp = ref_pos;
            leftmost_sym_tmp = ref_pos + 1;
            leftmost_off_tmp = 0;

            if (read_idx_base >= 0) {
                int need = min_seed_bc - seed_len_bases;
                if (need < 0) need = 0;

                int matched = ref_walk_pairwise_back_b4(
                    shortread_pattern, read_idx_base, ref_pos,
                    &mismatch_base, &ref_after_tmp, &new_read_idx,
                    need, &reached,
                    &leftmost_sym_tmp, &leftmost_off_tmp);

                seed_len_bases += matched;
                read_idx_base   = new_read_idx;
                ref_after       = ref_after_tmp;
            } else {
                ref_after = ref_pos;
            }
        }

        uint64_t diff = h - l;

        /* emit condition */
        if (seed_len_bases >= min_seed_bc && l < h
            && diff < (uint64_t)min_intv) {

            SMEM *s = &matchArray[*total_smem];
            s->rid = *rid;

            int start_idx = seed_end_base - (seed_len_bases - 1);
            if (start_idx < 0) start_idx = 0;
            int end_idx = seed_end_base;
            if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

            if (*chosen_strand == 0) {
                s->start_idx = start_idx;
                s->end_idx   = end_idx;
            } else {
                int f_start = (read_len - 1) - end_idx;
                int f_end   = (read_len - 1) - start_idx;
                s->start_idx = f_start;
                s->end_idx   = f_end;
            }

            s->unique_flag        = (diff == 1);
            s->unique_ref_pointer = (diff == 1) ? (ref_after + 1) : 0;
            s->unique_leftmost_twostep_pos = (diff == 1) ? leftmost_sym_tmp : 0;
            s->unique_leftmost_base_off    = (diff == 1) ? leftmost_off_tmp : 0;
            s->low_ptr_read    = l;
            s->low_ptr_rc_read = INT_MAX;
            s->smem_score      = (int)diff;
            s->seed_strand     = (uint8_t)(*chosen_strand);

#ifdef DEBUG_SANITY
            int computed_len = s->end_idx - s->start_idx + 1;
            if (computed_len != seed_len_bases) {
                fprintf(stderr,
                        "[SANITY-FAIL P3] rid=%d pivot_base=%d start=%d end=%d "
                        "expected_len=%d got_len=%d\n",
                        s->rid, pivot_base, s->start_idx, s->end_idx,
                        seed_len_bases, computed_len);
                exit(EXIT_FAILURE);
            }
#endif
            (*total_smem)++;
        }
        (void)pivot_base; 
    } /* for pivots PASS 4 */
} /* phaseII_routine */
