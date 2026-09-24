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
#include "phaseA_slot_state_compact.h"

__thread int short_read_len;

extern uint64_t counter_readunique_after1stjump;

static inline __attribute__((always_inline))
uint64_t compute_jumpN_from_base4(const uint8_t *pat4,
                                  int seed_end_base,
                                  int jump_len_nt)
{
    uint64_t addr = 0;
    int b = seed_end_base;
    for (int k = 0; k < jump_len_nt; ++k, --b) {
        if (b < 0 || pat4[b] >= 4) return UINT64_MAX;   /* N or out of range */
        addr |= ((uint64_t)pat4[b]) << (2 * k);
    }
    return addr;
}

static inline __attribute__((always_inline))
int jump_addr_is_valid(uint64_t addr) { return addr != UINT64_MAX; }

static inline __attribute__((always_inline))
int fm_step1_b4(const uint8_t *pat4,
                int *read_idx_base,
                uint64_t *l,
                uint64_t *h)
{
    int i = *read_idx_base;
    if (i < 0) return 0;

    uint8_t c = pat4[i];
    if (c >= 4) return 0;   /* N-HANDLING: stop extension at an ambiguous base */

    uint64_t L2 = *l;
    uint64_t H2 = *h;

    fm_index_mapping_backward_search(&L2, &H2, c);

    if (L2 >= H2) return 0;

    *l = L2;
    *h = H2;
    *read_idx_base = i - 1;
    return 1;
}

static inline __attribute__((always_inline))
int ref_walk_left_b4(const uint8_t *pat,
                     int read_idx_base,
                     uint64_t ref_pos_left_extend,
                     int *mismatch_base,
                     uint64_t *ref_after,
                     int *new_read_idx_base,
                     uint64_t *leftmost_ref_pos)
{
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
        if (r == 0) { i--; break; }
        r--;
        i--;
    }

    *ref_after = r;
    *new_read_idx_base = i;
    *leftmost_ref_pos = best_ref;
    return matched;
}

static inline __attribute__((always_inline))
int next_pivot_after_seed_base(int seed_end_base,
                               int seed_len_bases,
                               int emitted,
                               const uint8_t *pat4)
{
    int backoff_bases =
        emitted ? (seed_len_bases > 0 ? seed_len_bases : 1) : 1;
    int next_base = seed_end_base - backoff_bases - 1;

    /* N-HANDLING: skip backward over an entire run of N in one step, and also
       past a pivot whose jump-table window still covers an N, such a window
       can never produce a lookup, so re-probing each position wastes work. */
    for (;;) {
        if (next_base < 0) return -1;
        if (pat4[next_base] >= 4) { --next_base; continue; }
        int lo = next_base - (JT_LEN_NT - 1);
        if (lo < 0) lo = 0;
        int rightmost_n = -1;
        for (int p = next_base - 1; p >= lo; --p)
            if (pat4[p] >= 4) { rightmost_n = p; break; }
        if (rightmost_n < 0) break;
        next_base = rightmost_n - 1;
    }

    return next_base;
}

static inline __attribute__((always_inline))
void emit_phaseA_seed(SMEM *matchArray,
                      uint64_t *numberofSMEMs,
                      int *read_counter,
                      int read_len_bases,
                      int chosen_strand,
                      int seed_end_base,
                      int seed_len_bases,
                      uint64_t l,
                      uint64_t h,
                      int unique,
                      uint64_t unique_raw_start)
{
    if (*numberofSMEMs >= MATCH_ARRAY_CAPACITY) {
        fprintf(stderr,
                "[ERROR] PhaseA matchArray overflow: rid=%d cap=%d\n",
                *read_counter, MATCH_ARRAY_CAPACITY);
        exit(EXIT_FAILURE);
    }

    SMEM *s = &matchArray[*numberofSMEMs];

    int start_idx = seed_end_base - seed_len_bases + 1;
    int end_idx   = seed_end_base;
    if (start_idx < 0) start_idx = 0;
    if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

    s->rid = *read_counter;

    if (chosen_strand == 0) {
        s->start_idx = start_idx;
        s->end_idx   = end_idx;
    } else {
        int f_start = (read_len_bases - 1) - end_idx;
        int f_end   = (read_len_bases - 1) - start_idx;
        s->start_idx = f_start;
        s->end_idx   = f_end;
    }

    s->unique_flag = unique ? true : false;
    s->unique_ref_pointer = unique ? unique_raw_start : 0;
    s->low_ptr_read = l;
    s->low_ptr_rc_read = UINT64_MAX;
    s->smem_score = unique ? 1 : (int)(h - l);
    s->seed_strand = (uint8_t)chosen_strand;
    s->unique_leftmost_base_off = 0;

#ifdef DEBUG_SANITY
    int computed_len = s->end_idx - s->start_idx + 1;
    if (computed_len != seed_len_bases) {
        fprintf(stderr,
                "[SANITY-FAIL] rid=%d start=%d end=%d expected_len=%d got_len=%d strand=%d\n",
                s->rid, s->start_idx, s->end_idx,
                seed_len_bases, computed_len, chosen_strand);
        exit(EXIT_FAILURE);
    }
#endif

    (*numberofSMEMs)++;
}

/* =======================================================================
   PHASE I 
   ======================================================================= */
void phaseI_routine(
    SMEM *matchArray, uint64_t *numberofSMEMs, int *read_counter,
    uint8_t pat_f4[], uint8_t pat_rc4[], int *chosen_strand, int read_len)
{
    const int min_seed_a = g_min_seed_len_A;
    const int phase1_cap = g_phase1_cap;
    const int read_len_bases = read_len;
    const int E0_base = read_len_bases - 1;

    uint64_t l_f = 0, h_f = 0, diff_f = 0;
    uint64_t l_r = 0, h_r = 0, diff_r = 0;
    uint64_t jp_f = 0, jp_r = 0;

#ifdef ENABLE_F_RC_CHOICE
    {
        uint64_t addr_f = compute_jumpN_from_base4(pat_f4, E0_base, JT_LEN_NT);
        if (jump_addr_is_valid(addr_f)) {
            jp_f = jump_pointers[addr_f];
            extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
        } else { jp_f = 0; l_f = h_f = 0; diff_f = 0; }
    }
    {
        uint64_t addr_r = compute_jumpN_from_base4(pat_rc4, E0_base, JT_LEN_NT);
        if (jump_addr_is_valid(addr_r)) {
            jp_r = jump_pointers[addr_r];
            extract_jump_bounds(jp_r, &l_r, &h_r, &diff_r);
        } else { jp_r = 0; l_r = h_r = 0; diff_r = 0; }
    }
    const int use_rc = (diff_r > 0 && (diff_f == 0 || diff_r < diff_f));
    *chosen_strand = use_rc;
#else
    {
        uint64_t addr_f = compute_jumpN_from_base4(pat_f4, E0_base, JT_LEN_NT);
        if (jump_addr_is_valid(addr_f)) {
            jp_f = jump_pointers[addr_f];
            extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
        } else { jp_f = 0; l_f = h_f = 0; diff_f = 0; }
    }
    const int use_rc = 0;
    *chosen_strand = 0;
#endif

    uint8_t *pat   = use_rc ? pat_rc4 : pat_f4;
    uint64_t l0    = use_rc ? l_r   : l_f;
    uint64_t h0    = use_rc ? h_r   : h_f;
    uint64_t d0    = use_rc ? diff_r : diff_f;
    uint64_t jp0   = use_rc ? jp_r  : jp_f;

    int pivot_base = E0_base;
    int first_iter = 1;

    while (pivot_base >= JT_LEN_NT - 1) {
        if ((pivot_base + 1) < min_seed_a) break;

        uint64_t l = 0, h = 0, diff = 0, jp = 0;

        if (first_iter) {
            l = l0; h = h0; diff = d0; jp = jp0;
            first_iter = 0;
        } else {
            uint64_t addr = compute_jumpN_from_base4(pat, pivot_base, JT_LEN_NT);
            if (jump_addr_is_valid(addr)) {
                jp = jump_pointers[addr];
                extract_jump_bounds(jp, &l, &h, &diff);
            } else { jp = 0; l = h = 0; diff = 0; }
        }

        int seed_end_base  = pivot_base;
        int seed_len_bases = JT_LEN_NT;
        int read_idx_base  = pivot_base - JT_LEN_NT;

        if (diff == 0 || l >= h) {
            pivot_base = next_pivot_after_seed_base(seed_end_base, seed_len_bases, 0, pat);
            continue;
        }

        /* Case 1: unique immediately from jump table */
        if (((jp >> 63) & 1ULL) && diff == 1) {
#ifdef ENABLE_COUNTERS
            counter_readunique_after1stjump++;
#endif
            uint64_t ref_start = jp & 0x1FFFFFFFFULL;
            uint64_t leftmost_ref_start = ref_start;

            if (read_idx_base >= 0 && ref_start > 0) {
                int mismatch_base = -1;
                uint64_t ref_after = 0;
                int new_read_idx = read_idx_base;

                int matched = ref_walk_left_b4(
                    pat, read_idx_base, ref_start - 1,
                    &mismatch_base, &ref_after,
                    &new_read_idx, &leftmost_ref_start);

                seed_len_bases += matched;
                read_idx_base = new_read_idx;
            }

            int emitted = 0;
            if (seed_len_bases >= min_seed_a) {
                emitted = 1;
                emit_phaseA_seed(matchArray, numberofSMEMs, read_counter,
                                 read_len_bases, *chosen_strand,
                                 seed_end_base, seed_len_bases,
                                 l, h, 1, leftmost_ref_start);
            }

            pivot_base = next_pivot_after_seed_base(seed_end_base, seed_len_bases, emitted, pat);
            continue;
        }

        /* Case 2: non-unique interval: extend one base at a time */
        uint64_t unique_raw_start = 0;

        while (read_idx_base >= 0) {
            uint64_t L_before = l, H_before = h;
            int idx_before = read_idx_base;

            int ok = fm_step1_b4(pat, &read_idx_base, &l, &h);

            if (!ok) {
                l = L_before;
                h = H_before;
                read_idx_base = idx_before;
                break;
            }

            seed_len_bases++;
            diff = h - l;

            if (diff == 1) {
                uint64_t ref_start = reconstruct_ref_base_index(l);
                unique_raw_start = ref_start;

                if (read_idx_base >= 0 && ref_start > 0) {
                    int mismatch_base = -1;
                    uint64_t ref_after = 0;
                    int new_read_idx = read_idx_base;
                    uint64_t leftmost_ref_start = ref_start;

                    int matched = ref_walk_left_b4(
                        pat, read_idx_base, ref_start - 1,
                        &mismatch_base, &ref_after,
                        &new_read_idx, &leftmost_ref_start);

                    seed_len_bases += matched;
                    read_idx_base = new_read_idx;
                    unique_raw_start = leftmost_ref_start;
                }
                break;
            }
        }

        diff = h - l;
        if (diff == 1 && unique_raw_start == 0)
            unique_raw_start = reconstruct_ref_base_index(l);

        int interval = (int)diff;
        int emitted = 0;

        if (seed_len_bases >= min_seed_a && l < h
            && (phase1_cap == 0 || interval <= phase1_cap)) {
            emitted = 1;
            emit_phaseA_seed(matchArray, numberofSMEMs, read_counter,
                             read_len_bases, *chosen_strand,
                             seed_end_base, seed_len_bases,
                             l, h, (diff == 1), unique_raw_start);
        }

        pivot_base = next_pivot_after_seed_base(seed_end_base, seed_len_bases, emitted, pat);
    }
}

/* =======================================================================
   COROUTINE BATCH PHASE I: 1-step version
   Same design as 2-step coroutine but uses fm_step1_b4 (1 base/step)
   and ref_walk_left_b4 (single-base ref walk).
   All seed logic, exit conditions, emit conditions unchanged.
   ======================================================================= */

/* Per-pivot working state for one read in the coroutine scheduler */
typedef struct {
    const uint8_t *pat;
    int      read_len;
    int      chosen_strand;
    int      read_counter;

    int      pivot_base;
    int      first_iter;
    uint64_t l0, h0, d0, jp0;

    int      seed_end_base;
    int      seed_len_bases;
    int      read_idx_base;
    uint64_t l, h, diff, jp;
    uint64_t unique_raw_start;

    /* phase: 0=jump/classify  1=FM step loop  2=emit+advance */
    int      phase;
    int      fm_loop_done;

    SMEM    *matchArray;
    uint64_t *numberofSMEMs;

    int      min_seed_a;
    int      phase1_cap;
    int      read_len_bases;
} PhaseI1_PivotState;

/* -----------------------------------------------------------------------
   phaseI1_slot_init()
   Strand probe + first jump + prefetch cp_occ blocks
   ----------------------------------------------------------------------- */
static void phaseI1_slot_init(
    PhaseI1_PivotState *st,
    const uint8_t *f4,
    const uint8_t *rc4,
    int read_len,
    int read_counter,
    SMEM *matchArray,
    uint64_t *numberofSMEMs)
{
    st->read_len       = read_len;
    st->read_len_bases = read_len;
    st->read_counter   = read_counter;
    st->matchArray     = matchArray;
    st->numberofSMEMs  = numberofSMEMs;
    st->min_seed_a     = g_min_seed_len_A;
    st->phase1_cap     = g_phase1_cap;
    st->first_iter     = 1;
    st->phase          = 0;
    st->fm_loop_done   = 0;
    st->unique_raw_start = 0;

    const int E0_base = read_len - 1;

    uint64_t l_f = 0, h_f = 0, diff_f = 0, jp_f = 0;
    uint64_t l_r = 0, h_r = 0, diff_r = 0, jp_r = 0;

#ifdef ENABLE_F_RC_CHOICE
    {
        uint64_t addr_f = compute_jumpN_from_base4(f4, E0_base, JT_LEN_NT);
        if (jump_addr_is_valid(addr_f)) {
            jp_f = jump_pointers[addr_f];
            extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
        } else { jp_f = 0; l_f = h_f = 0; diff_f = 0; }
    }
    {
        uint64_t addr_r = compute_jumpN_from_base4(rc4, E0_base, JT_LEN_NT);
        if (jump_addr_is_valid(addr_r)) {
            jp_r = jump_pointers[addr_r];
            extract_jump_bounds(jp_r, &l_r, &h_r, &diff_r);
        } else { jp_r = 0; l_r = h_r = 0; diff_r = 0; }
    }
    const int use_rc = (diff_r > 0 && (diff_f == 0 || diff_r < diff_f));
    st->chosen_strand = use_rc;
    st->pat  = use_rc ? rc4 : f4;
    st->l0   = use_rc ? l_r  : l_f;
    st->h0   = use_rc ? h_r  : h_f;
    st->d0   = use_rc ? diff_r : diff_f;
    st->jp0  = use_rc ? jp_r : jp_f;
#else
    {
        uint64_t addr_f = compute_jumpN_from_base4(f4, E0_base, JT_LEN_NT);
        if (jump_addr_is_valid(addr_f)) {
            jp_f = jump_pointers[addr_f];
            extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
        } else { jp_f = 0; l_f = h_f = 0; diff_f = 0; }
    }
    st->chosen_strand = 0;
    st->pat  = f4;
    st->l0   = l_f;
    st->h0   = h_f;
    st->d0   = diff_f;
    st->jp0  = jp_f;
#endif

    st->pivot_base = E0_base;
    st->l = st->l0; st->h = st->h0;
    st->diff = st->d0; st->jp = st->jp0;

    /* prefetch first cp_occ blocks */
    uint64_t blk_l = st->l0 >> 5;
    uint64_t blk_h = st->h0 >> 5;
    __builtin_prefetch((const char *)&cp_occ[blk_l],      0, 1);
    __builtin_prefetch((const char *)&cp_occ[blk_l] + 64, 0, 1);
    if (blk_h != blk_l) {
        __builtin_prefetch((const char *)&cp_occ[blk_h],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[blk_h] + 64, 0, 1);
    }
}

/* -----------------------------------------------------------------------
   phaseI1_slot_step()
   Advance one slot by ONE fm_step1_b4() call (or one non-FM action).
   Prefetches next cp_occ block after each successful FM step.
   Returns 1 = still active, 0 = slot done.
   ----------------------------------------------------------------------- */
static int phaseI1_slot_step(PhaseI1_PivotState *st)
{
    const int min_seed_a     = st->min_seed_a;
    const int phase1_cap     = st->phase1_cap;
    const int read_len_bases = st->read_len_bases;
    const uint8_t *pat       = st->pat;

    /* ============================================================
       PHASE 0: jump for current pivot, classify unique vs non-unique
       ============================================================ */
    if (st->phase == 0) {
        int pivot_base = st->pivot_base;

        if (pivot_base < JT_LEN_NT - 1 || (pivot_base + 1) < min_seed_a) {
            st->pivot_base = -1; return 0;
        }

        if (st->first_iter) {
            st->l = st->l0; st->h = st->h0;
            st->diff = st->d0; st->jp = st->jp0;
            st->first_iter = 0;
        } else {
            uint64_t addr = compute_jumpN_from_base4(pat, pivot_base, JT_LEN_NT);
            if (jump_addr_is_valid(addr)) {
                st->jp = jump_pointers[addr];
                extract_jump_bounds(st->jp, &st->l, &st->h, &st->diff);
            } else { st->jp = 0; st->l = st->h = 0; st->diff = 0; }
        }

        st->seed_end_base   = pivot_base;
        st->seed_len_bases  = JT_LEN_NT;
        st->read_idx_base   = pivot_base - JT_LEN_NT;
        st->unique_raw_start = 0;
        st->fm_loop_done    = 0;

        if (st->diff == 0 || st->l >= st->h) {
            /* dead pivot: advance immediately */
            st->pivot_base = next_pivot_after_seed_base(
                st->seed_end_base, st->seed_len_bases, 0, pat);
            if (st->pivot_base < JT_LEN_NT - 1) { st->pivot_base = -1; return 0; }
            if (st->pivot_base + 1 < min_seed_a) { st->pivot_base = -1; return 0; }
            st->phase = 0;
            return 1;
        }

        /* Case 1: unique at jump: (jp>>63)&1 flag and diff==1 */
        if (((st->jp >> 63) & 1ULL) && st->diff == 1) {
#ifdef ENABLE_COUNTERS
            counter_readunique_after1stjump++;
#endif
            uint64_t ref_start = st->jp & 0x1FFFFFFFFULL;
            uint64_t leftmost_ref_start = ref_start;

            if (st->read_idx_base >= 0 && ref_start > 0) {
                int mismatch_base = -1;
                uint64_t ref_after = 0;
                int new_read_idx = st->read_idx_base;

                int matched = ref_walk_left_b4(
                    pat, st->read_idx_base, ref_start - 1,
                    &mismatch_base, &ref_after,
                    &new_read_idx, &leftmost_ref_start);

                st->seed_len_bases += matched;
                st->read_idx_base = new_read_idx;
            }

            int emitted = 0;
            if (st->seed_len_bases >= min_seed_a) {
                emitted = 1;
                emit_phaseA_seed(st->matchArray, st->numberofSMEMs,
                                 &st->read_counter, read_len_bases,
                                 st->chosen_strand, st->seed_end_base,
                                 st->seed_len_bases, st->l, st->h,
                                 1, leftmost_ref_start);
            }

            st->pivot_base = next_pivot_after_seed_base(
                st->seed_end_base, st->seed_len_bases, emitted, pat);

            if (st->pivot_base < JT_LEN_NT - 1) { st->pivot_base = -1; return 0; }
            if (st->pivot_base + 1 < min_seed_a) { st->pivot_base = -1; return 0; }

            /* prefetch next pivot's jump entry (skip if window touches N) */
            if (st->pivot_base >= JT_LEN_NT - 1) {
                uint64_t next_addr =
                    compute_jumpN_from_base4(pat, st->pivot_base, JT_LEN_NT);
                if (jump_addr_is_valid(next_addr))
                    __builtin_prefetch(&jump_pointers[next_addr], 0, 1);
            }
            st->phase = 0;
            return 1;
        }

        /* Case 2: non-unique: transition to FM step loop */
        st->phase = 1;

        /* prefetch cp_occ for first FM step */
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5] + 64, 0, 1);
        return 1;
    }

    /* ============================================================
       PHASE 1: one fm_step1_b4 step
       ============================================================ */
    if (st->phase == 1) {

        if (st->fm_loop_done) { st->phase = 2; return 1; }

        if (st->read_idx_base < 0) {
            st->fm_loop_done = 1; st->phase = 2; return 1;
        }

        uint64_t L_before = st->l, H_before = st->h;
        int      idx_before = st->read_idx_base;

        int ok = fm_step1_b4(pat, &st->read_idx_base, &st->l, &st->h);

        if (!ok) {
            /* rollback: identical to original */
            st->l = L_before;
            st->h = H_before;
            st->read_idx_base = idx_before;
            st->fm_loop_done = 1;
            st->phase = 2;
            return 1;
        }

        st->seed_len_bases++;   /* 1-step increments by 1, not 2 */
        st->diff = st->h - st->l;

        if (st->diff == 1) {
            /* became unique mid-extension */
            uint64_t ref_start = reconstruct_ref_base_index(st->l);
            st->unique_raw_start = ref_start;

            if (st->read_idx_base >= 0 && ref_start > 0) {
                int mismatch_base = -1;
                uint64_t ref_after = 0;
                int new_read_idx = st->read_idx_base;
                uint64_t leftmost_ref_start = ref_start;

                int matched = ref_walk_left_b4(
                    pat, st->read_idx_base, ref_start - 1,
                    &mismatch_base, &ref_after,
                    &new_read_idx, &leftmost_ref_start);

                st->seed_len_bases += matched;
                st->read_idx_base = new_read_idx;
                st->unique_raw_start = leftmost_ref_start;
            }

            st->fm_loop_done = 1;
            st->phase = 2;
            return 1;
        }

        if (st->read_idx_base < 0) {
            st->fm_loop_done = 1; st->phase = 2; return 1;
        }

        /* prefetch NEXT step's cp_occ blocks before yielding */
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5] + 64, 0, 1);
        return 1;
    }

    /* ============================================================
       PHASE 2: FM loop done: emit seed, advance pivot
       ============================================================ */
    if (st->phase == 2) {

        st->diff = st->h - st->l;

        if (st->diff == 1 && st->unique_raw_start == 0)
            st->unique_raw_start = reconstruct_ref_base_index(st->l);

        int interval = (int)st->diff;
        int emitted  = 0;

        if (st->seed_len_bases >= min_seed_a && st->l < st->h
            && (phase1_cap == 0 || interval <= phase1_cap)) {
            emitted = 1;
            emit_phaseA_seed(st->matchArray, st->numberofSMEMs,
                             &st->read_counter, read_len_bases,
                             st->chosen_strand, st->seed_end_base,
                             st->seed_len_bases, st->l, st->h,
                             (st->diff == 1), st->unique_raw_start);
        }

        st->pivot_base = next_pivot_after_seed_base(
            st->seed_end_base, st->seed_len_bases, emitted, pat);

        if (st->pivot_base < JT_LEN_NT - 1) { st->pivot_base = -1; return 0; }
        if (st->pivot_base + 1 < min_seed_a) { st->pivot_base = -1; return 0; }

        /* prefetch next pivot's jump entry (skip if window touches N) */
        if (st->pivot_base >= JT_LEN_NT - 1) {
            uint64_t next_addr =
                compute_jumpN_from_base4(pat, st->pivot_base, JT_LEN_NT);
            if (jump_addr_is_valid(next_addr))
                __builtin_prefetch(&jump_pointers[next_addr], 0, 1);
        }

        st->phase        = 0;
        st->fm_loop_done = 0;
        st->unique_raw_start = 0;
        return 1;
    }

    return 0;
}

/* -----------------------------------------------------------------------
   phaseI1_batch_interleaved()
   Public entry point: called by bridge instead of phaseI_routine().
   Round-robins RS_BATCH reads through phaseI1_slot_step() one step
   at a time, prefetching cp_occ after every FM step.
   ----------------------------------------------------------------------- */
void phaseI1_batch_interleaved(
    const uint8_t **f4,
    const uint8_t **rc4,
    const int      *read_lens,
    const int      *read_ctrs,
    SMEM          **matchArrays,
    uint64_t      **nSMEMs,
    int            *chosen_strands,
    int             nreads)
{
    PhaseI1_PivotState slots[RS_BATCH];

    for (int s = 0; s < nreads; ++s) {
        phaseI1_slot_init(&slots[s],
                          f4[s], rc4[s],
                          read_lens[s], read_ctrs[s],
                          matchArrays[s], nSMEMs[s]);
        chosen_strands[s] = slots[s].chosen_strand;
    }

    int active = nreads;
    int alive[RS_BATCH];
    for (int s = 0; s < nreads; ++s) alive[s] = 1;

    while (active > 0) {
        for (int s = 0; s < nreads; ++s) {
            if (!alive[s]) continue;
            if (!phaseI1_slot_step(&slots[s])) {
                alive[s] = 0;
                --active;
            }
        }
    }

    for (int s = 0; s < nreads; ++s)
        chosen_strands[s] = slots[s].chosen_strand;
}
