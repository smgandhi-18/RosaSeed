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
#include "get_reads.h"
#include "helper_functions.h"
#include "bwa.h"

__thread int short_read_len;

extern uint64_t counter_readunique_after1stjump;

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
    RIGHT=pat4[i].  Construct sym=(left<<2)|right.
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
    uint64_t *leftmost_sym,
    uint8_t  *leftmost_off)
{
    *mismatch_base = -1;
    uint64_t r = ref_pos_rightbase;
    int matched = 0;

    // Default: before any extra extension, the already-known unique anchor
    // starts at the LEFT base of symbol (r + 1).
    uint64_t best_sym = r + 1;
    uint8_t  best_off = 0;

    DEBUG_PRINTF("Entering the ref-read match 'for' loop:\n");

    int i = left_base_idx;
    for (; i >= 1; i -= 2)
    {
        uint8_t sym16 = REF_AT(r);

        DEBUG_PRINTF("ref sym:%d\n", sym16);

        uint8_t rbase = (sym16 & 0x3);
        uint8_t lbase = (sym16 >> 2) & 0x3;

        uint8_t right_read = pat[i];
        uint8_t left_read  = pat[i - 1];

        // Compare RIGHT base first
        if (rbase != right_read) {
            *mismatch_base  = (int)r;
            *ref_after      = r;
            *new_read_idx_base = i;
            *leftmost_sym   = best_sym;
            *leftmost_off   = best_off;
            DEBUG_PRINTF("right base matching failed!!\n");
            DEBUG_PRINTF("total bases matched so far : %d\n", matched);
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 1;   // right base of symbol r is now the leftmost matched base so far

        // Compare LEFT base of the same symbol
        if (lbase != left_read) {
            *mismatch_base  = (int)r;
            *ref_after      = r;
            *new_read_idx_base = i - 1;
            *leftmost_sym   = best_sym;
            *leftmost_off   = best_off;
            DEBUG_PRINTF("left base matching failed!!\n");
            DEBUG_PRINTF("total bases matched so far : %d\n", matched);
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 0;   // left base of symbol r is now the leftmost matched base so far

        --r;
    }

    if (i == -1) {
        *new_read_idx_base = i;
        *ref_after = r;
        *leftmost_sym = best_sym;
        *leftmost_off = best_off;
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
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 1;   // only right base matched
        *new_read_idx_base = -1;
    }

    *ref_after = r;
    *leftmost_sym = best_sym;
    *leftmost_off = best_off;

    DEBUG_PRINTF("Reached left end of the read!!\n");
    DEBUG_PRINTF("total bases matched so far : %d\n", matched);
    DEBUG_PRINTF("Current REF position: %ld\n", *ref_after);

    return matched;
}

/* Decide next pivot in *base* coordinates.
   - seed_end_base: base index of rightmost base of seed
   - seed_len_bases: total bases consumed for this seed
   - emitted: 1 if we actually emitted a seed; 0 otherwise

   If emitted: skip exactly the seed length -> disjoint seeds.
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

    int next_base = seed_end_base - backoff_bases - 1;
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

/* =======================================================================
   PHASE A
   ======================================================================= */
void phaseI_routine(
    SMEM *matchArray, uint64_t *numberofSMEMs, int *read_counter,
    uint8_t pat_f4[], uint8_t pat_rc4[], int *chosen_strand, int read_len)
{
    const int     min_seed_a   = g_min_seed_len_A;
    const int     phase1_cap   = g_phase1_cap;

    const int read_len_bases = read_len; 
    const int E0_base = read_len_bases - 1;

    uint64_t l_f = 0, h_f = 0, diff_f = 0;
    uint64_t l_r = 0, h_r = 0, diff_r = 0;
    uint8_t tail_f = 0, tail_r = 0;

#ifdef ENABLE_F_RC_CHOICE
    short_read_len = E0_base;
    {
        int seed_end_base = E0_base;
        uint64_t addr_f = compute_jumpN_from_base4(pat_f4, seed_end_base, JT_LEN_NT);
        tail_f = 0;
        uint64_t jp_f = jump_pointers[addr_f];
        extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
    }

    short_read_len = E0_base;
    {
        int seed_end_base = E0_base;
        uint64_t addr_r = compute_jumpN_from_base4(pat_rc4, seed_end_base, JT_LEN_NT);
        tail_r = 0;
        uint64_t jp_r = jump_pointers[addr_r];
        extract_jump_bounds(jp_r, &l_r, &h_r, &diff_r);
    }

    const int use_rc = (diff_r > 0 && (diff_f == 0 || diff_r < diff_f));
    *chosen_strand = use_rc;
    uint8_t *pat = use_rc ? pat_rc4 : pat_f4;

    uint64_t l0 = use_rc ? l_r : l_f;
    uint64_t h0 = use_rc ? h_r : h_f;
    uint64_t d0 = h0 - l0;
    uint8_t  t0 = use_rc ? tail_r : tail_f;
#else
    short_read_len = E0_base;
    uint64_t addr_f;
    {
        int seed_end_base = E0_base;
        addr_f = compute_jumpN_from_base4(pat_f4, seed_end_base, JT_LEN_NT);
        tail_f = 0;
        uint64_t jp_f = jump_pointers[addr_f];
        extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
    }
    const int use_rc = 0;
    *chosen_strand = 0;
    uint8_t *pat = pat_f4;
    uint64_t l0 = l_f;
    uint64_t h0 = h_f;
    uint64_t d0 = h0 - l0;
    uint8_t  t0 = tail_f;
#endif

    int pivot_base = E0_base;
    int first_iter = 1;
    uint64_t l = 0, h = 0;
    uint64_t diff = 0;
    uint8_t tail_b = 0;

    while (pivot_base >= JT_LEN_NT - 1) {
        if ((pivot_base + 1) < min_seed_a) {
            break;
        }

        if (first_iter) {
            l = l0; h = h0; diff = d0; tail_b = t0;
            first_iter = 0;
        } else {
            short_read_len = pivot_base;
            tail_b = 0;
            int seed_end_base = pivot_base;
            uint32_t addr = compute_jumpN_from_base4(pat, seed_end_base, JT_LEN_NT);
            uint64_t jp = jump_pointers[addr];
            extract_jump_bounds(jp, &l, &h, &diff);
        }

        int seed_end_base   = pivot_base;
        int seed_len_bases  = JT_LEN_NT;
        int read_idx_base   = pivot_base - JT_LEN_NT;
        uint64_t leftmost_sym_tmp = 0;
        uint8_t  leftmost_off_tmp = 0;

        if (diff == 1) {
#ifdef ENABLE_COUNTERS
            counter_readunique_after1stjump++;
#endif
            uint64_t ref_pos = l - 1;
            uint64_t ref_after = 0;
            int mismatch_base = -1, new_read_idx = read_idx_base;
            uint64_t ref_after_tmp = ref_pos;
            leftmost_sym_tmp = ref_pos + 1;
            leftmost_off_tmp = 0;

            if (read_idx_base >= 0) {
                int matched = ref_walk_pairwise_back_b4(
                    pat, read_idx_base, ref_pos,
                    &mismatch_base, &ref_after_tmp, &new_read_idx,
                    &leftmost_sym_tmp, &leftmost_off_tmp);

                seed_len_bases += matched;
                ref_after = ref_after_tmp;
                read_idx_base = new_read_idx;
            } else {
                ref_after = ref_pos;
            }

            int emitted = 0;
            if (seed_len_bases >= min_seed_a) {
                emitted = 1;
                SMEM *s = &matchArray[*numberofSMEMs];
                s->rid = *read_counter;

                int start_idx = seed_end_base - (seed_len_bases - 1);
                if (start_idx < 0) start_idx = 0;
                int end_idx   = seed_end_base;
                if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

                if(*chosen_strand==0){
                    s->start_idx = start_idx;
                    s->end_idx   = end_idx;
                }else{
                    int f_start = (read_len - 1) - end_idx;
                    int f_end   = (read_len - 1) - start_idx;
                    s->start_idx = f_start;
                    s->end_idx   = f_end;
                }
                s->unique_flag = true;
                s->unique_ref_pointer = ref_after + 1;
                s->unique_leftmost_twostep_pos = leftmost_sym_tmp;
                s->unique_leftmost_base_off    = leftmost_off_tmp;
                s->low_ptr_read = l;
                s->low_ptr_rc_read = INT_MAX;
                s->smem_score = 1;
                s->seed_strand = (uint8_t)(*chosen_strand);

                (*numberofSMEMs)++;
            } else {
                DEBUG_PRINTF("[phase1] unique seed too short: len=%d < MIN %d\n",
                             seed_len_bases, min_seed_a);
            }

            pivot_base = next_pivot_after_seed_base(seed_end_base,
                                                    seed_len_bases,
                                                    emitted);

            if (pivot_base < JT_LEN_NT - 1) break;
            if (pivot_base + 1 < min_seed_a){
                break;
            }
            continue;
        }

        uint64_t ref_after = 0;

        while (1)
        {
            int step_ok = 0;
            uint64_t L_before = l, H_before = h;
            int idx_before = read_idx_base;

            step_ok = fm_step2_b4(pat, &read_idx_base, &l, &h);

        #ifdef TRY_FOUR_CASES_ON_FAIL
            if (!step_ok && idx_before > 0) {
                uint64_t bestL = L_before, bestH = H_before;
                int rescue = fm_rescue_left_b4(pat, idx_before,
                                            L_before, H_before,
                                            &bestL, &bestH);
                if (rescue) {
                    l = bestL;
                    h = bestH;
                    read_idx_base = idx_before - 1;
                    seed_len_bases += 1;
                    step_ok = 1;

                    DEBUG_PRINTF("[phase1] rescue OK: L=%lu H=%lu read_idx=%d seed_len=%d\n",
                                l, h, read_idx_base, seed_len_bases);
                    break;                        
                }
            }
        #endif

            if (!step_ok) {
                l = L_before;
                h = H_before;
                read_idx_base = idx_before;
                break;
            }

            seed_len_bases += 2;
            diff = h - l;

            if (diff == 1) {
                uint64_t ref_pos = reconstruct_ref_base_index(l);
                int mismatch_base = -1, new_read_idx = read_idx_base;
                uint64_t ref_after_tmp = ref_pos;
                leftmost_sym_tmp = ref_pos + 1;
                leftmost_off_tmp = 0;
                
                int matched = ref_walk_pairwise_back_b4(
                    pat, read_idx_base, ref_pos,
                    &mismatch_base, &ref_after_tmp, &new_read_idx,
                    &leftmost_sym_tmp, &leftmost_off_tmp);

                read_idx_base = new_read_idx;
                seed_len_bases += matched;
                ref_after = ref_after_tmp;
                break;
            }

            if (read_idx_base < 0)
                break;
        }

        diff = h - l;
        if (diff == 1 && ref_after == 0) {
            ref_after = reconstruct_ref_base_index(l);
            leftmost_sym_tmp = ref_after + 1;
            leftmost_off_tmp = 0;
        }

        int emitted = 0;
        int interval = (int)(h - l);

        if (seed_len_bases >= min_seed_a && l < h
            && (phase1_cap == 0 || interval <= phase1_cap)
        ) {
            emitted = 1;
            SMEM *s = &matchArray[*numberofSMEMs];
            s->rid = *read_counter;

            int start_idx = seed_end_base - (seed_len_bases - 1);
            if (start_idx < 0) start_idx = 0;
            int end_idx   = seed_end_base;
            if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

            if(*chosen_strand==0){
                s->start_idx = start_idx;
                s->end_idx   = end_idx;
            }else{
                int f_start = (read_len - 1) - end_idx;
                int f_end   = (read_len - 1) - start_idx;
                s->start_idx = f_start;
                s->end_idx   = f_end;
            }
            s->unique_flag = (h - l == 1);
            s->unique_ref_pointer = (h - l == 1) ? (ref_after + 1) : 0;
            s->unique_leftmost_twostep_pos = (h - l == 1) ? leftmost_sym_tmp : 0;
            s->unique_leftmost_base_off    = (h - l == 1) ? leftmost_off_tmp : 0;
            s->low_ptr_read = l;
            s->low_ptr_rc_read = INT_MAX;
            s->smem_score = interval;
            s->seed_strand = (uint8_t)(*chosen_strand);

#ifdef DEBUG_SANITY
            int computed_len = s->end_idx - s->start_idx + 1;
            if (computed_len != seed_len_bases) {
                fprintf(stderr,
                        "[SANITY-FAIL] rid=%d pivot_base=%d start=%d end=%d "
                        "expected_len=%d got_len=%d\n",
                        s->rid, pivot_base, s->start_idx, s->end_idx,
                        seed_len_bases, computed_len);
                exit(EXIT_FAILURE);
            }
#endif
            (*numberofSMEMs)++;
            DEBUG_PRINTF("[phase1] EMIT seed: m=%d n=%d diff=%d l=%lu h=%lu\n",
                         s->start_idx, s->end_idx, interval, s->low_ptr_read, h);
        } else {
            DEBUG_PRINTF("[phase1] NO-EMIT seed_len=%d interval=%d L=%lu H=%lu\n",
                         seed_len_bases, interval, (unsigned long)l, (unsigned long)h);
        }

        pivot_base = next_pivot_after_seed_base(seed_end_base,
                                                seed_len_bases,
                                                emitted);
        if (pivot_base < JT_LEN_NT - 1) break;
        if (pivot_base + 1 < min_seed_a){
            break;
        }
    } /* while pivot_base */
}

/* =======================================================================
   COROUTINE BATCH PHASE I
   =======================================================================

   Design principle: hide BWT RAM latency by interleaving FM steps across
   RS_BATCH reads.  After every single fm_step2_b4() call we know the new
   l and h — we immediately prefetch cp_occ[l>>5] and cp_occ[h>>5] for
   THAT read, then move to the next read's step.  
   ======================================================================= */

#include "phaseI_slot_state.h"

/* -----------------------------------------------------------------------
   Per-pivot working state for one read in the coroutine interleaver.
   ----------------------------------------------------------------------- */
typedef struct {
    const uint8_t *pat;           /* chosen-strand base-4 read pointer    */
    int            read_len;
    int            chosen_strand;
    int            read_counter;

    /* --- pivot-level state --- */
    int      pivot_base;          
    int      first_iter;          /* 1 = still on first pivot of this read  */
    uint64_t l0, h0, d0;         

    /* --- seed-extension state (per pivot) --- */
    int      seed_end_base;
    int      seed_len_bases;
    int      read_idx_base;
    uint64_t l, h, diff;
    uint64_t ref_after;
    uint64_t leftmost_sym_tmp;
    uint8_t  leftmost_off_tmp;

    /* --- phase within this pivot --- */
    /*  0 = need to do jump / check unique
        1 = inside FM extension loop (non-unique path)
        2 = pivot fully processed, advance to next pivot   */
    int      phase;

    /* --- unique path needs ref_walk, handled immediately --- */
    /* non-unique: fm_loop_done signals the while(1) has exited */
    int      fm_loop_done;

    /* --- output --- */
    SMEM    *matchArray;
    uint64_t *numberofSMEMs;

    /* --- config (cached from globals) --- */
    int      min_seed_a;
    int      phase1_cap;
    int      read_len_bases;
} PhaseI_PivotState;

/* -----------------------------------------------------------------------
   phaseI_slot_init()
   Set up a slot for one read.  Does strand probe + first jump (same as
   phaseI_routine's step 0) and prefetches the first cp_occ block.
   ----------------------------------------------------------------------- */
static void phaseI_slot_init(
    PhaseI_PivotState *st,
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
    st->ref_after      = 0;
    st->leftmost_sym_tmp = 0;
    st->leftmost_off_tmp = 0;

    const int E0_base = read_len - 1;

    /* ---- strand probe + first jump — same logic as phaseI_routine ---- */
    uint64_t l_f = 0, h_f = 0, diff_f = 0;
    uint64_t l_r = 0, h_r = 0, diff_r = 0;

#ifdef ENABLE_F_RC_CHOICE
    {
        uint64_t addr_f = compute_jumpN_from_base4(f4, E0_base, JT_LEN_NT);
        uint64_t jp_f   = jump_pointers[addr_f];
        extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
    }
    {
        uint64_t addr_r = compute_jumpN_from_base4(rc4, E0_base, JT_LEN_NT);
        uint64_t jp_r   = jump_pointers[addr_r];
        extract_jump_bounds(jp_r, &l_r, &h_r, &diff_r);
    }
    const int use_rc = (diff_r > 0 && (diff_f == 0 || diff_r < diff_f));
    st->chosen_strand = use_rc;
    st->pat  = use_rc ? rc4 : f4;
    st->l0   = use_rc ? l_r : l_f;
    st->h0   = use_rc ? h_r : h_f;
    st->d0   = st->h0 - st->l0;
#else
    {
        uint64_t addr_f = compute_jumpN_from_base4(f4, E0_base, JT_LEN_NT);
        uint64_t jp_f   = jump_pointers[addr_f];
        extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
    }
    st->chosen_strand = 0;
    st->pat  = f4;
    st->l0   = l_f;
    st->h0   = h_f;
    st->d0   = h_f - l_f;
#endif

    st->pivot_base = E0_base;
    st->l = st->l0;
    st->h = st->h0;
    st->diff = st->d0;

    /* prefetch the first cp_occ block we'll need */
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
   phaseI_slot_emit_seed()
   Emit a completed seed into matchArray. 
   ----------------------------------------------------------------------- */
static inline void phaseI_slot_emit_seed(
    PhaseI_PivotState *st,
    int seed_end_base,
    int seed_len_bases,
    int interval,         /* (int)(h-l) — 0 means unique (use 1)   */
    int is_unique)
{
    const int read_len      = st->read_len;
    const int read_len_bases = st->read_len_bases;

    SMEM *s = &st->matchArray[*st->numberofSMEMs];
    s->rid = st->read_counter;

    int start_idx = seed_end_base - (seed_len_bases - 1);
    if (start_idx < 0) start_idx = 0;
    int end_idx = seed_end_base;
    if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

    if (st->chosen_strand == 0) {
        s->start_idx = start_idx;
        s->end_idx   = end_idx;
    } else {
        int f_start = (read_len - 1) - end_idx;
        int f_end   = (read_len - 1) - start_idx;
        s->start_idx = f_start;
        s->end_idx   = f_end;
    }

    if (is_unique) {
        s->unique_flag = true;
        s->unique_ref_pointer          = st->ref_after + 1;
        s->unique_leftmost_twostep_pos = st->leftmost_sym_tmp;
        s->unique_leftmost_base_off    = st->leftmost_off_tmp;
        s->smem_score  = 1;
    } else {
        s->unique_flag = (st->h - st->l == 1);
        s->unique_ref_pointer          = (st->h - st->l == 1) ? (st->ref_after + 1) : 0;
        s->unique_leftmost_twostep_pos = (st->h - st->l == 1) ? st->leftmost_sym_tmp : 0;
        s->unique_leftmost_base_off    = (st->h - st->l == 1) ? st->leftmost_off_tmp : 0;
        s->smem_score  = interval;
    }

    s->low_ptr_read    = st->l;
    s->low_ptr_rc_read = INT_MAX;
    s->seed_strand     = (uint8_t)(st->chosen_strand);

#ifdef DEBUG_SANITY
    int computed_len = s->end_idx - s->start_idx + 1;
    if (computed_len != seed_len_bases) {
        fprintf(stderr,
                "[SANITY-FAIL] rid=%d pivot_base=%d start=%d end=%d "
                "expected_len=%d got_len=%d\n",
                st->read_counter, st->pivot_base,
                s->start_idx, s->end_idx,
                seed_len_bases, computed_len);
        exit(EXIT_FAILURE);
    }
#endif

    (*st->numberofSMEMs)++;
}

/* -----------------------------------------------------------------------
   phaseI_slot_step()

   Advance one slot by exactly ONE fm_step2_b4() call (or complete one
   non-FM action: jump, unique-branch ref_walk, or pivot advance).

   Returns:
     1  — slot still has work to do, caller should continue scheduling it
     0  — slot is done (pivot_base went negative / pruned), caller removes it
   ----------------------------------------------------------------------- */
static int phaseI_slot_step(PhaseI_PivotState *st)
{
    const int min_seed_a  = st->min_seed_a;
    const int phase1_cap  = st->phase1_cap;
    const int read_len_bases = st->read_len_bases;
    const uint8_t *pat    = st->pat;

    /* ================================================================
       PHASE 0: do jump for current pivot, classify unique vs non-unique
       ================================================================ */
    if (st->phase == 0) {

        int pivot_base = st->pivot_base;

        /* prune check */
        if (pivot_base < JT_LEN_NT - 1 || (pivot_base + 1) < min_seed_a) {
            st->pivot_base = -1;
            return 0;  /* slot done */
        }

        /* jump */
        if (st->first_iter) {
            st->l    = st->l0;
            st->h    = st->h0;
            st->diff = st->d0;
            st->first_iter = 0;
        } else {
            uint64_t addr = compute_jumpN_from_base4(pat, pivot_base, JT_LEN_NT);
            uint64_t jp   = jump_pointers[addr];
            extract_jump_bounds(jp, &st->l, &st->h, &st->diff);
        }

        st->seed_end_base   = pivot_base;
        st->seed_len_bases  = JT_LEN_NT;
        st->read_idx_base   = pivot_base - JT_LEN_NT;
        st->ref_after       = 0;
        st->leftmost_sym_tmp = 0;
        st->leftmost_off_tmp = 0;
        st->fm_loop_done    = 0;

        /* ---- unique immediately at jump ---- */
        if (st->diff == 1) {
#ifdef ENABLE_COUNTERS
            counter_readunique_after1stjump++;
#endif
            /* ref_walk is cheap (sequential ref memory, well-prefetched
               by hardware), so we do it in full here rather than splitting
               it across scheduling rounds.                               */
            uint64_t ref_pos = st->l - 1;
            st->ref_after    = ref_pos;
            st->leftmost_sym_tmp = ref_pos + 1;
            st->leftmost_off_tmp = 0;

            if (st->read_idx_base >= 0) {
                int mismatch_base = -1, new_read_idx = st->read_idx_base;
                uint64_t ref_after_tmp = ref_pos;

                int matched = ref_walk_pairwise_back_b4(
                    (uint8_t *)pat, st->read_idx_base, ref_pos,
                    &mismatch_base, &ref_after_tmp, &new_read_idx,
                    &st->leftmost_sym_tmp, &st->leftmost_off_tmp);

                st->seed_len_bases += matched;
                st->ref_after       = ref_after_tmp;
                st->read_idx_base   = new_read_idx;
            }

            int emitted = 0;
            if (st->seed_len_bases >= min_seed_a) {
                emitted = 1;
                phaseI_slot_emit_seed(st,
                    st->seed_end_base, st->seed_len_bases,
                    1 /*interval — unused for unique*/, 1 /*is_unique*/);
            }

            st->pivot_base = next_pivot_after_seed_base(
                st->seed_end_base, st->seed_len_bases, emitted);

            /* check exit conditions — same as phaseI_routine */
            if (st->pivot_base < JT_LEN_NT - 1) { st->pivot_base = -1; return 0; }
            if (st->pivot_base + 1 < min_seed_a) { st->pivot_base = -1; return 0; }

            /* stay in phase 0 for next pivot */
            st->phase = 0;

            /* prefetch the jump entry for the NEXT pivot so it's warm */
            if (st->pivot_base >= JT_LEN_NT - 1) {
                uint64_t next_addr =
                    compute_jumpN_from_base4(pat, st->pivot_base, JT_LEN_NT);
                __builtin_prefetch(&jump_pointers[next_addr], 0, 1);
            }
            return 1;
        }

        /* non-unique: transition to FM extension phase */
        st->phase = 1;

        /* prefetch cp_occ for the FIRST fm_step2_b4 call we're about to make
           next time this slot is scheduled */
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5] + 64, 0, 1);
        return 1;
    }

    /* ================================================================
       PHASE 1: execute ONE fm_step2_b4 step of the non-unique path
       ================================================================ */
    if (st->phase == 1) {

        if (st->fm_loop_done) {
            /* FM loop already finished — go to emit phase */
            st->phase = 2;
            return 1;
        }

        uint64_t L_before = st->l, H_before = st->h;
        int      idx_before = st->read_idx_base;

        int step_ok = fm_step2_b4(pat, &st->read_idx_base, &st->l, &st->h);

#ifdef TRY_FOUR_CASES_ON_FAIL
        if (!step_ok && idx_before > 0) {
            uint64_t bestL = L_before, bestH = H_before;
            int rescue = fm_rescue_left_b4(pat, idx_before,
                                           L_before, H_before,
                                           &bestL, &bestH);
            if (rescue) {
                st->l = bestL;
                st->h = bestH;
                st->read_idx_base = idx_before - 1;
                st->seed_len_bases += 1;
                /* rescue exits the FM loop — same as original */
                st->fm_loop_done = 1;
                st->phase = 2;
                return 1;
            }
        }
#endif

        if (!step_ok) {
            /* rollback — identical to original */
            st->l = L_before;
            st->h = H_before;
            st->read_idx_base = idx_before;
            st->fm_loop_done = 1;
            st->phase = 2;
            return 1;
        }

        /* successful step */
        st->seed_len_bases += 2;
        st->diff = st->h - st->l;

        if (st->diff == 1) {
            /* became unique mid-extension — do ref_walk immediately
               (same logic as original, sequential memory access) */
            uint64_t ref_pos = reconstruct_ref_base_index(st->l);
            int mismatch_base = -1, new_read_idx = st->read_idx_base;
            uint64_t ref_after_tmp = ref_pos;
            st->leftmost_sym_tmp = ref_pos + 1;
            st->leftmost_off_tmp = 0;

            int matched = ref_walk_pairwise_back_b4(
                (uint8_t *)pat, st->read_idx_base, ref_pos,
                &mismatch_base, &ref_after_tmp, &new_read_idx,
                &st->leftmost_sym_tmp, &st->leftmost_off_tmp);

            st->read_idx_base   = new_read_idx;
            st->seed_len_bases += matched;
            st->ref_after       = ref_after_tmp;

            st->fm_loop_done = 1;
            st->phase = 2;
            return 1;
        }

        if (st->read_idx_base < 0) {
            st->fm_loop_done = 1;
            st->phase = 2;
            return 1;
        }

        /* more FM steps to do — prefetch the NEXT step's cp_occ blocks
           before yielding so the scheduler can work on other slots      */
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[st->h >> 5] + 64, 0, 1);

        /* stay in phase 1 for the next FM step */
        return 1;
    }

    /* ================================================================
       PHASE 2: FM loop done — reconstruct ref_after if needed, emit seed,
                advance to next pivot
       ================================================================ */
    if (st->phase == 2) {

        st->diff = st->h - st->l;

        /* reconstruct ref_after if unique but not yet set */
        if (st->diff == 1 && st->ref_after == 0) {
            st->ref_after        = reconstruct_ref_base_index(st->l);
            st->leftmost_sym_tmp = st->ref_after + 1;
            st->leftmost_off_tmp = 0;
        }

        int emitted  = 0;
        int interval = (int)(st->h - st->l);

        if (st->seed_len_bases >= min_seed_a
            && st->l < st->h
            && (phase1_cap == 0 || interval <= phase1_cap))
        {
            emitted = 1;
            phaseI_slot_emit_seed(st,
                st->seed_end_base, st->seed_len_bases,
                interval, 0 /*not unique path*/);
        }

        st->pivot_base = next_pivot_after_seed_base(
            st->seed_end_base, st->seed_len_bases, emitted);

        if (st->pivot_base < JT_LEN_NT - 1) { st->pivot_base = -1; return 0; }
        if (st->pivot_base + 1 < min_seed_a) { st->pivot_base = -1; return 0; }

        /* prefetch jump entry for NEXT pivot */
        if (st->pivot_base >= JT_LEN_NT - 1) {
            uint64_t next_addr =
                compute_jumpN_from_base4(pat, st->pivot_base, JT_LEN_NT);
            __builtin_prefetch(&jump_pointers[next_addr], 0, 1);
        }

        /* reset for next pivot */
        st->phase        = 0;
        st->fm_loop_done = 0;
        st->ref_after    = 0;
        return 1;
    }

    return 0; /* unreachable */
}

/* -----------------------------------------------------------------------
   phaseI_batch_interleaved()

   Entry point called by rosaseed_core_bridge_batched.cpp 
   ----------------------------------------------------------------------- */
void phaseI_batch_interleaved(
    const uint8_t **f4,
    const uint8_t **rc4,
    const int      *read_lens,
    const int      *read_ctrs,
    SMEM          **matchArrays,
    uint64_t      **nSMEMs,
    int            *chosen_strands,
    int             nreads)
{
    PhaseI_PivotState slots[RS_BATCH];

    /* initialise all slots */
    for (int s = 0; s < nreads; ++s) {
        phaseI_slot_init(&slots[s],
                         f4[s], rc4[s],
                         read_lens[s],
                         read_ctrs[s],
                         matchArrays[s],
                         nSMEMs[s]);
        chosen_strands[s] = slots[s].chosen_strand;
    }

    /* round-robin scheduler */
    int active = nreads;
    int alive[RS_BATCH];
    for (int s = 0; s < nreads; ++s) alive[s] = 1;

    while (active > 0) {
        for (int s = 0; s < nreads; ++s) {
            if (!alive[s]) continue;

            int still_going = phaseI_slot_step(&slots[s]);

            if (!still_going) {
                alive[s] = 0;
                --active;
            }
        }
    }

    /* write chosen_strand back from slots (init already did it but
       phaseI_slot_step never changes it, so this is a no-op safety copy) */
    for (int s = 0; s < nreads; ++s)
        chosen_strands[s] = slots[s].chosen_strand;
}

void phaseI_probe_first_pivot(PhaseI_SlotState *slot)
{
    const int read_len  = slot->read_len;
    const int E0_base   = read_len - 1;

    uint64_t l_f = 0, h_f = 0, diff_f = 0;
    uint64_t l_r = 0, h_r = 0, diff_r = 0;

#ifdef ENABLE_F_RC_CHOICE
    {
        uint64_t addr_f = compute_jumpN_from_base4(slot->f4, E0_base, JT_LEN_NT);
        uint64_t jp_f   = jump_pointers[addr_f];
        extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
    }
    {
        uint64_t addr_r = compute_jumpN_from_base4(slot->rc4, E0_base, JT_LEN_NT);
        uint64_t jp_r   = jump_pointers[addr_r];
        extract_jump_bounds(jp_r, &l_r, &h_r, &diff_r);
    }
    const int use_rc = (diff_r > 0 && (diff_f == 0 || diff_r < diff_f));
    slot->chosen_strand = use_rc;
    slot->l0 = use_rc ? l_r : l_f;
    slot->h0 = use_rc ? h_r : h_f;
    slot->d0 = slot->h0 - slot->l0;
#else
    {
        uint64_t addr_f = compute_jumpN_from_base4(slot->f4, E0_base, JT_LEN_NT);
        uint64_t jp_f   = jump_pointers[addr_f];
        extract_jump_bounds(jp_f, &l_f, &h_f, &diff_f);
    }
    slot->chosen_strand = 0;
    slot->l0 = l_f;
    slot->h0 = h_f;
    slot->d0 = h_f - l_f;
#endif

    const uint64_t blk_l = slot->l0 >> 5;
    const uint64_t blk_h = slot->h0 >> 5;

    __builtin_prefetch((const char *)&cp_occ[blk_l],     0, 1);
    __builtin_prefetch((const char *)&cp_occ[blk_l] + 64, 0, 1);

    if (blk_h != blk_l) {
        __builtin_prefetch((const char *)&cp_occ[blk_h],     0, 1);
        __builtin_prefetch((const char *)&cp_occ[blk_h] + 64, 0, 1);
    }
}
