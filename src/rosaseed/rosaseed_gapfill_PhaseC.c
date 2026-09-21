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

#include "gap_fill_phase.h"
#include "seeding_kernel_mt.h"
#include "file_dec.h"
#include "macros.h"
#include "bwa.h"
#include "helper_functions.h"

#ifdef DEBUG_GAPFILL
  #define GF_PRINT(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)
#else
  #define GF_PRINT(fmt, ...) do {} while(0)
#endif

#define DEBUG_PRINTF_GAPFILLPHASE GF_PRINT

static inline __attribute__((always_inline))
uint64_t compute_jumpN_from_base4(const uint8_t *pat4,
                                   int seed_end_base,
                                   int jump_len_nt)
{
    uint64_t addr = 0;
    int b = seed_end_base;
    for (int k = 0; k < jump_len_nt; ++k, --b)
        addr |= ((uint64_t)pat4[b]) << (2 * k);
    return addr;
}


static inline __attribute__((always_inline))
int fm_step2_b4(const uint8_t *pat4,
                int *read_idx_base,
                uint64_t *l,
                uint64_t *h)
{
    int i = *read_idx_base;
    if (i <= 0) return 0;

    base16_t sym = (base16_t)((pat4[i-1] << 2) | pat4[i]);
    uint64_t L2 = *l, H2 = *h;
    fm_index_mapping_backward_search(&L2, &H2, sym);
    if (L2 >= H2) return 0;

    *l = L2; *h = H2;
    *read_idx_base = i - 2;
    return 1;
}

int ref_walk_pairwise_back_b4(uint8_t *pat,
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
    *mismatch_base  = -1;
    *reached_target = 0;
    uint64_t r = ref_pos_rightbase;
    int matched = 0;

    uint64_t best_sym = r + 1;
    uint8_t  best_off = 0;

    int i = left_base_idx;
    for (; i >= 1; i -= 2) {
        uint8_t sym16 = REF_AT(r);
        uint8_t rbase = sym16 & 0x3;
        uint8_t lbase = (sym16 >> 2) & 0x3;

        // RIGHT base
        if (rbase != pat[i]) {
            *mismatch_base  = (int)r;
            *ref_after      = r;
            *new_read_idx_base = i;
            *leftmost_sym   = best_sym;
            *leftmost_off   = best_off;
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 1;

        if (need_at_least > 0 && matched >= need_at_least) {
            *reached_target = 1;
            *ref_after      = r;
            *new_read_idx_base = i - 1;
            *leftmost_sym   = best_sym;
            *leftmost_off   = best_off;
            return matched;
        }

        // LEFT base
        if (lbase != pat[i-1]) {
            *mismatch_base  = (int)r;
            *ref_after      = r;
            *new_read_idx_base = i - 1;
            *leftmost_sym   = best_sym;
            *leftmost_off   = best_off;
            return matched;
        }

        ++matched;
        best_sym = r;
        best_off = 0;

        if (need_at_least > 0 && matched >= need_at_least) {
            *reached_target = 1;
            *ref_after      = r;
            *new_read_idx_base = i - 2;
            *leftmost_sym   = best_sym;
            *leftmost_off   = best_off;
            return matched;
        }

        --r;
    }

    if (i == -1) {
        *new_read_idx_base = i;
        *ref_after = r;
        *leftmost_sym = best_sym;
        *leftmost_off = best_off;
        if (need_at_least > 0 && matched >= need_at_least) *reached_target = 1;
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
            if (need_at_least > 0 && matched >= need_at_least) *reached_target = 1;
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
    if (need_at_least > 0 && matched >= need_at_least) *reached_target = 1;
    return matched;
}

/* ============================================================
   run_single_pivot()
   pivot_base    : index into pat[] — chosen-strand coordinate
   active_strand : 0=forward pat_f4, 1=RC pat_rc4
   ============================================================ */
static inline __attribute__((always_inline))
int run_single_pivot(
    const uint8_t  *pat,
    int             pivot_base,
    SMEM           *matchArray,
    uint64_t       *total_smem,
    int             rid,
    int64_t         min_intv,
    int             active_strand,
    int             read_len_bases)
{
    const int min_seed_bc = g_min_seed_len_BC;  
    if (pivot_base < JT_LEN_NT - 1)            return 0;
    if (pivot_base + 1 < min_seed_bc) return 0;

    uint64_t addr = compute_jumpN_from_base4(pat, pivot_base, JT_LEN_NT);
    uint64_t jp   = jump_pointers[addr];

    uint64_t l, h, diff;
    extract_jump_bounds(jp, &l, &h, &diff);
    if (l >= h) return 0;

    GF_PRINT("          jt_interval=%-8lu", (unsigned long)diff);

    int      seed_end_base  = pivot_base;
    int      seed_len_bases = JT_LEN_NT;
    int      read_idx_base  = pivot_base - JT_LEN_NT;
    uint64_t ref_after      = 0;
    uint64_t leftmost_sym_tmp;
    uint8_t  leftmost_off_tmp;

    if (diff == 1) {
        uint64_t ref_pos = l - 1;
        if (read_idx_base >= 0) {
            int mismatch_base = -1, new_read_idx = read_idx_base, reached = 0;
            uint64_t ref_after_tmp = ref_pos;
            leftmost_sym_tmp = ref_pos + 1;
            leftmost_off_tmp = 0;

            int need = min_seed_bc - seed_len_bases;
            if (need < 0) need = 0;
            
            int matched = ref_walk_pairwise_back_b4(
                (uint8_t*)pat, read_idx_base, ref_pos,
                &mismatch_base, &ref_after_tmp, &new_read_idx,
                need, &reached,
                &leftmost_sym_tmp, &leftmost_off_tmp);

            seed_len_bases += matched;
            read_idx_base   = new_read_idx;
            ref_after       = ref_after_tmp;
        } else {
            ref_after = ref_pos;
        }
        diff = h - l;
    } else {
        while (1) {
            uint64_t L_before = l, H_before = h;
            int      idx_before = read_idx_base;

            int step_ok = fm_step2_b4(pat, &read_idx_base, &l, &h);

            if (!step_ok) {
                l = L_before; h = H_before; read_idx_base = idx_before;
                break;
            }

            seed_len_bases += 2;
            diff = h - l;

            if (diff < (uint64_t)min_intv &&
                seed_len_bases >= min_seed_bc)
                break;

            if (diff == 1) {
                uint64_t ref_pos = reconstruct_ref_base_index(l);
                
                int mismatch_base = -1, new_read_idx = read_idx_base, reached = 0;
                uint64_t ref_after_tmp = ref_pos;
                leftmost_sym_tmp = ref_pos + 1;
                leftmost_off_tmp = 0;
                
                int need = min_seed_bc - seed_len_bases;
                if (need < 0) need = 0;
                
                int matched = ref_walk_pairwise_back_b4(
                    (uint8_t*)pat, read_idx_base, ref_pos,
                    &mismatch_base, &ref_after_tmp, &new_read_idx,
                    need, &reached,
                    &leftmost_sym_tmp, &leftmost_off_tmp);

                seed_len_bases += matched;
                read_idx_base   = new_read_idx;
                ref_after       = ref_after_tmp;
                break;
            }

            if (read_idx_base < 0) break;
        }
        diff = h - l;
    }

    if (seed_len_bases  < min_seed_bc) return 0;
    if (l >= h)                                   return 0;
    if ((h - l) >= (uint64_t)min_intv)            return 0;

    /* If interval is 1 but ref_after was never set (rollback case),
    reconstruct it now from the SA entry.                          */
    if (diff == 1 && ref_after == 0) {
        ref_after = reconstruct_ref_base_index(l);
        leftmost_sym_tmp = ref_after + 1;
        leftmost_off_tmp = 0;
    }

    int start_idx = seed_end_base - (seed_len_bases - 1);
    if (start_idx < 0) start_idx = 0;
    int end_idx = seed_end_base;
    if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

    int fwd_start, fwd_end;
    if (active_strand == 0) {
        fwd_start = start_idx;
        fwd_end   = end_idx;
    } else {
        fwd_start = (read_len_bases - 1) - end_idx;
        fwd_end   = (read_len_bases - 1) - start_idx;
    }

    GF_PRINT(" → EMIT [%d..%d] s=%d\n", fwd_start, fwd_end, (int)(h - l));

    SMEM *s               = &matchArray[*total_smem];
    s->rid                = rid;
    s->start_idx          = fwd_start;
    s->end_idx            = fwd_end;
    s->unique_flag        = (h - l == 1);
    s->unique_ref_pointer = (h - l == 1) ? (ref_after + 1) : 0;
    s->unique_leftmost_twostep_pos = (h - l == 1) ? leftmost_sym_tmp : 0;
    s->unique_leftmost_base_off    = (h - l == 1) ? leftmost_off_tmp : 0;
    s->low_ptr_read       = l;
    s->low_ptr_rc_read    = INT_MAX;
    s->smem_score         = (int)(h - l);
    s->seed_strand        = (uint8_t)(active_strand);

    (*total_smem)++;
    return 1;
}

/* ============================================================
   gap_covered()
   -------------
   Returns 1 if any seed in matchArray[0..n_seeds) overlaps
   the forward-coordinate gap [gap_start..gap_end].
   A seed overlaps if: seed.start <= gap_end AND seed.end >= gap_start
   ============================================================ */
static inline int gap_covered(
    const SMEM *matchArray,
    uint64_t    n_seeds,
    int         gap_start,
    int         gap_end)
{
    for (uint64_t k = 0; k < n_seeds; ++k) {
        // if (matchArray[k].start_idx <= gap_end &&
        //     matchArray[k].end_idx   >= gap_start)
        if (matchArray[k].start_idx <= gap_start &&
            matchArray[k].end_idx >= gap_end)
            return 1;
    }
    return 0;
}

/* ============================================================
   run_gap_pivots()
   ----------------
   For a single gap [gap_start..gap_end] in forward coordinates:
     1. Build right-anchored pivot list for this gap
     2. Run on pat (active_strand)
     3. Return seeds emitted

   pivot_fwd positions are always in forward coordinates.
   strand_pivot = pivot_fwd         if active_strand=0
   strand_pivot = (read_len-1-pivot_fwd)   if active_strand=1
   ============================================================ */
/* ============================================================
   GapFill pivot coroutine state.
   Holds everything run_single_pivot() needs between s-steps
   so it can be suspended and resumed by the interleaver.
   ============================================================ */
typedef struct {
    int      strand_pivot;
    int      pivot_fwd;        /* original pivots_fwd[i] value */
    int      active;           /* 1 = still running in PASS 3 scheduler  */
    int      valid;            /* 1 = slot was successfully initialised
                                   in PASS 2 and should be considered
                                   for emission in PASS 4. Explicitly
                                   set: never read uninitialised.       */

    uint64_t l, h, diff;
    int      seed_end_base;
    int      seed_len_bases;
    int      read_idx_base;
    uint64_t ref_after;
    uint64_t leftmost_sym_tmp;
    uint8_t  leftmost_off_tmp;

    /* phase: 0=needs jump  1=FM extension loop  2=emit        */
    int      phase;
    int      unique_at_jump;   /* 1 = diff==1 after first jump  */
} GFPivotSlot;

/* One s-step for one GapFill pivot, with prefetch for the next step.
   Returns 1 = still running, 0 = done (emit or skip).
   Mirrors run_single_pivot()'s inner while(1) exactly:
     - same rollback on step failure
     - same diff==1 ref_walk path
     - same read_idx_base < 0 exit
     - same min_intv early-exit
   All emit logic happens in run_gap_pivots_coro after this returns 0. */
static inline int gf_pivot_step(
    GFPivotSlot    *sl,
    const uint8_t  *pat,
    int64_t         min_intv,
    int             min_seed_bc)
{
    /* Phase 0: initial jump (called once per pivot) */
    if (sl->phase == 0) {
        /* jump table already done in setup, go straight to s-step */
        if (sl->diff == 1) {
            /* unique at jump, ref_walk handles extension, no s-steps */
            sl->unique_at_jump = 1;
            sl->phase = 2;  /* go to emit */
            return 0;
        }
        sl->unique_at_jump = 0;
        sl->phase = 1;

        /* prefetch cp_occ for first FM step */
        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5] + 64, 0, 1);
        return 1;   
    }

    /* Phase 1: one fm_step2_b4 step */
    if (sl->phase == 1) {
        uint64_t L_before    = sl->l;
        uint64_t H_before    = sl->h;
        int      idx_before  = sl->read_idx_base;

        int step_ok = fm_step2_b4(pat, &sl->read_idx_base, &sl->l, &sl->h);

        if (!step_ok) {
            sl->l = L_before;
            sl->h = H_before;
            sl->read_idx_base = idx_before;
            sl->phase = 2;  /* go to emit check */
            return 0;
        }

        sl->seed_len_bases += 2;
        sl->diff = sl->h - sl->l;

        /* early-exit on specificity threshold*/
        if (sl->diff < (uint64_t)min_intv && sl->seed_len_bases >= min_seed_bc) {
            sl->phase = 2;
            return 0;
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
                (uint8_t *)pat, sl->read_idx_base, ref_pos,
                &mismatch_base, &ref_after_tmp, &new_read_idx,
                need, &reached,
                &sl->leftmost_sym_tmp, &sl->leftmost_off_tmp);

            sl->seed_len_bases += matched;
            sl->read_idx_base   = new_read_idx;
            sl->ref_after       = ref_after_tmp;
            sl->phase = 2;
            return 0;
        }

        /* read exhausted */
        if (sl->read_idx_base < 0) {
            sl->phase = 2;
            return 0;
        }

        /* more s-steps, prefetch next cp_occ before yielding */
        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5] + 64, 0, 1);
        return 1;   
    }

    return 0;
}

static int run_gap_pivots(
    const uint8_t  *pat,
    int             active_strand,
    int             gap_start,
    int             gap_end,
    SMEM           *matchArray,
    uint64_t       *total_smem,
    uint64_t        matchArray_capacity,
    int             rid,
    int64_t         min_intv,
    int             read_len_bases)
{
    const int min_seed_bc   = g_min_seed_len_BC;
    const int gap_left_step = g_gap_left_step;
    const int min_pivot     = min_seed_bc - 1;
    const int gap_size      = gap_end - gap_start + 1;

    const int step = (gap_size <= g_min_seed_len_BC)
                   ? (gap_left_step > 1 ? gap_left_step / 2 : 1)
                   : gap_left_step;

    const int p_min = (active_strand == 0)
                    ? (gap_start > min_pivot ? gap_start : min_pivot)
                    : gap_start;

    const int max_pivots = (read_len_bases / step) + 4;
    int pivots_fwd[max_pivots];
    int n_pivots = 0;

    int p = gap_end;
    while (p >= p_min && n_pivots < max_pivots) {
        pivots_fwd[n_pivots++] = p;
        p -= step;
    }

    if (n_pivots == 0) return 0;

    GFPivotSlot slots[max_pivots];
    int n_active = 0;

    for (int i = 0; i < n_pivots; ++i) {
        slots[i].active          = 0;
        slots[i].valid           = 0;
        slots[i].l               = 0;
        slots[i].h               = 0;
        slots[i].diff            = 0;
        slots[i].seed_len_bases  = 0;
        slots[i].read_idx_base   = 0;
        slots[i].ref_after       = 0;
        slots[i].leftmost_sym_tmp = 0;
        slots[i].leftmost_off_tmp = 0;
        slots[i].phase           = 0;
        slots[i].unique_at_jump  = 0;
    }

    /* PASS 1: prefetch jump table entries upfront */
    for (int i = 0; i < n_pivots; ++i) {
        int strand_pivot = (active_strand == 0)
                         ? pivots_fwd[i]
                         : (read_len_bases - 1) - pivots_fwd[i];
        if (strand_pivot < min_pivot) continue;   /* leave inactive */
        uint64_t addr = compute_jumpN_from_base4(pat, strand_pivot, JT_LEN_NT);
        __builtin_prefetch(&jump_pointers[addr], 0, 1);
        slots[i].active = 1;   /* candidate: still needs PASS 2 validation */
    }

    /* PASS 2: read warm jump entries, init slots, prefetch first cp_occ */
    for (int i = 0; i < n_pivots; ++i) {
        if (!slots[i].active) continue;

        int strand_pivot = (active_strand == 0)
                         ? pivots_fwd[i]
                         : (read_len_bases - 1) - pivots_fwd[i];

        uint64_t addr = compute_jumpN_from_base4(pat, strand_pivot, JT_LEN_NT);
        uint64_t jp   = jump_pointers[addr];

        uint64_t l, h, diff;
        extract_jump_bounds(jp, &l, &h, &diff);

        if (l >= h) { slots[i].active = 0; continue; }   /* leave valid=0 */

        slots[i].valid           = 1;   /* explicitly mark as a real,
                                            emittable candidate for PASS 4 */
        slots[i].strand_pivot    = strand_pivot;
        slots[i].pivot_fwd       = pivots_fwd[i];
        slots[i].l               = l;
        slots[i].h               = h;
        slots[i].diff            = diff;
        slots[i].seed_end_base   = strand_pivot;
        slots[i].seed_len_bases  = JT_LEN_NT;
        slots[i].read_idx_base   = strand_pivot - JT_LEN_NT;
        slots[i].ref_after       = 0;
        slots[i].leftmost_sym_tmp = 0;
        slots[i].leftmost_off_tmp = 0;
        slots[i].phase           = 0;
        slots[i].unique_at_jump  = 0;

        /* prefetch cp_occ for first s-step of this pivot */
        uint64_t blk_l = l >> 5, blk_h = h >> 5;
        __builtin_prefetch((const char *)&cp_occ[blk_l],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[blk_l] + 64, 0, 1);
        if (blk_h != blk_l) {
            __builtin_prefetch((const char *)&cp_occ[blk_h],      0, 1);
            __builtin_prefetch((const char *)&cp_occ[blk_h] + 64, 0, 1);
        }

        n_active++;
    }

    if (n_active == 0) return 0;

    /* PASS 3: round-robin interleaver */
    while (n_active > 0) {
        for (int i = 0; i < n_pivots; ++i) {
            if (!slots[i].active) continue;

            int still_running = gf_pivot_step(
                &slots[i], pat, min_intv, min_seed_bc);

            if (!still_running) {
                slots[i].active = 0;
                n_active--;
            }
        }
    }

    /* PASS 4: emit seeds from completed slots */
    int emitted = 0;
    for (int i = 0; i < n_pivots; ++i) {
        if (*total_smem >= matchArray_capacity - 1) break;

        GFPivotSlot *sl = &slots[i];

        /* Skip slots that were never successfully initialised in PASS 2
           (bad strand_pivot, or l>=h at jump). */
        if (!sl->valid) continue;

        /* handle unique-at-jump: ref_walk needed, done here like run_single_pivot */
        if (sl->unique_at_jump) {
            uint64_t ref_pos = sl->l - 1;
            if (sl->read_idx_base >= 0) {
                int mismatch_base = -1, new_read_idx = sl->read_idx_base, reached = 0;
                uint64_t ref_after_tmp = ref_pos;
                sl->leftmost_sym_tmp = ref_pos + 1;
                sl->leftmost_off_tmp = 0;

                int need = min_seed_bc - sl->seed_len_bases;
                if (need < 0) need = 0;

                int matched = ref_walk_pairwise_back_b4(
                    (uint8_t *)pat, sl->read_idx_base, ref_pos,
                    &mismatch_base, &ref_after_tmp, &new_read_idx,
                    need, &reached,
                    &sl->leftmost_sym_tmp, &sl->leftmost_off_tmp);

                sl->seed_len_bases += matched;
                sl->read_idx_base   = new_read_idx;
                sl->ref_after       = ref_after_tmp;
            } else {
                sl->ref_after = ref_pos;
            }
            sl->diff = sl->h - sl->l;
        }

        /* recompute diff in case it wasn't updated */
        sl->diff = sl->h - sl->l;

        if (sl->seed_len_bases < min_seed_bc) continue;
        if (sl->l >= sl->h)                   continue;
        if ((sl->h - sl->l) >= (uint64_t)min_intv) continue;

        if (sl->diff == 1 && sl->ref_after == 0) {
            sl->ref_after = reconstruct_ref_base_index(sl->l);
            sl->leftmost_sym_tmp = sl->ref_after + 1;
            sl->leftmost_off_tmp = 0;
        }

        int start_idx = sl->seed_end_base - (sl->seed_len_bases - 1);
        if (start_idx < 0) start_idx = 0;
        int end_idx = sl->seed_end_base;
        if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

        int fwd_start, fwd_end;
        if (active_strand == 0) {
            fwd_start = start_idx;
            fwd_end   = end_idx;
        } else {
            fwd_start = (read_len_bases - 1) - end_idx;
            fwd_end   = (read_len_bases - 1) - start_idx;
        }

        SMEM *s               = &matchArray[*total_smem];
        s->rid                = rid;
        s->start_idx          = fwd_start;
        s->end_idx            = fwd_end;
        s->unique_flag        = (sl->h - sl->l == 1);
        s->unique_ref_pointer = (sl->h - sl->l == 1) ? (sl->ref_after + 1) : 0;
        s->unique_leftmost_twostep_pos = (sl->h - sl->l == 1) ? sl->leftmost_sym_tmp : 0;
        s->unique_leftmost_base_off    = (sl->h - sl->l == 1) ? sl->leftmost_off_tmp : 0;
        s->low_ptr_read       = sl->l;
        s->low_ptr_rc_read    = INT_MAX;
        s->smem_score         = (int)(sl->h - sl->l);
        s->seed_strand        = (uint8_t)(active_strand);

        (*total_smem)++;
        emitted++;

#ifdef GAPFILL_EARLY_EXIT
        if (emitted > 0 &&
            gap_covered(matchArray, *total_smem, gap_start, gap_end))
            break;
#endif
    }

    return emitted;
}


/* ============================================================
   process_gap()
   -------------
   Handles ONE gap completely:
     Pass 1: run pivots on chosen strand
     Check:  is gap still uncovered?
     Pass 2: if still uncovered → run same pivot positions on
             opposite strand

   n_seeds_at_entry: number of seeds before this gap was processed
   ============================================================ */
static void process_gap(
    int             gap_start,
    int             gap_end,
    const uint8_t  *pat_chosen,
    int             chosen_strand,
    const uint8_t  *pat_opposite,
    int             opposite_strand,
    SMEM           *matchArray,
    uint64_t       *total_smem,
    uint64_t        matchArray_capacity,
    int             rid,
    int64_t         min_intv,
    int             read_len_bases)
{
    const int gap_thresh = g_gap_threshold;   
    int gap_size = gap_end - gap_start + 1;
    if (gap_size < gap_thresh) return;

    GF_PRINT("  GAP [%d..%d] size=%d\n", gap_start, gap_end, gap_size);

    /* chosen strand */
    uint64_t before_p1 = *total_smem;

    GF_PRINT("    [PASS1 strand=%d]\n", chosen_strand);
    run_gap_pivots(pat_chosen, chosen_strand,
                   gap_start, gap_end,
                   matchArray, total_smem, matchArray_capacity,
                   rid, min_intv, read_len_bases);

    uint64_t p1_added = *total_smem - before_p1;
    GF_PRINT("    pass1: %lu seed(s) added", (unsigned long)p1_added);

    int covered_after_p1 =
        gap_covered(matchArray, *total_smem, gap_start, gap_end);
    
    #ifdef GAPFILL_ALWAYS_RUN_BOTH_STRANDS
        int run_pass2 = 1;
    #else
        int run_pass2 = !covered_after_p1;
    #endif
    
    if (!run_pass2) {
        GF_PRINT(" —> gap COVERED, skipping PASS2\n");
        return;
    }
    
    GF_PRINT(covered_after_p1
             ? " —> gap COVERED, but PASS2 forced\n"
             : " —> gap still OPEN\n");
    
    GF_PRINT("    [PASS2 strand=%d]\n", opposite_strand);
    
    /* opposite strand, same gap */
    uint64_t before_p2 = *total_smem;
    
    run_gap_pivots(pat_opposite, opposite_strand,
                   gap_start, gap_end,
                   matchArray, total_smem, matchArray_capacity,
                   rid, min_intv, read_len_bases);
    
    GF_PRINT("    pass2: %lu seed(s) added\n",
             (unsigned long)(*total_smem - before_p2));
}


/* ============================================================
   gap_fill_phase()  —  main entry point
   ============================================================
   Processes each gap independently:
     - Left edge:     [0 .. first_seed.start - 1]
     - Internal gaps: [seed[i].end+1 .. seed[i+1].start-1]
     - Right edge:    [last_seed.end+1 .. 149]
     - No seeds:      [0 .. 149] treated as one left-edge gap
   ============================================================ */
void gap_fill_phase(
    const uint8_t  *pat_f4,
    const uint8_t  *pat_rc4,
    SMEM           *matchArray,
    uint64_t       *total_smem,
    uint64_t        matchArray_capacity,
    int             rid,
    int64_t         min_intv,
    int             chosen_strand,
    int             read_len)
{
    const int read_len_bases = read_len;   
    const uint64_t n_seeds   = *total_smem;

    GF_PRINT("=== gap_fill rid=%d  phase1_seeds=%lu  chosen_strand=%d ===\n",
             rid, (unsigned long)n_seeds, chosen_strand);

    /* Identify chosen and opposite strand arrays */
    const uint8_t *pat_chosen   = (chosen_strand == 0) ? pat_f4 : pat_rc4;
    const uint8_t *pat_opposite = (chosen_strand == 0) ? pat_rc4 : pat_f4;
    const int      opp_strand   = 1 - chosen_strand;

    /* ---- NO SEEDS: treat full read as left-edge gap ---- */
    if (n_seeds == 0) {
        process_gap(0, read_len_bases - 1,
                    pat_chosen,   chosen_strand,
                    pat_opposite, opp_strand,
                    matchArray, total_smem, matchArray_capacity,
                    rid, min_intv, read_len_bases);
        return;
    }

    /* ---- LEFT EDGE: [0 .. first_seed.start-1] ---- */
    {
        int gap_start = 0;
        int gap_end   = matchArray[0].start_idx - 1;

        if (gap_end >= 0) {
            process_gap(gap_start, gap_end,
                        pat_chosen,   chosen_strand,
                        pat_opposite, opp_strand,
                        matchArray, total_smem, matchArray_capacity,
                        rid, min_intv, read_len_bases);
        }
    }

    /* ---- INTERNAL GAPS ---- */
    for (uint64_t si = 0; si + 1 < n_seeds; ++si) {
        int gap_start = matchArray[si].end_idx + 1;
        int gap_end   = matchArray[si+1].start_idx - 1;

        process_gap(gap_start, gap_end,
                    pat_chosen,   chosen_strand,
                    pat_opposite, opp_strand,
                    matchArray, total_smem, matchArray_capacity,
                    rid, min_intv, read_len_bases);
    }

    /* ---- RIGHT EDGE: [last_seed.end+1 .. read_end] ---- */
    {
        int gap_start = (int)matchArray[n_seeds-1].end_idx + 1;
        int gap_end   = read_len_bases - 1;

        process_gap(gap_start, gap_end,
                    pat_chosen,   chosen_strand,
                    pat_opposite, opp_strand,
                    matchArray, total_smem, matchArray_capacity,
                    rid, min_intv, read_len_bases);
    }
}