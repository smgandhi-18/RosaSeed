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

1-step
*****************************************************************************************/
#include "gap_fill_phase_compact.h"
#include "seeding_kernel.h"
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

/* ============================================================
   run_single_pivot()
   UNCHANGED: kept for reference / any direct callers.
   The coroutine path below (run_gap_pivots) no longer calls this.
   ============================================================ */
static inline __attribute__((always_inline))
int run_single_pivot(
    const uint8_t  *pat,
    int             pivot_base,
    SMEM           *matchArray,
    uint64_t       *total_smem,
    uint64_t        matchArray_capacity,
    int             rid,
    int64_t         min_intv,
    int             active_strand,
    int             read_len_bases)
{
    const int min_seed_bc = g_min_seed_len_BC;

    if (pivot_base < JT_LEN_NT - 1) return 0;
    if (pivot_base + 1 < min_seed_bc) return 0;
    if (*total_smem >= matchArray_capacity) return 0;

    uint64_t addr = compute_jumpN_from_base4(pat, pivot_base, JT_LEN_NT);
    if (!jump_addr_is_valid(addr)) return 0;   /* pivot window touches N */
    uint64_t jp   = jump_pointers[addr];

    uint64_t l = 0, h = 0, diff = 0;
    extract_jump_bounds(jp, &l, &h, &diff);

    if (diff == 0 || l >= h) return 0;

    int seed_end_base  = pivot_base;
    int seed_len_bases = JT_LEN_NT;
    int read_idx_base  = pivot_base - JT_LEN_NT;

    uint64_t unique_raw_start = 0;

    GF_PRINT("          jt_interval=%-8lu", (unsigned long)diff);

    if (((jp >> 63) & 1ULL) && diff == 1) {
        uint64_t ref_start = jp & 0x1FFFFFFFFULL;
        unique_raw_start = ref_start;

        if (read_idx_base >= 0 && ref_start > 0) {
            int need = min_seed_bc - seed_len_bases;
            if (need < 0) need = 0;

            int reached = 0;
            int mismatch_base = -1;
            uint64_t ref_after = 0;
            int new_read_idx = read_idx_base;
            uint64_t leftmost_ref_start = ref_start;

            int matched = ref_walk_left_b4_until(
                pat,
                read_idx_base,
                ref_start - 1,
                need,
                &reached,
                &mismatch_base,
                &ref_after,
                &new_read_idx,
                &leftmost_ref_start);

            seed_len_bases += matched;
            read_idx_base = new_read_idx;
            unique_raw_start = leftmost_ref_start;
        }
    } else {
        while (read_idx_base >= 0) {
            uint64_t L_before = l;
            uint64_t H_before = h;
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

            if (diff < (uint64_t)min_intv &&
                seed_len_bases >= min_seed_bc) {
                if (diff == 1)
                    unique_raw_start = reconstruct_ref_base_index(l);
                break;
            }

            if (diff == 1) {
                uint64_t ref_start = reconstruct_ref_base_index(l);
                unique_raw_start = ref_start;

                if (read_idx_base >= 0 && ref_start > 0) {
                    int need = min_seed_bc - seed_len_bases;
                    if (need < 0) need = 0;

                    int reached = 0;
                    int mismatch_base = -1;
                    uint64_t ref_after = 0;
                    int new_read_idx = read_idx_base;
                    uint64_t leftmost_ref_start = ref_start;

                    int matched = ref_walk_left_b4_until(
                        pat,
                        read_idx_base,
                        ref_start - 1,
                        need,
                        &reached,
                        &mismatch_base,
                        &ref_after,
                        &new_read_idx,
                        &leftmost_ref_start);

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
    }

    diff = h - l;

    if (seed_len_bases < min_seed_bc) return 0;
    if (l >= h) return 0;
    if (diff >= (uint64_t)min_intv) return 0;

    int start_idx = seed_end_base - seed_len_bases + 1;
    int end_idx   = seed_end_base;

    if (start_idx < 0) start_idx = 0;
    if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

    int fwd_start, fwd_end;

    if (active_strand == 0) {
        fwd_start = start_idx;
        fwd_end   = end_idx;
    } else {
        fwd_start = (read_len_bases - 1) - end_idx;
        fwd_end   = (read_len_bases - 1) - start_idx;
    }

    GF_PRINT(" → EMIT [%d..%d] s=%d\n",
             fwd_start, fwd_end, (int)diff);

    SMEM *s = &matchArray[*total_smem];

    s->rid = rid;
    s->start_idx = fwd_start;
    s->end_idx = fwd_end;

    s->unique_flag = (diff == 1);
    s->unique_ref_pointer = (diff == 1) ? unique_raw_start : 0;

    s->low_ptr_read = l;
    s->low_ptr_rc_read = UINT64_MAX;

    s->smem_score = (diff == 1) ? 1 : (int)diff;
    s->seed_strand = (uint8_t)active_strand;

    s->unique_leftmost_base_off = 0;

    (*total_smem)++;
    return 1;
}

/* ============================================================
   gap_covered()
   ============================================================ */
static inline int gap_covered(
    const SMEM *matchArray,
    uint64_t    n_seeds,
    int         gap_start,
    int         gap_end)
{
    for (uint64_t k = 0; k < n_seeds; ++k) {
        if (matchArray[k].start_idx <= gap_start &&
            matchArray[k].end_idx >= gap_end)
            return 1;
    }
    return 0;
}

/* ============================================================
   run_gap_pivots()
   ----------------
   GapFill 1-step pivot coroutine.

   Same principle as the 2-step coroutines: interleave
   the FM extension one step (1 base) at a time across all pivots
   of this gap, prefetching each pivot's next cp_occ block before
   yielding to the next pivot. Hides BWT RAM latency for every FM
   step, not just the first. All exit/emit conditions identical to
   run_single_pivot() above.
   ============================================================ */

typedef struct {
    int      strand_pivot;
    int      pivot_fwd;
    int      active;   /* 1 = still running in PASS 3 scheduler */
    int      valid;    /* 1 = successfully initialised in PASS 2;
                           explicitly set, never inferred from
                           potentially-garbage field values.     */

    uint64_t l, h, diff;
    uint64_t jp;             /* raw jump pointer, needed for the
                                 (jp>>63)&1ULL unique-at-jump flag  */
    int      seed_end_base;
    int      seed_len_bases;
    int      read_idx_base;
    uint64_t unique_raw_start;

    /* phase: 0=classify  1=FM extension loop  2=emit */
    int      phase;
    int      unique_at_jump;
} GF1PivotSlot;

/* One FM step (1 base) for one GapFill 1-step pivot, with prefetch
   for the next step. Returns 1 = still running, 0 = done.
   Mirrors run_single_pivot()'s inner while(read_idx_base >= 0) loop
   exactly, including the min_intv early-exit's inline unique_raw_start
   reconstruction. All emit logic happens in PASS 4 after this returns 0. */
static inline int gf1_pivot_step(
    GF1PivotSlot   *sl,
    const uint8_t  *pat,
    int64_t         min_intv,
    int             min_seed_bc)
{
    /* Phase 0: classify unique-at-jump vs non-unique (jump already done in setup) */
    if (sl->phase == 0) {
        if (((sl->jp >> 63) & 1ULL) && sl->diff == 1) {
            sl->unique_at_jump = 1;
            sl->phase = 2;
            return 0;
        }
        sl->unique_at_jump = 0;
        sl->phase = 1;

        /* check read_idx_base >= 0, same guard as the original while() */
        if (sl->read_idx_base < 0) {
            sl->phase = 2;
            return 0;
        }

        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5] + 64, 0, 1);
        return 1;
    }

    /* Phase 1: one fm_step1_b4 step */
    if (sl->phase == 1) {
        uint64_t L_before   = sl->l;
        uint64_t H_before   = sl->h;
        int      idx_before = sl->read_idx_base;

        int ok = fm_step1_b4(pat, &sl->read_idx_base, &sl->l, &sl->h);

        if (!ok) {
            /* rollback: same as run_single_pivot */
            sl->l = L_before;
            sl->h = H_before;
            sl->read_idx_base = idx_before;
            sl->phase = 2;
            return 0;
        }

        sl->seed_len_bases++;
        sl->diff = sl->h - sl->l;

        /* early-exit on specificity threshold */
        if (sl->diff < (uint64_t)min_intv && sl->seed_len_bases >= min_seed_bc) {
            if (sl->diff == 1)
                sl->unique_raw_start = reconstruct_ref_base_index(sl->l);
            sl->phase = 2;
            return 0;
        }

        /* became unique mid-extension */
        if (sl->diff == 1) {
            uint64_t ref_start = reconstruct_ref_base_index(sl->l);
            sl->unique_raw_start = ref_start;

            if (sl->read_idx_base >= 0 && ref_start > 0) {
                int need = min_seed_bc - sl->seed_len_bases;
                if (need < 0) need = 0;

                int reached = 0;
                int mismatch_base = -1;
                uint64_t ref_after = 0;
                int new_read_idx = sl->read_idx_base;
                uint64_t leftmost_ref_start = ref_start;

                int matched = ref_walk_left_b4_until(
                    pat,
                    sl->read_idx_base,
                    ref_start - 1,
                    need,
                    &reached,
                    &mismatch_base,
                    &ref_after,
                    &new_read_idx,
                    &leftmost_ref_start);

                sl->seed_len_bases  += matched;
                sl->read_idx_base    = new_read_idx;
                sl->unique_raw_start = leftmost_ref_start;
            }

            sl->phase = 2;
            return 0;
        }

        /* read exhausted */
        if (sl->read_idx_base < 0) {
            sl->phase = 2;
            return 0;
        }

        /* more FM steps */
        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->l >> 5] + 64, 0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5],      0, 1);
        __builtin_prefetch((const char *)&cp_occ[sl->h >> 5] + 64, 0, 1);
        return 1;
    }

    /* Phase 2: done */
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

    /* ----------------------------------------------------------------
       COROUTINE INTERLEAVER for 1-step GapFill pivots.
       ---------------------------------------------------------------- */

    GF1PivotSlot slots[max_pivots];
    int n_active = 0;

    for (int i = 0; i < n_pivots; ++i) {
        slots[i].active           = 0;
        slots[i].valid            = 0;
        slots[i].l                = 0;
        slots[i].h                = 0;
        slots[i].diff             = 0;
        slots[i].jp               = 0;
        slots[i].seed_len_bases   = 0;
        slots[i].read_idx_base    = 0;
        slots[i].unique_raw_start = 0;
        slots[i].phase            = 0;
        slots[i].unique_at_jump   = 0;
    }

    /* PASS 1: prefetch jump_pointers[] entries for all pivots */
    for (int i = 0; i < n_pivots; ++i) {
        int strand_pivot = (active_strand == 0)
                          ? pivots_fwd[i]
                          : (read_len_bases - 1) - pivots_fwd[i];
        if (strand_pivot < min_pivot) continue;   /* leave inactive */
        uint64_t addr = compute_jumpN_from_base4(pat, strand_pivot, JT_LEN_NT);
        if (jump_addr_is_valid(addr))
            __builtin_prefetch(&jump_pointers[addr], 0, 1);
        slots[i].active = 1;   /* candidate: still needs PASS 2 validation */
    }

    /* PASS 2: read warm jump entries, init slots, prefetch first cp_occ */
    for (int i = 0; i < n_pivots; ++i) {
        if (!slots[i].active) continue;

        int strand_pivot = (active_strand == 0)
                          ? pivots_fwd[i]
                          : (read_len_bases - 1) - pivots_fwd[i];

        if (strand_pivot < JT_LEN_NT - 1) { slots[i].active = 0; continue; }

        uint64_t addr = compute_jumpN_from_base4(pat, strand_pivot, JT_LEN_NT);
        if (!jump_addr_is_valid(addr)) { slots[i].active = 0; continue; }
        uint64_t jp   = jump_pointers[addr];

        uint64_t l, h, diff;
        extract_jump_bounds(jp, &l, &h, &diff);

        /* same guard as run_single_pivot: diff==0 also bails */
        if (diff == 0 || l >= h) { slots[i].active = 0; continue; }   /* leave valid=0 */

        /* initialise slot state */
        slots[i].valid            = 1;
        slots[i].strand_pivot     = strand_pivot;
        slots[i].pivot_fwd        = pivots_fwd[i];
        slots[i].l                = l;
        slots[i].h                = h;
        slots[i].diff             = diff;
        slots[i].jp               = jp;
        slots[i].seed_end_base    = strand_pivot;
        slots[i].seed_len_bases   = JT_LEN_NT;
        slots[i].read_idx_base    = strand_pivot - JT_LEN_NT;
        slots[i].unique_raw_start = 0;
        slots[i].phase            = 0;
        slots[i].unique_at_jump   = 0;

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

    /* PASS 3: round-robin FM interleaver */
    while (n_active > 0) {
        for (int i = 0; i < n_pivots; ++i) {
            if (!slots[i].active) continue;

            int still_running = gf1_pivot_step(
                &slots[i], pat, min_intv, min_seed_bc);

            if (!still_running) {
                slots[i].active = 0;
                n_active--;
            }
        }
    }

    /* PASS 4: emit seeds from completed slots: same logic as
       run_single_pivot()'s tail end */
    int emitted = 0;
    for (int i = 0; i < n_pivots; ++i) {
        if (*total_smem >= matchArray_capacity - 1) break;

        GF1PivotSlot *sl = &slots[i];

        if (!sl->valid) continue;

        /* handle unique-at-jump: ref_walk needed, done here like
           run_single_pivot's ((jp>>63)&1ULL) branch */
        if (sl->unique_at_jump) {
            uint64_t ref_start = sl->jp & 0x1FFFFFFFFULL;
            sl->unique_raw_start = ref_start;

            if (sl->read_idx_base >= 0 && ref_start > 0) {
                int need = min_seed_bc - sl->seed_len_bases;
                if (need < 0) need = 0;

                int reached = 0;
                int mismatch_base = -1;
                uint64_t ref_after = 0;
                int new_read_idx = sl->read_idx_base;
                uint64_t leftmost_ref_start = ref_start;

                int matched = ref_walk_left_b4_until(
                    pat,
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
        } else {
            /* non-unique path: re-derive diff and the trailing
               unique_raw_start reconstruction, same as
               run_single_pivot's post-loop block */
            sl->diff = sl->h - sl->l;
            if (sl->diff == 1 && sl->unique_raw_start == 0)
                sl->unique_raw_start = reconstruct_ref_base_index(sl->l);
        }

        sl->diff = sl->h - sl->l;

        /* same emit checks as run_single_pivot */
        if (sl->seed_len_bases < min_seed_bc) continue;
        if (sl->l >= sl->h)                   continue;
        if (sl->diff >= (uint64_t)min_intv)   continue;

        int start_idx = sl->seed_end_base - sl->seed_len_bases + 1;
        int end_idx   = sl->seed_end_base;

        if (start_idx < 0) start_idx = 0;
        if (end_idx >= read_len_bases) end_idx = read_len_bases - 1;

        int fwd_start, fwd_end;
        if (active_strand == 0) {
            fwd_start = start_idx;
            fwd_end   = end_idx;
        } else {
            fwd_start = (read_len_bases - 1) - end_idx;
            fwd_end   = (read_len_bases - 1) - start_idx;
        }

        SMEM *s = &matchArray[*total_smem];

        s->rid       = rid;
        s->start_idx = fwd_start;
        s->end_idx   = fwd_end;

        s->unique_flag        = (sl->diff == 1);
        s->unique_ref_pointer = (sl->diff == 1) ? sl->unique_raw_start : 0;

        s->low_ptr_read    = sl->l;
        s->low_ptr_rc_read = UINT64_MAX;

        s->smem_score  = (sl->diff == 1) ? 1 : (int)sl->diff;
        s->seed_strand = (uint8_t)active_strand;

        s->unique_leftmost_base_off = 0;

        (*total_smem)++;
        emitted++;

#ifdef GAPFILL_EARLY_EXIT
        if (emitted > 0 &&
            gap_covered(matchArray, *total_smem, gap_start, gap_end)) {
            break;
        }
#endif
    }

    return emitted;
}

/* ============================================================
   process_gap()
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

    /* ---- Pass 1: chosen strand ---- */
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

#ifndef GAPFILL_ALWAYS_RUN_BOTH_STRANDS
    if (covered_after_p1) {
        GF_PRINT(", gap covered after PASS1\n");
        return;
    }
#endif

    GF_PRINT(", gap still OPEN or PASS2 forced\n");
    GF_PRINT("    [PASS2 strand=%d]\n", opposite_strand);

    uint64_t before_p2 = *total_smem;

    run_gap_pivots(pat_opposite, opposite_strand,
                    gap_start, gap_end,
                    matchArray, total_smem, matchArray_capacity,
                    rid, min_intv, read_len_bases);

    GF_PRINT("    pass2: %lu seed(s) added\n",
            (unsigned long)(*total_smem - before_p2));
}


/* ============================================================
   gap_fill_phase() 
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
    const int read_len_bases = read_len;   /* 150 */
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

    /* ---- RIGHT EDGE: [last_seed.end+1 .. 149] ---- */
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