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

#include "seeding_kernel_mt.h"
#include "file_dec.h"
#include "load_data.h"
#include "macros.h"
#include "profiling.h"
#include "read_init.h"
#include "helper_functions.h"

uint64_t increment_sa_compute = 0;

#define TWO_STEP_TOTAL_LEN (2ULL * rosaseed_L)

/* SA / BWT row space: [0, N).  SA values wrap modulo N.
   NOTE: this is NOT 2*rosaseed_L: for 2-step, 2L == N+1, and for
   1-step, 2L == N-1. */
#define BWT_ROW_SPACE    (BWT_SIZE_REFERENCE_SIZE)

/* BWA global coordinate space: [0, 2L), forward then reverse-complement. */
#define BWA_COORD_SPACE  (2ULL * rosaseed_L)

static int rosaseed_compare_smem(const void *a, const void *b)
{
    const SMEM *pa = (const SMEM *)a;
    const SMEM *pb = (const SMEM *)b;

    if (pa->start_idx < pb->start_idx) return -1;
    if (pa->start_idx > pb->start_idx) return  1;
    if (pa->end_idx   < pb->end_idx)   return -1;
    if (pa->end_idx   > pb->end_idx)   return  1;

    return 0;
}

void sort_smems_for_read(SMEM *matchArray, int64_t num_smem)
{
    qsort(matchArray, (size_t)num_smem, sizeof(SMEM), rosaseed_compare_smem);
}

/* --------------------------------------------------------------------------
   Compressed SA helpers
   -------------------------------------------------------------------------- */

#if SA_COMPRESSION_FACTOR_POWER > 0

static inline __attribute__((always_inline))
int bwt_symbol_at(uint64_t pos)
{
    uint64_t blk = pos >> 5;
    uint32_t off = pos & 31;
    uint32_t bit = 1u << (31 - off);

    for (int s = 0; s < ALPHABET_SIZE; ++s) {
        if (cp_occ[blk].one_hot[s] & bit)
            return s;
    }

    return -1;  /* sentinel */
}

static inline __attribute__((always_inline))
uint64_t LF_step_RS(uint64_t pos)
{
    if (pos == (uint64_t)sentinel_index)
        return 0;

    int sym = bwt_symbol_at(pos);
    if (sym < 0)
        return 0;

    uint64_t occ_excl;
    GET_OCC32(pos, sym, occ_excl);

    return c_vec[sym] + occ_excl;
}

static inline __attribute__((always_inline))
uint64_t get_sa_entry_compressed_RS(uint64_t pos)
{
    uint64_t steps = 0;

    while ((pos & SA_INDEX_AND_SAMPLE) != 0) {
        pos = LF_step_RS(pos);
        steps++;
    }

    uint64_t idx  = pos >> SA_COMPRESSION_FACTOR_POWER;
    uint64_t base = ((uint64_t)sa_ms_byte[idx] << 32) |
                     (uint64_t)sa_ls_word[idx];

    uint64_t val = base + steps;

    if (val >= BWT_ROW_SPACE) val -= BWT_ROW_SPACE;      /* was TWO_STEP_TOTAL_LEN */

    return val;
}

#endif

/* --------------------------------------------------------------------------
   SA reconstruction
   -------------------------------------------------------------------------- */

static inline __attribute__((always_inline))
uint64_t sa_rank_to_twostep_pos_fast(uint64_t rank)
{
#if SA_COMPRESSION_FACTOR_POWER == 0
    const uint64_t lo  = (uint64_t)sa_ls_word[rank];
    const uint64_t msb = (uint64_t)(sa_ms_byte[rank] & 1u);
    return (msb << 32) | lo;
#else
    return get_sa_entry_compressed_RS(rank);
#endif
}

uint64_t reconstruct_ref_base_index(uint64_t sa_rank_low)
{
    uint64_t sa = sa_rank_to_twostep_pos_fast(sa_rank_low);
    if (sa == 0)
        return BWT_SIZE_REFERENCE_SIZE - 1;
    return sa - 1;
}

void extract_jump_bounds(uint64_t jump_entry,
                         uint64_t *l,
                         uint64_t *h,
                         uint64_t *diff)
{
    *l = jump_entry & 0x1FFFFFFFFULL;
    *diff = ((jump_entry >> 63) & 1)
          ? 1
          : ((jump_entry >> 34) & 0x1FFFFFFF);
    *h = *l + *diff;
}

/* --------------------------------------------------------------------------
   Two-step coordinate map
   -------------------------------------------------------------------------- */
void init_two_step_mapinfo(TwoStepMapInfo *info, uint64_t L)
{
    info->L = L;

    if (L & 1) {
        info->N_even = L / 2;
        info->N_odd  = L / 2;
    } else {
        info->N_even = L / 2;
        info->N_odd  = (L / 2) - 1;
    }

    info->offset_fo  = info->N_even;
    info->offset_rce = info->N_even + info->N_odd;
    info->offset_rco = info->offset_rce + info->N_even;
}

static inline __attribute__((always_inline))
uint32_t smem_seed_len(const SMEM *s)
{
    return (uint32_t)(s->end_idx - s->start_idx + 1);
}

static inline __attribute__((always_inline))
uint64_t translate_two_step_to_bwa_global_fast(
    const TwoStepMapInfo *info,
    uint64_t p)
{
    if (p < info->offset_fo) {
        return p << 1;                                      /* FE */
    } else if (p < info->offset_rce) {
        return ((p - info->offset_fo) << 1) + 1;             /* FO */
    } else if (p < info->offset_rco) {
        return info->L + ((p - info->offset_rce) << 1);      /* RCE */
    } else {
        return info->L + ((p - info->offset_rco) << 1) + 1;  /* RCO */
    }
}

static inline __attribute__((always_inline))
uint64_t finalize_rbeg_from_raw_fast(
    const uint64_t twoL,
    const uint32_t slen,
    const int is_rc,
    uint64_t raw)
{
    if (slen == 0) return UINT64_MAX;
    if (raw >= twoL) return UINT64_MAX;
    if (raw + (uint64_t)slen > twoL) return UINT64_MAX;

    return is_rc ? (twoL - raw - (uint64_t)slen) : raw;
}

static inline __attribute__((always_inline))
uint64_t unique_seed_to_bwa_rbeg_fast(
    const TwoStepMapInfo *info,
    const SMEM *s)
{
    if (!s->unique_flag)
        return UINT64_MAX;

    const uint64_t twoL = 2ULL * info->L;
    const uint32_t slen = smem_seed_len(s);
    const int is_rc = (s->seed_strand == 1);

    uint64_t raw = translate_two_step_to_bwa_global_fast(
        info,
        s->unique_leftmost_twostep_pos);

    raw += (uint64_t)s->unique_leftmost_base_off;

    return finalize_rbeg_from_raw_fast(twoL, slen, is_rc, raw);
}

uint32_t materialize_seed_rbegs_base4(
    const TwoStepMapInfo *info,
    const SMEM *s,
    uint32_t max_occ,
    int64_t *rbeg_out)
{
    if (!info || !s || !rbeg_out || max_occ == 0)
        return 0;

    if (s->unique_flag) {
        uint64_t rbeg = unique_seed_to_bwa_rbeg_fast(info, s);
        if (rbeg == UINT64_MAX)
            return 0;

        rbeg_out[0] = (int64_t)rbeg;
        return 1;
    }

    const uint32_t width = (uint32_t)s->smem_score;
    if (width == 0)
        return 0;

    const uint64_t start_rank = (uint64_t)s->low_ptr_read;

    // fast old sampling
    const uint32_t step = (width > max_occ) ? (width / max_occ) : 1u;

    const uint64_t twoL = 2ULL * info->L;
    const uint32_t slen = smem_seed_len(s);
    const int is_rc = (s->seed_strand == 1);

    uint32_t out_count = 0;

    for (uint32_t k = 0; k < width && out_count < max_occ; k += step) {
        const uint64_t rank = start_rank + (uint64_t)k;

#if SA_COMPRESSION_FACTOR_POWER == 0
        __builtin_prefetch(&sa_ls_word[rank + 32], 0, 1);
        __builtin_prefetch(&sa_ms_byte[rank + 64], 0, 1);
#else
        __builtin_prefetch(&sa_ls_word[(rank >> SA_COMPRESSION_FACTOR_POWER) + 32], 0, 1);
        __builtin_prefetch(&sa_ms_byte[(rank >> SA_COMPRESSION_FACTOR_POWER) + 64], 0, 1);
#endif

        const uint64_t p = sa_rank_to_twostep_pos_fast(rank);

        const uint64_t raw =
            translate_two_step_to_bwa_global_fast(info, p);

        const uint64_t rbeg =
            finalize_rbeg_from_raw_fast(twoL, slen, is_rc, raw);

        if (rbeg == UINT64_MAX)
            continue;

        rbeg_out[out_count++] = (int64_t)rbeg;
    }

    return out_count;
}

uint32_t materialize_read_rbegs_base4(
    const TwoStepMapInfo *info,
    const SMEM *smems,
    uint32_t n_sm,
    uint32_t max_occ,
    int64_t *rbeg_out,
    uint32_t *seed_offsets,
    uint32_t *seed_counts)
{
    if (!info || !smems || !rbeg_out || !seed_offsets || !seed_counts)
        return 0;

    uint32_t total_hits = 0;

    for (uint32_t i = 0; i < n_sm; ++i) {
        seed_offsets[i] = total_hits;

        uint32_t n_hits = materialize_seed_rbegs_base4(
            info,
            &smems[i],
            max_occ,
            &rbeg_out[total_hits]);

        seed_counts[i] = n_hits;
        total_hits += n_hits;
    }

    return total_hits;
}

/* --------------------------------------------------------------------------
   Optional SAL touch benchmark helper
   -------------------------------------------------------------------------- */
#ifndef LIKELY
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#endif

uint64_t sal_touch_nonunique_seeds(
    const TwoStepMapInfo *info,
    const SMEM *smems,
    uint64_t n_sm)
{
    uint64_t checksum = 0;

    for (uint64_t i = 0; i < n_sm; ++i) {
        const uint32_t width = (uint32_t)smems[i].smem_score;

        if (LIKELY(width <= 1)) {
            uint64_t rbeg = unique_seed_to_bwa_rbeg_fast(info, &smems[i]);
            if (rbeg != UINT64_MAX)
                checksum += rbeg;
            continue;
        }

        uint64_t r = (uint64_t)smems[i].low_ptr_read;

        for (uint32_t k = 0; k < width; ++k, ++r) {
            __builtin_prefetch(&sa_ls_word[r + 32], 0, 1);
            __builtin_prefetch(&sa_ms_byte[r + 64], 0, 1);

            uint64_t p = sa_rank_to_twostep_pos_fast(r);
            uint64_t raw = translate_two_step_to_bwa_global_fast(info, p);

            checksum += raw;
        }
    }

    return checksum;
}

