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

uint64_t increment_sa_compute = 0;

extern int64_t sentinel_index;

int rosaseed_compare_smem(const void *a, const void *b)
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

#if SA_COMPRESSION_FACTOR_POWER > 0

static inline int bwt_symbol_at(uint64_t pos)
{
    uint64_t blk = pos >> 5;
    uint32_t off = pos & 31;
    uint32_t bit = 1u << (31 - off);

    for (int s = 0; s < ALPHABET_SIZE; ++s) {
        if (cp_occ[blk].one_hot[s] & bit) return s;
    }

    return -1;  // sentinel
}

static inline uint64_t LF_step(uint64_t pos)
{
    if (pos == (uint64_t)sentinel_index) return 0;

    int sym = bwt_symbol_at(pos);
    if (sym < 0) return 0;

    uint64_t occ_excl;
    GET_OCC32(pos, sym, occ_excl);

    return c_vec[sym] + occ_excl;
}

static inline uint64_t get_sa_entry_compressed(uint64_t pos)
{
    uint64_t steps = 0;

    while ((pos & SA_INDEX_AND_SAMPLE) != 0) {
        pos = LF_step(pos);
        steps++;
    }

    uint64_t idx = pos >> SA_COMPRESSION_FACTOR_POWER;
    uint64_t base = ((uint64_t)sa_ms_byte[idx] << 32) |
                    (uint64_t)sa_ls_word[idx];

    return base + steps;
}

#endif

uint64_t reconstruct_ref_base_index(uint64_t sa_rank_low)
{
#if SA_COMPRESSION_FACTOR_POWER == 0
    uint64_t sa = ((uint64_t)sa_ms_byte[sa_rank_low] << 32) |
                  (uint64_t)sa_ls_word[sa_rank_low];
    return sa;
#else
    return get_sa_entry_compressed(sa_rank_low);
#endif
}

void extract_jump_bounds(uint64_t jump_entry,
                         uint64_t *l,
                         uint64_t *h,
                         uint64_t *diff)
{
    if (jump_entry == 0) {
        *l = 0;
        *h = 0;
        *diff = 0;
        return;
    }

    if ((jump_entry >> 63) & 1ULL) {
        *l = jump_entry & 0x1FFFFFFFFULL;
        *h = *l + 1;
        *diff = 1;
        return;
    }

    *l = jump_entry & 0x3FFFFFFFFULL;
    *diff = jump_entry >> 34;
    *h = *l + *diff;
}

void init_one_step_mapinfo(OneStepMapInfo *info, uint64_t L)
{
    info->L = L;
    info->twoL = 2ULL * L;
    info->N = BWT_SIZE_REFERENCE_SIZE;
}

static inline uint32_t smem_seed_len(const SMEM *s)
{
    return (uint32_t)(s->end_idx - s->start_idx + 1);
}

static inline uint64_t finalize_1step_rbeg_from_raw(
    const OneStepMapInfo *info,
    uint64_t raw,
    const SMEM *s)
{
    const uint32_t slen = smem_seed_len(s);

    if (slen == 0) return UINT64_MAX;

    /*
     * Terminator position and boundary-crossing intervals are invalid
     * as biological seed coordinates.
     */
    if (raw >= info->twoL) return UINT64_MAX;
    if (raw + (uint64_t)slen > info->twoL) return UINT64_MAX;

    if (s->seed_strand == 1) {
        return info->twoL - raw - (uint64_t)slen;
    }

    return raw;
}

static inline uint64_t unique_seed_to_bwa_rbeg(
    const OneStepMapInfo *info,
    const SMEM *s)
{
    if (!s->unique_flag) return UINT64_MAX;

    /*
     * 1-step unique_ref_pointer is already raw SA start coordinate.
     */
    return finalize_1step_rbeg_from_raw(info, s->unique_ref_pointer, s);
}

uint32_t materialize_seed_rbegs_base4(
    const OneStepMapInfo *info,
    const SMEM *s,
    uint32_t max_occ,
    int64_t *rbeg_out)
{
    if (!info || !s || !rbeg_out || max_occ == 0) return 0;

    if (s->unique_flag) {
        uint64_t rbeg = unique_seed_to_bwa_rbeg(info, s);
        if (rbeg == UINT64_MAX) return 0;

        rbeg_out[0] = (int64_t)rbeg;
        return 1;
    }

    uint32_t width = (uint32_t)s->smem_score;
    if (width == 0) return 0;

    uint64_t start_rank = (uint64_t)s->low_ptr_read;

    uint32_t step = (width > max_occ) ? (width / max_occ) : 1;
    if (step == 0) step = 1;

    uint32_t out_count = 0;

    for (uint32_t k = 0; k < width && out_count < max_occ; k += step) {
        uint64_t rank = start_rank + k;
        uint64_t raw = reconstruct_ref_base_index(rank);

        uint64_t rbeg = finalize_1step_rbeg_from_raw(info, raw, s);
        if (rbeg == UINT64_MAX) continue;

        rbeg_out[out_count++] = (int64_t)rbeg;
    }

    return out_count;
}

uint32_t materialize_read_rbegs_base4(
    const OneStepMapInfo *info,
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

#ifndef LIKELY
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#endif

uint64_t sal_touch_nonunique_seeds(
    const OneStepMapInfo *info,
    const SMEM *smems,
    uint64_t n_sm)
{
    uint64_t checksum = 0;

    for (uint64_t i = 0; i < n_sm; ++i) {
        const uint32_t width = (uint32_t)smems[i].smem_score;

        if (LIKELY(width <= 1)) {
            uint64_t rbeg = unique_seed_to_bwa_rbeg(info, &smems[i]);
            if (rbeg != UINT64_MAX) checksum += rbeg;
            continue;
        }

        uint64_t r = (uint64_t)smems[i].low_ptr_read;

        for (uint32_t k = 0; k < width; ++k, ++r) {
            __builtin_prefetch(&sa_ls_word[r + 32], 0, 1);
            __builtin_prefetch(&sa_ms_byte[r + 128], 0, 1);

            uint64_t raw = reconstruct_ref_base_index(r);
            uint64_t rbeg = finalize_1step_rbeg_from_raw(info, raw, &smems[i]);

            if (rbeg != UINT64_MAX)
                checksum += rbeg;
        }
    }

    return checksum;
}