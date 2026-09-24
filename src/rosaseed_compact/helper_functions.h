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
#ifndef _HELPER_FUNCTIONS_H
#define _HELPER_FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <x86intrin.h>
#include <unistd.h>
#include <xmmintrin.h>
#include <assert.h>

#include "macros.h"
#include "file_dec.h"

int compare_smem(const void *a, const void *b);
void sort_smems_for_read(SMEM *matchArray, int64_t num_smem);

/*
 * Reconstruct SA[row].
 *
 * 1-step convention:
 *   SA[row] is the raw start coordinate in:
 *      forward_reference + reverse_complement_reference + terminator
 *
 * Do NOT subtract 1 here.
 */
uint64_t reconstruct_ref_base_index(uint64_t sa_rank_low);

void extract_jump_bounds(uint64_t jump_entry,
                         uint64_t *l,
                         uint64_t *h,
                         uint64_t *diff);

typedef struct {
    uint64_t L;      // original forward reference length
    uint64_t twoL;   // 2*L, valid biological coordinate space [0, twoL)
    uint64_t N;      // BWT length including terminator
} OneStepMapInfo;

void init_one_step_mapinfo(OneStepMapInfo *info, uint64_t L);

uint32_t materialize_seed_rbegs_base4(
    const OneStepMapInfo *info,
    const SMEM *s,
    uint32_t max_occ,
    int64_t *rbeg_out);

uint32_t materialize_read_rbegs_base4(
    const OneStepMapInfo *info,
    const SMEM *smems,
    uint32_t n_sm,
    uint32_t max_occ,
    int64_t *rbeg_out,
    uint32_t *seed_offsets,
    uint32_t *seed_counts);

uint64_t sal_touch_nonunique_seeds(
    const OneStepMapInfo *info,
    const SMEM *smems,
    uint64_t n_sm);

#endif