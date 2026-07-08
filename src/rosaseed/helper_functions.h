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

#ifndef _HELPER_FUNCTIONS_H
#define _HELPER_FUNCTIONS_H

// #define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "macros.h"
#include "file_dec.h"
#include <xmmintrin.h>
#include <assert.h>

//int compare_smem(const void *a, const void *b);
void sort_smems_for_read(SMEM *matchArray, int64_t num_smem);
uint64_t reconstruct_ref_base_index(uint64_t sa_rank_low);
void extract_jump_bounds(uint64_t jump_entry, uint64_t *l, uint64_t *h, uint64_t *diff);

typedef struct {
   uint64_t L;          // original forward reference length in bases
   uint64_t N_even;     // FE length in base16 symbols
   uint64_t N_odd;      // FO length in base16 symbols
   uint64_t offset_fo;  // start of FO
   uint64_t offset_rce; // start of RCE
   uint64_t offset_rco; // start of RCO
} TwoStepMapInfo;
void init_two_step_mapinfo(TwoStepMapInfo *info, uint64_t L);

uint64_t sal_touch_nonunique_seeds(
   const TwoStepMapInfo *info,
   const SMEM *smems,
   uint64_t n_sm
);

uint32_t materialize_seed_rbegs_base4(
    const TwoStepMapInfo *info,
    const SMEM *s,
    uint32_t max_occ,
    int64_t *rbeg_out);

uint32_t materialize_read_rbegs_base4(
    const TwoStepMapInfo *info,
    const SMEM *smems,
    uint32_t n_sm,
    uint32_t max_occ,
    int64_t *rbeg_out,
    uint32_t *seed_offsets,
    uint32_t *seed_counts);

void rosaseed_read_index_metadata(const char *index_dir);

#endif
