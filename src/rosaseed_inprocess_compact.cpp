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

/* RosaSeed-Compact implementation of the variant-independent interface in
   rosaseed_inprocess.h (the 2-step implementation is rosaseed_inprocess.cpp;
   the Makefile links exactly one of them). */
#include "rosaseed_inprocess.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rosaseed_core_bridge_compact.h"

int rosaseed_inprocess_init(const rosaseed_inprocess_config_t *cfg)
{
    rosaseed_core_config_t core_cfg;
    memset(&core_cfg, 0, sizeof(core_cfg));

    core_cfg.index_dir       = cfg->index_dir;
    core_cfg.read_len        = cfg->read_len;
    core_cfg.min_seed_len_A  = cfg->min_seed_len_A;
    core_cfg.min_seed_len_BC = cfg->min_seed_len_BC;
    core_cfg.min_intv        = cfg->min_intv;
    core_cfg.phase1_cap      = cfg->phase1_cap;
    core_cfg.gap_threshold   = cfg->gap_threshold;
    core_cfg.num_pivots_B    = cfg->num_pivots_B;
    core_cfg.gap_left_step   = cfg->gap_left_step;
    core_cfg.sa_max_occ      = cfg->sa_max_occ;

    return rosaseed_core_init_compact(&core_cfg);
}

int64_t rosaseed_inprocess_seed_batch(
    const bseq1_t *seqs,
    int64_t nseq,
    int64_t rid_start,
    rosaseed_hit_t *hitArray,
    int64_t max_hits
)
{
    const int64_t STACK_LIMIT = 2048;
    const char **seq_ptrs  = NULL;
    int        *rlens      = NULL;
    int        *rids_arr   = NULL;
    int64_t    *offsets    = NULL;

    const char *seq_ptrs_stk [STACK_LIMIT];
    int         rlens_stk    [STACK_LIMIT];
    int         rids_arr_stk [STACK_LIMIT];
    int64_t     offsets_stk  [STACK_LIMIT];

    if (nseq <= STACK_LIMIT) {
        seq_ptrs = seq_ptrs_stk;
        rlens    = rlens_stk;
        rids_arr = rids_arr_stk;
        offsets  = offsets_stk;
    } else {
        seq_ptrs = (const char **)malloc(nseq * sizeof(const char *));
        rlens    = (int *)        malloc(nseq * sizeof(int));
        rids_arr = (int *)        malloc(nseq * sizeof(int));
        offsets  = (int64_t *)    malloc(nseq * sizeof(int64_t));
        if (!seq_ptrs || !rlens || !rids_arr || !offsets) {
            fprintf(stderr, "[ERROR] RosaSeed-Compact batch alloc failed\n");
            exit(EXIT_FAILURE);
        }
    }

    for (int64_t i = 0; i < nseq; i++) {
        seq_ptrs [i] = seqs[i].seq;
        rlens    [i] = seqs[i].l_seq;
        rids_arr [i] = (int)(rid_start + i);
        offsets  [i] = 0;
    }

    int64_t total_hits = rosaseed_core_seed_batch_interleaved_compact(
        seq_ptrs, rlens, rids_arr, (int)nseq,
        (rosaseed_core_hit_t *)hitArray, max_hits, offsets);

    if (total_hits < 0) exit(EXIT_FAILURE);

    for (int64_t i = 0; i < nseq; i++) {
        int64_t hi_start = offsets[i];
        int64_t hi_end   = (i + 1 < nseq) ? offsets[i + 1] : total_hits;
        int     local_rid = (int)i;
        for (int64_t h = hi_start; h < hi_end; ++h)
            hitArray[h].rid = local_rid;
    }

    if (nseq > STACK_LIMIT) {
        free(seq_ptrs); free(rlens); free(rids_arr); free(offsets);
    }

    return total_hits;
}

void rosaseed_inprocess_destroy(void)
{
    rosaseed_core_destroy();
}
