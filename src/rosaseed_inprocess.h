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

#ifndef ROSASEED_INPROCESS_H
#define ROSASEED_INPROCESS_H

#include <stdint.h>
#include "bwamem.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char *index_dir;

    int read_len;

    int min_seed_len_A;      // -ms
    int min_seed_len_BC;     // derived as A + 1 i.e., -ms + 1

    int64_t min_intv;        // -mi
    int phase1_cap;          // -cap
    int gap_threshold;       // -gap

    int num_pivots_B;        // -pb
    int gap_left_step;       // -gs

    uint32_t sa_max_occ;     // -occ
} rosaseed_inprocess_config_t;

int rosaseed_inprocess_init(const rosaseed_inprocess_config_t *cfg);

int64_t rosaseed_inprocess_seed_batch(
    const bseq1_t *seqs,
    int64_t nseq,
    int64_t rid_start,
    rosaseed_hit_t *hitArray,
    int64_t max_hits
);

void rosaseed_inprocess_destroy(void);

#ifdef __cplusplus
}
#endif

#endif