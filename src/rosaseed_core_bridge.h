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

#ifndef ROSASEED_CORE_BRIDGE_H
#define ROSASEED_CORE_BRIDGE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char *index_dir;
    int read_len;
    int min_seed_len_A;
    int min_seed_len_BC;
    int64_t min_intv;
    int phase1_cap;
    int gap_threshold;
    int num_pivots_B;
    int gap_left_step;
    uint32_t sa_max_occ;
} rosaseed_core_config_t;

int rosaseed_core_init(const rosaseed_core_config_t *cfg);
void rosaseed_core_destroy(void);
void rosaseed_core_print_timing(double proc_freq);

typedef struct {
    int32_t rid;
    int32_t qbeg;
    int32_t len;
    int64_t rbeg;
    int32_t occ;
} rosaseed_core_hit_t;

int64_t rosaseed_core_seed_one_read(
    const char *seq,
    int read_len,
    int rid,
    rosaseed_core_hit_t *hits,
    int64_t max_hits
);

int64_t rosaseed_core_seed_batch_interleaved(
    const char **seqs,
    const int   *read_lens,
    const int   *rids,
    int          nreads,
    rosaseed_core_hit_t *hits,
    int64_t      max_hits,
    int64_t     *hit_offsets   /* can be NULL if you don't need per-read offsets */
);

#ifdef __cplusplus
}
#endif

#endif