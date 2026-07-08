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

#ifndef GAP_FILL_PHASE_H
#define GAP_FILL_PHASE_H

#include <stdint.h>
#include "macros.h"
#include "file_dec.h"
#include "bwa.h"
#include "helper_functions.h"

/* -----------------------------------------------------------------------
   Configuration
   ----------------------------------------------------------------------- */
#define GAP_THRESHOLD       5    // compile-time default only — use g_gap_threshold at runtime
#define GAP_LEFT_STEP       4    // compile-time default only — use g_gap_left_step at runtime

/* MAX_PIVOTS_PER_GAP is kept for reference; run_gap_pivots sizes its VLA at runtime
   using (read_len_bases / g_gap_left_step) + 4 to handle any step value correctly. */
#define MAX_PIVOTS_PER_GAP  ((READ_LEN / GAP_LEFT_STEP) + 4)
#define MAX_ADAPTIVE_PIVOTS 64   // unused in .c but harmless to keep
   
// Runtime-configurable gap parameters.
extern int g_gap_threshold;   // min gap size to fill   (default 5)
extern int g_gap_left_step;   // pivot stepping size    (default 4)

/*declarations */
void     fm_index_mapping_backward_search(uint64_t *low, uint64_t *high,
                                          base16_t nuc_in_read);
uint64_t reconstruct_ref_base_index(uint64_t sa_rank_low);
void     extract_jump_bounds(uint64_t jump_entry,
                             uint64_t *l, uint64_t *h, uint64_t *diff);

/* -----------------------------------------------------------------------
   API
   ----------------------------------------------------------------------- */
void gap_fill_phase(
    const uint8_t  *pat_f4,         /* always the FORWARD base-4 read   */
    const uint8_t  *pat_rc4,        /* always the RC base-4 read        */
    SMEM           *matchArray,
    uint64_t       *total_smem,
    uint64_t        matchArray_capacity,
    int             rid,
    int64_t         min_intv,
    int             chosen_strand,
    int             read_len); /* 0=fwd 1=RC from Phase1           */

#endif /* GAP_FILL_PHASE_H */
