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

/*************************************************************************************
 * phaseA_slot_state_compact.h
 *
 * Shared header for RosaSeed-Compact coroutine Phase I.
 * Included by:
 *   - src/rosaseed_compact/rosaseed_compact_phaseA.c  (defines phaseI1_batch_interleaved)
 *   - src/rosaseed_core_bridge_compact.cpp  (calls it)
 *************************************************************************************/
#ifndef PHASEA_SLOT_STATE_COMPACT_H
#define PHASEA_SLOT_STATE_COMPACT_H

#include <stdint.h>
#include "macros.h"   /* for SMEM typedef */

#ifndef RS_BATCH_MAX_READ_LEN
#  define RS_BATCH_MAX_READ_LEN 251
#endif

#ifndef RS_BATCH
#  define RS_BATCH 32
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * phaseI1_batch_interleaved()
 *
 * RosaSeed-Compact coroutine Phase I. Same scheduler design as the
 * 2-step version but uses fm_step1_b4 (1 base/step) and ref_walk_left_b4.
 * All logic, exit conditions and emit conditions identical to
 * phaseI_routine() in src/rosaseed_compact/rosaseed_compact_phaseA.c.
 */
void phaseI1_batch_interleaved(
    const uint8_t **f4,
    const uint8_t **rc4,
    const int      *read_lens,
    const int      *read_ctrs,
    SMEM          **matchArrays,
    uint64_t      **nSMEMs,
    int            *chosen_strands,
    int             nreads);

#ifdef __cplusplus
}
#endif

#endif /* PHASEA_SLOT_STATE_COMPACT_H */

/* Also declare the batch interleaved bridge function for completeness */
#ifdef __cplusplus
extern "C" {
#endif

/* Declared in src/rosaseed_core_bridge_compact.cpp */
/* int64_t rosaseed_core_seed_batch_interleaved_compact(...) -- see bridge header */

#ifdef __cplusplus
}
#endif
