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
#ifndef PHASEA_DEC_H
#define PHASEA_DEC_H

#include <stdint.h>
#include <stddef.h>

/* Maximum read length the slot arrays are sized for.
   Must be >= any actual read_len passed at runtime   */
#ifndef RS_BATCH_MAX_READ_LEN
#  define RS_BATCH_MAX_READ_LEN 251
#endif

/* Number of reads processed simultaneously in the interleaver  */
#ifndef RS_BATCH
#  define RS_BATCH 32
#endif

/*
 * Everything phaseI_routine uses as "local" state for one read.
 * The probe function fills chosen_strand / l0 / h0 / d0 and sets
 * done=0.  The bridge then calls phaseI_routine normally, the
 * probe result is NOT passed back in; phaseI_routine recomputes
 * the first jump internally (it is O(1) table lookup so the cost
 * is negligible compared to the FM extension that follows).
 *
 * The reason we keep the full f4/rc4 arrays here is so the bridge
 * can prefetch cp_occ blocks for MULTIPLE reads before processing
 * any of them, hiding cross-read RAM latency.
 */
typedef struct {
    /* base-4 encoded read arrays, filled before the probe */
    uint8_t  f4 [RS_BATCH_MAX_READ_LEN];
    uint8_t  rc4[RS_BATCH_MAX_READ_LEN];

    /* strand + first-pivot BWT interval, filled by probe */
    int      chosen_strand;   /* 0=fwd, 1=RC                  */
    uint64_t l0;              /* BWT low  after first jump     */
    uint64_t h0;              /* BWT high after first jump     */
    uint64_t d0;              /* h0 - l0                       */

    /* bookkeeping */
    int      read_len;
    int      global_rid;
    int      active;          /* 1 = slot still needs processing */
} PhaseI_SlotState;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * phaseI_probe_first_pivot()
 *
 * Does ONLY the strand-choice probe + first jump-table lookup for one read.
 * Writes chosen_strand, l0, h0, d0 into *slot.
 * Issues __builtin_prefetch for the cp_occ blocks those BWT positions
 * live in, so the hardware can start fetching them from RAM while the
 * caller goes on to probe other reads.
 *
 * This function does NOT touch matchArray or numberofSMEMs.
 * It does NOT perform any FM backward extension.
 * phaseI_routine() must still be called afterwards for the full extension.
 */
void phaseI_probe_first_pivot(PhaseI_SlotState *slot);

void phaseI_batch_interleaved(
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

#endif /* PHASEA_DEC_H */
