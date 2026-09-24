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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "profiling.h"
#include "macros.h"

uint64_t bwa_single_step_fn_access_count = 0;


void display_stats()
{

    double p1 = tprof[SMEMI] * 1.0 / proc_freq;
    double p2 = tprof[SMEMII] * 1.0 / proc_freq;
    double s1 = tprof[SORT1] * 1.0 / proc_freq;
    double gf = tprof[GAPFILL] * 1.0 / proc_freq;
    double s2 = tprof[SORT2] * 1.0 / proc_freq;
    double sal = tprof[SAL] * 1.0 / proc_freq;
    double seed_total_pI_gf = p1 + s1 + gf + s2;
    double seed_total_pI_pII = p1 + p2 + s2;

    fprintf(stderr, "\n==== RosaSeed timing ====\n");
    /* ---- Summary stats ---- */
    fprintf(stderr, "\n--- Timing breakdown (tprof entries) ---\n");
    fprintf(stderr, "  PhaseI     (Phase I)                      : %0.6lf\n", tprof[SMEMI]* 1.0 / proc_freq);

#ifdef ENABLE_PHASE_II
    fprintf(stderr, "  PhaseII    (Phase II)                     : %0.6lf\n", tprof[SMEMII]* 1.0 / proc_freq);
    fprintf(stderr, "  SORT       (Sort phase I and II seeds)    : %0.6lf\n", tprof[SORT2]* 1.0 / proc_freq);
#else
    fprintf(stderr, "  SORT1      (Sort PI seeds)                : %0.6lf\n", tprof[SORT1]* 1.0 / proc_freq);
    fprintf(stderr, "  Gap fill   (Gap fill phase)               : %0.6lf\n", tprof[GAPFILL]* 1.0 / proc_freq);
    fprintf(stderr, "  SORT2      (Sort after gapfill seeds)     : %0.6lf\n", tprof[SORT2]* 1.0 / proc_freq);
#endif
    
    fprintf(stderr, "  SAL        (SA lookup)                    : %0.6lf\n", tprof[SAL]* 1.0 / proc_freq);

#ifdef ENABLE_PHASE_II
        fprintf(stderr, "  Total seeding time (PhaseI + PhaseII + Sort2) = %0.6lf\n", seed_total_pI_pII);
#else
        fprintf(stderr, "  Total seeding time (Phase1+ Sort1 + GapFill + Sort2) = %0.6lf\n", seed_total_pI_gf);
#endif
    
}
