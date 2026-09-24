/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#include <stdio.h>
#include "macro.h"
#include <stdint.h>
#include <assert.h>
#include "profiling.h"

int find_opt(uint64_t *a, int len, uint64_t *max, uint64_t *min, double *avg)
{
    *max = 0;
    *min = 1e15;
    *avg = 0;

    int i=0;
    for (i=0; i<len; i++)
    {
        if (a[i] > *max) *max = a[i];
        if (a[i] < *min) *min = a[i];
        *avg += a[i];
    }
    *avg /= len;

    return 1;
}

int display_stats(int nthreads)
{
    uint64_t max, min;
    double avg;
    fprintf(stderr, "No. of OMP threads: %d\n", nthreads);
    fprintf(stderr, "Processor is running @%lf MHz\n", proc_freq*1.0/1e6);
    fprintf(stderr, "Runtime profile:\n");

    fprintf(stderr, "\n\tTime taken for main_mem function: %0.2lf sec\n\n",
            tprof[MEM][0]*1.0/proc_freq);

    fprintf(stderr, "\tIO times (sec) :\n");
    find_opt(tprof[READ_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tReading IO time (reads) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[SAM_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tWriting IO time (SAM) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[REF_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tReading IO time (Reference Genome) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[FMI], 1, &max, &min, &avg);
    fprintf(stderr, "\tIndex read time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    fprintf(stderr, "\n\tOverall time (sec) (Excluding Index reading time):\n");
    // find_opt(tprof[PROCESS], 1, &max, &min, &avg);
    fprintf(stderr, "\tPROCESS() (Total compute time + (read + SAM) IO time) : %0.2lf\n",
            tprof[PROCESS][0]*1.0/proc_freq);

    find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
    fprintf(stderr, "\tMEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    fprintf(stderr, "\n\t SAM Processing time (sec):\n");
    find_opt(tprof[WORKER20], 1, &max, &min, &avg);
    fprintf(stderr, "\t--WORKER_SAM avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
#if HIDE
    find_opt(tprof[SAM1], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[SAM2], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM2 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[ALIGN1], 1, &max, &min, &avg);
    fprintf(stderr, "\t\t\tWORKER_ALIGN1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[SAM3], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM3 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
    
    /* ------------------------------------------------------------------
       Kernel breakdown.

       The categories below are disjoint and add up to the total, so no
       component is counted twice:

         seed processing  seed search + seed sorting + SA (suffix array)
                          lookup.
                          RosaSeed: tprof[MEM_COLLECT] wraps the whole
                            rosaseed_inprocess_seed_batch() call, which already
                            performs Phase A/B/C, the seed sort and the SA
                            materialisation (tprof[MEM_SA] is unused, 0).
                          BWA-MEM2: tprof[MEM_COLLECT] is mem_collect_smem()
                            (SMEM search incl. sortSMEMs()) and tprof[MEM_SA]
                            is the separately timed SA lookup, which happens
                            inside mem_chain_seeds().
         chaining         chain construction from the seeds.
                          RosaSeed: mem_chain_preexpanded_hits(), timed whole
                            in tprof[MEM_SA_BLOCK] (includes the inline
                            pre-chain singleton filter, which has no timer of
                            its own).
                          BWA-MEM2: mem_chain_seeds() minus the SA lookup
                            already attributed to seed processing, i.e.
                            tprof[MEM_SA_BLOCK] - tprof[MEM_SA].
         chain filtering  mem_chain_flt() + mem_flt_chained_seeds().
         BSW              banded Smith-Waterman extension.
         unattributed     total kernel minus the above; bookkeeping residue.
       ------------------------------------------------------------------ */
    fprintf(stderr, "\n\tKernels' compute time (sec):\n");

    double t_total = 0.0, t_smem = 0.0, t_sa = 0.0, t_sa_block = 0.0;
    double t_chain_flt = 0.0, t_flt_seeds = 0.0, t_bsw = 0.0;

    find_opt(tprof[WORKER10], 1, &max, &min, &avg);
    t_total = avg*1.0/proc_freq;
    find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
    t_smem = avg*1.0/proc_freq;
    find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
    t_sa = avg*1.0/proc_freq;
    find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
    t_sa_block = avg*1.0/proc_freq;
    find_opt(tprof[MEM_CHAIN_FLT], nthreads, &max, &min, &avg);
    t_chain_flt = avg*1.0/proc_freq;
    find_opt(tprof[FLT_CHAINED_SEEDS], nthreads, &max, &min, &avg);
    t_flt_seeds = avg*1.0/proc_freq;
    find_opt(tprof[MEM_ALN2], nthreads, &max, &min, &avg);
    t_bsw = avg*1.0/proc_freq;

#ifdef ROSASEED_INPROCESS
    const double t_seed_proc = t_smem;              /* SA lookup is inside it */
    const double t_chaining  = t_sa_block;          /* mem_chain_preexpanded_hits */
#else
    const double t_seed_proc = t_smem + t_sa;       /* SMEM search + SA lookup */
    const double t_chaining  = t_sa_block - t_sa;   /* chain build, SA removed */
#endif
    const double t_sum = t_seed_proc + t_chaining + t_chain_flt + t_flt_seeds + t_bsw;

    fprintf(stderr, "\tTotal kernel (seeding+SA+chaining+BSW) time: %0.2lf\n", t_total);

#ifdef ROSASEED_INPROCESS
    fprintf(stderr, "\t\tSeed processing (RosaSeed seeding + sort + SA lookup): %0.2lf\n",
            t_seed_proc);
    fprintf(stderr, "\t\t\t(per-phase split printed by the RosaSeed bridge below)\n");
    fprintf(stderr, "\t\tChaining (mem_chain_preexpanded_hits, incl. pre-chain filter): %0.2lf\n",
            t_chaining);
#else
    fprintf(stderr, "\t\tSeed processing (SMEM search + sort + SA lookup): %0.2lf\n",
            t_seed_proc);
    fprintf(stderr, "\t\t\tSMEM search + sort (mem_collect_smem): %0.2lf\n", t_smem);
    fprintf(stderr, "\t\t\tSA lookup (get_sa_entries*): %0.2lf\n", t_sa);
    fprintf(stderr, "\t\tChaining (mem_chain_seeds minus SA lookup): %0.2lf\n", t_chaining);
#endif
    fprintf(stderr, "\t\tChain filtering (mem_chain_flt): %0.2lf\n", t_chain_flt);
    fprintf(stderr, "\t\tChained-seed filtering (mem_flt_chained_seeds): %0.2lf\n", t_flt_seeds);
    fprintf(stderr, "\t\tBSW extension: %0.2lf\n", t_bsw);
    fprintf(stderr, "\t\tUnattributed (total - above): %0.2lf\n", t_total - t_sum);

    #if HIDE
    int agg1 = 0, agg2 = 0, agg3 = 0;
    for (int i=0; i<nthreads; i++) {
        agg1 += tprof[PE11][i];
        agg2 += tprof[PE12][i];
        agg3 += tprof[PE13][i];
    }
    if (agg1 != agg3) 
        fprintf(stderr, "There is a discrepancy re-allocs, plz rectify!!\n");

    if(agg2 > 0)
    {
        fprintf(stderr, "\n\tTotal re-allocs: %d out of total requests: %d, Rate: %0.2f\n",
                agg1, agg2, agg1*1.0/agg2);
    }

    double res, max_ = 0, min_=1e10;
    for (int i=0; i<nthreads; i++) {
        double val = (tprof[ALIGN1][i]*1.0) / tprof[MEM_CHAIN][i];
        res += val;
        if (max_ < val) max_ = val;
        if (min_ > val) min_ = val;
    }
    fprintf(stderr, "\tAvg. FM-index traversal per get_sa_entry(): avg: %lf, max: %lf, min: %lf\n",
            res/nthreads, max_, min_);

    int64_t tot_inst1 = 0, tot_inst2 = 0;
    for (int i=0; i<nthreads; i++) {
        tot_inst1 += tprof[SAM1][i];
        tot_inst2 += tprof[SAM2][i];
    }
    
    fprintf(stderr, "\ttot_inst1: %ld, tot_inst2: %ld, over %d threads\n",
            tot_inst1, tot_inst2, nthreads);
    #endif
    
#if HIDE
    fprintf(stderr, "\n BSW Perf.:\n");
    find_opt(tprof[MEM_ALN2_B], 1, &max, &min, &avg);
    fprintf(stderr, "\tLeft 16-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_D], 1, &max, &min, &avg);
    fprintf(stderr, "\tLeft 8-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_C], 1, &max, &min, &avg);
    fprintf(stderr, "\tRight 16-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_E], 1, &max, &min, &avg);
    fprintf(stderr, "\tRight 8-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif

#if HIDE
    fprintf(stderr, "\nSTATSV\n");
    fprintf(stderr, "%0.2lf\n", tprof[PROCESS][0]*1.0/proc_freq);
    find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[WORKER10], 1, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    double val = 0;
    find_opt(tprof[MEM_ALN2_UP], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[CLEFT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[CRIGHT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[MEM_ALN2_DOWN], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    fprintf(stderr, "%0.2lf\n", val*1.0/proc_freq);
    fprintf(stderr, "%0.2lf\n", (tprof[REF_IO][0] + tprof[FMI][0])*1.0/proc_freq);
    fprintf(stderr, "%0.2lf\n", tprof[READ_IO][0]*1.0/proc_freq);
#endif
    // printf("\tMemory usage (GB):\n");
    // printf("\tAvg: %0.2lf, Peak: %0.2lf\n", tprof[PE21][0]*1.0/1e9, tprof[PE22][0]*1.0/1e9);
        
    return 1;
}

