#include "rosaseed_inprocess.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rosaseed_core_bridge.h"

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

    return rosaseed_core_init(&core_cfg);
}

int64_t rosaseed_inprocess_seed_batch(
    const bseq1_t *seqs,
    int64_t nseq,
    int64_t rid_start,
    rosaseed_hit_t *hitArray,
    int64_t max_hits
)
{
    /* ---- 2-step path: BATCH INTERLEAVED ---- */
    const int64_t STACK_LIMIT = 2048;

    const char **seq_ptrs  = NULL;
    int        *rlens      = NULL;
    int        *rids_arr   = NULL;
    int64_t    *offsets    = NULL;   
    /* stack fallback arrays */
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
            fprintf(stderr,
                "[ERROR] rosaseed_inprocess_seed_batch: alloc failed\n");
            exit(EXIT_FAILURE);
        }
    }

    for (int64_t i = 0; i < nseq; i++) {
        seq_ptrs [i] = seqs[i].seq;
        rlens    [i] = seqs[i].l_seq;
        rids_arr [i] = (int)(rid_start + i);
        offsets  [i] = 0;
    }

    int64_t total_hits = rosaseed_core_seed_batch_interleaved(
        seq_ptrs,
        rlens,
        rids_arr,
        (int)nseq,
        (rosaseed_core_hit_t *)hitArray,
        max_hits,
        offsets);

    if (total_hits < 0) exit(EXIT_FAILURE);

    /*
     * offsets[i] is the index of the first hit for read i in hitArray[].
     * seed_counts for read i are the hits between offsets[i] and
     * offsets[i+1] (or total_hits for the last read).
     */
    for (int64_t i = 0; i < nseq; i++) {
        int64_t hi_start = offsets[i];
        int64_t hi_end   = (i + 1 < nseq) ? offsets[i + 1] : total_hits;
        int     local_rid = (int)i;
        for (int64_t h = hi_start; h < hi_end; ++h)
            hitArray[h].rid = local_rid;
    }

    if (nseq > STACK_LIMIT) {
        free(seq_ptrs);
        free(rlens);
        free(rids_arr);
        free(offsets);
    }

    return total_hits;
}

void rosaseed_inprocess_destroy(void)
{
    rosaseed_core_destroy();
}
