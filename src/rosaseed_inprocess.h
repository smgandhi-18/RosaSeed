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