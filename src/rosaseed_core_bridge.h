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