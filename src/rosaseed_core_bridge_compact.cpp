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

#include "rosaseed_core_bridge_compact.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>

#include "rosaseed_compact/file_dec.h"
#include "rosaseed_compact/load_data.h"
#include "rosaseed_compact/helper_functions.h"
#include "rosaseed_compact/macros.h"
#include "rosaseed_compact/seeding_kernel.h"
#include "rosaseed_compact/gap_fill_phase_compact.h"
#include "rosaseed_compact/phaseA_slot_state_compact.h"

static thread_local SMEM *tl_batch_matchArrays_compact[RS_BATCH];
static thread_local int   tl_batch_bufs_ready_compact = 0;

static void ensure_batch_match_arrays_compact(void)
{
    if (tl_batch_bufs_ready_compact) return;
    for (int s = 0; s < RS_BATCH; ++s) {
        tl_batch_matchArrays_compact[s] =
            (SMEM *)_mm_malloc(MATCH_ARRAY_CAPACITY * sizeof(SMEM), 64);
        if (!tl_batch_matchArrays_compact[s]) {
            fprintf(stderr,
                "[ERROR] RosaSeed 1-step batch matchArray[%d] alloc failed\n", s);
            exit(EXIT_FAILURE);
        }
    }
    tl_batch_bufs_ready_compact = 1;
}

static thread_local SMEM *tl_matchArray = NULL;
static thread_local int64_t *tl_rbeg_out = NULL;
static thread_local uint32_t *tl_seed_offsets = NULL;
static thread_local uint32_t *tl_seed_counts = NULL;
static thread_local int tl_buffers_initialized = 0;

static thread_local uint8_t *tl_loc_f4 = NULL;
static thread_local uint8_t *tl_loc_rc4 = NULL;
static thread_local int tl_read_buf_capacity = 0;







//     tl_buffers_initialized = 1;
// }

static thread_local int64_t tl_rbeg_capacity = 0;

static void rosaseed_core_ensure_thread_buffers(int64_t max_hits_for_call, int32_t read_len)
{
    if (!tl_matchArray) {
        tl_matchArray = (SMEM *)_mm_malloc(MATCH_ARRAY_CAPACITY * sizeof(SMEM), 64);
        tl_seed_offsets = (uint32_t *)malloc(MATCH_ARRAY_CAPACITY * sizeof(uint32_t));
        tl_seed_counts  = (uint32_t *)malloc(MATCH_ARRAY_CAPACITY * sizeof(uint32_t));
    }

    // int64_t needed = MATCH_ARRAY_CAPACITY * g_sa_max_occ;
    // if (needed > max_hits_for_call) needed = max_hits_for_call;
    int64_t needed =
    (int64_t)MATCH_ARRAY_CAPACITY * (int64_t)g_sa_max_occ;

    if (tl_rbeg_capacity < needed) {
        int64_t *newbuf = (int64_t *)realloc(tl_rbeg_out, needed * sizeof(int64_t));
        if (!newbuf) {
            fprintf(stderr, "[ERROR] failed to allocate tl_rbeg_out\n");
            exit(EXIT_FAILURE);
        }
        tl_rbeg_out = newbuf;
        tl_rbeg_capacity = needed;
    }

    if (!tl_matchArray || !tl_seed_offsets || !tl_seed_counts || !tl_rbeg_out) {
        fprintf(stderr, "[ERROR] RosaSeed thread-local buffer allocation failed\n");
        exit(EXIT_FAILURE);
    }

    if (tl_read_buf_capacity < read_len) {
        uint8_t *new_f4 = (uint8_t *)realloc(tl_loc_f4, read_len * sizeof(uint8_t));
        uint8_t *new_rc4 = (uint8_t *)realloc(tl_loc_rc4, read_len * sizeof(uint8_t));
    
        if (!new_f4 || !new_rc4) {
            fprintf(stderr, "[ERROR] failed to allocate read buffers\n");
            exit(EXIT_FAILURE);
        }
    
        tl_loc_f4 = new_f4;
        tl_loc_rc4 = new_rc4;
        tl_read_buf_capacity = read_len;
    }
}


static void rosaseed_core_free_thread_buffers(void)
{
    if (tl_matchArray) {
        _mm_free(tl_matchArray);
        tl_matchArray = NULL;
    }

    if (tl_rbeg_out) {
        free(tl_rbeg_out);
        tl_rbeg_out = NULL;
    }

    if (tl_seed_offsets) {
        free(tl_seed_offsets);
        tl_seed_offsets = NULL;
    }

    if (tl_seed_counts) {
        free(tl_seed_counts);
        tl_seed_counts = NULL;
    }

    if (tl_loc_f4) {
        free(tl_loc_f4);
        tl_loc_f4 = NULL;
    }

    if (tl_loc_rc4) {
        free(tl_loc_rc4);
        tl_loc_rc4 = NULL;
    }

    tl_rbeg_capacity = 0;
    tl_read_buf_capacity = 0;
    tl_buffers_initialized = 0;
}

/* RosaSeed globals expected by RosaSeed files */
const char *REFERENCE_BIN_FILE    = NULL;
const char *JUMP_TABLE_FILE       = NULL;
const char *READS_FILE            = NULL;
int         NUM_OF_READS_IN_FILE  = 2147483647;

const char *SA_LS_WORD_BIN_FILE   = NULL;
const char *SA_MSB_BIT_BIN_FILE   = NULL;
const char *OCC_BIN_FILE          = NULL;

static std::atomic<uint64_t> rs_t_phase1{0};
static std::atomic<uint64_t> rs_t_phase2_or_gapfill{0};
static std::atomic<uint64_t> rs_t_sort{0};
static std::atomic<uint64_t> rs_t_sal{0};
static std::atomic<uint64_t> rs_t_total_core{0};

int g_min_seed_len_A  = 19;
int g_min_seed_len_BC = 20;
int64_t g_min_intv    = 20;
int g_phase1_cap      = 5000;
int g_gap_threshold   = 5;
int g_num_pivots_B    = 3;
int g_gap_left_step   = 4;
int g_sa_max_occ      = 500;
int g_is_fastq = 0;

static int g_rosaseed_core_initialized = 0;
static OneStepMapInfo g_rosaseed_mapinfo;

static char path_ref4_bin[1024];
static char path_jt[1024];
static char path_sa_bin[1024];
static char path_msb_bin[1024];
static char path_occ_bin[1024];

/*
 * rosaseed_to_b4(), 1-step version. Same N-handling contract as the 2-step
 * bridge: an ambiguous base becomes the sentinel 4 and seed extension stops
 * there. See rosaseed_core_bridge.cpp for the full rationale.
 */
static inline uint8_t rosaseed_to_b4(unsigned char c)
{
    if (c <= 3) return c;   // BWA already encoded A/C/G/T as 0/1/2/3
    if (c == 4) return 4;   // BWA's pre-encoded N (nst_nt4_table), sentinel

    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;   // raw 'N' or garbage, sentinel
    }
}

static void rosaseed_seq_to_b4_arrays(const char *seq,
                                      int read_len,
                                      uint8_t *fwd,
                                      uint8_t *rc)
{
    for (int i = 0; i < read_len; ++i) {
        fwd[i] = rosaseed_to_b4((unsigned char)seq[i]);
    }

    for (int i = 0; i < read_len; ++i) {
        uint8_t f = fwd[read_len - 1 - i];
        rc[i] = (f == 4) ? 4 : (uint8_t)(f ^ 0x3u);   /* N stays N under RC */
    }
}

int64_t rosaseed_core_seed_one_read_compact(
    const char *seq,
    int read_len,
    int rid,
    rosaseed_core_hit_t *hits,
    int64_t max_hits)
{

    if (!g_rosaseed_core_initialized) {
        fprintf(stderr, "[ERROR] RosaSeed core not initialized before Phase I\n");
        return -1;
    }


    if (read_len <= 0) return -1;

    rosaseed_core_ensure_thread_buffers(max_hits, read_len);

    SMEM *matchArray = tl_matchArray;
    int64_t *rbeg_out = tl_rbeg_out;
    uint32_t *seed_offsets = tl_seed_offsets;
    uint32_t *seed_counts = tl_seed_counts;

    // uint8_t loc_f4[READ_LENGTH_BASE4]  __attribute__((aligned(64)));
    // uint8_t loc_rc4[READ_LENGTH_BASE4] __attribute__((aligned(64)));

    uint8_t *loc_f4 = tl_loc_f4;
    uint8_t *loc_rc4 = tl_loc_rc4;



    rosaseed_seq_to_b4_arrays(seq, read_len, loc_f4, loc_rc4);



    uint64_t numberofSMEMs = 0;
    int chosen_strand = 0;
    int read_ctr = rid;

    /* Phase I */
    uint64_t t0 = __rdtsc();
    phaseI_routine(matchArray,
                   &numberofSMEMs,
                   &read_ctr,
                   loc_f4,
                   loc_rc4,
                   &chosen_strand,
                   read_len);
    rs_t_phase1.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
    
    uint64_t totalsmems = numberofSMEMs;
    uint64_t smems_before = totalsmems;
    
    #ifdef ENABLE_PHASE_II
    
        uint64_t noofSMEMspII = 0;
    
    #  ifdef ENABLE_F_RC_CHOICE
        if (chosen_strand == 0) {
            t0 = __rdtsc();
            phaseII_routine(loc_f4,
                            matchArray + numberofSMEMs,
                            &noofSMEMspII,
                            &read_ctr,
                            g_min_intv,
                            &chosen_strand,
                            read_len);
            rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
        } else {
            t0 = __rdtsc();
            phaseII_routine(loc_rc4,
                            matchArray + numberofSMEMs,
                            &noofSMEMspII,
                            &read_ctr,
                            g_min_intv,
                            &chosen_strand,
                            read_len);
            rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
        }
    #  else
        t0 = __rdtsc();
        phaseII_routine(loc_f4,
                        matchArray + numberofSMEMs,
                        &noofSMEMspII,
                        &read_ctr,
                        g_min_intv,
                        &chosen_strand,
                        read_len);
        rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
    #  endif
    
        totalsmems += noofSMEMspII;
    
        if (totalsmems > 0) {
            t0 = __rdtsc();
            sort_smems_for_read(matchArray, totalsmems);
            rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
        }
    
    #else
    
        if (totalsmems > 0) {
            t0 = __rdtsc();
            sort_smems_for_read(matchArray, totalsmems);
            rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
        }
    
        smems_before = totalsmems;
        
        t0 = __rdtsc();
        gap_fill_phase(loc_f4,
                        loc_rc4,
                        matchArray,
                        &totalsmems,
                        (uint64_t)MATCH_ARRAY_CAPACITY,
                        read_ctr,
                        g_min_intv,
                        chosen_strand,
                        read_len);
        rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
    
        if (totalsmems > smems_before) {
            t0 = __rdtsc();
            sort_smems_for_read(matchArray, totalsmems);
            rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
        }
    
    #endif

#ifdef DEBUG_ROSASEED_INPROCESS_MODE
    #  ifdef ENABLE_PHASE_II
        fprintf(stderr,
                "[RS_MODE] rid=%d mode=PhaseI+PhaseII phaseI=%lu final=%lu added=%lu\n",
                rid,
                (unsigned long)numberofSMEMs,
                (unsigned long)totalsmems,
                (unsigned long)(totalsmems - numberofSMEMs));
    #  else
        fprintf(stderr,
                "[RS_MODE] rid=%d mode=PhaseI+GapFill phaseI=%lu final=%lu added=%lu\n",
                rid,
                (unsigned long)numberofSMEMs,
                (unsigned long)totalsmems,
                (unsigned long)(totalsmems - numberofSMEMs));
    #  endif
#endif

    memset(seed_offsets, 0, MATCH_ARRAY_CAPACITY * sizeof(uint32_t));
    memset(seed_counts, 0, MATCH_ARRAY_CAPACITY * sizeof(uint32_t));

    t0 = __rdtsc();
    materialize_read_rbegs_base4(
        &g_rosaseed_mapinfo,
        matchArray,
        totalsmems,
        g_sa_max_occ,
        rbeg_out,
        seed_offsets,
        seed_counts
    );
    rs_t_sal.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);

    int64_t out_count = 0;

    for (uint64_t i = 0; i < totalsmems; ++i) {
        int qbeg = matchArray[i].start_idx;
        int len  = matchArray[i].end_idx - matchArray[i].start_idx + 1;
    
        for (uint32_t h = 0; h < seed_counts[i]; ++h) {
            if (out_count >= max_hits) {
                fprintf(stderr,
                        "[ERROR] RosaSeed core hit overflow rid=%d max_hits=%ld\n",
                        rid, (long)max_hits);
                return -1;
            }
            hits[out_count].rid  = rid;
            hits[out_count].qbeg = qbeg;
            hits[out_count].len  = len;
            hits[out_count].rbeg = rbeg_out[seed_offsets[i] + h];
            hits[out_count].occ  = seed_counts[i];
            out_count++;
        }
    }
    
    
    return out_count;
}

int rosaseed_core_init_compact(const rosaseed_core_config_t *cfg)
{
    if (g_rosaseed_core_initialized) return 0;

    if (cfg == NULL || cfg->index_dir == NULL) {
        fprintf(stderr, "[ERROR] RosaSeed core init requires index_dir\n");
        return -1;
    }

    g_min_seed_len_A  = cfg->min_seed_len_A;
    g_min_seed_len_BC = cfg->min_seed_len_BC;
    g_min_intv        = cfg->min_intv;
    g_phase1_cap      = cfg->phase1_cap;
    g_gap_threshold   = cfg->gap_threshold;
    g_num_pivots_B    = cfg->num_pivots_B;
    g_gap_left_step   = cfg->gap_left_step;
    g_sa_max_occ      = (int)cfg->sa_max_occ;

    snprintf(path_ref4_bin, sizeof(path_ref4_bin),
            "%s/ref4_packed.bin",
            cfg->index_dir);

    int cf = 1 << SA_COMPRESSION_FACTOR_POWER;

    snprintf(path_sa_bin, sizeof(path_sa_bin),
             "%s/sa_ls_word_cf%d.bin",
             cfg->index_dir, cf);

    snprintf(path_msb_bin, sizeof(path_msb_bin),
             "%s/sa_ms_byte_cf%d.bin",
             cfg->index_dir, cf);

    snprintf(path_occ_bin, sizeof(path_occ_bin),
             "%s/cp_occ_compact.bin",
             cfg->index_dir);

#if defined(LOAD_JTABLE_15nt)
    snprintf(path_jt, sizeof(path_jt),
             "%s/jumptable_15nt.bin",
             cfg->index_dir);
#elif defined(LOAD_JTABLE_16nt)
    snprintf(path_jt, sizeof(path_jt),
             "%s/jumptable_16nt.bin",
             cfg->index_dir);
#else
    snprintf(path_jt, sizeof(path_jt),
             "%s/jumptable_14nt.bin",
             cfg->index_dir);
#endif

    REFERENCE_BIN_FILE = path_ref4_bin;
    JUMP_TABLE_FILE       = path_jt;
    SA_LS_WORD_BIN_FILE   = path_sa_bin;
    SA_MSB_BIT_BIN_FILE   = path_msb_bin;
    OCC_BIN_FILE          = path_occ_bin;

    fprintf(stderr, "[INFO] RosaSeed core init\n");
    fprintf(stderr, "[INFO] RosaSeed index dir: %s\n", cfg->index_dir);
    fprintf(stderr, "[INFO] RosaSeed ref4 bin: %s\n", REFERENCE_BIN_FILE);
    fprintf(stderr, "[INFO] RosaSeed OCC     : %s\n", OCC_BIN_FILE);
    fprintf(stderr, "[INFO] RosaSeed SA LS   : %s\n", SA_LS_WORD_BIN_FILE);
    fprintf(stderr, "[INFO] RosaSeed SA MS   : %s\n", SA_MSB_BIT_BIN_FILE);
    fprintf(stderr, "[INFO] RosaSeed JT      : %s\n", JUMP_TABLE_FILE);

    fprintf(stderr,
        "[RS_DEBUG] SA_COMPRESSION_FACTOR_POWER=%d CF=%d\n",
        SA_COMPRESSION_FACTOR_POWER,
        1 << SA_COMPRESSION_FACTOR_POWER);

    rosaseed_compact_read_index_metadata();   /* sets BWT_SIZE_REFERENCE_SIZE, rosaseed_L */
    allocate_memory();
    load_bwt_data_structures();
    init_one_step_mapinfo(&g_rosaseed_mapinfo, rosaseed_L);

    g_rosaseed_core_initialized = 1;

    fprintf(stderr, "Phase 1 cap : %d\n", g_phase1_cap);
    fprintf(stderr, "[INFO] RosaSeed core index loaded successfully\n");
    return 0;
}


void rosaseed_core_destroy(void)
{
    if (!g_rosaseed_core_initialized) return;

    rosaseed_core_free_thread_buffers();

    /* free per-slot batch matchArrays */
    if (tl_batch_bufs_ready_compact) {
        for (int s = 0; s < RS_BATCH; ++s) {
            if (tl_batch_matchArrays_compact[s]) {
                _mm_free(tl_batch_matchArrays_compact[s]);
                tl_batch_matchArrays_compact[s] = NULL;
            }
        }
        tl_batch_bufs_ready_compact = 0;
    }

    free_memory();
    g_rosaseed_core_initialized = 0;
}

/* -----------------------------------------------------------------------
   rosaseed_core_seed_batch_interleaved_compact()
   -----------------------------------------------------------------------
   Batch coroutine entry point for the 1-step path.
   Mirrors rosaseed_core_seed_batch_interleaved() from the 2-step bridge
   exactly, but calls phaseI1_batch_interleaved() and uses 1-step types.
   ----------------------------------------------------------------------- */
int64_t rosaseed_core_seed_one_read_compact(
    const char *seq,
    int read_len,
    int rid,
    rosaseed_core_hit_t *hits,
    int64_t max_hits);

int64_t rosaseed_core_seed_batch_interleaved_compact(
    const char **seqs,
    const int   *read_lens,
    const int   *rids,
    int          nreads,
    rosaseed_core_hit_t *hits,
    int64_t      max_hits,
    int64_t     *hit_offsets   /* can be NULL */
);   /* forward-decl; defined below */

int64_t rosaseed_core_seed_batch_interleaved_compact(
    const char **seqs,
    const int   *read_lens,
    const int   *rids,
    int          nreads,
    rosaseed_core_hit_t *hits,
    int64_t      max_hits,
    int64_t     *hit_offsets)
{
    if (!g_rosaseed_core_initialized) {
        fprintf(stderr, "[ERROR] RosaSeed 1-step core not initialized\n");
        return -1;
    }
    if (nreads <= 0) return 0;

    ensure_batch_match_arrays_compact();

    int max_rlen = 0;
    for (int i = 0; i < nreads; ++i)
        if (read_lens[i] > max_rlen) max_rlen = read_lens[i];

    if (max_rlen > RS_BATCH_MAX_READ_LEN) {
        fprintf(stderr,
            "[ERROR] 1-step read_len=%d > RS_BATCH_MAX_READ_LEN=%d\n",
            max_rlen, RS_BATCH_MAX_READ_LEN);
        return -1;
    }

    rosaseed_core_ensure_thread_buffers(max_hits, max_rlen);

    int64_t  *rbeg_out     = tl_rbeg_out;
    uint32_t *seed_offsets = tl_seed_offsets;
    uint32_t *seed_counts  = tl_seed_counts;

    /* per-slot encode buffers on stack */
    uint8_t f4_bufs [RS_BATCH][RS_BATCH_MAX_READ_LEN];
    uint8_t rc4_bufs[RS_BATCH][RS_BATCH_MAX_READ_LEN];

    int64_t total_hits_out = 0;
    int     base           = 0;

    while (base < nreads) {
        int batch_n = nreads - base;
        if (batch_n > RS_BATCH) batch_n = RS_BATCH;

        /* encode all reads in this micro-batch */
        const uint8_t *f4_ptrs [RS_BATCH];
        const uint8_t *rc4_ptrs[RS_BATCH];
        int             rlens_b [RS_BATCH];
        int             rctrs_b [RS_BATCH];
        SMEM           *ma_ptrs [RS_BATCH];
        uint64_t        nsmems_b[RS_BATCH];
        uint64_t       *nsmems_ptrs[RS_BATCH];
        int             chosen_strands[RS_BATCH];

        for (int s = 0; s < batch_n; ++s) {
            int idx  = base + s;
            int rlen = read_lens[idx];

            rosaseed_seq_to_b4_arrays(seqs[idx], rlen,
                                      f4_bufs[s], rc4_bufs[s]);

            f4_ptrs [s]    = f4_bufs[s];
            rc4_ptrs[s]    = rc4_bufs[s];
            rlens_b [s]    = rlen;
            rctrs_b [s]    = rids[idx];
            ma_ptrs [s]    = tl_batch_matchArrays_compact[s];
            nsmems_b[s]    = 0;
            nsmems_ptrs[s] = &nsmems_b[s];
            chosen_strands[s] = 0;
        }

        /* coroutine Phase I */
        uint64_t t0 = __rdtsc();
        phaseI1_batch_interleaved(
            f4_ptrs, rc4_ptrs,
            rlens_b, rctrs_b,
            ma_ptrs, nsmems_ptrs,
            chosen_strands,
            batch_n);
        rs_t_phase1.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);

        /* post-Phase-I + materialise per read */
        for (int s = 0; s < batch_n; ++s) {
            int     idx    = base + s;
            int     rlen   = read_lens[idx];
            int     rid    = rids[idx];
            SMEM   *ma     = tl_batch_matchArrays_compact[s];
            uint64_t numberofSMEMs = nsmems_b[s];
            int      chosen_strand = chosen_strands[s];
            int      read_ctr = rid;

            uint64_t totalsmems  = numberofSMEMs;
            uint64_t smems_before = totalsmems;

            uint8_t *lf4  = f4_bufs[s];
            uint8_t *lrc4 = rc4_bufs[s];

#ifdef ENABLE_PHASE_II
            /* Suggestion 2: cross-read Phase II prefetch (1-step).
               Before running phaseII_routine for read s, prefetch cp_occ
               for reads s+1..s+pf_dist's first Phase II pivot.
               Uses JT_LEN_NT (1-step) for jump address computation. */
            {
                const int pf_dist = 3;
                for (int ahead = 1; ahead <= pf_dist && s + ahead < batch_n; ++ahead) {
                    const uint8_t *pf_pat =
                        (chosen_strands[s + ahead] == 0)
                        ? f4_bufs[s + ahead]
                        : rc4_bufs[s + ahead];
                    int pf_rlen = read_lens[base + s + ahead];

                    int pf_pivot = pf_rlen - 1;   /* pivots[0] = rightmost base */

                    if (pf_pivot >= JT_LEN_NT - 1) {
                        uint64_t pf_addr = 0;
                        int pf_valid = 1;
                        for (int k = 0; k < JT_LEN_NT; ++k) {
                            uint8_t nt = pf_pat[pf_pivot - k];
                            if (nt >= 4) { pf_valid = 0; break; }
                            pf_addr |= ((uint64_t)nt) << (2 * k);
                        }

                        if (pf_valid) {
                            uint64_t pf_jp = jump_pointers[pf_addr];
                            uint64_t pf_l = 0, pf_h = 0, pf_diff = 0;
                            extract_jump_bounds(pf_jp, &pf_l, &pf_h, &pf_diff);

                            if (pf_l < pf_h) {
                                __builtin_prefetch((const char *)&cp_occ[pf_l >> 5],      0, 1);
                                __builtin_prefetch((const char *)&cp_occ[pf_l >> 5] + 64, 0, 1);
                                if ((pf_h >> 5) != (pf_l >> 5)) {
                                    __builtin_prefetch((const char *)&cp_occ[pf_h >> 5],      0, 1);
                                    __builtin_prefetch((const char *)&cp_occ[pf_h >> 5] + 64, 0, 1);
                                }
                            }
                        }
                        /* pf_valid==0: N in the prefetch window, nothing
                           valid to prefetch, skip silently (this is only
                           a speculative prefetch hint, not correctness-critical) */
                    }
                }
            }

            uint64_t noofSMEMspII = 0;
#  ifdef ENABLE_F_RC_CHOICE
            if (chosen_strand == 0) {
                t0 = __rdtsc();
                phaseII_routine(lf4, ma + numberofSMEMs, &noofSMEMspII,
                                &read_ctr, g_min_intv, &chosen_strand, rlen);
                rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                    std::memory_order_relaxed);
            } else {
                t0 = __rdtsc();
                phaseII_routine(lrc4, ma + numberofSMEMs, &noofSMEMspII,
                                &read_ctr, g_min_intv, &chosen_strand, rlen);
                rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                    std::memory_order_relaxed);
            }
#  else
            t0 = __rdtsc();
            phaseII_routine(lf4, ma + numberofSMEMs, &noofSMEMspII,
                            &read_ctr, g_min_intv, &chosen_strand, rlen);
            rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                std::memory_order_relaxed);
#  endif
            totalsmems += noofSMEMspII;
            if (totalsmems > 0) {
                t0 = __rdtsc();
                sort_smems_for_read(ma, totalsmems);
                rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
            }
#else
            if (totalsmems > 0) {
                t0 = __rdtsc();
                sort_smems_for_read(ma, totalsmems);
                rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
            }
            smems_before = totalsmems;
            t0 = __rdtsc();
            gap_fill_phase(lf4, lrc4, ma, &totalsmems,
                           (uint64_t)MATCH_ARRAY_CAPACITY,
                           read_ctr, g_min_intv, chosen_strand, rlen);
            rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                std::memory_order_relaxed);
            if (totalsmems > smems_before) {
                t0 = __rdtsc();
                sort_smems_for_read(ma, totalsmems);
                rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
            }
#endif

            /* materialise */
            if (hit_offsets) hit_offsets[idx] = total_hits_out;

            memset(seed_offsets, 0, totalsmems * sizeof(uint32_t));
            memset(seed_counts,  0, totalsmems * sizeof(uint32_t));

            t0 = __rdtsc();
            materialize_read_rbegs_base4(
                &g_rosaseed_mapinfo, ma, totalsmems,
                g_sa_max_occ, rbeg_out, seed_offsets, seed_counts);
            rs_t_sal.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);

            for (uint64_t i = 0; i < totalsmems; ++i) {
                int qbeg = ma[i].start_idx;
                int len  = ma[i].end_idx - ma[i].start_idx + 1;
                for (uint32_t h = 0; h < seed_counts[i]; ++h) {
                    if (total_hits_out >= max_hits) {
                        fprintf(stderr,
                            "[ERROR] 1-step hit overflow rid=%d\n", rid);
                        return -1;
                    }
                    hits[total_hits_out].rid  = rid;
                    hits[total_hits_out].qbeg = qbeg;
                    hits[total_hits_out].len  = len;
                    hits[total_hits_out].rbeg = rbeg_out[seed_offsets[i] + h];
                    hits[total_hits_out].occ  = seed_counts[i];
                    total_hits_out++;
                }
            }
        }

        base += batch_n;
    }

    return total_hits_out;
}

void rosaseed_core_print_timing(double proc_freq)
{
    fprintf(stderr, "Phase I        : %.3f sec\n",
        rs_t_phase1.load(std::memory_order_relaxed) / proc_freq);

    fprintf(stderr, "PhaseII/GapFill: %.3f sec\n",
            rs_t_phase2_or_gapfill.load(std::memory_order_relaxed) / proc_freq);

    fprintf(stderr, "Sort           : %.3f sec\n",
            rs_t_sort.load(std::memory_order_relaxed) / proc_freq);

    fprintf(stderr, "SAL/materialize: %.3f sec\n",
            rs_t_sal.load(std::memory_order_relaxed) / proc_freq);
}