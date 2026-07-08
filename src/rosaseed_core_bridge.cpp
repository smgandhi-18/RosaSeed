#include "rosaseed_core_bridge.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>

#include "rosaseed/file_dec.h"
#include "rosaseed/read_init.h"
#include "rosaseed/load_data.h"
#include "rosaseed/helper_functions.h"
#include "rosaseed/macros.h"
#include "rosaseed/seeding_kernel_mt.h"
#include "rosaseed/gap_fill_phase.h"

#include "rosaseed/phaseA_dec.h"

static thread_local SMEM     *tl_matchArray        = NULL;
static thread_local int64_t  *tl_rbeg_out          = NULL;
static thread_local uint32_t *tl_seed_offsets      = NULL;
static thread_local uint32_t *tl_seed_counts       = NULL;
static thread_local int       tl_buffers_initialized = 0;

static thread_local uint8_t  *tl_loc_f4            = NULL;
static thread_local uint8_t  *tl_loc_rc4           = NULL;
static thread_local int       tl_read_buf_capacity  = 0;
static thread_local int64_t   tl_rbeg_capacity      = 0;

/* -----------------------------------------------------------------------
   Per-slot match arrays for the batch interleaver.
   Each slot gets its own MATCH_ARRAY_CAPACITY-sized SMEM buffer so that
   the RS_BATCH reads can be processed independently without interfering.
   These are thread-local so multi-threaded BWA-MEM2 workers don't share.
   ----------------------------------------------------------------------- */
static thread_local SMEM     *tl_batch_matchArrays[RS_BATCH];
static thread_local int       tl_batch_bufs_ready = 0;

static void ensure_batch_match_arrays(void)
{
    if (tl_batch_bufs_ready) return;
    for (int s = 0; s < RS_BATCH; ++s) {
        tl_batch_matchArrays[s] =
            (SMEM *)_mm_malloc(MATCH_ARRAY_CAPACITY * sizeof(SMEM), 64);
        if (!tl_batch_matchArrays[s]) {
            fprintf(stderr,
                "[ERROR] RosaSeed batch matchArray[%d] alloc failed\n", s);
            exit(EXIT_FAILURE);
        }
    }
    tl_batch_bufs_ready = 1;
}

/* -----------------------------------------------------------------------
   RosaSeed globals
   ----------------------------------------------------------------------- */
const char *REFERENCE_GENOME_FILE = NULL;
const char *C_VEC_FILE            = NULL;
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

int     g_min_seed_len_A  = 19;
int     g_min_seed_len_BC = 20;
int64_t g_min_intv        = 20;
int     g_phase1_cap      = 5000;
int     g_gap_threshold   = 5;
int     g_num_pivots_B    = 3;
int     g_gap_left_step   = 4;
int     g_sa_max_occ      = 500;
int     g_is_fastq        = 0;

static int            g_rosaseed_core_initialized = 0;
static TwoStepMapInfo g_rosaseed_mapinfo;

static char path_ref16_bin[1024];
static char path_cvec[1024];
static char path_jt[1024];
static char path_sa_bin[1024];
static char path_msb_bin[1024];
static char path_occ_bin[1024];

/* -----------------------------------------------------------------------
   Helpers 
   ----------------------------------------------------------------------- */
static inline uint8_t rosaseed_to_b4(unsigned char c)
{
    if (c <= 3) return c;
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 0;
    }
}

static void rosaseed_seq_to_b4_arrays(const char *seq,
                                      int read_len,
                                      uint8_t *fwd,
                                      uint8_t *rc)
{
    for (int i = 0; i < read_len; ++i)
        fwd[i] = rosaseed_to_b4((unsigned char)seq[i]);
    for (int i = 0; i < read_len; ++i)
        rc[i] = fwd[read_len - 1 - i] ^ 0x3u;
}

static void rosaseed_core_ensure_thread_buffers(int64_t max_hits_for_call,
                                                int32_t read_len)
{
    if (!tl_matchArray) {
        tl_matchArray   = (SMEM *)_mm_malloc(
                              MATCH_ARRAY_CAPACITY * sizeof(SMEM), 64);
        tl_seed_offsets = (uint32_t *)malloc(
                              MATCH_ARRAY_CAPACITY * sizeof(uint32_t));
        tl_seed_counts  = (uint32_t *)malloc(
                              MATCH_ARRAY_CAPACITY * sizeof(uint32_t));
    }

    int64_t needed = (int64_t)MATCH_ARRAY_CAPACITY * g_sa_max_occ;
    if (needed > max_hits_for_call) needed = max_hits_for_call;

    if (tl_rbeg_capacity < needed) {
        int64_t *newbuf =
            (int64_t *)realloc(tl_rbeg_out, needed * sizeof(int64_t));
        if (!newbuf) {
            fprintf(stderr, "[ERROR] failed to allocate tl_rbeg_out\n");
            exit(EXIT_FAILURE);
        }
        tl_rbeg_out      = newbuf;
        tl_rbeg_capacity = needed;
    }

    if (!tl_matchArray || !tl_seed_offsets ||
        !tl_seed_counts || !tl_rbeg_out) {
        fprintf(stderr,
            "[ERROR] RosaSeed thread-local buffer allocation failed\n");
        exit(EXIT_FAILURE);
    }

    if (tl_read_buf_capacity < read_len) {
        uint8_t *new_f4  = (uint8_t *)realloc(tl_loc_f4,
                               read_len * sizeof(uint8_t));
        uint8_t *new_rc4 = (uint8_t *)realloc(tl_loc_rc4,
                               read_len * sizeof(uint8_t));
        if (!new_f4 || !new_rc4) {
            fprintf(stderr,
                "[ERROR] failed to allocate read buffers\n");
            exit(EXIT_FAILURE);
        }
        tl_loc_f4  = new_f4;
        tl_loc_rc4 = new_rc4;
        tl_read_buf_capacity = read_len;
    }
}

/* -----------------------------------------------------------------------
   Helper: run everything AFTER Phase I for one read.
   Extracted so both the serial path and the batch path can call it
   without duplicating the #ifdef ENABLE_PHASE_II / gap_fill logic.
   matchArray must already contain the Phase I seeds (numberofSMEMs).
   On return, totalsmems_out holds the final count.
   ----------------------------------------------------------------------- */
static void run_post_phaseI(
    SMEM        *matchArray,
    uint64_t     numberofSMEMs,
    uint64_t    *totalsmems_out,
    int         *chosen_strand,
    int          read_ctr,
    int          read_len,
    uint8_t     *loc_f4,
    uint8_t     *loc_rc4)
{
    uint64_t totalsmems = numberofSMEMs;
    uint64_t smems_before;
    uint64_t t0;

#ifdef ENABLE_PHASE_II

    uint64_t noofSMEMspII = 0;

#  ifdef ENABLE_F_RC_CHOICE
    if (*chosen_strand == 0) {
        t0 = __rdtsc();
        phaseII_routine(loc_f4,
                        matchArray + numberofSMEMs,
                        &noofSMEMspII,
                        &read_ctr,
                        g_min_intv,
                        chosen_strand,
                        read_len);
        rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                                         std::memory_order_relaxed);
    } else {
        t0 = __rdtsc();
        phaseII_routine(loc_rc4,
                        matchArray + numberofSMEMs,
                        &noofSMEMspII,
                        &read_ctr,
                        g_min_intv,
                        chosen_strand,
                        read_len);
        rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                                         std::memory_order_relaxed);
    }
#  else
    t0 = __rdtsc();
    phaseII_routine(loc_f4,
                    matchArray + numberofSMEMs,
                    &noofSMEMspII,
                    &read_ctr,
                    g_min_intv,
                    chosen_strand,
                    read_len);
    rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                                     std::memory_order_relaxed);
#  endif

    totalsmems += noofSMEMspII;

    if (totalsmems > 0) {
        t0 = __rdtsc();
        sort_smems_for_read(matchArray, totalsmems);
        rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
    }

#else  /* gap-fill path */

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
                   *chosen_strand,
                   read_len);
    rs_t_phase2_or_gapfill.fetch_add(__rdtsc() - t0,
                                     std::memory_order_relaxed);

    if (totalsmems > smems_before) {
        t0 = __rdtsc();
        sort_smems_for_read(matchArray, totalsmems);
        rs_t_sort.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);
    }

#endif

    *totalsmems_out = totalsmems;
}

/* -----------------------------------------------------------------------
   Helper: materialise hits from matchArray into hits[] array.
   Returns number of hits written, or -1 on overflow.
   ----------------------------------------------------------------------- */
static int64_t materialise_hits(
    SMEM           *matchArray,
    uint64_t        totalsmems,
    int64_t        *rbeg_out,
    uint32_t       *seed_offsets,
    uint32_t       *seed_counts,
    rosaseed_core_hit_t *hits,
    int64_t         max_hits,
    int             rid)
{
    memset(seed_offsets, 0, totalsmems * sizeof(uint32_t));
    memset(seed_counts,  0, totalsmems * sizeof(uint32_t));

    uint64_t t0 = __rdtsc();
    materialize_read_rbegs_base4(
        &g_rosaseed_mapinfo,
        matchArray,
        (uint32_t)totalsmems,
        (uint32_t)g_sa_max_occ,
        rbeg_out,
        seed_offsets,
        seed_counts);
    rs_t_sal.fetch_add(__rdtsc() - t0, std::memory_order_relaxed);

    int64_t out_count = 0;
    for (uint64_t i = 0; i < totalsmems; ++i) {
        int qbeg = matchArray[i].start_idx;
        int len  = matchArray[i].end_idx - matchArray[i].start_idx + 1;
        for (uint32_t h = 0; h < seed_counts[i]; ++h) {
            if (out_count >= max_hits) {
                fprintf(stderr,
                    "[ERROR] RosaSeed core hit overflow rid=%d "
                    "max_hits=%ld\n", rid, (long)max_hits);
                return -1;
            }
            hits[out_count].qbeg = qbeg;
            hits[out_count].len  = len;
            hits[out_count].rbeg = rbeg_out[seed_offsets[i] + h];
            hits[out_count].occ  = seed_counts[i];
            out_count++;
        }
    }
    return out_count;
}

/* -----------------------------------------------------------------------
   rosaseed_core_seed_one_read()
   ----------------------------------------------------------------------- */
int64_t rosaseed_core_seed_one_read(
    const char *seq,
    int read_len,
    int rid,
    rosaseed_core_hit_t *hits,
    int64_t max_hits)
{
    if (!g_rosaseed_core_initialized) {
        fprintf(stderr,
            "[ERROR] RosaSeed core not initialized before Phase I\n");
        return -1;
    }
    if (read_len <= 0) return -1;

    rosaseed_core_ensure_thread_buffers(max_hits, read_len);

    SMEM     *matchArray   = tl_matchArray;
    int64_t  *rbeg_out     = tl_rbeg_out;
    uint32_t *seed_offsets = tl_seed_offsets;
    uint32_t *seed_counts  = tl_seed_counts;
    uint8_t  *loc_f4       = tl_loc_f4;
    uint8_t  *loc_rc4      = tl_loc_rc4;

    rosaseed_seq_to_b4_arrays(seq, read_len, loc_f4, loc_rc4);

    uint64_t numberofSMEMs = 0;
    int      chosen_strand = 0;
    int      read_ctr      = rid;

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

    /* Phase II / gap-fill / sort — via shared helper */
    uint64_t totalsmems = 0;
    run_post_phaseI(matchArray, numberofSMEMs, &totalsmems,
                    &chosen_strand, read_ctr, read_len,
                    loc_f4, loc_rc4);

#ifdef DEBUG_ROSASEED_INPROCESS_MODE
#  ifdef ENABLE_PHASE_II
    fprintf(stderr,
        "[RS_MODE] rid=%d mode=PhaseI+PhaseII phaseI=%lu final=%lu "
        "added=%lu\n", rid,
        (unsigned long)numberofSMEMs,
        (unsigned long)totalsmems,
        (unsigned long)(totalsmems - numberofSMEMs));
#  else
    fprintf(stderr,
        "[RS_MODE] rid=%d mode=PhaseI+GapFill phaseI=%lu final=%lu "
        "added=%lu\n", rid,
        (unsigned long)numberofSMEMs,
        (unsigned long)totalsmems,
        (unsigned long)(totalsmems - numberofSMEMs));
#  endif
#endif

    return materialise_hits(matchArray, totalsmems,
                            rbeg_out, seed_offsets, seed_counts,
                            hits, max_hits, rid);
}

/* -----------------------------------------------------------------------
   rosaseed_core_seed_batch_interleaved()
   -----------------------------------------------------------------------
   Processes up to RS_BATCH reads simultaneously using a coroutine-style
   Phase I interleaver.
   --------------------------------------------- */
int64_t rosaseed_core_seed_batch_interleaved(
    const char **seqs,          /* array of nreads C-strings            */
    const int   *read_lens,     /* length of each read                  */
    const int   *rids,          /* global rid for each read             */
    int          nreads,        /* number of reads in this call         */
    rosaseed_core_hit_t *hits,  /* output: all hits concatenated        */
    int64_t      max_hits,      /* capacity of hits[]                   */
    int64_t     *hit_offsets)   /* out: hits[hit_offsets[i]] = read i   */
{
    if (!g_rosaseed_core_initialized) {
        fprintf(stderr,
            "[ERROR] RosaSeed core not initialized\n");
        return -1;
    }
    if (nreads <= 0) return 0;

    ensure_batch_match_arrays();

    int max_rlen = 0;
    for (int i = 0; i < nreads; ++i)
        if (read_lens[i] > max_rlen) max_rlen = read_lens[i];

    /* Guard: RS_BATCH_MAX_READ_LEN must cover max_rlen */
    if (max_rlen > RS_BATCH_MAX_READ_LEN) {
        fprintf(stderr,
            "[ERROR] read_len=%d exceeds RS_BATCH_MAX_READ_LEN=%d; "
            "rebuild with larger RS_BATCH_MAX_READ_LEN\n",
            max_rlen, RS_BATCH_MAX_READ_LEN);
        return -1;
    }

    rosaseed_core_ensure_thread_buffers(max_hits, max_rlen);

    int64_t  *rbeg_out     = tl_rbeg_out;
    uint32_t *seed_offsets = tl_seed_offsets;
    uint32_t *seed_counts  = tl_seed_counts;

    int64_t total_hits_out = 0;
    int     base           = 0;   /* index of first read in current micro-batch */

    while (base < nreads) {

        int batch_n = nreads - base;
        if (batch_n > RS_BATCH) batch_n = RS_BATCH;

        PhaseI_SlotState slots[RS_BATCH];

        for (int s = 0; s < batch_n; ++s) {
            int idx  = base + s;
            int rlen = read_lens[idx];
            slots[s].read_len   = rlen;
            slots[s].global_rid = rids[idx];
            slots[s].active     = 1;
            rosaseed_seq_to_b4_arrays(seqs[idx], rlen,
                                      slots[s].f4, slots[s].rc4);
        }

        /* ============================================================
           COROUTINE PHASE I
           ============================================================ */
        const uint8_t  *f4_ptrs [RS_BATCH];
        const uint8_t  *rc4_ptrs[RS_BATCH];
        int             rlens_b [RS_BATCH];
        int             rctrs_b [RS_BATCH];
        SMEM           *ma_ptrs [RS_BATCH];
        uint64_t        nsmems_b[RS_BATCH];
        uint64_t       *nsmems_ptrs[RS_BATCH];
        int             chosen_strands[RS_BATCH];

        for (int s = 0; s < batch_n; ++s) {
            f4_ptrs [s]    = slots[s].f4;
            rc4_ptrs[s]    = slots[s].rc4;
            rlens_b [s]    = slots[s].read_len;
            rctrs_b [s]    = slots[s].global_rid;
            ma_ptrs [s]    = tl_batch_matchArrays[s];
            nsmems_b[s]    = 0;
            nsmems_ptrs[s] = &nsmems_b[s];
            chosen_strands[s] = 0;
        }

        uint64_t t0_phase1 = __rdtsc();
        phaseI_batch_interleaved(
            f4_ptrs, rc4_ptrs,
            rlens_b, rctrs_b,
            ma_ptrs, nsmems_ptrs,
            chosen_strands,
            batch_n);
        rs_t_phase1.fetch_add(__rdtsc() - t0_phase1,
                              std::memory_order_relaxed);

        /* ============================================================
           POST-PHASE-I + MATERIALISE: serial per-read 
           Before running Phase II for read s, prefetch the cp_occ
           blocks for reads s+1, s+2, s+3's Phase II first pivots.
           Only active when ENABLE_PHASE_II is defined.
           ============================================================ */
        for (int s = 0; s < batch_n; ++s) {
            int     idx    = base + s;
            int     rlen   = slots[s].read_len;
            int     rid    = slots[s].global_rid;
            SMEM   *ma     = tl_batch_matchArrays[s];
            uint64_t numberofSMEMs = nsmems_b[s];
            int      chosen_strand = chosen_strands[s];
            int      read_ctr      = rid;

#ifdef ENABLE_PHASE_II
            {
                const int pf_dist = 3;
                const int num_piv = g_num_pivots_B;
                for (int ahead = 1; ahead <= pf_dist && s + ahead < batch_n; ++ahead) {
                    const uint8_t *pf_pat =
                        (chosen_strands[s + ahead] == 0)
                        ? slots[s + ahead].f4
                        : slots[s + ahead].rc4;
                    int pf_rlen = slots[s + ahead].read_len;
                    (void)num_piv;

                    /* pivots[0] = read_len_bases - 1 (rightmost base) */
                    int pf_pivot = pf_rlen - 1;

                    if (pf_pivot >= JT_LEN_NT - 1) {
                        uint64_t pf_addr = 0;
                        for (int k = 0; k < JT_LEN_NT; ++k)
                            pf_addr |= ((uint64_t)pf_pat[pf_pivot - k]) << (2 * k);

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
                }
            }
#endif /* ENABLE_PHASE_II */

            /* ------ Post-Phase-I ------ */
            uint64_t totalsmems = 0;
            run_post_phaseI(ma, numberofSMEMs, &totalsmems,
                            &chosen_strand, read_ctr, rlen,
                            slots[s].f4, slots[s].rc4);

            /* ------ Materialise ------ */
            if (hit_offsets) hit_offsets[idx] = total_hits_out;

            int64_t nh = materialise_hits(
                ma, totalsmems,
                rbeg_out, seed_offsets, seed_counts,
                hits + total_hits_out,
                max_hits  - total_hits_out,
                rid);

            if (nh < 0) return -1;
            total_hits_out += nh;
        }

        base += batch_n;
    }

    return total_hits_out;
}

/* -----------------------------------------------------------------------
   rosaseed_core_init()
   ----------------------------------------------------------------------- */
int rosaseed_core_init(const rosaseed_core_config_t *cfg)
{
    if (g_rosaseed_core_initialized) return 0;

    if (cfg == NULL || cfg->index_dir == NULL) {
        fprintf(stderr,
            "[ERROR] RosaSeed core init requires index_dir\n");
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

    snprintf(path_ref16_bin, sizeof(path_ref16_bin),
             "%s/ref16_packed.bin", cfg->index_dir);
    snprintf(path_cvec, sizeof(path_cvec),
             "%s/c_vector.txt",
             cfg->index_dir);

    int cf = 1 << SA_COMPRESSION_FACTOR_POWER;
    snprintf(path_sa_bin,  sizeof(path_sa_bin),
             "%s/sa_ls_word_cf%d.bin",  cfg->index_dir, cf);
    snprintf(path_msb_bin, sizeof(path_msb_bin),
             "%s/sa_ms_byte_cf%d.bin",  cfg->index_dir, cf);
    snprintf(path_occ_bin, sizeof(path_occ_bin),
             "%s/cp_occ_full.bin",       cfg->index_dir);

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

    REFERENCE_GENOME_FILE = path_ref16_bin;
    C_VEC_FILE            = path_cvec;
    JUMP_TABLE_FILE       = path_jt;
    SA_LS_WORD_BIN_FILE   = path_sa_bin;
    SA_MSB_BIT_BIN_FILE   = path_msb_bin;
    OCC_BIN_FILE          = path_occ_bin;

    fprintf(stderr, "[INFO] RosaSeed core init\n");
    fprintf(stderr, "[INFO] RosaSeed index dir: %s\n", cfg->index_dir);
    fprintf(stderr,
        "[RS_DEBUG] SA_COMPRESSION_FACTOR_POWER=%d CF=%d\n",
        SA_COMPRESSION_FACTOR_POWER,
        1 << SA_COMPRESSION_FACTOR_POWER);

    init_char2base16();
    rosaseed_read_index_metadata(cfg->index_dir);  
    allocate_memory();
    load_bwt_data_structures();
    init_two_step_mapinfo(&g_rosaseed_mapinfo, rosaseed_L);

    g_rosaseed_core_initialized = 1;

    fprintf(stderr, "Phase 1 cap : %d\n", g_phase1_cap);
    fprintf(stderr, "[INFO] RosaSeed core index loaded successfully\n");
    return 0;
}

/* -----------------------------------------------------------------------
   rosaseed_core_destroy()
   ----------------------------------------------------------------------- */
void rosaseed_core_destroy(void)
{
    if (!g_rosaseed_core_initialized) return;

    /* free per-slot match arrays if allocated */
    if (tl_batch_bufs_ready) {
        for (int s = 0; s < RS_BATCH; ++s) {
            if (tl_batch_matchArrays[s]) {
                _mm_free(tl_batch_matchArrays[s]);
                tl_batch_matchArrays[s] = NULL;
            }
        }
        tl_batch_bufs_ready = 0;
    }

    free_memory();
    g_rosaseed_core_initialized = 0;
}

/* -----------------------------------------------------------------------
   rosaseed_core_print_timing()
   ----------------------------------------------------------------------- */
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
