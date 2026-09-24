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
#ifndef _MACROS_H
#define _MACROS_H

#include <stdint.h>
#include <stdbool.h>

// #define ENABLE_F_RC_CHOICE
#define MATCH_ARRAY_CAPACITY  512

#ifndef PHASE1_MAX_INTERVAL
    #define PHASE1_MAX_INTERVAL  5000   // default: primary paper config
#endif
// ==============================
//  CONFIGURABLE PARAMETERS 
// ==============================

static const int LINE_LENGTH_FASTQ = 1000;
static const int READ_LENGTH_FASTQ = 4;
static const int MINIMUM_SEEDLENGTH    = 19;   // compile-time default only
static const int MINIMUM_SEEDLENGTH_P3 = 20;   // compile-time default only

// Runtime-configurable equivalents (set from argv; g_min_seed_len_BC is always
// g_min_seed_len_A + 1, one user-facing parameter, the other is derived).
extern int g_min_seed_len_A;    // Phase A  minimum seed length  (default 19)
extern int g_min_seed_len_BC;   // Phase B/C minimum seed length (default 20)

static const int64_t min_inv_threshold = 20;   // compile-time default only

// Runtime-configurable Phase B/C specificity threshold and Phase A interval cap.
extern int64_t g_min_intv;     // Phase B/C early-stop / emission threshold  (default 20)
extern int     g_phase1_cap;   // Phase A interval cap (0 = no cap)          (default 5000)
extern int     g_num_pivots_B; // Phase B number of fixed pivots             (default 3)
extern int     g_sa_max_occ;   // max SA hits per seed for materialization / sampling

#ifndef READ_LEN
#  define READ_LEN 150
#endif
static const int READ_LENGTH_BASE4 = READ_LEN;

/*
 * One-step RosaSeed reference model:
 *
 *   indexed text = forward_reference + reverse_complement_reference + terminator
 *
 * L is the original forward reference length.
 * BWA/RosaSeed coordinate space is approximately 2*L plus the terminator.
 */
extern uint64_t rosaseed_L;   // = (BWT_SIZE_REFERENCE_SIZE - 1) / 2, set from the index header


// ==============================
//  DEBUGGING MACROS
// ==============================

#ifdef DEBUG_MODE
#define DEBUG_PRINTF(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINTF(fmt, ...) do {} while (0)
#endif

#ifdef DEBUG_MODE_P2
#define DEBUG_PRINTF_P2(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINTF_P2(fmt, ...) do {} while (0)
#endif

// ==============================
//  RETURN CODES
// ==============================

#define NOT_FOUND 0
#define FOUND     1

typedef struct smem_struct
{
    uint64_t low_ptr_read;        // FM interval lower bound for forward active strand
    uint64_t low_ptr_rc_read;     // FM interval lower bound for RC active strand

    uint64_t unique_ref_pointer;  // raw SA/reference start position for unique seed

    int start_idx;                // forward read coordinate
    int end_idx;                  // forward read coordinate
    int smem_score;
    int rid;

    bool unique_flag;

    uint8_t seed_strand;          // 0 = forward read, 1 = RC read

    /*
     * 1-step does not need two-step symbol offset.
     * Kept as padding/compat field for now; should remain 0.
     */
    uint8_t unique_leftmost_base_off;

} SMEM;

// ==============================
//  FEATURE FLAGS
// ==============================
// #define PRINT_SMEM

#define MAX_PROF_ENTRIES 256
#define MAX_PATTERN_SIZE 250


// ==============================
//  JUMP TABLE SIZE
//  Set via compiler flag:
//    -DLOAD_JTABLE_14nt
//    -DLOAD_JTABLE_15nt
//    -DLOAD_JTABLE_16nt
// ==============================

#if defined(LOAD_JTABLE_15nt)
    #define JT_LEN_NT 15
#elif defined(LOAD_JTABLE_16nt)
    #define JT_LEN_NT 16
#else
    #define JT_LEN_NT 14
#endif

#define KMER_BASES JT_LEN_NT

// ==============================
//  SUFFIX ARRAY COMPRESSION
//  Set via -DSA_COMPRESSION_FACTOR_POWER=N
//  0 = CF=1
//  1 = CF=2
//  2 = CF=4
//  3 = CF=8
// ==============================

#ifndef SA_COMPRESSION_FACTOR_POWER
#define SA_COMPRESSION_FACTOR_POWER 0
#endif

#define SA_INDEX_AND_SAMPLE ((1u << SA_COMPRESSION_FACTOR_POWER) - 1)

#define SA_SAMPLED_SIZE \
    ((BWT_SIZE_REFERENCE_SIZE + SA_INDEX_AND_SAMPLE) >> SA_COMPRESSION_FACTOR_POWER)

#endif // _MACROS_H