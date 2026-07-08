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

#ifndef _MACROS_H
#define _MACROS_H

#include <stdint.h>
#include <stdbool.h>

#define ENABLE_F_RC_CHOICE
#define MATCH_ARRAY_CAPACITY  512

/********************************************************************************
 * Runtime-configurable equivalents (set from argv; g_min_seed_len_BC is always
 * g_min_seed_len_A + 1 — one user-facing parameter, the other is derived).
 * Runtime-configurable Phase B/C specificity threshold and Phase A interval cap.
 * Both set from argv; compile-time values above are kept as documentation defaults.
 ********************************************************************************/

extern int     g_min_seed_len_A;    // Phase A  minimum seed length  (default 19)
extern int     g_min_seed_len_BC;   // Phase B/C minimum seed length (default 20)
extern int64_t g_min_intv;     // Phase B/C early-stop / emission threshold  (default 20)
extern int     g_phase1_cap;   // Phase A interval cap (0 = no cap)          (default 2000)
extern int     g_num_pivots_B; // Phase B number of fixed pivots             (default 3)
extern int     g_sa_max_occ;   // max SA hits per seed for materialization / sampling

#ifndef READ_LEN
#  define READ_LEN 302
#endif

extern uint64_t rosaseed_L;   // = (BWT_SIZE_REFERENCE_SIZE + 1) / 2

// ==============================
//  DEBUGGING MACROS
// ==============================
#ifdef DEBUG_MODE
#define DEBUG_PRINTF(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINTF(fmt, ...) \
    do                         \
    {                          \
    } while (0) // No-op in release mode
#endif

#ifdef DEBUG_MODE_P2
#define DEBUG_PRINTF_P2(fmt, ...) printf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINTF_P2(fmt, ...) \
    do                            \
    {                             \
    } while (0) // No-op in release mode
#endif

// ==============================
//  RETURN CODES
// ==============================
#define NOT_FOUND 0
#define FOUND 1

// ==============================
//  STRUCT DEFINITIONS
// ==============================
typedef struct smem_struct
{
    uint64_t low_ptr_read;
    uint64_t low_ptr_rc_read;
    uint64_t unique_ref_pointer;
    uint64_t unique_leftmost_twostep_pos;
    int start_idx;
    int end_idx;
    int smem_score;
    int rid;
    bool unique_flag;
    uint8_t seed_strand;    // 0 = forward active strand, 1 = RC active strand
    uint8_t  unique_leftmost_base_off;   // 0 = left base of symbol, 1 = right base
} SMEM;

// ==============================
//  FEATURE FLAGS
// ==============================
#define MAX_PROF_ENTRIES 256
#define MAX_PATTERN_SIZE 250

/**********************************************************
 * JUMP TABLE SIZE — auto-derived from LOAD_JTABLE_Xnt flag
 * Set via compiler -D flag: -DLOAD_JTABLE_14nt / _15nt / _16nt
 **********************************************************/

#if defined(LOAD_JTABLE_15nt)
    #define JT_LEN_NT       15
#elif defined(LOAD_JTABLE_16nt)
    #define JT_LEN_NT       16
#else
    #define JT_LEN_NT       14   // default — also covers fallback in file_dec_new.h
#endif

#define JT_PAIRS        (JT_LEN_NT / 2)
#define JT_TAIL_SINGLE  (JT_LEN_NT & 1)   // 1 for 15nt (odd), 0 for 14/16nt
#define KMER_BASES      JT_LEN_NT

/***************************************
 * SUFFIX ARRAY COMPRESSION
 * Set via -DSA_COMPRESSION_FACTOR_POWER=N
 * 0 = CF=1 (uncompressed, default)
 * 1 = CF=2  (sample every 2nd entry)
 * 2 = CF=4  (sample every 4th entry)
 * 3 = CF=8  (sample every 8th entry)
 ***************************************/
#ifndef SA_COMPRESSION_FACTOR_POWER
    #define SA_COMPRESSION_FACTOR_POWER 0   // default: uncompressed
#endif

// SA_INDEX_AND_SAMPLE: bitmask to test if a row is sampled
// CF=1 → 0 (every row sampled), CF=2 → 1, CF=4 → 3, CF=8 → 7
#define SA_INDEX_AND_SAMPLE  ((1u << SA_COMPRESSION_FACTOR_POWER) - 1)

// Number of sampled SA entries
#define SA_SAMPLED_SIZE  \
    ((BWT_SIZE_REFERENCE_SIZE + SA_INDEX_AND_SAMPLE) >> SA_COMPRESSION_FACTOR_POWER)
#endif // _MACROS_H