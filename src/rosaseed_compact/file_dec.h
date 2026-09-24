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

#ifndef FILE_DECLARATIONS_H
#define FILE_DECLARATIONS_H

// #define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>
#include <unistd.h>

extern const char *REFERENCE_BIN_FILE;
extern const char *C_VEC_FILE;
extern const char *BWT_FILE;
extern const char *JUMP_TABLE_FILE;
extern const char *READS_FILE;
extern int         NUM_OF_READS_IN_FILE;

// Binary index files
extern const char *SA_LS_WORD_BIN_FILE;   // sa_ls_word_cfX.bin
extern const char *SA_MSB_BIT_BIN_FILE;   // sa_ms_byte_cfX.bin
extern const char *OCC_BIN_FILE;          // cp_occ.bin

extern uint64_t BWT_SIZE_REFERENCE_SIZE;   // set at startup from the cp_occ_compact.bin header

/* Index format limits (same as 2-step): SA values use 33 bits, so the BWT
   length N = 2L + 1 must fit in 33 bits, i.e. the genome must be below 2^32 bp. */
#define RS_MAX_BWT_LEN   8589934592ULL   /* 2^33 */
#define RS_MAX_GENOME_BP 4294967296ULL   /* 2^32 */

// ==============================
//  JUMP TABLE CONFIGURATION 
// ==============================

// #define LOAD_JTABLE_16nt
#ifdef LOAD_JTABLE_14nt
    #define JUMP_TABLE_ENTRIES 268435456ULL
#endif    
#ifdef LOAD_JTABLE_15nt
    #define JUMP_TABLE_ENTRIES 1073741824ULL 
#endif
#ifdef LOAD_JTABLE_16nt
    #define JUMP_TABLE_ENTRIES 4294967296ULL 
#endif

// Real builds must pass -DLOAD_JTABLE_Xnt; a runtime check catches a missing flag.
#ifndef JUMP_TABLE_ENTRIES
    #define JUMP_TABLE_ENTRIES 268435456ULL 
#endif

extern uint64_t *jump_pointers;  // Dynamic allocation

// ==============================
// FILE POINTERS 
// ==============================
 extern FILE *reads_file_fp;
 extern FILE *bwt_to_ref_fp;
 extern FILE *jump_table_fp;
 extern FILE *ref_file_fp;
 extern FILE *org_bwt_to_ref_fp;
 extern FILE *c_vec_fp;
 extern FILE *bwt_fp;
 extern FILE *pivot_log_fp;
 extern FILE *read_output_fp; 


// ==============================
// DYNAMIC MEMORY ALLOCATION 
// ==============================
extern uint64_t *c_vec;
extern uint64_t *dummy;
extern uint64_t **occ_mat;
extern uint32_t *sa_ls_word;   // sampled SA lower 32 bits
extern uint8_t  *sa_ms_byte;   // sampled SA bit 32 (MSB)
#define sa_msb_bit sa_ms_byte  // backward compat alias

// ==============================
//  FUNCTION PROTOTYPES 
// ==============================
void allocate_memory();   // Dynamically allocate memory
void free_memory();       // Free allocated memory
void print_fm_memory_report(void);

// ===== One-step / base-4 FM-index config =====
// Alphabet: A=0, C=1, G=2, T=3. Terminator '$'/0x00 is not counted.
#define ALPHABET_SIZE 4      // symbols: A,C,G,T
#define OCC_INTERVAL  32

// One checkpoint block (32 bytes, 0.5 cache lines):
//   - 4 * uint32_t cumulative counts  (16 B)
//   - 4 * uint32_t one-hot bitmaps    (16 B)
typedef struct {
    uint32_t cp_count[ALPHABET_SIZE];  // cumulative counts at row k*OCC_INTERVAL - cumulative A/C/G/T counts before block
    uint32_t one_hot[ALPHABET_SIZE];   // bit i set if BWT char at (k*OCC_INTERVAL + i) == base
} cp_occ32_t;

extern cp_occ32_t *cp_occ;     // [num_blocks]
extern uint32_t   *mask32;     // [33] prefix masks for 32-bit windows

// Fast 32-bit popcount
#if defined(__clang__) || defined(__GNUC__)
  #define POPCOUNT32 __builtin_popcount
#else
  #include <immintrin.h>
  static inline int POPCOUNT32(uint32_t x){ return _mm_popcnt_u32(x); }
#endif

// Prefix masks: mask32[o] has the top o bits set (MSB-first). mask32[0] == 0.
static inline void build_masks32(void) {
    for (int i = 0; i <= 32; ++i) {
        mask32[i] = (i == 0) ? 0u : (~0u << (32 - i));
    }
}

// Occ query at absolute position `pos` (0..N), symbol `sym` (0..15).
// out = cp_count(sym, block) + popcount( one_hot(sym) & prefix_mask(offset) )
#define GET_OCC32(pos, sym, out) do {                             \
    uint64_t _blk = (pos) >> 5;           /* /32 */               \
    uint32_t _off = (pos) & 31;           /* %32 */               \
    uint32_t _acc = cp_occ[_blk].cp_count[(sym)];                 \
    uint32_t prefix_mask = _off ? (~0u << (32 - _off)) : 0u; \
    uint32_t _bits= cp_occ[_blk].one_hot[(sym)] & prefix_mask;   \
    (out) = (uint64_t)_acc + (uint64_t)POPCOUNT32(_bits);         \
} while(0)
    // if (_blk == TARGET_BLK) {                                        \
    //     printf("\n======================================================\n"); \
    //     printf("--- OCC DEBUGGER HIT (Sym %d) ---\n", (int)(sym));   \
    //     printf("Target Bucket (Block): %lu\n", (unsigned long)_blk); \
    //     printf("Absolute Position (pos): %lu\n", (unsigned long)(pos)); \
    //     printf("Offset (k): %u\n", _off);                            \
    //     printf("------------------------------------------------------\n"); \
    //     printf("1. Checkpoint Count (_acc): %u\n", _acc);             \
    //     printf("2. One-Hot Bitmask (one_hot[sym]): 0x%08x\n", cp_occ[_blk].one_hot[(sym)]); \
    //     printf("3. Prefix Mask (mask32[_off]): 0x%08x\n", mask32[_off]); \
    //     printf("4. Masked Bits (_bits) [AND Result]: 0x%08x\n", _bits); \
    //     printf("5. Prefix Popcount (Counted Bits): %u\n", POPCOUNT32(_bits)); \
    //     printf("6. FINAL Occ Value (out): %lu\n", (unsigned long)(out)); \
    //     printf("======================================================\n"); \
    // }                                                                \
    
typedef long int heap_data_t;

// --- Packed base4 reference for 1-step RosaSeed ---
// Four DNA bases per byte:
//
//   bits 7-6 : base i
//   bits 5-4 : base i+1
//   bits 3-2 : base i+2
//   bits 1-0 : base i+3
//
// Encoding:
//   A = 0 = 00
//   C = 1 = 01
//   G = 2 = 10
//   T = 3 = 11
//
// This matches the 1-step jump-table/base4 convention.

#define PACK_REF4 1

#if PACK_REF4

extern uint8_t *ref4_packed;

static inline __attribute__((always_inline))
uint8_t ref4_get(const uint8_t *buf, uint64_t i)
{
    uint8_t x = buf[i >> 2];                 // 4 bases per byte
    uint8_t shift = (uint8_t)(6 - 2 * (i & 3));
    return (uint8_t)((x >> shift) & 0x03);
}

static inline __attribute__((always_inline))
void ref4_set(uint8_t *buf, uint64_t i, uint8_t v)
{
    uint64_t b = i >> 2;
    uint8_t shift = (uint8_t)(6 - 2 * (i & 3));
    uint8_t mask = (uint8_t)(0x03u << shift);

    buf[b] = (uint8_t)((buf[b] & ~mask) | ((v & 0x03u) << shift));
}

#define REF_AT(pos) ref4_get(ref4_packed, (uint64_t)(pos))

#else

extern uint8_t *ref_4;

static inline __attribute__((always_inline))
uint8_t REF_AT(uint64_t pos)
{
    return ref_4[pos];
}

#endif


#endif // FILE_DECLERATIONS_H