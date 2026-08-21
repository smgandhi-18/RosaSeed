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
#include <unistd.h>

extern const char *REFERENCE_GENOME_FILE;
extern const char *C_VEC_FILE;
extern const char *BWT_FILE;
extern const char *JUMP_TABLE_FILE;
extern const char *READS_FILE;
extern int         NUM_OF_READS_IN_FILE;

extern const char *SA_LS_WORD_BIN_FILE;   // sa_ls_word.bin
extern const char *SA_MSB_BIT_BIN_FILE;   // sa_msb_bit.bin
extern const char *OCC_BIN_FILE;          // checkpoint_occ.bin

// and does not change per run, only per reference build
extern uint64_t BWT_SIZE_REFERENCE_SIZE;   // set at startup from cp_occ header
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

// Real builds must pass -DLOAD_JTABLE_Xnt — runtime check catches missing flag.
#ifndef JUMP_TABLE_ENTRIES
    #define JUMP_TABLE_ENTRIES 268435456ULL 
#endif

extern uint64_t *jump_pointers;  

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
// ENUM DEFINITIONS 
// ==============================

typedef enum {
    SYM_0=0, SYM_1, SYM_2, SYM_3,
    SYM_4, SYM_5, SYM_6, SYM_7,
    SYM_8, SYM_9, SYM_A, SYM_B,
    SYM_C, SYM_D, SYM_E, SYM_F,
    SYM_DOLLAR = 16
} base16_t;

// ==============================
//  FUNCTION PROTOTYPES 
// ==============================
void allocate_memory();   // Dynamically allocate memory
void free_memory();       // Free allocated memory
void print_fm_memory_report(void);

// ===== Two-step / base-16 FM-index config =====
#define ALPHABET_SIZE 16      // symbols: '0'..'9','A'..'F' (no '$' in the 16)
#define OCC_INTERVAL   32     // one checkpoint + 32-bit bitmaps fits 2 cache lines
// #define TARGET_BLK    108307284

// One checkpoint block (128 bytes, 2 cache lines):
//   - 16 * uint32_t cumulative counts  (64 B)
//   - 16 * uint32_t one-hot bitmaps    (64 B)
typedef struct __attribute__((aligned(64))) {
    uint32_t cp_count[ALPHABET_SIZE];  // cumulative counts at row k*OCC_INTERVAL
    uint32_t one_hot[ALPHABET_SIZE];   // bit i set if BWT char at (k*OCC_INTERVAL + i) == symbol
} cp_occ32_t;

extern cp_occ32_t *cp_occ;     // [num_blocks]
extern uint32_t   *mask32;     // [33] prefix masks for 32-bit windows

typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t occ_interval;
    uint32_t alphabet_size;
    uint32_t reserved;
    uint64_t bwt_len_non_dollar;
    uint64_t bwt_len_total_with_dollar;
    uint64_t num_blocks;
    int64_t sentinel_index;
} rs_occ_full_header_t;

// Fast 32-bit popcount
#if defined(__clang__) || defined(__GNUC__)
  #define POPCOUNT32 __builtin_popcount
#else
  #include <immintrin.h>
  static inline int POPCOUNT32(uint32_t x){ return _mm_popcnt_u32(x); }
#endif

// Map ASCII to base-16 symbol id [0..15]; return -1 for '$' or anything else.
static inline int symbol_index16(int ch) {
    if (ch >= '0' && ch <= '9') return ch - '0';
    ch &= ~0x20; // toupper without locale
    if (ch >= 'A' && ch <= 'F') return 10 + (ch - 'A');
    return -1; // '$' or others
}

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

// --- Packed base16 reference (optional) ---
#define PACK_REF16 1   // set 0 to fall back to old byte-per-symbol array

#if PACK_REF16
// Two base16 symbols per byte: even index -> low nibble, odd index -> high nibble
extern uint8_t *ref16_packed;

static inline __attribute__((always_inline))
uint8_t ref16_get(const uint8_t *buf, uint64_t i) {
    uint8_t x = buf[i >> 1];
    return (i & 1) ? (uint8_t)((x >> 4) & 0x0F) : (uint8_t)(x & 0x0F);
}

static inline __attribute__((always_inline))
void ref16_set(uint8_t *buf, uint64_t i, uint8_t v) {
    uint64_t b = i >> 1;
    if (i & 1) {          // odd -> high nibble
        buf[b] = (uint8_t)((buf[b] & 0x0F) | ((v & 0x0F) << 4));
    } else {              // even -> low nibble
        buf[b] = (uint8_t)((buf[b] & 0xF0) | (v & 0x0F));
    }
}

#define REF_AT(pos) ref16_get(ref16_packed, (uint64_t)(pos))
#else
extern base16_t *reference_genome;
// #define REF_AT(pos) (reference_genome[(pos)])
static inline uint8_t REF_AT(uint64_t pos) {
    return reference_genome[pos];
}
#endif


#endif // FILE_DECLERATIONS_H
