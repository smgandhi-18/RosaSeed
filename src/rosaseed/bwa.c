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

#include "bwa.h"
#include "file_dec.h"
#include "macros.h"
#include "helper_functions.h"
#include "load_data.h"
#include "profiling.h"
#include <xmmintrin.h>  
#include <immintrin.h>  

extern __thread int short_read_len;

// Right/left base (0..3) from a base16 symbol (0..15)
static inline uint8_t b4_right(uint8_t s16) { return s16 & 0x3; }
static inline uint8_t b4_left (uint8_t s16) { return (s16 >> 2) & 0x3; }

static inline uint64_t addr16_from_pairs_global(const base16_t *pat, int npairs) {
    uint64_t addr = 0;
    for (int p = 0; p < npairs; ++p) {
        addr |= ((uint64_t)pat[short_read_len]) << (4 * p);
        --short_read_len;
    }
    return addr;
}

// Wrapper function for FM-index backward search
void fm_index_mapping_backward_search(uint64_t *low, uint64_t *high, base16_t nuc_in_read) {
    if (*low >= *high) return;
    uint64_t occ_lo, occ_hi;
    GET_OCC32(*low,  nuc_in_read, occ_lo);
    GET_OCC32(*high, nuc_in_read, occ_hi);
    uint64_t base = c_vec[nuc_in_read];
    *low  = base + occ_lo;
    *high = base + occ_hi;   
}
