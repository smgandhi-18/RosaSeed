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

#include "read_init.h"
#include "macros.h"
#include "file_dec.h"
 
char base16_chars[ALPHABET_SIZE] =
    { '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F' };
 
base16_t char2base16[128];
 
// Reverse-complement mapping for base-16 two-step encoding
// (AA..TT → 0..F, from corrected table)
const uint8_t rc16[16] = {
    /*0*/ 0xF, /*1*/ 0xB, /*2*/ 0x7, /*3*/ 0x3,
    /*4*/ 0xE, /*5*/ 0xA, /*6*/ 0x6, /*7*/ 0x2,
    /*8*/ 0xD, /*9*/ 0x9, /*A*/ 0x5, /*B*/ 0x1,
    /*C*/ 0xC, /*D*/ 0x8, /*E*/ 0x4, /*F*/ 0x0
};
 
void init_char2base16() {
    for (int i = 0; i < 128; i++) char2base16[i] = SYM_DOLLAR; // default
    for (int i = 0; i < 10; i++) char2base16['0'+i] = (base16_t)i;
    for (int i = 0; i < 6;  i++) char2base16['A'+i] = (base16_t)(10+i);
    char2base16['$'] = SYM_DOLLAR;
}
 
extern FILE *reads_file_fp;
 
/* Convert DNA char → base-4 (0..3) */
static inline uint8_t nt2b4(char c)
{
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:
            fprintf(stderr, "Invalid nucleotide '%c' in read.\n", c);
            exit(EXIT_FAILURE);
    }
}
 
/* Reverse complement in base-4 (0..3) */
static inline uint8_t b4_rc(uint8_t b)
{
    /* A(0) ↔ T(3), C(1) ↔ G(2) */
    return (b ^ 0x3);
}
 