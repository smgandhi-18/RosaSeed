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

#ifndef _SEEDING_KERNEL_H
#define _SEEDING_KERNEL_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>
#include <unistd.h>
#include "macros.h"
#include "file_dec.h"
#include <xmmintrin.h>
#include <limits.h>

#include <stddef.h>

extern __thread int short_read_len;

void phaseI_routine(
   SMEM *matchArray, uint64_t *numberofSMEMs, int *read_counter,
   uint8_t pat_f4[], uint8_t pat_rc4[], int *chosen_strand, int read_len);                           

void phaseII_routine(
   uint8_t shortread_pattern[],  /* base-4 read, chosen strand */
   SMEM *matchArray,
   uint64_t *total_smem,
   int *rid,
   int64_t min_intv,
   int *chosen_strand,
   int read_len);   

void extract_jump_bounds(uint64_t jump_entry, uint64_t* l, uint64_t* h, uint64_t* diff);

#endif // _SEEDING_KERNEL_H
