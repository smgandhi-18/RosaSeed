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
#ifndef _PROFILE_H
#define _PROFILE_H

// #define _GNU_SOURCE

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
#include "seeding_kernel.h"

void display_stats();

extern uint64_t proc_freq, tprof[MAX_PROF_ENTRIES];

/*** Runtime profiling macros ***/

#define LOAD_INDEX   0
#define SMEMI        1
#define SMEMII       2
#define GAPFILL      3
#define SORT1        4
#define SORT2        5
#define SAL          6

extern uint64_t bwa_single_step_fn_access_count;  // Declare it as extern
extern uint64_t phase2_counter;
extern uint64_t num_calls;

extern uint64_t time_singlestep;
extern uint64_t time_compute_jump;
extern uint64_t time_heap_find;
extern uint64_t time_ref_match;

#endif //_PROFILE_H
