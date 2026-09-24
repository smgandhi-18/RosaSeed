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

#include "file_dec.h"
#include "mm_malloc.h"
#include "macros.h"

extern void unmap_all_index_regions(void);

FILE *reads_file_fp = NULL;
FILE *bwt_to_ref_fp = NULL;
FILE *jump_table_fp = NULL;
FILE *ref_file_fp = NULL;
FILE *occ_mat_fp = NULL;
FILE *read_output_fp = NULL;

uint32_t *sa_ls_word = NULL;
uint8_t  *sa_ms_byte = NULL;
uint64_t  *c_vec = NULL;
uint64_t *dummy = NULL;
uint64_t **occ_mat = NULL;
uint64_t *one_hot_mask_array = NULL;
cp_occ32_t *cp_occ = NULL;
uint32_t   *mask32 = NULL;
uint64_t BWT_SIZE_REFERENCE_SIZE = 0;   // set by rosaseed_compact_read_index_metadata()
uint64_t rosaseed_L = 0;
uint64_t *jump_pointers = NULL;

#if PACK_REF4
uint8_t *ref4_packed = NULL;
#else
uint8_t *ref_4 = NULL;
#endif

static inline double bytes_to_mib(uint64_t b)
{
    return (double)b / (1024.0 * 1024.0);
}

void allocate_memory()
{
    /* Jump table */
    if (posix_memalign((void**)&jump_pointers, 64,
                        JUMP_TABLE_ENTRIES * sizeof(uint64_t)) != 0) {
        fprintf(stderr, "Memory alignment failed for jump_pointers\n");
        exit(EXIT_FAILURE);
    }

    /* SA arrays */
    const size_t sa_n = (size_t)SA_SAMPLED_SIZE;
    if (posix_memalign((void**)&sa_ls_word, 64,
                        sa_n * sizeof(uint32_t)) != 0) {
        fprintf(stderr, "Memory alignment failed for sa_ls_word\n");
        exit(EXIT_FAILURE);
    }
    if (posix_memalign((void**)&sa_ms_byte, 64,
                        sa_n * sizeof(uint8_t)) != 0) {
        fprintf(stderr, "Memory alignment failed for sa_ms_byte\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "SA arrays (CF=%d):\n", 1 << SA_COMPRESSION_FACTOR_POWER);
    fprintf(stderr, "  sa_ls_word: %.3f GB\n",
            (sa_n * sizeof(uint32_t)) / 1073741824.0);
    fprintf(stderr, "  sa_ms_byte: %.3f GB\n",
            (sa_n * sizeof(uint8_t)) / 1073741824.0);
    fprintf(stderr, "  TOTAL SA  : %.3f GB\n\n",
            (sa_n * (sizeof(uint32_t) + sizeof(uint8_t))) / 1073741824.0);

    /* Reference storage */
#if PACK_REF4
    {
        size_t packed_bytes = ((size_t)BWT_SIZE_REFERENCE_SIZE + 3u) >> 2;
        if (posix_memalign((void**)&ref4_packed, 64, packed_bytes) != 0) {
            fprintf(stderr, "Memory alignment failed for ref4_packed\n");
            exit(EXIT_FAILURE);
        }
    }
#else
    ref_4 = (uint8_t *)_mm_malloc(
                (size_t)BWT_SIZE_REFERENCE_SIZE * sizeof(uint8_t), 64);
    if (!ref_4) {
        fprintf(stderr, "Memory allocation failed for ref_4\n");
        exit(EXIT_FAILURE);
    }
    /* NO memset */
#endif

    /* Dummy */
    dummy = (uint64_t *)malloc(4 * sizeof(uint64_t));
    if (!dummy) {
        fprintf(stderr, "Memory allocation failed for dummy\n");
        exit(EXIT_FAILURE);
    }

    /* C-vector */
    if (posix_memalign((void**)&c_vec, 64,
                        (ALPHABET_SIZE + 1) * sizeof(uint64_t)) != 0) {
        fprintf(stderr, "Memory alignment failed for c_vec\n");
        exit(EXIT_FAILURE);
    }

    /* cp_occ */
    const uint64_t num_blocks = BWT_SIZE_REFERENCE_SIZE / OCC_INTERVAL + 1;
    cp_occ = (cp_occ32_t*)_mm_malloc(num_blocks * sizeof(cp_occ32_t), 64);
    if (!cp_occ) { perror("cp_occ alloc"); exit(EXIT_FAILURE); }
    /* NO memset */

    /* Prefix masks, small, stays _mm_malloc'd */
    mask32 = (uint32_t*)_mm_malloc(33 * sizeof(uint32_t), 64);
    if (!mask32) { perror("mask32 alloc"); exit(EXIT_FAILURE); }
    build_masks32();

    print_fm_memory_report();
}

void free_memory()
{
    /* sa_ls_word, sa_ms_byte, jump_pointers, ref4_packed/ref_4, cp_occ
       are all mmap regions after load_bwt_data_structures() */
    unmap_all_index_regions();
    sa_ls_word    = NULL;
    sa_ms_byte    = NULL;
    jump_pointers = NULL;
    cp_occ        = NULL;
#if PACK_REF4
    ref4_packed   = NULL;
#else
    ref_4         = NULL;
#endif

    /* These three are still normal allocations */
    free(dummy);   dummy = NULL;
    free(c_vec);   c_vec = NULL;
    if (mask32) { _mm_free(mask32); mask32 = NULL; }
}

void print_fm_memory_report(void)
{
    const uint64_t N = (uint64_t)BWT_SIZE_REFERENCE_SIZE;
    fprintf(stderr, "\n==== 1-step FM-index / seeding memory report ====\n");
    {
        uint64_t bytes = (uint64_t)JUMP_TABLE_ENTRIES * sizeof(uint64_t);
        fprintf(stderr, "jump_pointers: %llu * %zu = %.2f MiB\n",
                (unsigned long long)JUMP_TABLE_ENTRIES,
                sizeof(uint64_t), bytes_to_mib(bytes));
    }
    {
        const uint64_t sa_n = (uint64_t)SA_SAMPLED_SIZE;
        uint64_t bytes = sa_n * sizeof(*sa_ls_word) +
                         sa_n * sizeof(*sa_ms_byte);
        fprintf(stderr, "SA split arrays (CF=%d): %llu sampled * (%zu+%zu) = %.2f MiB\n",
                1 << SA_COMPRESSION_FACTOR_POWER,
                (unsigned long long)sa_n,
                sizeof(*sa_ls_word), sizeof(*sa_ms_byte),
                bytes_to_mib(bytes));
    }
#if PACK_REF4
    {
        uint64_t packed_bytes = (N + 3u) >> 2;
        fprintf(stderr, "ref4_packed: ceil(%llu/4) = %llu bytes = %.2f MiB\n",
                (unsigned long long)N,
                (unsigned long long)packed_bytes,
                bytes_to_mib(packed_bytes));
    }
#else
    {
        uint64_t bytes = N * sizeof(uint8_t);
        fprintf(stderr, "ref_4: %llu * %zu = %.2f MiB\n",
                (unsigned long long)N, sizeof(uint8_t), bytes_to_mib(bytes));
    }
#endif
    {
        uint64_t bytes = (uint64_t)(ALPHABET_SIZE + 1) * sizeof(uint64_t);
        fprintf(stderr, "c_vec: %d * %zu = %.6f MiB\n",
                ALPHABET_SIZE + 1, sizeof(uint64_t), bytes_to_mib(bytes));
    }
    {
        uint64_t bytes = 4ull * sizeof(uint64_t);
        fprintf(stderr, "dummy: 4 * %zu = %.6f MiB\n",
                sizeof(uint64_t), bytes_to_mib(bytes));
    }
    {
        uint64_t num_blocks =
            (N + (uint64_t)OCC_INTERVAL - 1) / (uint64_t)OCC_INTERVAL;
        uint64_t bytes = num_blocks * sizeof(cp_occ32_t);
        fprintf(stderr, "cp_occ: %llu blocks * %zu = %.2f MiB\n",
                (unsigned long long)num_blocks,
                sizeof(cp_occ32_t), bytes_to_mib(bytes));
    }
    {
        uint64_t bytes = 33ull * sizeof(uint32_t);
        fprintf(stderr, "mask32: 33 * %zu = %.6f MiB\n",
                sizeof(uint32_t), bytes_to_mib(bytes));
    }
    fprintf(stderr, "=================================================\n\n");
}
