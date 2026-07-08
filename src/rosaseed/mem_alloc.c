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

#include "file_dec.h"
#include "mm_malloc.h"
#include "macros.h"

FILE *reads_file_fp = NULL;
FILE *bwt_to_ref_fp = NULL;
FILE *jump_table_fp = NULL;
FILE *ref_file_fp = NULL;
FILE *org_bwt_to_ref_fp = NULL;
FILE *c_vec_fp = NULL;
FILE *occ_mat_fp = NULL;
FILE *bwt_fp = NULL;
FILE *read_output_fp = NULL;
FILE *jump_masks_fp = NULL;

uint32_t *sa_ls_word = NULL;
uint8_t  *sa_ms_byte = NULL;   // was sa_msb_bit
uint64_t *c_vec = NULL;
uint64_t *dummy = NULL;
uint64_t **occ_mat = NULL;
uint64_t *one_hot_mask_array = NULL;
cp_occ32_t *cp_occ   = NULL;
uint32_t   *mask32   = NULL;
uint64_t *jump_pointers = NULL;

#if PACK_REF16
uint8_t  *ref16_packed = NULL;
#else
base16_t *reference_genome = NULL;
#endif

uint64_t BWT_SIZE_REFERENCE_SIZE = 0;   // initialized by rosaseed_read_index_metadata()
uint64_t rosaseed_L = 0;

extern void unmap_all_index_regions(void);

void allocate_memory() {

    if (posix_memalign((void**)&jump_pointers, 64, JUMP_TABLE_ENTRIES * sizeof(uint64_t)) != 0) {
        fprintf(stderr, "Memory alignment failed for jump_pointers\n");
        exit(EXIT_FAILURE);
    }

    // SA arrays — size depends on compression factor
    // CF=1: full N entries.  CF=2: N/2.  CF=4: N/4.  CF=8: N/8.
    const size_t sa_n = (size_t)SA_SAMPLED_SIZE;
    int rc=0;

    rc = posix_memalign((void**)&sa_ls_word, 64, sa_n * sizeof(uint32_t));
    if (rc) { perror("posix_memalign sa_ls_word"); exit(EXIT_FAILURE); }

    rc = posix_memalign((void**)&sa_ms_byte, 64, sa_n * sizeof(uint8_t));
    if (rc) { perror("posix_memalign sa_ms_byte"); exit(EXIT_FAILURE); }

    {
        size_t sa_ls_bytes = sa_n * sizeof(uint32_t);
        size_t sa_ms_bytes = sa_n * sizeof(uint8_t);
        fprintf(stderr, "SA arrays (CF=%d):\n", 1 << SA_COMPRESSION_FACTOR_POWER);
        fprintf(stderr, "  sa_ls_word: %.3f GB\n", sa_ls_bytes / 1073741824.0);
        fprintf(stderr, "  sa_ms_byte: %.3f GB\n", sa_ms_bytes / 1073741824.0);
        fprintf(stderr, "  TOTAL SA  : %.3f GB\n\n",
                (sa_ls_bytes + sa_ms_bytes) / 1073741824.0);
    }

    #if PACK_REF16
    size_t packed_bytes = (BWT_SIZE_REFERENCE_SIZE + 1) >> 1;  // ceil(N/2)
    if (posix_memalign((void**)&ref16_packed, 64, packed_bytes) != 0) {
        fprintf(stderr, "Memory alignment failed for ref16_packed\n");
        exit(EXIT_FAILURE);
    }
    #else
        reference_genome = (base16_t *)_mm_malloc(BWT_SIZE_REFERENCE_SIZE * sizeof(base16_t), 64);
    #endif

    dummy = (uint64_t *)malloc(4 * sizeof(uint64_t *));

    // Allocate aligned memory for c_vec 
    if (posix_memalign((void**)&c_vec, 64, (ALPHABET_SIZE+1) * sizeof(uint64_t)) != 0) {
        fprintf(stderr, "Memory alignment failed for c_vec\n");
        exit(EXIT_FAILURE);
    }

    const uint64_t num_blocks = (BWT_SIZE_REFERENCE_SIZE + OCC_INTERVAL) / OCC_INTERVAL; // ceil((N)/32)
    cp_occ = (cp_occ32_t*)_mm_malloc(num_blocks * sizeof(cp_occ32_t), 64);
    if (!cp_occ) { perror("cp_occ alloc"); exit(1); }

    mask32 = (uint32_t*)_mm_malloc((33) * sizeof(uint32_t), 64);
    if (!mask32) { perror("mask32 alloc"); exit(1); }
    build_masks32();
    
    if (cp_occ == NULL) {
        fprintf(stderr, "Memory allocation failed for cp_occ\n");
        exit(EXIT_FAILURE);
    }

    print_fm_memory_report();

}

void free_memory() {
    /* sa_ls_word, sa_ms_byte, jump_pointers, ref16_packed, cp_occ
       are now mmap regions — must use munmap */
    unmap_all_index_regions();
    sa_ls_word    = NULL;
    sa_ms_byte    = NULL;
    jump_pointers = NULL;
    cp_occ        = NULL;
    #if PACK_REF16
    ref16_packed  = NULL;
    #else
    free(reference_genome);
    reference_genome = NULL;
    #endif

    free(dummy);
    free(c_vec);
    if (mask32) { _mm_free(mask32); mask32 = NULL; }
}

static inline double bytes_to_mib(uint64_t b){
    return (double)b / (1024.0 * 1024.0);
}

void print_fm_memory_report(void){
    const uint64_t N = (uint64_t)BWT_SIZE_REFERENCE_SIZE;

    fprintf(stderr, "\n==== FM-index / seeding memory report (allocated bytes) ====\n");

    /* jump table */
    {
        uint64_t bytes = (uint64_t)JUMP_TABLE_ENTRIES * (uint64_t)sizeof(uint64_t);
        fprintf(stderr, "jump_pointers: %llu * %zu = %.2f MiB\n",
                (unsigned long long)JUMP_TABLE_ENTRIES, sizeof(uint64_t), bytes_to_mib(bytes));
    }

    /* SA / BWT-to-ref */
    {
        // Actual allocated size depends on compression factor
        const uint64_t sa_n   = BWT_SIZE_REFERENCE_SIZE >> SA_COMPRESSION_FACTOR_POWER;
        uint64_t bytes = sa_n * sizeof(*sa_ls_word) + sa_n * sizeof(*sa_ms_byte);
        fprintf(stderr, "SA split arrays (CF=%d): %llu sampled * (%zu+%zu) = %.2f MiB\n",
                1 << SA_COMPRESSION_FACTOR_POWER,
                (unsigned long long)sa_n,
                sizeof(*sa_ls_word), sizeof(*sa_ms_byte),
                bytes_to_mib(bytes));
    }

    /* reference storage (packed vs unpacked) */
#if PACK_REF16
    {
        uint64_t packed_bytes = (N + 1) >> 1;  /* ceil(N/2) bytes */
        fprintf(stderr, "ref16_packed: ceil(%llu/2) = %llu bytes = %.2f MiB\n",
                (unsigned long long)N,
                (unsigned long long)packed_bytes,
                bytes_to_mib(packed_bytes));
    }
#else
    {
        uint64_t bytes = N * (uint64_t)sizeof(base16_t);
        fprintf(stderr, "reference_genome: %llu * %zu = %.2f MiB\n",
                (unsigned long long)N, sizeof(base16_t), bytes_to_mib(bytes));
    }
#endif

    /* c_vec */
    {
        uint64_t bytes = (uint64_t)(ALPHABET_SIZE + 1) * (uint64_t)sizeof(uint64_t);
        fprintf(stderr, "c_vec: %d * %zu = %.6f MiB\n",
                (ALPHABET_SIZE + 1), sizeof(uint64_t), bytes_to_mib(bytes));
    }

    /* dummy */
    {
        uint64_t bytes = 4ull * (uint64_t)sizeof(uint64_t);
        fprintf(stderr, "dummy: 4 * %zu = %.6f MiB\n",
                sizeof(uint64_t), bytes_to_mib(bytes));
    }

    /* occ checkpoints */
    {
        uint64_t num_blocks = (N + (uint64_t)OCC_INTERVAL - 1) / (uint64_t)OCC_INTERVAL;
        uint64_t bytes = num_blocks * (uint64_t)sizeof(cp_occ32_t);
        fprintf(stderr, "cp_occ: %llu blocks * %zu = %.2f MiB\n",
                (unsigned long long)num_blocks, sizeof(cp_occ32_t), bytes_to_mib(bytes));
    }

    /* prefix masks */
    {
        uint64_t bytes = 33ull * (uint64_t)sizeof(uint32_t);
        fprintf(stderr, "mask32: 33 * %zu = %.6f MiB\n",
                sizeof(uint32_t), bytes_to_mib(bytes));
    }

    fprintf(stderr, "===========================================================\n\n");
}