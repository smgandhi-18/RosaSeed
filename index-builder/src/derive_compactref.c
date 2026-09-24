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

/*
 * derive_compactref.c
 *
 * Builds the RosaSeed-Compact indexed text from a clean reference:
 *
 *   input : single-line ACGT file from preprocess_genome (optional trailing newline)
 *   output: forward + reverse_complement(forward) + "\n"   (one line, for gsufsort)
 *
 * gsufsort appends the terminator, so the BWT length is N = 2L + 1.
 *
 * Usage: derive_compactref <clean.txt> <compact_ref.txt>
 */
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define RS_MAX_GENOME_BP 4294967296ULL   /* 2^32: RosaSeed index format limit */
#define OUT_BUF (64u * 1024u * 1024u)

int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <clean.txt> <compact_ref.txt>\n", argv[0]);
        return EXIT_FAILURE;
    }
    FILE *in = fopen(argv[1], "rb");
    if (!in) { perror(argv[1]); return EXIT_FAILURE; }
    if (fseeko(in, 0, SEEK_END) != 0) { perror(argv[1]); return EXIT_FAILURE; }
    off_t sz = ftello(in);
    if (sz < 0 || fseeko(in, 0, SEEK_SET) != 0) { perror(argv[1]); return EXIT_FAILURE; }

    char *seq = (char *)malloc((size_t)sz + 1);
    if (!seq) { fprintf(stderr, "ERROR: cannot allocate %lld bytes\n", (long long)sz); return EXIT_FAILURE; }
    if (fread(seq, 1, (size_t)sz, in) != (size_t)sz) { fprintf(stderr, "ERROR: read failed: %s\n", argv[1]); return EXIT_FAILURE; }
    fclose(in);

    uint64_t L = (uint64_t)sz;
    if (L > 0 && seq[L - 1] == '\n') L--;
    if (L == 0) { fprintf(stderr, "ERROR: %s contains no sequence\n", argv[1]); return EXIT_FAILURE; }
    if (L >= RS_MAX_GENOME_BP) {
        fprintf(stderr, "ERROR: genome is %llu bp; must be below %llu bp (see README, 'Genome size limits')\n",
                (unsigned long long)L, (unsigned long long)RS_MAX_GENOME_BP);
        return EXIT_FAILURE;
    }
    for (uint64_t i = 0; i < L; i++) {
        char c = seq[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            fprintf(stderr, "ERROR: non-ACGT byte 0x%02x at position %llu (run preprocess_genome first)\n",
                    (unsigned char)c, (unsigned long long)i);
            return EXIT_FAILURE;
        }
    }

    FILE *out = fopen(argv[2], "wb");
    if (!out) { perror(argv[2]); return EXIT_FAILURE; }
    if (fwrite(seq, 1, (size_t)L, out) != (size_t)L) { fprintf(stderr, "ERROR: write failed\n"); return EXIT_FAILURE; }

    char *buf = (char *)malloc(OUT_BUF);
    if (!buf) { fprintf(stderr, "ERROR: cannot allocate output buffer\n"); return EXIT_FAILURE; }
    size_t n = 0;
    for (uint64_t i = L; i-- > 0; ) {
        char c = seq[i];
        buf[n++] = (c == 'A') ? 'T' : (c == 'C') ? 'G' : (c == 'G') ? 'C' : 'A';
        if (n == OUT_BUF) {
            if (fwrite(buf, 1, n, out) != n) { fprintf(stderr, "ERROR: write failed\n"); return EXIT_FAILURE; }
            n = 0;
        }
    }
    buf[n++] = '\n';
    if (fwrite(buf, 1, n, out) != n || fclose(out) != 0) { fprintf(stderr, "ERROR: write failed\n"); return EXIT_FAILURE; }
    free(buf); free(seq);

    fprintf(stderr, "[derive_compactref] L = %llu bp  ->  text = 2L = %llu chars (+ newline)\n",
            (unsigned long long)L, (unsigned long long)(2 * L));
    return EXIT_SUCCESS;
}
