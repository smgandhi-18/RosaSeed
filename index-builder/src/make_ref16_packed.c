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
 * make_ref16_packed.c
 *
 * Packs a base-16 ASCII reference file into a binary nibble-packed file.
 * Replaces the slow fgetc-based version with fread (64 MB chunks).
 *
 * Input:  ASCII base-16 reference (.txt, single-line or multi-line)
 *         Valid chars: 0-9, A-F, a-f (case-insensitive)
 *         Whitespace, '>', '$' are skipped silently
 *
 * Output: binary file
 *   [0..7]   uint64_t n_symbols    (8-byte little-endian symbol count)
 *   [8..]    uint8_t  packed[]     (ceil(n_symbols/2) bytes, nibble-packed)
 *
 * Packing (matches load_data.c):
 *   byte b = packed[i>>1]
 *   even index i → low  nibble: (byte & 0x0F)
 *   odd  index i → high nibble: (byte >> 4) & 0x0F
 *   So: packed[b] = sym[2b] | (sym[2b+1] << 4)
 *
 * Usage:
 *   ./make_ref16_packed <input_base16.txt> <output_ref16_packed.bin>
 *
 * Example:
 *   ./make_ref16_packed genome_2step_ref.txt ref16_packed.bin
 */

#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#define READ_BUF  (64 * 1024 * 1024)   /* 64 MB read buffer */

/* lookup table: ASCII → nibble value (0-15), or 255 = skip */
static uint8_t lut[256];

static void init_lut(void) {
    memset(lut, 255, sizeof(lut));
    for (int i = 0; i <= 9; i++) lut['0'+i] = (uint8_t)i;
    for (int i = 0; i <= 5; i++) {
        lut['A'+i] = (uint8_t)(10+i);
        lut['a'+i] = (uint8_t)(10+i);
    }
}

static inline void ref16_set(uint8_t *buf, uint64_t i, uint8_t v) {
    uint64_t b = i >> 1;
    if (i & 1)
        buf[b] = (uint8_t)((buf[b] & 0x0F) | (uint8_t)((v & 0x0F) << 4));
    else
        buf[b] = (uint8_t)((buf[b] & 0xF0) | (v & 0x0F));
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr,
            "\nUsage: %s <input_base16.txt> <output_ref16_packed.bin>\n\n"
            "  input_base16.txt      ASCII 2-step reference (single or multi-line)\n"
            "  output_ref16_packed.bin  nibble-packed binary output\n\n"
            "Example:\n"
            "  %s genome_2step_ref.txt ref16_packed.bin\n\n",
            argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    init_lut();

    uint8_t *rbuf = (uint8_t *)malloc(READ_BUF);
    if (!rbuf) { fprintf(stderr,"malloc read buf failed\n"); return 1; }

    /* ── pass 1: count valid base-16 symbols ── */
    FILE *fin = fopen(argv[1], "rb");
    if (!fin) { perror(argv[1]); return 1; }

    fprintf(stderr, "[make_ref16_packed] Pass 1: counting symbols...\n");

    uint64_t n_symbols = 0;
    size_t nread;
    clock_t t0 = clock();

    while ((nread = fread(rbuf, 1, READ_BUF, fin)) > 0) {
        for (size_t i = 0; i < nread; i++)
            if (lut[rbuf[i]] != 255) n_symbols++;
    }
    fclose(fin);

    double t1 = (double)(clock()-t0)/CLOCKS_PER_SEC;
    fprintf(stderr, "  symbols found : %llu  (%.1fs)\n",
            (unsigned long long)n_symbols, t1);

    /* ── allocate packed buffer ── */
    uint64_t n_bytes = (n_symbols + 1) >> 1;   
    fprintf(stderr, "  packed bytes  : %llu  (%.2f GB)\n",
            (unsigned long long)n_bytes, (double)n_bytes / 1e9);

    uint8_t *packed = (uint8_t *)calloc((size_t)n_bytes, 1);
    if (!packed) {
        fprintf(stderr,"calloc failed for %llu bytes\n",
                (unsigned long long)n_bytes);
        return 1;
    }

    /* ── pass 2: pack symbols into nibbles ── */
    fin = fopen(argv[1], "rb");
    if (!fin) { perror(argv[1]); free(packed); return 1; }

    fprintf(stderr, "[make_ref16_packed] Pass 2: packing symbols...\n");

    uint64_t idx = 0;
    clock_t t2 = clock();

    while ((nread = fread(rbuf, 1, READ_BUF, fin)) > 0) {
        for (size_t i = 0; i < nread; i++) {
            uint8_t v = lut[rbuf[i]];
            if (v != 255) ref16_set(packed, idx++, v);
        }
        if (idx % 500000000ULL < (uint64_t)nread) {
            double el = (double)(clock()-t2)/CLOCKS_PER_SEC;
            fprintf(stderr, "\r  %llu / %llu symbols  (%.0fs)",
                    (unsigned long long)idx,
                    (unsigned long long)n_symbols, el);
        }
    }
    fclose(fin);

    double t3 = (double)(clock()-t2)/CLOCKS_PER_SEC;
    fprintf(stderr, "\r  %llu symbols packed  (%.1fs)\n",
            (unsigned long long)idx, t3);

    if (idx != n_symbols) {
        fprintf(stderr,"[ERROR] symbol count mismatch: pass1=%llu pass2=%llu\n",
                (unsigned long long)n_symbols,
                (unsigned long long)idx);
        free(packed); return 1;
    }

    FILE *out = fopen(argv[2], "wb");
    if (!out) { perror(argv[2]); free(packed); return 1; }

    if (fwrite(&n_symbols, sizeof(uint64_t), 1, out) != 1) {
        fprintf(stderr,"Write header failed\n"); free(packed); fclose(out); return 1;
    }
    if (fwrite(packed, 1, (size_t)n_bytes, out) != (size_t)n_bytes) {
        fprintf(stderr,"Write data failed\n"); free(packed); fclose(out); return 1;
    }
    fclose(out);
    free(packed);
    free(rbuf);

    double total = (double)(clock()-t0)/CLOCKS_PER_SEC;
    fprintf(stderr, "\n[make_ref16_packed] done in %.1fs\n", total);
    fprintf(stderr, "  n_symbols     : %llu\n",  (unsigned long long)n_symbols);
    fprintf(stderr, "  packed bytes  : %llu\n",  (unsigned long long)n_bytes);
    fprintf(stderr, "  output size   : %llu bytes (header + data)\n",
            (unsigned long long)(sizeof(uint64_t) + n_bytes));
    fprintf(stderr, "  output file   : %s\n", argv[2]);

    fprintf(stderr, "\n  Sanity checks:\n");
    fprintf(stderr, "  %s  n_symbols = %llu  (expected 6234551000 for T2T 2-step)\n",
            n_symbols == 6234551000ULL ? "[OK]  " : "[NOTE]",
            (unsigned long long)n_symbols);
    fprintf(stderr, "  %s  packed_bytes = ceil(n/2) = %llu\n",
            n_bytes == (n_symbols+1)/2 ? "[OK]  " : "[WARN]",
            (unsigned long long)n_bytes);

    return EXIT_SUCCESS;
}
