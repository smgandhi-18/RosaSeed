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

#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>

#define RS1_REF4_MAGIC 0x5253315245463401ULL
#define RS1_REF4_VERSION 1u

typedef struct {
    uint64_t magic;
    uint32_t version;
    uint32_t bases_per_byte;   /* 4 */
    uint64_t n_bases;          /* number of A/C/G/T bases written */
    uint64_t packed_bytes;     /* ceil(n_bases / 4) */
} rs1_ref4_header_t;

static inline int nt_to_base4(int ch)
{
    switch (toupper((unsigned char)ch)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

static inline void ref4_set(uint8_t *buf, uint64_t i, uint8_t v)
{
    uint64_t b = i >> 2;                         /* byte index */
    uint8_t shift = (uint8_t)(6 - 2 * (i & 3));  /* 7-6, 5-4, 3-2, 1-0 */
    uint8_t mask = (uint8_t)(0x03u << shift);

    buf[b] = (uint8_t)((buf[b] & ~mask) | ((v & 0x03u) << shift));
}

static uint64_t count_valid_bases(const char *input_path,
                                  uint64_t *headers,
                                  uint64_t *whitespace,
                                  uint64_t *dollars,
                                  uint64_t *invalid)
{
    FILE *fp = fopen(input_path, "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open input %s: %s\n",
                input_path, strerror(errno));
        exit(EXIT_FAILURE);
    }

    uint64_t n = 0;
    int ch;

    *headers = 0;
    *whitespace = 0;
    *dollars = 0;
    *invalid = 0;

    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '>') {
            (*headers)++;
            while ((ch = fgetc(fp)) != EOF && ch != '\n') {
                /* skip FASTA header */
            }
            continue;
        }

        if (isspace((unsigned char)ch)) {
            (*whitespace)++;
            continue;
        }

        if (ch == '$' || ch == 0) {
            (*dollars)++;
            continue;
        }

        if (nt_to_base4(ch) >= 0) {
            n++;
        } else {
            (*invalid)++;
        }
    }

    fclose(fp);
    return n;
}

static void pack_reference(const char *input_path,
                           uint8_t *packed,
                           uint64_t expected_bases)
{
    FILE *fp = fopen(input_path, "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot reopen input %s: %s\n",
                input_path, strerror(errno));
        exit(EXIT_FAILURE);
    }

    uint64_t idx = 0;
    int ch;

    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '>') {
            while ((ch = fgetc(fp)) != EOF && ch != '\n') {
                /* skip FASTA header */
            }
            continue;
        }

        if (isspace((unsigned char)ch)) continue;
        if (ch == '$' || ch == 0) continue;

        int v = nt_to_base4(ch);
        if (v < 0) continue;

        if (idx >= expected_bases) {
            fprintf(stderr, "ERROR: base count exceeded expected count\n");
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        ref4_set(packed, idx, (uint8_t)v);
        idx++;
    }

    fclose(fp);

    if (idx != expected_bases) {
        fprintf(stderr,
                "ERROR: packed base count mismatch: expected=%llu got=%llu\n",
                (unsigned long long)expected_bases,
                (unsigned long long)idx);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv)
{
    if (argc != 3 && argc != 4) {
        fprintf(stderr,
                "Usage:\n"
                "  %s <input_ref_ascii.fa/txt> <output_ref4.bin> [expected_bases]\n\n"
                "Example:\n"
                "  %s compact_ref.txt ref4_packed.bin\n",
                argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    const char *input_path = argv[1];
    const char *output_path = argv[2];

    uint64_t headers = 0, whitespace = 0, dollars = 0, invalid = 0;
    uint64_t n_bases = count_valid_bases(input_path,
                                         &headers,
                                         &whitespace,
                                         &dollars,
                                         &invalid);

    if (argc == 4) {
        uint64_t expected = strtoull(argv[3], NULL, 10);
        if (n_bases != expected) {
            fprintf(stderr,
                    "ERROR: expected_bases mismatch: counted=%llu expected=%llu\n",
                    (unsigned long long)n_bases,
                    (unsigned long long)expected);
            return EXIT_FAILURE;
        }
    }

    uint64_t packed_bytes = (n_bases + 3ULL) >> 2;

    uint8_t *packed = (uint8_t *)calloc((size_t)packed_bytes, 1);
    if (!packed) {
        fprintf(stderr,
                "ERROR: calloc failed for %llu bytes\n",
                (unsigned long long)packed_bytes);
        return EXIT_FAILURE;
    }

    pack_reference(input_path, packed, n_bases);

    FILE *out = fopen(output_path, "wb");
    if (!out) {
        fprintf(stderr, "ERROR: cannot open output %s: %s\n",
                output_path, strerror(errno));
        free(packed);
        return EXIT_FAILURE;
    }

    rs1_ref4_header_t hdr;
    hdr.magic = RS1_REF4_MAGIC;
    hdr.version = RS1_REF4_VERSION;
    hdr.bases_per_byte = 4;
    hdr.n_bases = n_bases;
    hdr.packed_bytes = packed_bytes;

    if (fwrite(&hdr, sizeof(hdr), 1, out) != 1) {
        fprintf(stderr, "ERROR: failed to write header\n");
        fclose(out);
        free(packed);
        return EXIT_FAILURE;
    }

    if (fwrite(packed, 1, (size_t)packed_bytes, out) != packed_bytes) {
        fprintf(stderr, "ERROR: failed to write packed reference\n");
        fclose(out);
        free(packed);
        return EXIT_FAILURE;
    }

    fclose(out);
    free(packed);

    fprintf(stderr, "Reference base4 packing complete\n");
    fprintf(stderr, "  Input           : %s\n", input_path);
    fprintf(stderr, "  Output          : %s\n", output_path);
    fprintf(stderr, "  Bases packed    : %llu\n", (unsigned long long)n_bases);
    fprintf(stderr, "  Packed bytes    : %llu\n", (unsigned long long)packed_bytes);
    fprintf(stderr, "  Headers skipped : %llu\n", (unsigned long long)headers);
    fprintf(stderr, "  Whitespace skip : %llu\n", (unsigned long long)whitespace);
    fprintf(stderr, "  Dollars skipped : %llu\n", (unsigned long long)dollars);
    fprintf(stderr, "  Invalid skipped : %llu\n", (unsigned long long)invalid);

    return EXIT_SUCCESS;
}
