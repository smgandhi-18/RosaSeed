/*
 * preprocess_genome.c
 *
 * Preprocesses any FASTA file into a clean single-sequence file
 * ready for derive_2stepref (2-step RosaSeed pipeline).
 *
 * What it does:
 *   1. Reads any FASTA — single or multi-chromosome, any line width
 *   2. Strips all FASTA headers (lines starting with '>')
 *   3. Converts lowercase to uppercase
 *   4. Replaces ALL non-ACGT characters (N, R, Y, S, W, K, M, B, D, H, V,
 *      and anything else) with 'A'
 *      Rationale: preserves sequence length and coordinates.
 *      Silent dropping (as derive_2stepref does internally) would shift
 *      all coordinates downstream of any N — replaced-with-A does not.
 *   5. Writes output as a single-line plain text file (no FASTA header)
 *      OR as a clean single-header FASTA — controlled by -f flag
 *
 * Usage:
 *   ./preprocess_genome <input.fa> <output.txt> [options]
 *
 * Options:
 *   -f    Write output as FASTA (adds >genome header line)
 *         Default: plain text, no header (ready for derive_2stepref)
 *
 * Examples:
 *   ./preprocess_genome hg38.fa hg38_clean.txt
 *   ./preprocess_genome T2T_CHM13.fa T2T_clean.txt
 *   ./preprocess_genome GRCh38_chroms.fa GRCh38_clean.txt
 *
 * Output stats printed to stderr:
 *   - Total bases read
 *   - Characters replaced (N count, other ambiguity count)
 *   - Chromosomes / headers found
 *   - Output length
 */

#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>

#define READ_BUF_SIZE  (64  * 1024 * 1024)   /* 64 MB read buffer  */
#define WRITE_BUF_SIZE (64  * 1024 * 1024)   /* 64 MB write buffer */

static void usage(const char *prog) {
    fprintf(stderr,
        "\nUsage: %s <input.fa> <output.txt> [-f]\n\n"
        "  input.fa     any FASTA file (single or multi-chromosome)\n"
        "  output.txt   clean output (single line, no header by default)\n"
        "  -f           write output as FASTA with >genome header\n\n"
        "Non-ACGT characters are replaced with 'A' (not dropped).\n"
        "Lowercase is converted to uppercase.\n\n"
        "Examples:\n"
        "  %s hg38.fa           hg38_clean.txt\n"
        "  %s T2T_CHM13.fa      T2T_clean.txt\n"
        "  %s multi_chrom.fa    genome_clean.txt\n\n",
        prog, prog, prog, prog);
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
    if (argc < 3 || argc > 4) usage(argv[0]);

    const char *in_path  = argv[1];
    const char *out_path = argv[2];
    int write_fasta = 0;
    if (argc == 4) {
        if (strcmp(argv[3], "-f") == 0) write_fasta = 1;
        else { fprintf(stderr, "Unknown option: %s\n", argv[3]); usage(argv[0]); }
    }

    FILE *fin = fopen(in_path, "rb");
    if (!fin) {
        fprintf(stderr, "Cannot open input '%s': %s\n", in_path, strerror(errno));
        return EXIT_FAILURE;
    }
    FILE *fout = fopen(out_path, "wb");
    if (!fout) {
        fprintf(stderr, "Cannot open output '%s': %s\n", out_path, strerror(errno));
        fclose(fin);
        return EXIT_FAILURE;
    }

    uint8_t *rbuf = (uint8_t *)malloc(READ_BUF_SIZE);
    uint8_t *wbuf = (uint8_t *)malloc(WRITE_BUF_SIZE);
    if (!rbuf || !wbuf) { fprintf(stderr, "malloc failed\n"); return EXIT_FAILURE; }

    /* write FASTA header if requested */
    if (write_fasta) {
        fprintf(fout, ">genome\n");
    }

    /* ── stats ── */
    uint64_t n_headers    = 0;
    uint64_t n_bases_read = 0;   /* total ACGT + non-ACGT (before replacement) */
    uint64_t n_replaced_N = 0;   /* N characters replaced with A */
    uint64_t n_replaced_other = 0; /* other ambiguity codes replaced */
    uint64_t n_out        = 0;   /* bases written */

    /* ── lookup table ──
     * lut[byte] = output byte to write, or 0 = skip (whitespace/header) */
    uint8_t lut[256];
    memset(lut, 0, sizeof(lut));

    /* ACGT → themselves (uppercase) */
    lut['A'] = lut['a'] = 'A';
    lut['C'] = lut['c'] = 'C';
    lut['G'] = lut['g'] = 'G';
    lut['T'] = lut['t'] = 'T';

    /* N → A (replacement, not drop — preserves coordinates) */
    lut['N'] = lut['n'] = 'A';

    /* IUPAC ambiguity codes → A */
    const char *ambig = "RYSWKMBDHVrysWkmbdhv";
    for (int i = 0; ambig[i]; i++) lut[(uint8_t)ambig[i]] = 'A';

    /* ── state machine ──
     * in_header: currently consuming a header line, skip until newline */
    int in_header   = 0;
    size_t w_pos    = 0;        /* write buffer fill level */

    clock_t t0 = clock();
    size_t nread;

    /* print banner */
    fprintf(stderr, "[preprocess_genome]\n");
    fprintf(stderr, "  Input  : %s\n", in_path);
    fprintf(stderr, "  Output : %s\n", out_path);
    fprintf(stderr, "  Format : %s\n", write_fasta ? "FASTA" : "plain text");
    fprintf(stderr, "  Processing...\n\n");

    while ((nread = fread(rbuf, 1, READ_BUF_SIZE, fin)) > 0) {
        for (size_t i = 0; i < nread; i++) {
            uint8_t c = rbuf[i];

            /* ── header line handling ── */
            if (c == '>') {
                in_header = 1;
                n_headers++;
                continue;
            }
            if (in_header) {
                if (c == '\n') in_header = 0;
                continue;
            }

            /* ── whitespace: skip silently ── */
            if (c == '\n' || c == '\r' || c == ' ' || c == '\t') continue;

            /* ── classify and replace ── */
            uint8_t out_c = lut[c];

            if (out_c == 0) {
                /* completely unknown character — replace with A */
                out_c = 'A';
                n_replaced_other++;
            } else if (c == 'N' || c == 'n') {
                n_replaced_N++;
            } else if (out_c == 'A' &&
                       c != 'A' && c != 'a') {
                /* ambiguity code replaced */
                n_replaced_other++;
            }

            n_bases_read++;

            /* ── write to buffer ── */
            wbuf[w_pos++] = out_c;
            if (w_pos == WRITE_BUF_SIZE) {
                fwrite(wbuf, 1, w_pos, fout);
                n_out += w_pos;
                w_pos = 0;

                /* progress every ~500M bases */
                if (n_bases_read % 500000000ULL < WRITE_BUF_SIZE) {
                    double el = (double)(clock()-t0)/CLOCKS_PER_SEC;
                    fprintf(stderr, "\r  %.2f Gbp read  (%.0fs)  "
                                    "N_replaced=%llu  other_replaced=%llu",
                            (double)n_bases_read/1e9, el,
                            (unsigned long long)n_replaced_N,
                            (unsigned long long)n_replaced_other);
                }
            }
        }
    }

    /* flush remaining */
    if (w_pos > 0) {
        fwrite(wbuf, 1, w_pos, fout);
        n_out += w_pos;
    }

    /* trailing newline */
    fputc('\n', fout);

    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;

    fclose(fin); fclose(fout);
    free(rbuf);  free(wbuf);

    /* ── summary ── */
    uint64_t n_acgt_original = n_bases_read - n_replaced_N - n_replaced_other;

    fprintf(stderr, "\r                                                              \r");
    fprintf(stderr, "[preprocess_genome] done in %.1fs\n\n", elapsed);
    fprintf(stderr, "  ── Input stats ──────────────────────────────\n");
    fprintf(stderr, "  Chromosomes / headers found : %llu\n",
            (unsigned long long)n_headers);
    fprintf(stderr, "  Total characters processed  : %llu\n",
            (unsigned long long)n_bases_read);
    fprintf(stderr, "  Pure ACGT (kept as-is)      : %llu  (%.2f%%)\n",
            (unsigned long long)n_acgt_original,
            100.0 * n_acgt_original / (double)n_bases_read);
    fprintf(stderr, "  N replaced with A           : %llu  (%.2f%%)\n",
            (unsigned long long)n_replaced_N,
            100.0 * n_replaced_N / (double)n_bases_read);
    fprintf(stderr, "  Other ambiguity → A         : %llu  (%.2f%%)\n",
            (unsigned long long)n_replaced_other,
            100.0 * n_replaced_other / (double)n_bases_read);
    fprintf(stderr, "\n  ── Output stats ─────────────────────────────\n");
    fprintf(stderr, "  Output bases written        : %llu\n",
            (unsigned long long)n_out);
    fprintf(stderr, "  Output file                 : %s\n", out_path);
    fprintf(stderr, "\n");

    if (n_replaced_N > 0 || n_replaced_other > 0) {
        fprintf(stderr,
            "  [NOTE] Non-ACGT characters were replaced with 'A'.\n"
            "         Sequence length is preserved — coordinates are intact.\n"
            "         Reads CAN align to replaced regions — if this is\n"
            "         undesirable, remove those chromosomes from the input.\n\n");
    } else {
        fprintf(stderr,
            "  [OK] Input was pure ACGT — no replacements needed.\n\n");
    }

    return EXIT_SUCCESS;
}
