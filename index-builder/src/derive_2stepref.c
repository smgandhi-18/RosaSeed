#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_REF_LEN      7000000000ULL    /* allocation cap */
#define RS_MAX_GENOME_BP 4294967296ULL    /* 2^32 — real format limit */
#define WRITE_BUF_SIZE (1 << 20)    // 1 MB write buffer
#define LINE_LEN 80                 // FASTA line width

// ---------- utility functions ----------
int nt_to_int(char nt) {
    switch (toupper(nt)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

char rev_comp_nt(char nt) {
    switch (toupper(nt)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

char encode_pair(char n1, char n2) {
    int v1 = nt_to_int(n1);
    int v2 = nt_to_int(n2);
    if (v1 < 0 || v2 < 0) return 'N';
    int code = v1 * 4 + v2;
    return (code < 10) ? ('0' + code) : ('A' + (code - 10));
}

// ---------- Base-4 Reverse Complement Generation (NEW) ----------
// Generates the true reverse complement of the original Base-4 (ACGT) reference.
void generate_base4_reverse_complement(const char *ref_in, char *ref_rc) {
    size_t len = strlen(ref_in);
    // Iterate backwards through the input reference
    for (size_t i = 0; i < len; i++) {
        // Map the base from the end of ref_in to the start of ref_rc
        char original_base = ref_in[len - 1 - i];
        ref_rc[i] = rev_comp_nt(original_base);
    }
    ref_rc[len] = '\0';
}

// ---------- main two-step reference generation ----------
// Encodes the Base-4 sequence (ref) into two Base-16 sequences (ref_even and ref_odd).
void build_two_step_references(const char *ref, char *ref_even, char *ref_odd) {
    size_t len = strlen(ref);
    size_t even_idx = 0, odd_idx = 0;

    // Even indices (0, 2, 4, ...)
    for (size_t i = 0; i + 1 < len; i += 2)
        ref_even[even_idx++] = encode_pair(ref[i], ref[i+1]);
    ref_even[even_idx] = '\0';

    // Odd indices (1, 3, 5, ...)
    for (size_t i = 1; i + 1 < len; i += 2)
        ref_odd[odd_idx++] = encode_pair(ref[i], ref[i+1]);
    ref_odd[odd_idx] = '\0';
}



// Helper function to print the start and end of a segment
void print_segment_info(const char *name, const char *seq) {
    size_t len = strlen(seq);
    size_t print_len = (len < 5) ? len : 5;
    
    printf("\n--- %s ---\n", name);
    printf("Length: %zu symbols (Base-16)\n", len);

    // Print first 5 symbols
    printf("Start: ");
    for (size_t i = 0; i < print_len; i++) {
        putchar(seq[i]);
    }
    printf("\n");

    // Print last 5 symbols
    if (len > 5) {
        printf("End:   ");
        for (size_t i = len - 5; i < len; i++) {
            putchar(seq[i]);
        }
        printf("\n");
    } else if (len > 0) {
        printf("End:   (Sequence too short to show distinct last 5)\n");
    } else {
        printf("End:   (Empty sequence)\n");
    }
}

// NOTE: The previous build_reverse_complement function is removed as it RC'd the Base-16 string, which is now replaced by the more robust method in main().

// ---------- FASTA reading ----------
char* read_fasta(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening FASTA file");
        exit(1);
    }

    char *seq = malloc(MAX_REF_LEN);
    if (!seq) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    size_t len = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') continue; // skip header
        for (int i = 0; line[i]; i++) {
            char c = line[i];
            // Only include standard ACGT characters
            if (!isspace(c) && strchr("ACGTacgt", c))
                seq[len++] = toupper(c);
        }
        if (len >= MAX_REF_LEN - 1) {
            fprintf(stderr, "Reference exceeds MAX_REF_LEN limit.\n");
            exit(1);
        }
    }
    seq[len] = '\0';
    fclose(fp);
    return seq;
}

// ---------- FASTA writing (optimized, buffered) ----------
void write_fasta(const char *filename, const char *header, const char *seq) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error writing FASTA file");
        exit(1);
    }

    fprintf(fp, ">%s\n", header);
    size_t seq_len = strlen(seq);
    size_t count = 0;
    char *buf = malloc(WRITE_BUF_SIZE);
    if (!buf) {
        fprintf(stderr, "Buffer allocation failed\n");
        exit(1);
    }

    for (size_t i = 0; i < seq_len; i++) {
        buf[count++] = seq[i];
        if ((i + 1) % LINE_LEN == 0) buf[count++] = '\n';
        if (count >= WRITE_BUF_SIZE - LINE_LEN) {
            fwrite(buf, 1, count, fp);
            count = 0;
        }
    }

    if (count > 0) {
        // Write remaining sequence
        if (count > 0) fwrite(buf, 1, count, fp);
        // Ensure final newline if the last line wasn't exactly LINE_LEN
        if (seq_len % LINE_LEN != 0) fputc('\n', fp);
    }
    
    free(buf);
    fclose(fp);
}

// ---------- main driver ----------
int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <reference.fasta>\n", argv[0]);
        return 1;
    }

    // Step 1: Read original reference (Base-4)
    char *ref = read_fasta(argv[1]);
    size_t ref_len = strlen(ref);
    printf("Reference length: %zu bases\n", ref_len);
    
    if (ref_len > RS_MAX_GENOME_BP) {
        fprintf(stderr,
            "\n[FATAL] Reference is %zu bp, exceeding the RosaSeed index format\n"
            "        limit of %llu bp (~4.29 Gbp).  BWT positions are stored in\n"
            "        33-bit fields; a larger reference would be silently truncated.\n"
            "        See README, 'Genome size limits'.\n",
            ref_len, (unsigned long long)RS_MAX_GENOME_BP);
        exit(EXIT_FAILURE);
    }

    // Step 2: Allocate memory
    size_t half_len_base16 = ref_len / 2 + 2;
    size_t rc_base4_len = ref_len + 1;

    char *ref_rc_base4 = malloc(rc_base4_len); // New buffer for Base-4 RC sequence
    char *ref_even = malloc(half_len_base16);
    char *ref_odd  = malloc(half_len_base16);
    char *rc_even  = malloc(half_len_base16);
    char *rc_odd   = malloc(half_len_base16);

    if (!ref_rc_base4 || !ref_even || !ref_odd || !rc_even || !rc_odd) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // --- PHASE I: FORWARD SEGMENTS ---
    // 3A. Build two-step forward references (FE, FO)
    build_two_step_references(ref, ref_even, ref_odd);

    // --- PHASE II: REVERSE COMPLEMENT SEGMENTS (RCE, RCO) ---
    // 3B. Generate the Base-4 Reverse Complement sequence
    generate_base4_reverse_complement(ref, ref_rc_base4);

    // 3C. Build two-step references from the Base-4 RC sequence
    // This correctly generates the Base-16 RC Even and RC Odd segments
    build_two_step_references(ref_rc_base4, rc_even, rc_odd);

    // --- Inspection (NEW) ---
    printf("\n\n=============== Two-Step Reference Segment Details ==============\n");
    print_segment_info("Forward Even (FE)", ref_even);
    print_segment_info("Forward Odd (FO)", ref_odd);
    print_segment_info("RC Even (RCE) - Mapped from Base-4 RC", rc_even);
    print_segment_info("RC Odd (RCO) - Mapped from Base-4 RC", rc_odd);
    printf("=================================================================\n");
    

    // Step 4: Concatenate all four segments (FE + FO + RCE + RCO)
    size_t len_even = strlen(ref_even);
    size_t len_odd = strlen(ref_odd);
    size_t len_rc_even = strlen(rc_even);
    size_t len_rc_odd = strlen(rc_odd);

    fprintf(stderr, "FE length: %ld\nFO length: %ld\nRCE length: %ld\nRCO length: %ld\n", 
            len_even, len_odd, len_rc_even, len_rc_odd);
    size_t total_len = len_even + len_odd + len_rc_even + len_rc_odd;

    char *final_ref = malloc(total_len + 1);
    if (!final_ref) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    char *dst = final_ref;

	// Correct concatenation to match actual file:
	memcpy(dst, ref_even,  len_even);  dst += len_even;
	memcpy(dst, ref_odd,   len_odd);   dst += len_odd;
	memcpy(dst, rc_even,   len_rc_even); dst += len_rc_even;  // rc_even THIRD
	memcpy(dst, rc_odd,    len_rc_odd);  dst += len_rc_odd;   // rc_odd FOURTH
    
    *dst = '\0';

    // Step 5: Write final combined reference (buffered)
    write_fasta("GaplessGenomeT2T_2step_ref_verify.fasta", "two_step_reference", final_ref);
    printf("Wrote concatenated two-step reference to GaplessGenomeT2T_2step_ref_verify.fasta (length = %zu)\n", total_len);

    // Cleanup
    free(ref); free(ref_rc_base4);
    free(ref_even); free(ref_odd); 
    free(rc_even); free(rc_odd);
    free(final_ref);
    return 0;
}
