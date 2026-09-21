#!/usr/bin/env python3
import pysam
import sys

def compute_nm_from_cigar(read):
    """Approximate NM tag from CIGAR if NM tag missing."""
    if read is None or read.cigartuples is None:
        return None
    nm = 0
    for op, length in read.cigartuples:
        if op in (1, 2, 8):  # I or D or mismatch
            nm += length
    return nm


def edit_burden(read):
    """Total burden: I + D + X (structural edit burden)."""
    if read is None or read.cigartuples is None:
        return None
    burden = 0
    for op, length in read.cigartuples:
        if op in (1, 2, 8):
            burden += length
    return burden


def relaxed_structural_correct(t_read, a_read, tolerance=0.05):
    if t_read is None or a_read is None or a_read.is_unmapped:
        return False

    et = edit_burden(t_read)
    ea = edit_burden(a_read)
    if et is None or ea is None:
        return False

    read_len = t_read.query_length
    if read_len is None or read_len <= 0:
        return False

    return (abs(et - ea) / read_len) <= tolerance


def scc_sequence_consistent_correct(a_read, threshold=0.10):
    if a_read is None or a_read.is_unmapped:
        return False

    read_len = a_read.query_length
    if read_len is None or read_len <= 0:
        return False

    try:
        nm = a_read.get_tag("NM")
    except KeyError:
        nm = compute_nm_from_cigar(a_read)

    if nm is None:
        return False

    error_rate = nm / read_len
    return error_rate <= threshold

def open_alignment(filename):
    if filename.endswith(".bam"):
        return pysam.AlignmentFile(filename, "rb")
    else:
        return pysam.AlignmentFile(filename, "r")
# ----------------------------------------------------------------------
# Main Evaluation
# ----------------------------------------------------------------------

def evaluate(reference_sam, test_sam, tol_strict=5, tol_standard=50, ignore_refname=False):
    """
    reference_sam: proxy gold/reference aligner SAM (e.g. BWA-MEM2)
    test_sam: aligner under evaluation (e.g. Rosaseed)
    """

    truth = open_alignment(reference_sam)
    align = open_alignment(test_sam)

    # Key must stay (QNAME, mate): paired mates share a QNAME, so keying on the
    # name alone drops one mate per pair and can compare truth R2 against test R1.
    truth_dict = {
        (r.query_name, r.is_read1): r
        for r in truth.fetch(until_eof=True)
        if (not r.is_secondary) and (not r.is_supplementary)
    }
    align_dict = {
        (r.query_name, r.is_read1): r
        for r in align.fetch(until_eof=True)
        if (not r.is_secondary) and (not r.is_supplementary)
    }

    all_keys = set(truth_dict.keys()) | set(align_dict.keys())

    # Strict (±5 bp)
    TP = FP = FN = TN = 0

    # Standard (±50 bp)
    TP50 = FP50 = FN50 = TN50 = 0

    # Structural
    TP_struct_relaxed = FP_struct_relaxed = 0
    struct_considered = 0

    # SCC
    TP_scc = FP_scc = 0
    scc_considered = 0

    # Paper-friendly metrics
    correct_std_50 = 0
    correct_strict_5 = 0
    not_mapped = 0

    unmapped_aligner = 0
    total = 0

    FN_truth_recs = []
    FN_align_recs = []

    # Reads mapped by Rosaseed but unmapped/missing in BWA-MEM2
    test_only_mapped_ids = []
    test_only_truth_recs = []
    test_only_align_recs = []

    # New: keep FP read IDs only (lightweight)
    standard_fp_ids = []
    strict_fp_ids = []

    for key in all_keys:
        qname = key[0]
        t_read = truth_dict.get(key)   # BWA-MEM2 in your current usage
        a_read = align_dict.get(key)   # Rosaseed in your current usage
        total += 1

        # Label used in the dumped read-ID lists: name/1, name/2 when paired,
        # plain name for single-end.
        _rec = t_read if t_read is not None else a_read
        read_id = (f"{qname}/{1 if _rec.is_read1 else 2}"
                   if (_rec is not None and _rec.is_paired) else qname)

        # Treat missing record like unmapped for bookkeeping
        t_missing_or_unmapped = (t_read is None) or t_read.is_unmapped
        a_missing_or_unmapped = (a_read is None) or a_read.is_unmapped

        if a_missing_or_unmapped:
            not_mapped += 1
            unmapped_aligner += 1

        # Case 1: both missing/unmapped
        if t_missing_or_unmapped and a_missing_or_unmapped:
            TN += 1
            TN50 += 1
            continue

        # Case 2: reference mapped, Rosaseed missing/unmapped
        elif (not t_missing_or_unmapped) and a_missing_or_unmapped:
            FN += 1
            FN50 += 1
            if a_read is not None and a_read.is_unmapped:
                FN_align_recs.append(a_read)
            if t_read is not None:
                FN_truth_recs.append(t_read)
            continue

        # Case 3: reference missing/unmapped, Rosaseed mapped
        elif t_missing_or_unmapped and (not a_missing_or_unmapped):
            FP += 1
            FP50 += 1
            strict_fp_ids.append(read_id)
            standard_fp_ids.append(read_id)

            test_only_mapped_ids.append(read_id)
            if t_read is not None:
                test_only_truth_recs.append(t_read)
            if a_read is not None:
                test_only_align_recs.append(a_read)

            # SCC (mapped alignments only)
            if (a_read is not None) and (not a_read.is_unmapped):
                scc_considered += 1
                if scc_sequence_consistent_correct(a_read):
                    TP_scc += 1
                else:
                    FP_scc += 1

            continue

        # Case 4: both mapped -> compare positions/strand/reference
        else:
            same_ref = ignore_refname or (t_read.reference_name == a_read.reference_name)
            same_strand = (t_read.is_reverse == a_read.is_reverse)
            pos_diff = abs(t_read.reference_start - a_read.reference_start)

            if same_ref and same_strand and (pos_diff <= tol_standard):
                correct_std_50 += 1
            if same_ref and same_strand and (pos_diff <= tol_strict):
                correct_strict_5 += 1

            if same_ref and same_strand and (pos_diff <= tol_strict):
                TP += 1
            else:
                FP += 1
                strict_fp_ids.append(read_id)

            if same_ref and same_strand and (pos_diff <= tol_standard):
                TP50 += 1
            else:
                FP50 += 1
                standard_fp_ids.append(read_id)

        # Structural accuracy (mapped alignments only)
        if (a_read is not None) and (not a_read.is_unmapped) and (t_read is not None):
            struct_considered += 1
            if relaxed_structural_correct(t_read, a_read):
                TP_struct_relaxed += 1
            else:
                FP_struct_relaxed += 1

        # SCC (mapped alignments only)
        if (a_read is not None) and (not a_read.is_unmapped):
            scc_considered += 1
            if scc_sequence_consistent_correct(a_read):
                TP_scc += 1
            else:
                FP_scc += 1

    def metrics(TP, FP, FN, TN):
        precision = TP / (TP + FP) if (TP + FP) > 0 else 0
        recall = TP / (TP + FN) if (TP + FN) > 0 else 0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0
        accuracy = (TP + TN) / (TP + FP + FN + TN) if (TP + FP + FN + TN) > 0 else 0
        return precision, recall, f1, accuracy

    prec, rec, f1, acc = metrics(TP, FP, FN, TN)
    prec50, rec50, f150, acc50 = metrics(TP50, FP50, FN50, TN50)

    print("=== Evaluation ===")
    print(f"Compared reads (union of reference/test primary reads): {total}")
    print(f"Ignore refname: {ignore_refname}")
    print(f"Standard tolerance (bp): {tol_standard}")
    print(f"Strict tolerance (bp):   {tol_strict}")

    print("\n--- Paper-friendly metrics (recommended) ---")
    if total > 0:
        print(f"Correctly mapped (±{tol_standard} bp, same strand): {correct_std_50}  ({correct_std_50/total:.5f})")
        print(f"Strict accuracy  (±{tol_strict} bp, same strand): {correct_strict_5}  ({correct_strict_5/total:.5f})")
        print(f"Not mapped fraction (unmapped OR missing record):   {not_mapped}  ({not_mapped/total:.5f})")

    print(f"\n--- Standard accuracy (±{tol_standard} bp; TP/FP/FN/TN bookkeeping) ---")
    print(f"TP50: {TP50}, FP50: {FP50}, FN50: {FN50}, TN50: {TN50}")
    print(f"Standard Precision: {prec50:.5f}")
    print(f"Standard Recall:    {rec50:.5f}")
    print(f"Standard F1 Score:  {f150:.5f}")
    print(f"Standard Accuracy:  {acc50:.5f}")

    print(f"\n--- STRICT accuracy (±{tol_strict} bp; TP/FP/FN/TN bookkeeping) ---")
    print(f"TP: {TP}, FP: {FP}, FN: {FN}, TN: {TN}")
    print(f"Strict Precision: {prec:.5f}")
    print(f"Strict Recall:    {rec:.5f}")
    print(f"Strict F1 Score:  {f1:.5f}")
    print(f"Strict Accuracy:  {acc:.5f}")

    print("\n--- Relaxed structural accuracy (reference-guided; mapped-alignments only) ---")
    print(f"TP_struct_relaxed: {TP_struct_relaxed}")
    print(f"FP_struct_relaxed: {FP_struct_relaxed}")
    print(f"Structural considered (test aligned): {struct_considered}")
    if struct_considered > 0:
        print(f"Relaxed Structural accuracy: {TP_struct_relaxed/struct_considered:.5f}")
    else:
        print("Relaxed Structural accuracy: N/A")

    print("\n--- SCA: Sequence-consistent accuracy (test alignments only) ---")
    print(f"TP_sca: {TP_scc}")
    print(f"FP_sca: {FP_scc}")
    print(f"SCA considered (test aligned): {scc_considered}")
    if scc_considered > 0:
        print(f"SCA accuracy: {TP_scc/scc_considered:.5f}")
    else:
        print("SCA accuracy: N/A")

    if FN_truth_recs:
        with pysam.AlignmentFile("FN_truth_records.sam", "w", header=truth.header) as out1:
            for r in FN_truth_recs:
                out1.write(r)

        if FN_align_recs:
            with pysam.AlignmentFile("FN_align_records.sam", "w", header=align.header) as out2:
                for r in FN_align_recs:
                    out2.write(r)

        print(f"\nWrote {len(FN_truth_recs)} FN truth/reference records.")
        if FN_align_recs:
            print(f"Wrote {len(FN_align_recs)} FN align/test records (unmapped cases only).")

    if test_only_mapped_ids:
        with open("test_only_mapped_read_ids.txt", "w") as f:
            for rid in test_only_mapped_ids:
                f.write(rid + "\n")

        if test_only_truth_recs:
            with pysam.AlignmentFile("test_only_reference_records.sam", "w", header=truth.header) as out1:
                for r in test_only_truth_recs:
                    out1.write(r)

        with pysam.AlignmentFile("test_only_align_records.sam", "w", header=align.header) as out2:
            for r in test_only_align_recs:
                out2.write(r)

        print(f"\nTest aligner mapped {len(test_only_mapped_ids)} reads that reference left unmapped or missing.")
        print("Wrote test_only_mapped_read_ids.txt")
        if test_only_truth_recs:
            print("Wrote test_only_reference_records.sam")
        print("Wrote test_only_align_records.sam")

    if standard_fp_ids:
        with open("standard_fp_all_read_ids.txt", "w") as f:
            for rid in standard_fp_ids:
                f.write(rid + "\n")
        print(f"Wrote {len(standard_fp_ids)} standard FP read IDs to standard_fp_all_read_ids.txt")

    if strict_fp_ids:
        with open("strict_fp_all_read_ids.txt", "w") as f:
            for rid in strict_fp_ids:
                f.write(rid + "\n")
        print(f"Wrote {len(strict_fp_ids)} strict FP read IDs to strict_fp_all_read_ids.txt")

    print("\n--- Unmapped summary ---")
    print(f"Unmapped by test aligner (includes missing records): {unmapped_aligner}")
    if total > 0:
        print(f"Not-mapped fraction (same as above):               {not_mapped/total:.5f}")

    truth.close()
    align.close()


# ----------------------------------------------------------------------
# ENTRY
# ----------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python evaluate_alignments.py reference.sam test.sam [tol_strict] [--tol50=50] [--ignore-refname]")
        sys.exit(1)

    reference_sam = sys.argv[1]
    test_sam = sys.argv[2]

    tol_strict = 5
    tol_standard = 50
    ignore_refname = False

    for arg in sys.argv[3:]:
        if arg.startswith("--tol50="):
            tol_standard = int(arg.split("=", 1)[1])
        elif arg == "--ignore-refname":
            ignore_refname = True
        elif arg.isdigit():
            tol_strict = int(arg)

    evaluate(
        reference_sam,
        test_sam,
        tol_strict=tol_strict,
        tol_standard=tol_standard,
        ignore_refname=ignore_refname
    )
