#!/usr/bin/env python3
# =============================================================================
#                            The MIT License
#
#    RosaSeed (RosaSeed: Faster and Accurate Short Read Alignment Using a Configurable Seeding Strategy),
#    Copyright (C) 2026  University of Alberta, Gandhi Shyama.
#
#    Permission is hereby granted, free of charge, to any person obtaining
#    a copy of this software and associated documentation files (the
#    "Software"), to deal in the Software without restriction, including
#    without limitation the rights to use, copy, modify, merge, publish,
#    distribute, sublicense, and/or sell copies of the Software, and to
#    permit persons to whom the Software is furnished to do so, subject to
#    the following conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.
#
#    Contacts: Shyama Gandhi <smgandhi@ualberta.ca>
# =============================================================================
"""
filter_n_reads_paired.py — Remove PAIRED FASTQ reads if EITHER mate contains
ambiguous / non-ACGT bases.

For paired-end data, a read pair is only usable if BOTH mates are valid. If the
R1 mate OR the R2 mate contains a base outside {A,C,G,T}, the ENTIRE pair is
dropped from both output files. This keeps R1 and R2 synchronized (the i-th
record in the filtered R1 always corresponds to the i-th record in the filtered
R2). Dropping only one mate would desynchronize the mates and corrupt pairing.

Handles plain and gzip-compressed FASTQ. R1 and R2 must have their records in the
same order (standard for Illumina output) and contain the same number of reads.

Usage:
    python3 filter_n_reads_paired.py R1.fastq.gz R2.fastq.gz \\
            R1.filtered.fastq.gz R2.filtered.fastq.gz --log run.log

    python3 filter_n_reads_paired.py in_1.fq in_2.fq out_1.fq out_2.fq \\
            --strict-uppercase
"""

import argparse
import gzip
import sys
import time
from pathlib import Path


def open_maybe_gzip(path, mode):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def build_allowed(allow_lowercase):
    allowed = set("ACGT")
    if allow_lowercase:
        allowed |= set("acgt")
    return allowed


def _read_record(fh):
    """Read one 4-line FASTQ record. Returns (header, seq, plus, qual) or None at EOF."""
    header = fh.readline()
    if not header:
        return None
    seq = fh.readline()
    plus = fh.readline()
    qual = fh.readline()
    if not qual:
        raise ValueError("Truncated FASTQ record (incomplete 4-line block).")
    return header, seq, plus, qual


def filter_paired(r1_in, r2_in, r1_out, r2_out, log_path=None,
                  allow_lowercase=True, report_every=1_000_000):
    allowed = build_allowed(allow_lowercase)

    n_pairs = 0
    n_kept = 0
    n_dropped = 0
    dropped_because_r1 = 0   # pair dropped due to bad base in R1 (R2 may be clean)
    dropped_because_r2 = 0   # pair dropped due to bad base in R2 (R1 may be clean)
    dropped_because_both = 0
    dropped_base_counts = {}

    t0 = time.time()

    f1i = open_maybe_gzip(r1_in, "r")
    f2i = open_maybe_gzip(r2_in, "r")
    f1o = open_maybe_gzip(r1_out, "w")
    f2o = open_maybe_gzip(r2_out, "w")

    try:
        while True:
            rec1 = _read_record(f1i)
            rec2 = _read_record(f2i)

            # Both files should end together. If one ends before the other, the
            # inputs are mismatched — stop loudly rather than silently truncate.
            if rec1 is None and rec2 is None:
                break
            if (rec1 is None) != (rec2 is None):
                raise ValueError(
                    f"R1 and R2 have different numbers of reads "
                    f"(desynchronized near pair {n_pairs + 1}). "
                    f"Inputs must be matched mates in the same order."
                )

            n_pairs += 1
            seq1 = rec1[1].strip()
            seq2 = rec2[1].strip()

            bad1 = set(seq1) - allowed
            bad2 = set(seq2) - allowed

            if bad1 or bad2:
                # Drop the WHOLE pair — write nothing to either output.
                n_dropped += 1
                if bad1 and bad2:
                    dropped_because_both += 1
                elif bad1:
                    dropped_because_r1 += 1
                else:
                    dropped_because_r2 += 1
                for b in bad1:
                    dropped_base_counts[b] = dropped_base_counts.get(b, 0) + seq1.count(b)
                for b in bad2:
                    dropped_base_counts[b] = dropped_base_counts.get(b, 0) + seq2.count(b)
            else:
                # Both mates clean — write both, keeping the pair together.
                f1o.writelines(rec1)
                f2o.writelines(rec2)
                n_kept += 1

            if report_every and n_pairs % report_every == 0:
                elapsed = time.time() - t0
                rate = n_pairs / elapsed if elapsed > 0 else 0
                print(f"  ... {n_pairs:,} pairs processed "
                      f"({n_kept:,} kept, {n_dropped:,} dropped) "
                      f"[{rate:,.0f} pairs/s]", file=sys.stderr)
    finally:
        f1i.close(); f2i.close(); f1o.close(); f2o.close()

    elapsed = time.time() - t0
    pct_dropped = (100.0 * n_dropped / n_pairs) if n_pairs else 0.0

    stats = {
        "r1_in": str(r1_in), "r2_in": str(r2_in),
        "r1_out": str(r1_out), "r2_out": str(r2_out),
        "total_pairs": n_pairs,
        "kept_pairs": n_kept,
        "dropped_pairs": n_dropped,
        "pct_dropped": pct_dropped,
        "dropped_because_r1_only": dropped_because_r1,
        "dropped_because_r2_only": dropped_because_r2,
        "dropped_because_both": dropped_because_both,
        "dropped_base_counts": dropped_base_counts,
        "elapsed_sec": elapsed,
        "allowed_alphabet": "".join(sorted(allowed)),
    }
    _emit_log(stats, log_path)
    return stats


def _emit_log(stats, log_path):
    L = []
    L.append("=" * 66)
    L.append("PAIRED-END FASTQ ambiguous-base filtering — summary")
    L.append("=" * 66)
    L.append(f"R1 input          : {stats['r1_in']}")
    L.append(f"R2 input          : {stats['r2_in']}")
    L.append(f"R1 output         : {stats['r1_out']}")
    L.append(f"R2 output         : {stats['r2_out']}")
    L.append(f"Allowed alphabet  : {{{stats['allowed_alphabet']}}}")
    L.append(f"Total pairs       : {stats['total_pairs']:,}")
    L.append(f"Pairs kept        : {stats['kept_pairs']:,}")
    L.append(f"Pairs dropped     : {stats['dropped_pairs']:,} "
             f"({stats['pct_dropped']:.4f}%)")
    L.append(f"    dropped (R1 bad only) : {stats['dropped_because_r1_only']:,}")
    L.append(f"    dropped (R2 bad only) : {stats['dropped_because_r2_only']:,}")
    L.append(f"    dropped (both bad)    : {stats['dropped_because_both']:,}")
    if stats["dropped_base_counts"]:
        L.append("Dropped-base tally (base: occurrences across dropped mates):")
        for b, c in sorted(stats["dropped_base_counts"].items(), key=lambda kv: -kv[1]):
            L.append(f"    {b!r}: {c:,}")
    else:
        L.append("No pairs contained ambiguous bases; nothing dropped.")
    L.append(f"Elapsed           : {stats['elapsed_sec']:.1f} s")
    L.append("=" * 66)
    text = "\n".join(L)
    print(text, file=sys.stderr)
    if log_path:
        with open(log_path, "w") as lf:
            lf.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Remove paired FASTQ reads if EITHER mate has ambiguous bases."
    )
    ap.add_argument("r1_in", help="Input R1 FASTQ (.fastq/.fq, optionally .gz)")
    ap.add_argument("r2_in", help="Input R2 FASTQ (.fastq/.fq, optionally .gz)")
    ap.add_argument("r1_out", help="Output R1 FASTQ (.gz -> gzip)")
    ap.add_argument("r2_out", help="Output R2 FASTQ (.gz -> gzip)")
    ap.add_argument("--log", default=None, help="Optional summary log path.")
    ap.add_argument("--allow-lowercase", action="store_true", default=True,
                    help="Treat a,c,g,t as valid (default: on).")
    ap.add_argument("--strict-uppercase", dest="allow_lowercase",
                    action="store_false",
                    help="Only A,C,G,T valid; lowercase bases cause a drop.")
    args = ap.parse_args()

    for p in (args.r1_in, args.r2_in):
        if not Path(p).exists():
            sys.exit(f"ERROR: input file not found: {p}")

    print(f"Filtering pair:\n  {args.r1_in}\n  {args.r2_in}", file=sys.stderr)
    stats = filter_paired(
        args.r1_in, args.r2_in, args.r1_out, args.r2_out,
        log_path=args.log, allow_lowercase=args.allow_lowercase,
    )
    if stats["kept_pairs"] == 0:
        sys.exit("ERROR: no pairs kept — check inputs and alphabet settings.")


if __name__ == "__main__":
    main()
