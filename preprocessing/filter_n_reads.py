#!/usr/bin/env python3
"""
filter_n_reads.py — Remove FASTQ reads containing ambiguous / non-ACGT bases.

Used to preprocess sequencing datasets for the RosaSeed evaluation: any read whose
sequence contains a base outside {A, C, G, T} (upper- or lower-case) is discarded,
so that the 2-bit-encoded seeding pipeline sees only canonical bases. The same
filtered read set is used for every aligner benchmarked, ensuring identical input.

Handles plain and gzip-compressed FASTQ (.fastq / .fq / .fastq.gz / .fq.gz),
single-end or one mate at a time. Writes a filtered FASTQ and a log with counts.

Usage:
    python3 filter_n_reads.py  input.fastq.gz  output.fastq.gz
    python3 filter_n_reads.py  input.fastq     output.fastq     --log run.log
    python3 filter_n_reads.py  in.fq.gz  out.fq.gz  --allow-lowercase

Notes:
  * By default the accepted alphabet is {A,C,G,T,a,c,g,t}; a read is dropped if ANY
    base falls outside it (this includes 'N' and any IUPAC ambiguity code).
  * Output compression is chosen from the output filename (.gz -> gzip).
  * The script never modifies reads; it only keeps or drops whole records.
"""

import argparse
import gzip
import sys
import time
from pathlib import Path


def open_maybe_gzip(path, mode):
    """Open a path transparently as gzip if it ends in .gz, else plain."""
    path = str(path)
    if path.endswith(".gz"):
        # text mode so we iterate lines as str; gzip handles the (de)compression
        return gzip.open(path, mode + "t")
    return open(path, mode)


def build_allowed(allow_lowercase):
    allowed = set("ACGT")
    if allow_lowercase:
        allowed |= set("acgt")
    return allowed


def filter_fastq(in_path, out_path, log_path=None, allow_lowercase=True,
                 report_every=1_000_000):
    """
    Stream a FASTQ file, dropping any 4-line record whose sequence line contains
    a base outside the allowed alphabet. Returns a dict of counts.
    """
    allowed = build_allowed(allow_lowercase)
    # Precompute the set of *disallowed* detection via translation for speed:
    # a read is bad if set(seq) - allowed is non-empty.

    n_total = 0
    n_kept = 0
    n_dropped = 0
    dropped_base_counts = {}   # e.g. {'N': 12345, 'R': 3, ...}

    t0 = time.time()

    fin = open_maybe_gzip(in_path, "r")
    fout = open_maybe_gzip(out_path, "w")

    try:
        while True:
            header = fin.readline()
            if not header:
                break  # EOF
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            if not qual:
                # Truncated final record — FASTQ must come in 4-line groups.
                raise ValueError(
                    f"Truncated FASTQ record near read {n_total + 1} "
                    f"(incomplete 4-line block). File may be corrupt."
                )

            n_total += 1
            seq_stripped = seq.strip()

            # Find any bases not in the allowed alphabet.
            bad = set(seq_stripped) - allowed
            if bad:
                n_dropped += 1
                for b in bad:
                    dropped_base_counts[b] = (
                        dropped_base_counts.get(b, 0)
                        + seq_stripped.count(b)
                    )
            else:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)
                n_kept += 1

            if report_every and n_total % report_every == 0:
                elapsed = time.time() - t0
                rate = n_total / elapsed if elapsed > 0 else 0
                print(f"  ... {n_total:,} reads processed "
                      f"({n_kept:,} kept, {n_dropped:,} dropped) "
                      f"[{rate:,.0f} reads/s]", file=sys.stderr)
    finally:
        fin.close()
        fout.close()

    elapsed = time.time() - t0
    pct_dropped = (100.0 * n_dropped / n_total) if n_total else 0.0

    stats = {
        "input": str(in_path),
        "output": str(out_path),
        "total_reads": n_total,
        "kept_reads": n_kept,
        "dropped_reads": n_dropped,
        "pct_dropped": pct_dropped,
        "dropped_base_counts": dropped_base_counts,
        "elapsed_sec": elapsed,
        "allowed_alphabet": "".join(sorted(allowed)),
    }

    _emit_log(stats, log_path)
    return stats


def _emit_log(stats, log_path):
    """Write a human-readable summary to stderr and, if given, to a log file."""
    lines = []
    lines.append("=" * 60)
    lines.append("FASTQ ambiguous-base filtering — summary")
    lines.append("=" * 60)
    lines.append(f"Input file        : {stats['input']}")
    lines.append(f"Output file       : {stats['output']}")
    lines.append(f"Allowed alphabet  : {{{stats['allowed_alphabet']}}}")
    lines.append(f"Total reads       : {stats['total_reads']:,}")
    lines.append(f"Reads kept        : {stats['kept_reads']:,}")
    lines.append(f"Reads dropped     : {stats['dropped_reads']:,} "
                 f"({stats['pct_dropped']:.4f}%)")
    if stats["dropped_base_counts"]:
        lines.append("Dropped-base tally (base: total occurrences in dropped reads):")
        for b, c in sorted(stats["dropped_base_counts"].items(),
                           key=lambda kv: -kv[1]):
            lines.append(f"    {b!r}: {c:,}")
    else:
        lines.append("No reads contained ambiguous bases; nothing dropped.")
    lines.append(f"Elapsed           : {stats['elapsed_sec']:.1f} s")
    lines.append("=" * 60)

    text = "\n".join(lines)
    print(text, file=sys.stderr)

    if log_path:
        with open(log_path, "w") as lf:
            lf.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Remove FASTQ reads containing ambiguous / non-ACGT bases."
    )
    ap.add_argument("input", help="Input FASTQ (.fastq/.fq, optionally .gz)")
    ap.add_argument("output", help="Output FASTQ (.gz extension -> gzip output)")
    ap.add_argument("--log", default=None,
                    help="Optional path to write the summary log.")
    ap.add_argument("--allow-lowercase", action="store_true", default=True,
                    help="Treat a,c,g,t as valid (default: on).")
    ap.add_argument("--strict-uppercase", dest="allow_lowercase",
                    action="store_false",
                    help="Only A,C,G,T are valid; lowercase bases are dropped.")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"ERROR: input file not found: {args.input}")

    print(f"Filtering {args.input} -> {args.output}", file=sys.stderr)
    stats = filter_fastq(
        args.input, args.output,
        log_path=args.log,
        allow_lowercase=args.allow_lowercase,
    )
    # Non-zero exit if nothing was written (likely a problem worth noticing).
    if stats["kept_reads"] == 0:
        sys.exit("ERROR: no reads were kept — check input and alphabet settings.")


if __name__ == "__main__":
    main()
