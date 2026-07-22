#!/usr/bin/env python3
"""
filter_n_reads_fast.py — Faster single-end FASTQ ambiguous-base filter.

Same behavior as filter_n_reads.py (drops any read with a base outside
{A,C,G,T[,acgt]}), but faster:
  * Uses `pigz`/`gzip` as external decompression/compression subprocesses
    (multi-threaded pigz is much faster than Python's gzip module).
  * Uses a fast bytes translation-table membership check instead of building a
    Python set per read.
  * Reads/writes in large blocks.

Falls back to Python gzip automatically if pigz/gzip binaries are unavailable.

Usage:
    python3 filter_n_reads_fast.py in.fastq.gz out.fastq.gz --log run.log
    python3 filter_n_reads_fast.py in.fq out.fq --strict-uppercase
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path


def _decompress_proc(path):
    """Return a file-like stdout stream decompressing `path`, or a plain file."""
    p = str(path)
    if p.endswith(".gz"):
        exe = shutil.which("pigz") or shutil.which("gzip")
        if exe:
            proc = subprocess.Popen([exe, "-dc", p], stdout=subprocess.PIPE)
            return proc.stdout, proc
        import gzip
        return gzip.open(p, "rb"), None
    return open(p, "rb"), None


def _compress_proc(path):
    """Return a file-like stdin stream compressing to `path`, or a plain file."""
    p = str(path)
    if p.endswith(".gz"):
        exe = shutil.which("pigz") or shutil.which("gzip")
        if exe:
            fout = open(p, "wb")
            proc = subprocess.Popen([exe, "-c"], stdin=subprocess.PIPE, stdout=fout)
            return proc.stdin, (proc, fout)
        import gzip
        return gzip.open(p, "wb"), None
    return open(p, "wb"), None


def build_bad_table(allow_lowercase):
    """
    Build a 256-byte translation table where allowed bases map to 0 and everything
    else maps to 1. Then a read is 'bad' iff its translated bytes contain a 1 —
    i.e. iff seq.translate(table) has any nonzero byte. We detect that cheaply by
    translating and checking against an all-zero bytes object of the same length.
    Simpler & robust: use a 'delete everything allowed' approach.
    """
    allowed = b"ACGT" + (b"acgt" if allow_lowercase else b"")
    # deletion table: after removing all allowed bytes, a clean read becomes b"".
    return allowed


def filter_fast(in_path, out_path, log_path=None, allow_lowercase=True,
                report_every=1_000_000):
    allowed = build_bad_table(allow_lowercase)
    # A read line is "clean" iff removing all allowed bytes leaves nothing but
    # the trailing newline. We strip newline first.
    delete_allowed = bytes(b for b in range(256) if b in allowed)

    n_total = n_kept = n_dropped = 0
    dropped_base_counts = {}
    t0 = time.time()

    fin, pin = _decompress_proc(in_path)
    fout, pout = _compress_proc(out_path)

    try:
        readline = fin.readline
        write = fout.write
        while True:
            header = readline()
            if not header:
                break
            seq = readline()
            plus = readline()
            qual = readline()
            if not qual:
                raise ValueError(
                    f"Truncated FASTQ record near read {n_total + 1}."
                )
            n_total += 1

            seq_body = seq.rstrip(b"\n")
            # Fast check: delete all allowed bytes; if anything remains, it's bad.
            leftover = seq_body.translate(None, delete_allowed)
            if leftover:
                n_dropped += 1
                for b in set(leftover):
                    ch = chr(b)
                    dropped_base_counts[ch] = (
                        dropped_base_counts.get(ch, 0) + leftover.count(b)
                    )
            else:
                write(header); write(seq); write(plus); write(qual)
                n_kept += 1

            if report_every and n_total % report_every == 0:
                el = time.time() - t0
                rate = n_total / el if el > 0 else 0
                print(f"  ... {n_total:,} reads ({n_kept:,} kept, "
                      f"{n_dropped:,} dropped) [{rate:,.0f} reads/s]",
                      file=sys.stderr)
    finally:
        try: fin.close()
        except Exception: pass
        try: fout.close()
        except Exception: pass
        # let subprocesses flush/finish
        if pin:
            pin.wait()
        if pout:
            proc, fh = pout
            proc.wait(); fh.close()

    el = time.time() - t0
    pct = (100.0 * n_dropped / n_total) if n_total else 0.0
    stats = {
        "input": str(in_path), "output": str(out_path),
        "total_reads": n_total, "kept_reads": n_kept,
        "dropped_reads": n_dropped, "pct_dropped": pct,
        "dropped_base_counts": dropped_base_counts,
        "elapsed_sec": el,
        "allowed_alphabet": allowed.decode(),
    }
    _emit_log(stats, log_path)
    return stats


def _emit_log(stats, log_path):
    L = ["=" * 60,
         "FASTQ ambiguous-base filtering — summary (fast)",
         "=" * 60,
         f"Input file        : {stats['input']}",
         f"Output file       : {stats['output']}",
         f"Allowed alphabet  : {{{stats['allowed_alphabet']}}}",
         f"Total reads       : {stats['total_reads']:,}",
         f"Reads kept        : {stats['kept_reads']:,}",
         f"Reads dropped     : {stats['dropped_reads']:,} ({stats['pct_dropped']:.4f}%)"]
    if stats["dropped_base_counts"]:
        L.append("Dropped-base tally:")
        for b, c in sorted(stats["dropped_base_counts"].items(), key=lambda kv: -kv[1]):
            L.append(f"    {b!r}: {c:,}")
    else:
        L.append("No reads contained ambiguous bases; nothing dropped.")
    L.append(f"Elapsed           : {stats['elapsed_sec']:.1f} s")
    L.append("=" * 60)
    text = "\n".join(L)
    print(text, file=sys.stderr)
    if log_path:
        with open(log_path, "w") as lf:
            lf.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Faster single-end FASTQ ambiguous-base filter (pigz-accelerated)."
    )
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--log", default=None)
    ap.add_argument("--allow-lowercase", action="store_true", default=True)
    ap.add_argument("--strict-uppercase", dest="allow_lowercase",
                    action="store_false")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"ERROR: input not found: {args.input}")

    print(f"Filtering (fast) {args.input} -> {args.output}", file=sys.stderr)
    if not (shutil.which("pigz") or shutil.which("gzip")):
        print("  (note: pigz/gzip not found; using Python gzip — slower)",
              file=sys.stderr)
    stats = filter_fast(args.input, args.output, log_path=args.log,
                        allow_lowercase=args.allow_lowercase)
    if stats["kept_reads"] == 0:
        sys.exit("ERROR: no reads kept — check inputs.")


if __name__ == "__main__":
    main()
