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
filter_n_reads_paired_fast.py — Fast PAIRED-END FASTQ ambiguous-base filter.

Same behavior as filter_n_reads_paired.py: a read PAIR is dropped from BOTH
output files if EITHER mate contains a base outside {A,C,G,T[,acgt]}, keeping R1
and R2 synchronized. Accelerated the same way as filter_n_reads_fast.py:
  * pigz/gzip subprocesses for (de)compression (much faster than Python gzip)
  * fast bytes-translation membership check instead of a per-read Python set

Falls back to Python gzip if pigz/gzip binaries are unavailable.

Usage:
    python3 filter_n_reads_paired_fast.py R1.fastq.gz R2.fastq.gz \\
            R1.filtered.fastq.gz R2.filtered.fastq.gz --log run.log

    python3 filter_n_reads_paired_fast.py in_1.fq in_2.fq out_1.fq out_2.fq \\
            --strict-uppercase
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path


def _decompress_proc(path):
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


def _read_record(readline, idx):
    header = readline()
    if not header:
        return None
    seq = readline()
    plus = readline()
    qual = readline()
    if not qual:
        raise ValueError(f"Truncated FASTQ record near read {idx}.")
    return header, seq, plus, qual


def filter_paired_fast(r1_in, r2_in, r1_out, r2_out, log_path=None,
                       allow_lowercase=True, report_every=1_000_000):
    allowed = b"ACGT" + (b"acgt" if allow_lowercase else b"")
    delete_allowed = bytes(b for b in range(256) if b in allowed)

    n_pairs = n_kept = n_dropped = 0
    d_r1 = d_r2 = d_both = 0
    dropped_base_counts = {}
    t0 = time.time()

    f1i, p1i = _decompress_proc(r1_in)
    f2i, p2i = _decompress_proc(r2_in)
    f1o, p1o = _compress_proc(r1_out)
    f2o, p2o = _compress_proc(r2_out)

    rl1, rl2 = f1i.readline, f2i.readline
    w1, w2 = f1o.write, f2o.write

    try:
        while True:
            rec1 = _read_record(rl1, n_pairs + 1)
            rec2 = _read_record(rl2, n_pairs + 1)

            if rec1 is None and rec2 is None:
                break
            if (rec1 is None) != (rec2 is None):
                raise ValueError(
                    f"R1 and R2 have different read counts "
                    f"(desynchronized near pair {n_pairs + 1}). "
                    f"Inputs must be matched mates in the same order."
                )

            n_pairs += 1
            leftover1 = rec1[1].rstrip(b"\n").translate(None, delete_allowed)
            leftover2 = rec2[1].rstrip(b"\n").translate(None, delete_allowed)

            if leftover1 or leftover2:
                n_dropped += 1
                if leftover1 and leftover2:
                    d_both += 1
                elif leftover1:
                    d_r1 += 1
                else:
                    d_r2 += 1
                for lo in (leftover1, leftover2):
                    for b in set(lo):
                        ch = chr(b)
                        dropped_base_counts[ch] = (
                            dropped_base_counts.get(ch, 0) + lo.count(b)
                        )
            else:
                w1(rec1[0]); w1(rec1[1]); w1(rec1[2]); w1(rec1[3])
                w2(rec2[0]); w2(rec2[1]); w2(rec2[2]); w2(rec2[3])
                n_kept += 1

            if report_every and n_pairs % report_every == 0:
                el = time.time() - t0
                rate = n_pairs / el if el > 0 else 0
                print(f"  ... {n_pairs:,} pairs ({n_kept:,} kept, "
                      f"{n_dropped:,} dropped) [{rate:,.0f} pairs/s]",
                      file=sys.stderr)
    finally:
        for fh in (f1i, f2i):
            try: fh.close()
            except Exception: pass
        for fh in (f1o, f2o):
            try: fh.close()
            except Exception: pass
        for pin in (p1i, p2i):
            if pin: pin.wait()
        for pout in (p1o, p2o):
            if pout:
                proc, fh = pout
                proc.wait(); fh.close()

    el = time.time() - t0
    pct = (100.0 * n_dropped / n_pairs) if n_pairs else 0.0
    stats = {
        "r1_in": str(r1_in), "r2_in": str(r2_in),
        "r1_out": str(r1_out), "r2_out": str(r2_out),
        "total_pairs": n_pairs, "kept_pairs": n_kept,
        "dropped_pairs": n_dropped, "pct_dropped": pct,
        "dropped_because_r1_only": d_r1,
        "dropped_because_r2_only": d_r2,
        "dropped_because_both": d_both,
        "dropped_base_counts": dropped_base_counts,
        "elapsed_sec": el,
        "allowed_alphabet": allowed.decode(),
    }
    _emit_log(stats, log_path)
    return stats


def _emit_log(stats, log_path):
    L = ["=" * 66,
         "PAIRED-END FASTQ ambiguous-base filtering — summary (fast)",
         "=" * 66,
         f"R1 input          : {stats['r1_in']}",
         f"R2 input          : {stats['r2_in']}",
         f"R1 output         : {stats['r1_out']}",
         f"R2 output         : {stats['r2_out']}",
         f"Allowed alphabet  : {{{stats['allowed_alphabet']}}}",
         f"Total pairs       : {stats['total_pairs']:,}",
         f"Pairs kept        : {stats['kept_pairs']:,}",
         f"Pairs dropped     : {stats['dropped_pairs']:,} ({stats['pct_dropped']:.4f}%)",
         f"    dropped (R1 bad only) : {stats['dropped_because_r1_only']:,}",
         f"    dropped (R2 bad only) : {stats['dropped_because_r2_only']:,}",
         f"    dropped (both bad)    : {stats['dropped_because_both']:,}"]
    if stats["dropped_base_counts"]:
        L.append("Dropped-base tally (across dropped mates):")
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
        description="Fast paired-end FASTQ ambiguous-base filter (pigz-accelerated)."
    )
    ap.add_argument("r1_in")
    ap.add_argument("r2_in")
    ap.add_argument("r1_out")
    ap.add_argument("r2_out")
    ap.add_argument("--log", default=None)
    ap.add_argument("--allow-lowercase", action="store_true", default=True)
    ap.add_argument("--strict-uppercase", dest="allow_lowercase",
                    action="store_false")
    args = ap.parse_args()

    for p in (args.r1_in, args.r2_in):
        if not Path(p).exists():
            sys.exit(f"ERROR: input not found: {p}")

    print(f"Filtering pair (fast):\n  {args.r1_in}\n  {args.r2_in}", file=sys.stderr)
    if not (shutil.which("pigz") or shutil.which("gzip")):
        print("  (note: pigz/gzip not found; using Python gzip — slower)",
              file=sys.stderr)
    stats = filter_paired_fast(
        args.r1_in, args.r2_in, args.r1_out, args.r2_out,
        log_path=args.log, allow_lowercase=args.allow_lowercase,
    )
    if stats["kept_pairs"] == 0:
        sys.exit("ERROR: no pairs kept — check inputs.")


if __name__ == "__main__":
    main()
