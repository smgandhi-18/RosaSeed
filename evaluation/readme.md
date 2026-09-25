# Evaluation

`evaluate_alignments.py` compares one alignment against a reference alignment
and reports mapping accuracy. The reference is either another aligner's output
(for real reads, where no truth exists) or a simulator's truth SAM (for
simulated reads, where it is absolute truth).

SAM and BAM inputs are both accepted, through `pysam`. Only primary alignments
are compared; secondary and supplementary records are skipped. Reads are keyed
on (read name, mate), so paired-end files are handled correctly.

## Usage

```bash
# real reads, scored against a reference aligner
python3 evaluate_alignments.py bwa.bam rosaseed.bam --stream

# simulated reads, scored against simulator truth
python3 evaluate_alignments.py truth.sam rosaseed.bam --stream
```

The first argument is always the reference, the second the alignment under
evaluation.

| Option | Default | Meaning |
|---|---|---|
| `--stream` | off | Compare the two files in lockstep instead of loading them into memory. Requires both to be in the same read order, which is the normal case because aligner output follows the input FASTQ. Constant memory: 21 MB instead of 20 GB on a 10M-read set, and about twice as fast. It stops with an error if the order differs, so it cannot silently mispair records. |
| `--tol50=N` | 50 | Positional tolerance in bp for standard accuracy. |
| `<integer>` | 5 | A bare integer argument sets the strict tolerance in bp. |
| `--ignore-refname` | off | Compare positions without requiring the same reference sequence name. Useful when the two files use different naming conventions (`chr1` vs `NC_060925.1`). |
| `--dump` | off | Write the diagnostic read lists and records (false negatives, false positives, test-only mappings). |
| `--outdir=DIR` | `.` | Where `--dump` writes those files. |

Without `--stream` both files are held in dictionaries, which costs roughly
2 GB per million reads; use it for anything above a few million reads.

## Metrics

Let a\* be the reference alignment of a read and a the alignment under
evaluation.

**Standard accuracy** counts a read as agreeing when a\* and a give the same
reference name and strand and their leftmost positions differ by at most the
tolerance (50 bp by default), **or** when both leave the read unmapped. It is
reported as (TP + TN) / total.

Counting "both unmapped" as agreement is deliberate: against a reference
aligner this metric measures agreement, and two aligners declining the same
unmappable read do agree. Against simulator truth the distinction never arises,
because truth maps every read. The script also prints
`Correctly mapped (+/-50 bp, same strand)`, which is the stricter reading that
excludes those reads; the two differ by exactly the reference aligner's
unmapped fraction.

**Strict accuracy** is the same comparison at a tighter tolerance (5 bp by
default).

**Structural accuracy** compares edit burden rather than position:
|E(a\*) - E(a)| / N <= 0.05, where E(r) = I + D + X and N is the read length.
That quantity is exactly the NM tag, which the script reads when present and
otherwise derives from the CIGAR. Two alignments to different copies of a
repeat with the same edit burden count as structurally equivalent, which
separates genuine misalignment from an equally plausible placement. Evaluated
over reads the test aligner mapped.

**Sequence-consistent accuracy (SCA)** asks whether an alignment is plausible
in absolute terms: NM(a) / N <= 0.10, i.e. at most 15 edits on a 150 bp read.
This catches a placement that is positionally close to the reference but of
poor quality. Evaluated over reads the test aligner mapped.

Structural accuracy and SCA answer different questions: structural is relative
to the reference alignment, SCA is absolute against the genome. An alignment
can pass one and fail the other in either direction.

**Unmapped fraction** is the share of reads the test aligner left unmapped or
omitted entirely.

## Output

All metrics are printed to stdout, including the TP/FP/FN/TN counts behind
them. With `--dump`, these files are written to `--outdir`:

| File | Contents |
|---|---|
| `FN_truth_records.sam` | reference records for reads the test aligner failed to map |
| `FN_align_records.sam` | the corresponding unmapped test records |
| `standard_fp_all_read_ids.txt` | reads disagreeing at the standard tolerance |
| `strict_fp_all_read_ids.txt` | reads disagreeing at the strict tolerance |
| `test_only_mapped_read_ids.txt` | reads the test aligner mapped and the reference did not |
| `test_only_*_records.sam` | the records behind that list |
