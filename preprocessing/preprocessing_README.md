# Read preprocessing: ambiguous-base (N) filtering

These scripts remove reads containing ambiguous or non-canonical bases: any
base outside `{A, C, G, T}`, which includes `N` and all IUPAC ambiguity codes.

**Filtering is optional.** RosaSeed and RosaSeed-Compact both handle ambiguous
bases directly: seed extension stops at an `N` and resumes past it, as the
BWA-MEM2 SMEM search does, so unfiltered reads align correctly. Reads of
differing lengths in one file are also fine; no trimming is required. See
"Read handling" in the top-level README.

The scripts exist for benchmarking. The **same filtered read set is used for
every aligner compared** (RosaSeed, BWA-MEM2, Minimap2, Bowtie2, ERT, ERT2,
minibwa, Strobealign), so every tool receives byte-identical input and no result
can be attributed to one aligner handling ambiguous bases differently from
another. They are published so the exact read sets reported in the paper can be
reproduced from the public raw data.

## Scripts

| Script | Mode | Dependencies | Use for |
|--------|------|--------------|---------|
| `filter_n_reads.py` | single-end | Python 3 only | reference implementation, single-end |
| `filter_n_reads_paired.py` | paired-end | Python 3 only | reference implementation, paired-end |
| `filter_n_reads_fast.py` | single-end | Python 3 + `pigz`/`gzip` | fast processing, single-end |
| `filter_n_reads_paired_fast.py` | paired-end | Python 3 + `pigz`/`gzip` | fast processing, paired-end |

The `*_fast.py` variants use `pigz` (or `gzip`) as external (de)compression
subprocesses and a byte-translation membership test; they produce **identical
output** to the pure-Python reference scripts but run roughly 80–90x faster on
gzip-compressed input. The pure-Python scripts have no external dependencies and
are the canonical reference; the fast scripts are a convenience for large
datasets.

## Filtering rule

- **Single-end:** a read is discarded if its sequence contains any base outside
  the accepted alphabet.
- **Paired-end:** a read *pair* is discarded from **both** output files if
  **either** mate contains a disallowed base. Dropping only one mate would
  desynchronize R1 and R2; keeping pairs intact preserves correct mate pairing.

By default the accepted alphabet is `{A, C, G, T, a, c, g, t}`: lowercase
soft-masked bases are kept, which is what `--allow-lowercase` selects and is the
default. Pass `--strict-uppercase` to accept only `{A, C, G, T}`, which drops
any read containing a soft-masked base.

## Usage

Single-end:

```bash
python3 filter_n_reads.py  input.fastq.gz  output.filtered.fastq.gz  --log run.log
# fast:
python3 filter_n_reads_fast.py  input.fastq.gz  output.filtered.fastq.gz  --log run.log
```

Paired-end:

```bash
python3 filter_n_reads_paired.py \
    R1.fastq.gz R2.fastq.gz \
    R1.filtered.fastq.gz R2.filtered.fastq.gz \
    --log run.log
# fast:
python3 filter_n_reads_paired_fast.py \
    R1.fastq.gz R2.fastq.gz \
    R1.filtered.fastq.gz R2.filtered.fastq.gz \
    --log run.log
```

Input and output may be plain FASTQ (`.fastq`/`.fq`) or gzip-compressed
(`.fastq.gz`/`.fq.gz`); output compression is chosen from the output filename.

## Output

Each run writes the filtered FASTQ file(s) and a summary log reporting the total
number of reads/pairs, the number kept and dropped, the drop percentage, and a
tally of the ambiguous bases removed. For paired runs the log additionally
reports how many pairs were dropped because R1 only, R2 only, or both mates
contained a disallowed base.

## Reproducibility note

The reads counts reported in the paper reflect the datasets after this filtering.
On typical Illumina data the fraction removed is small (for example, on one
HG002 2×250 bp lane chunk of 4,000,000 pairs, 4,795 pairs (0.12%) were removed),
so the filtering does not materially alter the datasets and is applied uniformly
across all aligners.
