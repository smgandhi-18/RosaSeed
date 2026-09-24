# RosaSeed

[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Paper](https://img.shields.io/badge/Paper-bioRxiv-red.svg)](https://doi.org/...)

![Platform](https://img.shields.io/badge/Platform-Linux-blue.svg)
![C](https://img.shields.io/badge/C-C11-blue.svg)
![C++](https://img.shields.io/badge/C%2B%2B-C%2B%2B14-blue.svg)
![Python](https://img.shields.io/badge/Python-3.x-yellow.svg)

**RosaSeed** is a configurable, fast, and accurate seeding engine for short-read DNA
sequence alignment. It replaces the BWA-MEM2 seeding kernel and feeds
candidate seeds directly into the existing BWA-MEM2 chaining and alignment
extension pipeline, producing standard SAM output compatible with downstream
variant analysis tools.

> See [README-bwamem2.md](README-bwamem2.md) for the original BWA-MEM2 documentation.

---

## Repository structure

```
RosaSeed/
├── src/
│   ├── rosaseed/                   RosaSeed (2-step) seeding kernel
│   │   ├── rosaseed_phaseA.c       Phase A: jump table + backward search
│   │   ├── rosaseed_phaseB.c       Phase B: fixed-pivot supplementary seeding
│   │   ├── rosaseed_gapfill_PhaseC.c  Phase C: gap-directed seeding
│   │   ├── load_data_mmap.c        index loading (mmap)
│   │   ├── mem_alloc.c, helper_functions.c, bwa.c
│   │   └── file_dec.h, macros.h, ...
│   ├── rosaseed_compact/           RosaSeed-Compact (1-step) seeding kernel
│   │   ├── rosaseed_compact_phaseA.c
│   │   ├── rosaseed_compact_phaseB.c
│   │   ├── rosaseed_compact_gapfill_PhaseC.c
│   │   ├── load_data_compact_mmap.c, mem_alloc_compact.c
│   │   └── file_dec.h, macros.h, ...
│   ├── rosaseed_core_bridge.cpp    2-step kernel to BWA-MEM2 bridge
│   ├── rosaseed_core_bridge_compact.cpp   compact kernel bridge
│   ├── rosaseed_inprocess.cpp      seeding entry point used by bwamem.cpp
│   ├── rosaseed_inprocess_compact.cpp
│   └── [bwa-mem2 source files]     chaining, extension, SAM output
├── index-builder/                  FM-index build pipelines
│   ├── build_2step_pipeline.sh     2-step (base-16) index
│   ├── build_compact_pipeline.sh   compact (radix-4) index
│   ├── rss_monitor.py              RAM tracker used by the pipelines
│   ├── src/                        index-builder tools (C)
│   └── README.md
├── evaluation/                     alignment accuracy evaluation
│   ├── evaluate_alignments.py
│   └── readme.md
├── preprocessing/                  optional FASTQ preprocessing scripts
├── test/                           inherited BWA-MEM2 unit tests
├── ext/safestringlib/              inherited BWA-MEM2 dependency
├── images/
├── Makefile
├── LICENSE
├── NEWS.md
├── README.md
└── README-bwamem2.md
```

---

## Requirements

| Requirement | Notes |
|---|---|
| g++ with C++14 support | Tested with GCC 15.2; the sources also compile as C++11. `src/rosaseed/*.c` are compiled with `g++` (`-fpermissive`), so a C++ compiler is required, not just a C compiler |
| GNU Make | Standard build |
| zlib development headers | Linked with `-lz` (`zlib1g-dev` on Ubuntu) |
| x86-64 CPU, SSE4.1 or newer | AVX2 (or AVX-512) is strongly recommended for alignment speed. `arch=native` builds for the CPU it is compiled on and does not itself require AVX2 |
| Linux | Uses `mmap` and pthreads |
| RAM: ~50 GB (alignment) | ~49.6 GB combined peak at alignment time with the recommended configuration. The [Memory-efficient](#memory-efficient) build (`-DSA_COMPRESSION_FACTOR_POWER=3`) cuts index memory to ~35 GB and suits 48 GB systems, see [Genome size limits](#genome-size-limits) |
| RAM: ~64 GB (index building) | Only needed to build indexes: ~56 GB for the RosaSeed index, ~62 GB for `bwa-mem2 index` (human genome) |

**Install on Ubuntu/Debian:**

```bash
sudo apt install gcc g++ make git zlib1g-dev
```

**Index builder additional requirements:**

```bash
sudo apt install python3 python3-numpy
```

---

## Quick start

### 1. Clone

```bash
git clone https://github.com/smgandhi-18/RosaSeed.git
cd RosaSeed
```

### 2. Build RosaSeed

```bash
make clean
make arch=native CXX=g++ ROSASEED=1 \
  CPPFLAGS_EXTRA=" \
    -DLOAD_JTABLE_14nt \
    -DENABLE_F_RC_CHOICE \
    -DSA_COMPRESSION_FACTOR_POWER=1 \
    -DGAPFILL_ALWAYS_RUN_BOTH_STRANDS \
    -DGAPFILL_EARLY_EXIT \
    -DROSASEED_PRECHAIN_SINGLETON_SUPPRESS \
    -DROSASEED_PRECHAIN_TRIGGER=200 \
    -DROSASEED_PRECHAIN_WEAK_LEN=60"
```

### 3. Build the BWA-MEM2 index

RosaSeed replaces the seeding stage and modifies the chaining stage with
pruning heuristics. The alignment extension (Smith-Waterman) pipeline still
uses BWA-MEM2's own index files. Generate them with:

```bash
./bwa-mem2 index /path/to/genome.fna
```

This produces the following files alongside the FASTA:
```
genome.fna.0123
genome.fna.amb
genome.fna.ann
genome.fna.bwt.2bit.64
genome.fna.pac
```

BWA-MEM2 automatically locates these files by appending extensions to the
FASTA path supplied at alignment time - no extra flags needed.

### 4. Build the RosaSeed index

```bash
cd index-builder
chmod +x build_2step_pipeline.sh
./build_2step_pipeline.sh /path/to/genome.fna
cd ..
```

See [index-builder/README.md](index-builder/README.md) for full options.
Index files are written to `index-builder/index/<genome_name>/` by default
(~42.7 GB for the T2T human genome with 14 and 15-mer jump tables).
Build time: ~44 minutes, peak RAM ~56 GB.

### 5. Run alignment

```bash
./bwa-mem2 mem \
    -t 20 \
    -k 19 \
    --rs-index index-builder/index/<genome_name>/ \
    --rs-cap 2000 \
    /path/to/genome.fna \
    reads_R1.fastq.gz \
    > output.sam
```

> **Note:** The FASTA path supplied to `./bwa-mem2 mem` must be the same
> one used in Step 3. BWA-MEM2 looks for its index files (`.0123`, `.pac`
> etc.) in the same directory as the FASTA.

---

### 6. Evaluate accuracy (optional)

`evaluation/evaluate_alignments.py` compares an alignment against a reference
alignment (for real reads, typically BWA-MEM2) or against simulator truth, and
reports standard, structural and sequence-consistent accuracy plus the unmapped
fraction.

```bash
python3 evaluation/evaluate_alignments.py reference.bam test.bam --stream
```

| Option | Meaning |
|---|---|
| `--stream` | both files are in the same read order (the usual case: aligner output follows the FASTQ). Constant memory, required for 100M+ read sets; without it both files are loaded into RAM |
| `--tol50=N` | positional tolerance for standard accuracy (default 50 bp) |
| `--ignore-refname` | compare positions without requiring the same reference name |
| `--dump --outdir=DIR` | write the diagnostic false-negative / false-positive read lists and records (off by default) |

Both SAM and BAM inputs are accepted. Only primary alignments are considered;
secondary and supplementary records are skipped.

---

## Read handling

| Input property | Behaviour |
|---|---|
| Ambiguous bases (`N`) in reads | Supported. Seed extension stops at an `N` and resumes past it, the same way the BWA-MEM2 SMEM search does. No preprocessing is required; the scripts in `preprocessing/` are optional. |
| Variable read lengths in one FASTQ | Supported. Reads do not need to be trimmed to a uniform length. |
| Maximum read length | 251 bp for both RosaSeed and RosaSeed-Compact (`RS_BATCH_MAX_READ_LEN`, which sizes the per-read slot arrays). A longer read makes the seeder exit with an error naming the limit, rather than truncating silently. Raise it at compile time with `-DRS_BATCH_MAX_READ_LEN=<n>`; the cost is a few KB per thread. Standard 100/150/250 bp Illumina reads need no change. |
| Non-ACGT bases in the *reference* | Replaced with `A` by the index builder (`preprocess_genome`), preserving genome coordinates. |

---

## Index summary

RosaSeed requires two separate indexes:

| Index | Built by | Location | Size (T2T human) |
|---|---|---|---|
| BWA-MEM2 index | `./bwa-mem2 index genome.fna` | Same directory as FASTA | ~6.9 GB |
| RosaSeed index | `index-builder/build_2step_pipeline.sh` | `index-builder/index/<name>/` | ~42.7 GB |
| RosaSeed-Compact index | `index-builder/build_compact_pipeline.sh` | `index-builder/index/<name>_compact/` | see [RosaSeed-Compact](#rosaseed-compact) |

**Combined peak memory at alignment time: ~49.61 GB**

---

## RosaSeed-Compact

RosaSeed-Compact pairs conventional 1-base (radix-4) FM-index traversal with
8x suffix-array compression. It is intended for machines where memory is the
main constraint: in the preprint benchmark it needs 25.86 GB peak (about 48%
less than RosaSeed) at 44.34 µs/read, with 98.61% standard accuracy.

It is selected at build time and uses its own index. The compact index builder
takes any genome FASTA as downloaded, like the 2-step one (non-ACGT bases are
replaced with A), and builds every file needed for any jump-table size and SA
compression factor:

```bash
make clean
make arch=native CXX=g++ ROSASEED=1 ROSASEED_1STEP=1 \
  CPPFLAGS_EXTRA=" \
    -DLOAD_JTABLE_15nt \
    -DENABLE_F_RC_CHOICE \
    -DSA_COMPRESSION_FACTOR_POWER=3 \
    -DGAPFILL_EARLY_EXIT \
    -DROSASEED_PRECHAIN_SINGLETON_SUPPRESS \
    -DROSASEED_PRECHAIN_TRIGGER=200 \
    -DROSASEED_PRECHAIN_WEAK_LEN=60 \
    -DROSASEED_PRECHAIN_USE_ABUNDANCE \
    -DROSASEED_PRECHAIN_ABUNDANCE=500"

cd index-builder && ./build_compact_pipeline.sh /path/to/genome.fna && cd ..

./bwa-mem2 mem -t 20 -k 19 --rs-index index-builder/index/<genome_name>_compact/ \
    --rs-cap 5000 /path/to/genome.fna reads.fastq.gz > output.sam
```

This is the recommended RosaSeed-Compact configuration: a 15-nt jump table, 8x SA
compression (`-DSA_COMPRESSION_FACTOR_POWER=3`), gap fill with early exit but
without `-DGAPFILL_ALWAYS_RUN_BOTH_STRANDS`, the abundance-aware pre-chain
filter (SSF+A, threshold 500) and a Phase A interval cap of 5000.
RosaSeed-Compact and RosaSeed are different seeders, so their SAM output is not
expected to be identical.

The BWA-MEM2 index (Step 3 of Quick start) is still required. The jump-table
flag must match a table built with `-j`, and `-DSA_COMPRESSION_FACTOR_POWER`
selects which `sa_*_cf*.bin` pair is loaded, as for RosaSeed. The runtime
flags (`--rs-*`, `-k`, `-t`) are the same as for RosaSeed, including the 251 bp
read-length limit (see [Read handling](#read-handling)) and the 2^32 bp genome
limit.

> **Note:** `make clean` is required before switching between RosaSeed
> (2-step) and RosaSeed-Compact builds: object files are not rebuilt just
> because the flags changed. Also, `ROSASEED_1STEP=1` only takes effect when
> passed together with `ROSASEED=1`; on its own it selects nothing.

---

## Genome size limits

The RosaSeed index format stores BWT positions in 33-bit fields, giving a hard
limit of **4,294,967,296 bp (2³², ≈4.29 Gbp)** for the reference genome.

| Genome | Size | Field used | Supported |
|---|---|---|---|
| Human T2T-CHM13v2 | 3.12 Gbp | 72.6% | ✅ |
| Human GRCh38 | 3.10 Gbp | 72.2% | ✅ |
| Mouse GRCm39 | 2.70 Gbp | 62.9% | ✅ |
| Maize B73 | 2.30 Gbp | 53.6% | ✅ |
| Wheat IWGSC | 17.0 Gbp | — | ❌ |
| Axolotl | 32.0 Gbp | — | ❌ |

Two further limits are enforced at index-build time:

- **Per-symbol occurrence count** must fit in 32 bits. This is composition
  dependent and roughly 22 Gbp for a human-like base distribution, so it is not
  the binding constraint in practice.
- **No k-mer may occur more than 536,870,911 times.** The jump-table builder
  reports the observed maximum for each table, so the headroom is visible on
  every run.

All limits are checked explicitly. The pipeline exits with an error naming the
offending value rather than producing incorrect coordinates, and the aligner
re-checks at startup so an index built elsewhere cannot be used silently.

---

## Flags reference

RosaSeed is controlled through two mechanisms:
- **Build-time flags** (`-D` flags via `CPPFLAGS_EXTRA`): select algorithms and data structures at compile time
- **Runtime flags** (`--rs-*` and standard BWA-MEM2 flags): tune thresholds at run time

---

### Build-time flags (`CPPFLAGS_EXTRA`)

#### Jump table size

Controls which pre-computed k-mer jump table is loaded into memory.
The jump table maps a k-nucleotide pattern directly to its BWT interval,
allowing Phase A seeding to begin without an FM-index backward search.
Exactly **one** of these must be set per build.

| Flag | k-mer | Entries | Memory |
|------|-------|---------|--------|
| `-DLOAD_JTABLE_14nt` | 14 nt (7 base-16 pairs) | 268M | 2 GB |
| `-DLOAD_JTABLE_15nt` | 15 nt (7 pairs + 1 NT) | 1.07B | 8 GB |
| `-DLOAD_JTABLE_16nt` | 16 nt (8 base-16 pairs) | 4.29B | 32 GB |

The jump table must be built to match the flag:
`index-builder/build_2step_pipeline.sh -j 14` (default), `-j 15`, or `-j 16`.

---

#### SA Compression Factor (CF)

Controls how densely the suffix array (SA) is sampled on disk and in memory.
Higher compression = less memory, slower non-unique seed position lookup
(requires LF-mapping walk to reach the nearest sampled entry).

| Flag | CF | SA files used | Memory (T2T) | Notes |
|---|---|---|---|---|
| `-DSA_COMPRESSION_FACTOR_POWER=0` | 1 | `sa_*_cf1.bin` | ~29 GB | Full SA. Direct O(1) lookup. |
| `-DSA_COMPRESSION_FACTOR_POWER=1` | 2 | `sa_*_cf2.bin` | ~15 GB | Recommended. At most 1 LF step. |
| `-DSA_COMPRESSION_FACTOR_POWER=2` | 4 | `sa_*_cf4.bin` | ~7.3 GB | At most 3 LF steps. |
| `-DSA_COMPRESSION_FACTOR_POWER=3` | 8 | `sa_*_cf8.bin` | ~3.6 GB | Lowest memory. At most 7 LF steps. |

All CF files are generated automatically by the index builder regardless of
this setting. Only the matching pair is loaded at runtime.

---

#### Phase A: strand selection

| Flag | Description |
|---|---|
| `-DENABLE_F_RC_CHOICE` | At the very first pivot of every read, look up the jump table for **both** the forward and reverse-complement read pattern and select whichever gives the smaller (more specific) BWT interval. Costs one extra jump table lookup per pivot but improves seed specificity and is recommended. Omitting this flag always seeds on the forward strand. |

---

#### Supplementary seeding strategy selection

RosaSeed has two modes for finding seeds missed by Phase A. Only one is active per build.

| Flag | Mode | Description |
|---|---|---|
| *(not set)* | **Phase C** (default) | Gap-aware adaptive recovery: sorts Phase A seeds by read position, finds uncovered gaps, and runs targeted searches inside those gaps. |
| `-DENABLE_PHASE_II` | **Phase B** | Fixed-pivot supplementary seeding: places `g_num_pivots_B` evenly-spaced pivots across the read and runs full backward searches from each. Simpler but less targeted than Phase C. |

---

#### Phase C: gap-fill configuration

These flags apply when `ENABLE_PHASE_II` is **not**  set (i.e., the default Phase C path). 
They are ignored when `ENABLE_PHASE_II` is enabled.

| Flag | Description |
|---|---|
| `-DGAPFILL_ALWAYS_RUN_BOTH_STRANDS` | Inside each detected coverage gap, run seed searches on **both** forward and RC strands, not just the strand chosen by Phase A. Increases recall in repetitive or ambiguous regions at a modest compute cost. |
| `-DGAPFILL_EARLY_EXIT` | Stop gap-filling as soon as a seed of sufficient length and specificity is found covering the gap. Without this flag, all pivot positions in the gap are tried exhaustively. Recommended: improves speed with negligible accuracy impact. |

---

#### Pre-chain filtering

After Phase A and the selected supplementary seeding stage (Phase B or
Phase C), RosaSeed converts the generated hits into BWA-MEM2-compatible
seed records and passes them to chain construction. The optional
pre-chain singleton suppression filter prevents selected weak candidates
from initiating new singleton chains, thereby reducing downstream
chaining and alignment-extension work.

Two filtering modes are supported:

* **SSF** suppresses weak singleton-chain candidates based on the number of
  chains already created for the read and the candidate seed length.
* **SSF+A** additionally considers seed abundance, retaining short,
  relatively specific candidates while preferentially suppressing short,
  highly repetitive candidates.

| Flag                                     | Description                                                                                                                                                                                                                                                                 |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-DROSASEED_PRECHAIN_SINGLETON_SUPPRESS` | Enables the pre-chain singleton suppression filter (SSF). The filter is evaluated only when a candidate seed cannot be merged into an existing chain and would otherwise initiate a new singleton chain.                                                                    |
| `-DROSASEED_PRECHAIN_TRIGGER=N`          | Existing-chain-count threshold above which singleton suppression becomes active. Suppression is applied when the number of chains already created for the read is strictly greater than `N`. Default: `200`. Lower values activate suppression earlier.                     |
| `-DROSASEED_PRECHAIN_WEAK_LEN=N`         | Seed-length threshold used to classify weak singleton candidates. Candidates shorter than `N` nt may be suppressed when the other filtering conditions are satisfied. Default: `60`.                                                                                        |
| `-DROSASEED_PRECHAIN_USE_ABUNDANCE`      | Enables the abundance-aware variant (SSF+A). This option is meaningful only when `ROSASEED_PRECHAIN_SINGLETON_SUPPRESS` is enabled. A weak singleton candidate is suppressed only when its occurrence count is also greater than or equal to `ROSASEED_PRECHAIN_ABUNDANCE`. |
| `-DROSASEED_PRECHAIN_ABUNDANCE=N`        | Seed-abundance threshold used by SSF+A. Weak singleton candidates with occurrence count `>= N` are treated as highly repetitive and suppressed. Default: `500`.                                                                                                             |

---

#### Read batching and length

| Flag | Default | Description |
|---|---|---|
| `-DRS_BATCH_MAX_READ_LEN=N` | 251 | Maximum read length the seeder accepts, in both variants. It sizes the fixed per-read arrays in the batch slot state, so it is an allocation bound: a longer read makes the seeder stop with an error naming the limit. Raising it costs roughly `2 * RS_BATCH * N` bytes per thread. |
| `-DRS_BATCH=N` | 32 | Number of reads the Phase A interleaver processes at once. Larger batches give the prefetcher more independent memory accesses to overlap, at the cost of slot-state memory. |

---

### Runtime flags

These are passed on the command line at alignment time and can be tuned
without recompilation.

#### RosaSeed-specific runtime flags

| Flag | Default | Description |
|---|---|---|
| `--rs-index <dir>` | *(required)* | Path to the RosaSeed index directory produced by `build_2step_pipeline.sh`. Must contain `cp_occ_full.bin`, `c_vector.txt`, `ref16_packed.bin`, the SA split files, and the jump table. A RosaSeed-Compact build (`ROSASEED_1STEP=1`) instead expects a directory produced by `build_compact_pipeline.sh`, containing `cp_occ_compact.bin`, `ref4_packed.bin`, the SA split files, and the jump table. See [RosaSeed-Compact](#rosaseed-compact). |
| `--rs-cap <int>` | 2000 | Phase A SA interval cap. Seeds whose BWT interval width exceeds this value are skipped in Phase A (too repetitive to be useful). Lower values run faster but may miss seeds in repetitive regions. Recommended: 2000 for standard use, 200 for miniRosaSeed. |

#### Inherited BWA-MEM2 runtime flags (relevant to RosaSeed)

| Flag | Default | Description |
|---|---|---|
| `-t <int>` | 1 | Number of alignment threads. RosaSeed is fully multi-threaded - each thread runs its own Phase A/C pipeline independently. |
| `-k <int>` | 19 | Minimum seed length in nucleotides. Seeds shorter than this are not emitted by Phase A. Phase B/C use `k+1` as their threshold. |

---

## Recommended configurations

### Standard (recommended default)

Best for most use cases. Balances accuracy and speed.

```bash
make arch=native CXX=g++ ROSASEED=1 \
  CPPFLAGS_EXTRA=" \
    -DLOAD_JTABLE_14nt \
    -DENABLE_F_RC_CHOICE \
    -DSA_COMPRESSION_FACTOR_POWER=1 \
    -DGAPFILL_ALWAYS_RUN_BOTH_STRANDS \
    -DGAPFILL_EARLY_EXIT \
    -DROSASEED_PRECHAIN_SINGLETON_SUPPRESS \
    -DROSASEED_PRECHAIN_TRIGGER=200 \
    -DROSASEED_PRECHAIN_WEAK_LEN=60"

./bwa-mem2 mem -t 20 -k 19 --rs-index index-builder/index/<name>/ \
    --rs-cap 2000 genome.fna reads.fastq.gz > output.sam
```

### Memory-efficient

Reduces index memory using 8× SA compression.
Suitable for systems with 48 GB RAM.

```bash
make arch=native CXX=g++ ROSASEED=1 \
  CPPFLAGS_EXTRA=" \
    -DLOAD_JTABLE_14nt \
    -DENABLE_F_RC_CHOICE \
    -DSA_COMPRESSION_FACTOR_POWER=3 \
    -DGAPFILL_ALWAYS_RUN_BOTH_STRANDS \
    -DGAPFILL_EARLY_EXIT \
    -DROSASEED_PRECHAIN_SINGLETON_SUPPRESS \
    -DROSASEED_PRECHAIN_TRIGGER=200 \
    -DROSASEED_PRECHAIN_WEAK_LEN=60"
```

### miniRosaSeed

Speed-optimised configuration, benchmarked against minibwa. A loose interval
cap combined with an aggressive pre-chain singleton filter: the cap lets enough
seeds through for accuracy, and the filter removes the repetitive singletons
that would otherwise cost chaining and extension time.

```bash
make arch=native CXX=g++ ROSASEED=1 \
  CPPFLAGS_EXTRA=" \
    -DLOAD_JTABLE_14nt \
    -DENABLE_F_RC_CHOICE \
    -DSA_COMPRESSION_FACTOR_POWER=1 \
    -DGAPFILL_ALWAYS_RUN_BOTH_STRANDS \
    -DGAPFILL_EARLY_EXIT \
    -DROSASEED_PRECHAIN_SINGLETON_SUPPRESS \
    -DROSASEED_PRECHAIN_TRIGGER=25 \
    -DROSASEED_PRECHAIN_WEAK_LEN=80"

./bwa-mem2 mem -t 1 -k 19 --rs-index index-builder/index/<name>/ \
    --rs-cap 200 genome.fna reads.fastq.gz > output.sam
```

This configuration replaces the earlier one (`--rs-cap 50`,
`-DROSASEED_PRECHAIN_TRIGGER=50`, `-DROSASEED_PRECHAIN_WEAK_LEN=60`) used in the
preprint. On the 20M-read ERR dataset, single-threaded, it is 2.8% faster
(12.96 vs 13.33 µs/read, mean of 3 runs) while reaching 97.345% standard
accuracy instead of 96.947% and leaving 1.62% of reads unmapped instead of
2.09%. `--rs-cap 100` with the same filter settings is marginally faster again
(12.79 µs/read) at 97.057% accuracy. Neither configuration uses the
abundance-aware filter.

### Compact

Lowest memory (25.86 GB peak in the preprint benchmark): 15-nt jump table,
`-DSA_COMPRESSION_FACTOR_POWER=3`, SSF+A (`-DROSASEED_PRECHAIN_USE_ABUNDANCE`,
`-DROSASEED_PRECHAIN_ABUNDANCE=500`), gap fill with `-DGAPFILL_EARLY_EXIT` only
(no `-DGAPFILL_ALWAYS_RUN_BOTH_STRANDS`) and `--rs-cap 5000`. Full build and run
commands are in [RosaSeed-Compact](#rosaseed-compact).

---
## Key results

Compared against BWA-MEM2, Minimap2, Bowtie2, and the Enumerated Radix
Tree (ERT) algorithm on simulated and real short-read datasets:

| Comparison | Seed processing speedup | Total alignment speedup |
|---|---|---|
| vs BWA-MEM2 (SA cf=8) | **15.0×** | **4.04×** |
| vs ERT | 4.0× | 2.48× |
| vs ERT2 | 1.2× | 1.44× |
| vs Minimap2 | 4.8× | 2.48× |
| vs Bowtie2 | 33.2× | 10.8× |

Peak memory: **49.61 GB**, ~25% less than ERT/ERT2 (~66.3 GB each).
This comprises the RosaSeed 2-step FM-index (~42.7 GB) and the BWA-MEM2
index files (~6.9 GB) required by the downstream chaining and alignment
extension pipeline.

On a separate AMD Zen3 workstation, the speed-optimised
**miniRosaSeed** configuration achieves a **2.08× end-to-end speedup**
over minibwa at single-thread execution while attaining
**2.581 percentage points higher standard accuracy**, remaining faster
and more accurate simultaneously up to ~26 threads. Those figures are from
the preprint, measured with the earlier miniRosaSeed settings; the
configuration shipped here (see [miniRosaSeed](#minirosaseed)) is faster and
more accurate than those settings on the ERR dataset, so both margins widen.

RosaSeed exposes multiple runtime-memory operating configurations
spanning memory-efficient and high-performance designs, enabling users
to select operating points appropriate for diverse computational
environments while maintaining alignment accuracy comparable to BWA-MEM2.

---

### Downstream pipeline (BWA-MEM2)

After seeding, RosaSeed converts seed positions from s-step reference
coordinates back to BWA-MEM2's global coordinate space and hands them
to the standard BWA-MEM2 chaining and alignment pipeline:

- **Chaining** (`mem_chain_preexpanded_hits`): uses `bns` (chromosome table from BWA-MEM2 index) to assign seeds to chromosomes and group them into chains
- **Chain filtering** (`mem_flt_chained_seeds`): uses `pac` (packed reference sequence from BWA-MEM2 index) for Smith-Waterman scoring of short seeds
- **Alignment extension** (`mem_chain2aln`, `bwa_gen_cigar2`): uses `pac` and `bns` to produce CIGAR strings
- **SAM output** (`mem_reg2sam`, `mem_aln2sam`): uses `bns->anns[].name` for chromosome names

This is why both the RosaSeed index (`--rs-index`) and the BWA-MEM2 index
(`.0123`, `.pac` files) are required.

---

## Citation

If you use RosaSeed in your research, please cite:

> Gandhi Shyama, Cockburn Bruce (2026). RosaSeed: A Configurable Seeding Framework 
> for Fast and Accurate Short-Read Alignment, 2026, 
> University of Alberta. URL: 

Please also cite the original BWA-MEM2:

> Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient
> Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.
> IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019.
> URL: https://ieeexplore.ieee.org/document/8820962 

If you use gsufsort for index building, please cite:

> Louza, F.A., Telles, G.P., Gog, S., Prezza, N., Rosone, G.
> gsufsort: constructing suffix arrays, LCP arrays and BWTs for string collections.
> Algorithms Mol Biol 15, 18 (2020).
> URL: https://link.springer.com/article/10.1186/s13015-020-00177-y

---

## Acknowledgements

RosaSeed is built on [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
by Vasimuddin Md, Sanchit Misra, Heng Li, and Chirag Jain.
RosaSeed uses [gsufsort](https://github.com/felipelouza/gsufsort)
by Louza et al. to construct SA and BWT during index construction.
RosaSeed bundles [safestringlib](https://github.com/intel/safestringlib)
by Intel Corporation, inherited from BWA-MEM2, for bounds-checked string
and memory operations.

---

## License

RosaSeed: MIT License.
Copyright (c) 2026 Shyama Gandhi, Bruce F. Cockburn, University of Alberta.

BWA-MEM2 components: MIT License.
Copyright (c) 2019 Vasimuddin Md, Sanchit Misra, Heng Li, Chirag Jain.

safestringlib (bundled in `ext/safestringlib/`): MIT License.
Copyright (c) 2014-2018 Intel Corporation; Copyright (c) 2012, 2013 Cisco Systems.

gsufsort: GPL-3.0 License.
Copyright (c) 2020 Louza, F.A., Telles, G.P., Gog, S., Prezza, N., Rosone, G.

See [LICENSE](LICENSE) for full text.
