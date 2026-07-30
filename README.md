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
├── preprocessing/           FASTQ preprocessing scripts
├── evaluation/              Accuracy evaluation scripts
│   ├── evaluate_alignments.py
│   └── README.md
├── src/
│   ├── rosaseed/            RosaSeed seeding algorithm
│   │   ├── rosaseed_phaseA.c
│   │   ├── rosaseed_phaseB.c
│   │   ├── rosaseed_gapfill_PhaseC.c
│   │   ├── rosaseed_core_bridge.cpp
│   │   ├── load_data_mmap.c
│   │   ├── mem_alloc.c
│   │   ├── helper_functions.c
│   │   ├── bwa.c
│   │   ├── file_dec.h
│   │   ├── macros.h
│   │   └── ...
│   └── [bwa-mem2 source files]
├── index-builder/           RosaSeed FM-index build pipeline
│   ├── build_2step_pipeline.sh
│   ├── rss_monitor.py
│   ├── src/
│   └── README.md
├── Makefile
├── LICENSE
└── README.md
└── README-bwamem2.md
```

---

## Requirements

| Requirement | Notes |
|---|---|
| GCC / g++ ≥ 7 | C++14 support required |
| GNU Make | Standard build |
| x86-64 with AVX2 | Required for `arch=native` |
| RAM ≥ 64 GB | ~50 GB peak at alignment time |

**Index builder additional requirements:**

```bash
sudo apt install gcc make git python3 python3-pip
pip install numpy
```

RAM: ~64 GB minimum (peak ~56 GB during RosaSeed index construction for the
T2T gapless human genome). The BWA-MEM2 index build (`bwa-mem2 index`)
requires ~62 GB RAM for the human genome.

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
> one used in Step 3 — BWA-MEM2 looks for its index files (`.0123`, `.pac`
> etc.) in the same directory as the FASTA.

---

## Index summary

RosaSeed requires two separate indexes:

| Index | Built by | Location | Size (T2T human) |
|---|---|---|---|
| BWA-MEM2 index | `./bwa-mem2 index genome.fna` | Same directory as FASTA | ~6.9 GB |
| RosaSeed index | `index-builder/build_2step_pipeline.sh` | `index-builder/index/<name>/` | ~42.7 GB |

**Combined peak memory at alignment time: ~49.61 GB**

---

## Flags reference

RosaSeed is controlled through two mechanisms:
- **Build-time flags** (`-D` flags via `CPPFLAGS_EXTRA`) — select algorithms and data structures at compile time
- **Runtime flags** (`--rs-*` and standard BWA-MEM2 flags) — tune thresholds at run time

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

### Runtime flags

These are passed on the command line at alignment time and can be tuned
without recompilation.

#### RosaSeed-specific runtime flags

| Flag | Default | Description |
|---|---|---|
| `--rs-index <dir>` | *(required)* | Path to the RosaSeed index directory produced by `build_2step_pipeline.sh`. Must contain `cp_occ_full.bin`, `c_vector.txt`, `ref16_packed.bin`, the SA split files, and the jump table. |
| `--rs-cap <int>` | 2000 | Phase A SA interval cap. Seeds whose BWT interval width exceeds this value are skipped in Phase A (too repetitive to be useful). Lower values run faster but may miss seeds in repetitive regions. Recommended: 2000 for standard use, 50 for miniRosaSeed. |

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

Reduces index memory to ~35 GB by using 8× SA compression.
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

Reduced-aggressiveness configuration benchmarked against minibwa.
Achieves 2.08× end-to-end speedup over minibwa at single thread
with 2.581 percentage points higher standard accuracy.

```bash
make arch=native CXX=g++ ROSASEED=1 \
  CPPFLAGS_EXTRA=" \
    -DLOAD_JTABLE_14nt \
    -DENABLE_F_RC_CHOICE \
    -DSA_COMPRESSION_FACTOR_POWER=1 \
    -DGAPFILL_ALWAYS_RUN_BOTH_STRANDS \
    -DGAPFILL_EARLY_EXIT \
    -DROSASEED_PRECHAIN_SINGLETON_SUPPRESS \
    -DROSASEED_PRECHAIN_TRIGGER=50 \
    -DROSASEED_PRECHAIN_WEAK_LEN=60"

./bwa-mem2 mem -t 1 -k 19 --rs-index index-builder/index/<name>/ \
    --rs-cap 50 genome.fna reads.fastq.gz > output.sam
```

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

Peak memory: **49.61 GB** — ~25% less than ERT/ERT2 (~66.3 GB each).
This comprises the RosaSeed 2-step FM-index (~42.7 GB) and the BWA-MEM2
index files (~6.9 GB) required by the downstream chaining and alignment
extension pipeline.

On a separate AMD Zen3 workstation, the reduced-aggressiveness
**miniRosaSeed** configuration achieves a **2.08× end-to-end speedup**
over minibwa at single-thread execution while attaining
**2.581 percentage points higher standard accuracy**, remaining faster
and more accurate simultaneously up to ~26 threads.

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
> University of Alberta.

Please also cite the original BWA-MEM2:

> Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient
> Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.
> IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019.
> doi:10.1109/IPDPS.2019.00041

If you use gsufsort for index building, please cite:

> Louza, F.A., Telles, G.P., Gog, S., Prezza, N., Rosone, G.
> gsufsort: constructing suffix arrays, LCP arrays and BWTs for string collections.
> Algorithms Mol Biol 15, 18 (2020). (https://link.springer.com/article/10.1186/s13015-020-00177-y)

---

## Acknowledgements

RosaSeed is built on [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
by Vasimuddin Md, Sanchit Misra, Heng Li, and Chirag Jain.
RosaSeed uses [gsufsort](https://github.com/felipelouza/gsufsort)
by Louza et al. to construct SA and BWT during index construction.

---

## License

RosaSeed: MIT License.
Copyright (c) 2026 Shyama Gandhi, University of Alberta.

BWA-MEM2 components: MIT License.
Copyright (c) 2019 Vasimuddin Md, Sanchit Misra, Heng Li, Chirag Jain.

gsufsort: GPL-3.0 License.
Copyright (c) 2020 Louza, F.A., Telles, G.P., Gog, S., Prezza, N., Rosone, G.

See [LICENSE](LICENSE) for full text.
