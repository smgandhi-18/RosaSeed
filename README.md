# RosaSeed

**RosaSeed** is a fast, configurable seeding algorithm for short-read DNA
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
│   ├── rosaseed/            RosaSeed seeding algorithm
│   │   ├── rosaseed_phaseA.c
│   │   ├── rosaseed_phaseB.c
│   │   ├── rosaseed_gapfill_PhaseC.c
│   │   ├── rosaseed_core_bridge.cpp
│   │   ├── load_data_mmap.c
│   │   ├── mem_alloc.c
│   │   ├── helper_functions.c
│   │   ├── bwa.c
│   │   ├── file_dec_new.h
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

RosaSeed replaces the seeding stage and rewrites the chaining stage with pruning heuristics but
alignment extension (BSW) pipeline still uses BWA-MEM2's own index files.
Generate them with:

```bash
./bwa-mem2 index /path/to/genome.fna
```

This produces the following files alongside the FASTA (~6.9 GB total):
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

Build the RosaSeed index for your reference genome (~44 min, ~56 GB RAM
peak for human T2T genome):

```bash
cd index-builder
chmod +x build_2step_pipeline.sh
./build_2step_pipeline.sh /path/to/genome.fna
cd ..
```

See [index-builder/README.md](index-builder/README.md) for full options.
Index files are written to `index-builder/index/<genome_name>/` by default
(~42.7 GB for the T2T human genome with 14 and 15-mer jump tables).

### 5. Run alignment

```bash
./bwa-mem2 mem \
    -t 1 \
    -k 19 \
    --rs-index index-builder/index/<genome_name>/ \
    --rs-cap 2000 \
    /path/to/genome.fna \
    reads_R1.fastq.gz \
    > output.sam
```

> **Note:** The FASTA path supplied to `./bwa-mem2 mem` must be the same
> one used in Step 3 : BWA-MEM2 looks for its index files (`.0123`, `.pac`
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

RAM: ~64 GB minimum (peak ~56 GB during RosaSeed index construction for T2T 
gapless human genome). The BWA-MEM2 index build (`bwa-mem2 index`) requires
~62 GB RAM for the human genome.

---

## Build flags reference

These `-D` flags are passed via `CPPFLAGS_EXTRA` at build time.

### Jump table size

| Flag | Entries | Size | Description |
|---|---|---|---|
| `-DLOAD_JTABLE_14nt` | 268M | 2 GiB | 14-mer table (default) |
| `-DLOAD_JTABLE_15nt` | 1.07B | 8 GiB | 15-mer table |
| `-DLOAD_JTABLE_16nt` | 4.29B | 32 GiB | 16-mer table |

Use exactly one per build.

### SA compression factor

| Flag | CF | SA size | Description |
|---|---|---|---|
| `-DSA_COMPRESSION_FACTOR_POWER=0` | 1 | ~29 GB | Full SA |
| `-DSA_COMPRESSION_FACTOR_POWER=1` | 2 | ~15 GB | Half SA |
| `-DSA_COMPRESSION_FACTOR_POWER=2` | 4 | ~7.3 GB | Quarter SA |
| `-DSA_COMPRESSION_FACTOR_POWER=3` | 8 | ~3.6 GB | Eighth SA  |

### Other flags

| Flag | Description |
|---|---|
| `-DENABLE_F_RC_CHOICE` | Try both strands at first pivot, pick smallest one |
| `-DGAPFILL_ALWAYS_RUN_BOTH_STRANDS` | Run gap fill on both strands |
| `-DGAPFILL_EARLY_EXIT` | Exit gap fill early when gap is covered |
| `-DROSASEED_PRECHAIN_SINGLETON_SUPPRESS` | Suppress singleton seeds before chaining |
| `-DROSASEED_PRECHAIN_TRIGGER=N` | Apply pre-chain filter when seed count > N |
| `-DROSASEED_PRECHAIN_WEAK_LEN=N` | Seeds shorter than N bp treated as weak |

## Runtime flags

| Flag | Description |
|---|---|
| `-t <int>` | Number of threads |
| `-k <int>` | Minimum seed length (recommended: 19) |
| `--rs-index <dir>` | Path to RosaSeed index directory |
| `--rs-cap <int>` | Maximum seed interval size to extend (recommended: 2000) |

---

### Recommended configurations

**Recommended default configuration:**
```bash
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
Use the alignment command as:
```
./bwa-mem2 mem \
    -t 1 \
    -k 19 \
    --rs-index index-builder/index/<genome_name>/ \
    --rs-cap 2000 \
    /path/to/genome.fna \
    reads_R1.fastq.gz \
    > output.sam
```

**Memory-efficient:**
```bash
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

**miniRosaSeed (reduced aggressiveness, benchmarked vs minibwa):**
```bash
CPPFLAGS_EXTRA=" \
  -DLOAD_JTABLE_14nt \
  -DENABLE_F_RC_CHOICE \
  -DSA_COMPRESSION_FACTOR_POWER=1 \
  -DGAPFILL_ALWAYS_RUN_BOTH_STRANDS \
  -DGAPFILL_EARLY_EXIT \
  -DROSASEED_PRECHAIN_SINGLETON_SUPPRESS \
  -DROSASEED_PRECHAIN_TRIGGER=50 \
  -DROSASEED_PRECHAIN_WEAK_LEN=60"
```
Use the alignment command as:
```
./bwa-mem2 mem \
    -t 1 \
    -k 19 \
    --rs-index index-builder/index/<genome_name>/ \
    --rs-cap 50 \
    /path/to/genome.fna \
    reads_R1.fastq.gz \
    > output.sam
```
---

## Key results

Compared against BWA-MEM2, Minimap2, Bowtie2, and the Enumerated Radix
Tree (ERT) algorithm on simulated and real short-read datasets:

| Comparison | Seed processing speedup | Total alignment speedup |
|---|---|---|
| vs BWA-MEM2 | **15.0×** | **3.84×** |
| vs ERT | 4.0× | 2.48× |
| vs ERT2 | 1.2× | 1.44× |
| vs Minimap2 | 4.8× | 2.48× |
| vs Bowtie2 | 33.2× | 10.8× |

Peak memory: **49.61 GB** which is ~25% less than ERT/ERT2 (~66.3 GB each).
This comprises the RosaSeed 2-step FM-index (~42.7 GB) and the BWA-MEM2
index files (~6.9 GB) required by the downstream chaining and alignment
extension pipeline.

On a separate AMD Zen3 workstation, the reduced-aggressiveness
**miniRosaSeed** configuration achieves a **2.08× end-to-end speedup**
over minibwa at single-thread execution while attaining
**2.581 percentage points higher standard accuracy**, remaining faster
and more accurate simultaneously up to ~28 threads.

RosaSeed exposes multiple runtime-memory operating configurations
spanning memory-efficient and high-performance designs, enabling users
to select operating points appropriate for diverse computational
environments while maintaining alignment accuracy comparable to BWA-MEM2.

---

## Citation

If you use RosaSeed in your research, please cite:

> Gandhi Shyama, Cockburn Bruce (2026). RosaSeed: A Fast and Configurable Seeding
> Algorithm for Short-Read DNA Sequence Alignment.
> University of Alberta.

Please also cite the original BWA-MEM2:

> Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient
> Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.
> IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019.
> doi:https://doi.org/10.1109/IPDPS.2019.00041

If you use gsufsort for index building, please cite:

> Louza, F.A., Telles, G.P., Gog, S., Prezza, N., Rosone, G.. gsufsort: 
> constructing suffix arrays, LCP arrays and BWTs for string collections. 
> Algorithms Mol Biol 15, 18 (2020). https://doi.org/10.1186/s13015-020-00177-y
---

## Acknowledgements

RosaSeed is built on [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
by Vasimuddin Md, Sanchit Misra, Heng Li, and Chirag Jain.
RosaSeed uses [gsufsort](https://github.com/felipelouza/gsufsort)
by Louza et al., to construct SA and BWT during index construction.

---

## License

RosaSeed: MIT License.
Copyright (c) 2026 Shyama Gandhi, University of Alberta.

BWA-MEM2 components: MIT License.
Copyright (c) 2019 Vasimuddin Md, Sanchit Misra, Heng Li, Chirag Jain.

See [LICENSE](LICENSE) for full text.

gsufsort: GPL-3.0 License
Copyright (c) 2020 Louza, F.A., Telles, G.P., Gog, S., Prezza, N., Rosone, G.
