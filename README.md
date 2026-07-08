# RosaSeed

**RosaSeed** is a fast, configurable seeding algorithm for short-read DNA
sequence alignment. It replaces the BWA-MEM2 seeding kernel and feeds
candidate seeds directly into the existing BWA-MEM2 chaining and alignment
extension pipeline, producing standard SAM output compatible with downstream
variant analysis tools.

> See [README-bwamem2.md](README-bwamem2.md) for the original BWA-MEM2 documentation.

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

Peak memory: **49.61 GB** — ~25% less than ERT/ERT2 (~66.3 GB each).

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

## Repository structure

```
RosaSeed/
├── src/
│   ├── rosaseed/            RosaSeed seeding algorithm
│   │   ├── rosaseed_phaseI_coroutine.c
│   │   ├── rosaseed_phaseII_coro.c
│   │   ├── rosaseed_gapfillphase_pf.c
│   │   ├── rosaseed_core_bridge_coroutine_p2.cpp
│   │   ├── load_data_mmap.c
│   │   ├── mem_alloc.c
│   │   ├── helper_functions.c
│   │   ├── bwa.c
│   │   ├── file_dec_new.h
│   │   ├── macros.h
│   │   └── ...
│   └── [bwa-mem2 source files]
├── index-builder/           FM-index build pipeline
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

### 2. Build the index

Before running alignment, build the 2-step FM-index for your reference genome:

```bash
cd index-builder
chmod +x build_2step_pipeline.sh
./build_2step_pipeline.sh /path/to/genome.fna
cd ..
```

See [index-builder/README.md](index-builder/README.md) for full options.
Index files are written to `index-builder/index/<genome_name>/` by default.

### 3. Build RosaSeed

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

### 4. Run alignment

```bash
./bwa-mem2 mem \
    -t 20 \
    -k 19 \
    --rs-index /path/to/index-builder/index/genome/ \
    --rs-cap 2000 \
    genome.fna \
    reads_R1.fastq.gz \
    > output.sam
```

---

## Requirements

| Requirement | Notes |
|---|---|
| GCC / g++ ≥ 7 | C++14 support required |
| GNU Make | Standard build |
| x86-64 with AVX2 | Required for `arch=native` |
| RAM ≥ 64 GB | Index load: ~50 GB peak |

**Index builder requirements:**

```bash
sudo apt install gcc make git python3 python3-pip
pip install numpy
```

RAM: ~64 GB minimum (peak ~56 GB during index construction for human genome).

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
| `-DSA_COMPRESSION_FACTOR_POWER=0` | 1 | ~30 GB | Full SA (highest accuracy) |
| `-DSA_COMPRESSION_FACTOR_POWER=1` | 2 | ~15 GB | Half SA (recommended) |
| `-DSA_COMPRESSION_FACTOR_POWER=2` | 4 | ~7.3 GB | Quarter SA |
| `-DSA_COMPRESSION_FACTOR_POWER=3` | 8 | ~3.6 GB | Eighth SA (lowest memory) |

### Other flags

| Flag | Description |
|---|---|
| `-DENABLE_F_RC_CHOICE` | Try both strands at first pivot, pick better one |
| `-DGAPFILL_ALWAYS_RUN_BOTH_STRANDS` | Run gap fill on both strands |
| `-DGAPFILL_EARLY_EXIT` | Exit gap fill early when good seed found |
| `-DROSASEED_PRECHAIN_SINGLETON_SUPPRESS` | Suppress singleton seeds before chaining |
| `-DROSASEED_PRECHAIN_TRIGGER=N` | Apply pre-chain filter when seed count > N |
| `-DROSASEED_PRECHAIN_WEAK_LEN=N` | Seeds shorter than N bp treated as weak |

### Recommended configurations

**High-performance (recommended default):**
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

**Memory-efficient:**
```bash
CPPFLAGS_EXTRA=" \
  -DLOAD_JTABLE_14nt \
  -DENABLE_F_RC_CHOICE \
  -DSA_COMPRESSION_FACTOR_POWER=3"
```

**miniRosaSeed (reduced aggressiveness, benchmarked vs minibwa):**
```bash
CPPFLAGS_EXTRA=" \
  -DLOAD_JTABLE_14nt \
  -DENABLE_F_RC_CHOICE \
  -DSA_COMPRESSION_FACTOR_POWER=1"
```

---

## Runtime flags

| Flag | Description |
|---|---|
| `-t <int>` | Number of threads |
| `-k <int>` | Minimum seed length (recommended: 19) |
| `--rs-index <dir>` | Path to RosaSeed index directory |
| `--rs-cap <int>` | Maximum seed interval size to extend (recommended: 2000) |

---

## Citation

If you use RosaSeed in your research, please cite:

> Gandhi, S. et al. (2026). RosaSeed: A Fast and Configurable Seeding
> Algorithm for Short-Read DNA Sequence Alignment.
> University of Alberta.

Please also cite the original BWA-MEM2:

> Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. Efficient Architecture-Aware Acceleration
> of BWA-MEM for Multicore Systems. IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. 10.1109/IPDPS.2019.00041 

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
