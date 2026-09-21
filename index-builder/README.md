# 2-step RosaSeed FM-index Build Pipeline

Builds a complete FM-index from any FASTA genome file for the 2-step RosaSeed.

## Repository structure

```
rosaseed-index/
├── build_2step_pipeline.sh   ← run this
├── rss_monitor.py            ← background RAM tracker (auto-used by script)
├── src/
│   ├── preprocess_genome.c       FASTA → clean single-line ACGT
│   ├── derive_2stepref.c         ACGT → base-16 2-step reference (FE+FO+RCE+RCO)
│   ├── convert_sa_to_bin_CF.c    binary SA → compressed SA files (CF=1,2,4,8)
│   ├── make_cp_occ_2step.c       binary BWT → OCC checkpoints + bitvectors
│   ├── make_ref16_packed.c       base-16 ASCII ref → nibble-packed binary
│   └── make_jumptable_2step.c    cp_occ + SA → k-mer jump tables
├── tools/                        auto-populated on first run
│   └── gsufsort/                 auto-cloned from GitHub if not found
└── README.md
```

All C programs are compiled automatically on first run.
gsufsort is downloaded and compiled automatically if not already present.

## Requirements

```bash
sudo apt install gcc make git python3 python3-numpy
```

## Quick start

```bash
git clone <this-repo> rosaseed-index
cd rosaseed-index
chmod +x build_2step_pipeline.sh

# Build index: 14+15-mer jump tables by default
./build_2step_pipeline.sh /path/to/genome.fna
```

Index files are written to `./index/<genome_basename>/` by default.

## Usage

```
./build_2step_pipeline.sh [OPTIONS] <genome.fa>

Options:
  -o <dir>   Output directory  (default: ./index/<genome_basename>/)
  -g <path>  Path to gsufsort-64 binary  (default: auto-install to tools/)
  -k         Keep original .8.sa file  (default: delete after CF files)
  -j <k>     Jump table k-mer size:
               14      14-mer only  ( 2 GiB output, ~10 min)
               15      15-mer only  ( 8 GiB output, ~25 min)
               16      16-mer only  (32 GiB output, ~20 min)
               all     14 + 15 + 16
               14,15   default: recommended for most read lengths
  -h         Help
```

## Examples

```bash
# T2T CHM13: default 14+15-mer tables, auto output directory
./build_2step_pipeline.sh GCF_009914755.1_T2T-CHM13v2.0_genomic.fna

# Custom output directory, all three jump tables
./build_2step_pipeline.sh -o /data/T2T_index -j all T2T_CHM13.fna

# GRCh38: 15-mer only, keep original SA file
./build_2step_pipeline.sh -k -j 15 -o /scratch/hg38_index hg38.fa
```

## Output files

| File | Size (T2T CHM13) | Description |
|------|-----------------|-------------|
| `*_clean.txt` | 3.0 GB | Clean ACGT genome (single line, N→A) |
| `*_2step_ref.txt` | 5.9 GB | Base-16 2-step reference (FE+FO+RCE+RCO) |
| `*.bwt` | 5.9 GB | BWT (binary, trimmed) |
| `cp_occ_full.bin` | 24 GB | OCC checkpoints + bitvectors |
| `c_vector.txt` | tiny | C[] prefix sums (16 tab-separated values) |
| `ref16_packed.bin` | 3.0 GB | Nibble-packed base-16 reference |
| `sa_ls_word_cf1.bin` | 24 GB | SA lower 32 bits (CF=1, all rows) |
| `sa_ms_byte_cf1.bin` | 5.9 GB | SA bits 32-39 (CF=1) |
| `sa_ls_word_cf2.bin` | 12 GB | SA CF=2 |
| `sa_ms_byte_cf2.bin` | 3.0 GB | |
| `sa_ls_word_cf4.bin` | 5.9 GB | SA CF=4 |
| `sa_ms_byte_cf4.bin` | 1.5 GB | |
| `sa_ls_word_cf8.bin` | 3.0 GB | SA CF=8 |
| `sa_ms_byte_cf8.bin` | 744 MB | |
| `jumptable_14nt.bin` | 2.0 GiB | 16^7 jump pointers |
| `jumptable_15nt.bin` | 8.0 GiB | 4^15 jump pointers |
| `jumptable_16nt.bin` | 32 GiB | 4^16 jump pointers (optional, -j 16 or -j all) |

## Resource requirements (T2T CHM13 v2.0, 3.1 Gbp, observed on Threadripper P620)

| Step | RAM peak | Time |
|------|----------|------|
| preprocess + derive_2stepref | ~7 GB | ~1 min |
| gsufsort (SA + BWT) | **~56 GB** ← absolute peak | ~25 min |
| trim + verify | < 1 GB | ~5 min |
| convert_sa_to_bin_CF | < 1 GB | ~2 min |
| make_cp_occ_2step | < 1 GB | ~1 min |
| make_ref16_packed | ~3 GB | < 1 min |
| make_jumptable (14+15-mer) | ~54 GB | ~5 min |
| **Total (14+15-mer default)** | **~56 GB peak** | **~41 min** |
| Total (all three tables) | ~54 GB | ~61 min |

> **Minimum system RAM: 64 GB.**
> The hard ceiling is gsufsort at ~56 GB for the T2T human genome (~12× the
> 5.9 GB single-line reference file). The jump table step peaks at ~54 GB
> (cp_occ 25 GB + SA CF=1 30 GB). Both fit within 64 GB with the OS overhead.
>
> The pipeline prints an observed peak RSS at the end of each run.
> The RSS monitor samples every 1 second; gsufsort's true peak may be slightly
> higher than the reported value, use your system resource monitor for
> independent confirmation.

## Pipeline steps

```
[ 1/15 ] Check dependencies         gcc, make, git, python3, numpy
[ 2/15 ] Compile tools              src/*.c → binaries (only if source is newer)
[ 3/15 ] Locate gsufsort-64         auto-clone + compile into tools/ if needed
[ 4/15 ] Preprocess FASTA           strip headers, join chromosomes, N/IUPAC→A
[ 5/15 ] Build 2-step reference     FE + FO + RCE + RCO base-16 encoding
[ 6/15 ] Validate reference         hex chars only, correct segment lengths
[ 7/15 ] Single-line conversion     gsufsort requires single-line input
[ 8/15 ] Run gsufsort               SA + BWT construction  (~25 min, ~56 GB RAM)
[ 9/15 ] Trim first entry           remove null terminator in-place from SA + BWT
[10/15 ] Verify BWT                 base-16 char counts, 1 terminator, no junk
[11/15 ] Compressed SA files        CF=1,2,4,8: 5-byte split per entry
[12/15 ] Delete original SA         frees ~50 GB  (use -k to keep)
[13/15 ] cp_occ_full.bin            OCC checkpoint + bitvector file (~25 GB)
[14/15 ] ref16_packed.bin           nibble-packed base-16 reference (~3 GB)
[15/15 ] Jump table(s)              14+15-mer default; -j 16 or -j all for more
```

## What the output looks like

```
══════════════════════════════════════════════════════════════
  Index build complete
══════════════════════════════════════════════════════════════

  Output: ./index/GCF_009914755.1_T2T-CHM13v2.0_genomic/

  ✓  *_clean.txt                             3.0G
  ✓  *_2step_ref.txt                         5.9G
  ✓  *.bwt                                   5.9G
  ✓  cp_occ_full.bin                         24G
  ✓  c_vector.txt                            165 bytes
  ✓  ref16_packed.bin                        3.0G
  ✓  sa_ls_word_cf1.bin                      24G
  ...
  ✓  jumptable_14nt.bin                      2.0G
  ✓  jumptable_15nt.bin                      8.0G
  –  jumptable_16nt.bin                      (not generated)

══════════════════════════════════════════════════════════════
  Resource summary
══════════════════════════════════════════════════════════════
  Wall time : 0h 40m 55s
  Peak RSS  : 52.4 GB  (during: jumptable_15nt)
```

## Non-ACGT character handling

N and all IUPAC ambiguity codes (R,Y,S,W,K,M,B,D,H,V) are replaced with `A`.  
Sequence length is preserved, genome coordinates remain intact.  
Reads may align to replaced regions; exclude those chromosomes from the input FASTA if needed.

## The 2-step encoding

The 2-step reference encodes pairs of consecutive nucleotides as a single
base-16 symbol (0-9, A-F), reducing the alphabet from 4 to 16 symbols.

```
NT pair  → base-16 symbol
AA=0  AC=1  AG=2  AT=3
CA=4  CC=5  CG=6  CT=7
GA=8  GC=9  GG=A  GT=B
TA=C  TC=D  TG=E  TT=F
```

Four segments are built from genome G of length L:
```
FE  (Forward Even)  : pairs at positions (0,1), (2,3), ..., floor(L/2) symbols
FO  (Forward Odd)   : pairs at positions (1,2), (3,4), ..., floor((L-1)/2) symbols
RCE (RC Even)       : same on reverse complement
RCO (RC Odd)        : same on reverse complement
```

Concatenated as FE+FO+RCE+RCO, this allows any 14-nt (7 pairs), 15-nt
(7 pairs + 1 NT init), or 16-nt (8 pairs) seed to be looked up with a
single index query.

## References

- gsufsort: https://github.com/felipelouza/gsufsort
- Louza et al., "gsufsort: constructing suffix arrays, LCP arrays and BWTs
  for string collections", Algorithms for Molecular Biology, 2020.
