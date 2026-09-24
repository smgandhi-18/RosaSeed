#!/usr/bin/env bash
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
#
# build_compact_pipeline.sh
#
# Complete FM-index build pipeline for RosaSeed-Compact (1-step, radix-4).
# Takes any genome FASTA as downloaded (e.g. from NCBI): headers are stripped,
# contigs joined, and N / IUPAC ambiguity codes replaced with A, exactly as in
# the 2-step pipeline. The indexed text is forward + reverse complement.
#
# Steps:
#   1.  Check dependencies
#   2.  Compile all tools  (from src/)
#   3.  Locate or clone gsufsort-64  (into tools/)
#   4.  Preprocess FASTA  → clean single-line ACGT .txt
#   5.  Build forward + reverse-complement text  (single line)
#   6.  Run gsufsort-64  → .8.sa  .bwt
#   7.  Trim first entry from SA and BWT  (in-place)
#   8.  Verify BWT character counts
#   9.  Build compressed SA files  (CF = 1, 2, 4, 8)
#  10.  Delete original .8.sa  (-k to keep)
#  11.  Build cp_occ_compact.bin
#  12.  Build ref4_packed.bin
#  13.  Build jump table(s)
#
# Usage:
#   ./build_compact_pipeline.sh [OPTIONS] <genome.fa>
#
# Options:
#   -o <dir>   Output directory  (default: ./index/<genome_basename>_compact/)
#   -g <path>  Path to gsufsort-64 binary  (default: auto-install to tools/)
#   -k         Keep original .8.sa file  (default: delete after CF files)
#   -j <k>     Jump table k-mer size:
#                14      14-mer only  ( 2 GiB)
#                15      15-mer only  ( 8 GiB)
#                16      16-mer only  (32 GiB)
#                all     14 + 15 + 16
#                14,15   default (same as 2-step)
#                14,16  15,16
#   -h         Help
#
# Examples:
#   ./build_compact_pipeline.sh genome.fna
#   ./build_compact_pipeline.sh -o /data/index_compact -j 14,15 genome.fna
#
# Requirements: gcc  make  git  python3  numpy
# =============================================================================

set -euo pipefail

# ── colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; BLUE='\033[0;34m'
YELLOW='\033[1;33m'; BOLD='\033[1m'; RESET='\033[0m'
info()  { echo -e "${BLUE}[INFO]${RESET}  $*"; }
ok()    { echo -e "${GREEN}[ OK ]${RESET}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
die()   { echo -e "${RED}[FAIL]${RESET}  $*" >&2; exit 1; }
step()  { echo ""; echo -e "${BOLD}── $* ──${RESET}"; }
hr()    { echo "══════════════════════════════════════════════════════════════"; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/src"
TOOLS_DIR="$SCRIPT_DIR/tools"

# ── defaults ──────────────────────────────────────────────────────────────────
OUT_DIR=""
GSUFSORT=""
KEEP_SA=0
KMER_OPT="14,15"

usage() {
    sed -n '/^# Usage/,/^# Requirements/p' "$0" | sed 's/^# \{0,2\}//'
    exit 0
}

while getopts ":o:g:j:kh" opt; do
    case $opt in
        o) OUT_DIR="$OPTARG" ;;
        g) GSUFSORT="$OPTARG" ;;
        j) KMER_OPT="$OPTARG" ;;
        k) KEEP_SA=1 ;;
        h) usage ;;
        :) die "Option -$OPTARG requires an argument." ;;
       \?) die "Unknown option: -$OPTARG  (use -h for help)" ;;
    esac
done
shift $((OPTIND - 1))
[[ $# -lt 1 ]] && { warn "No genome FASTA supplied."; usage; }

INPUT_FA="$(realpath "$1")"
[[ -f "$INPUT_FA" ]] || die "Input file not found: $INPUT_FA"

INPUT_BASE="$(basename "$INPUT_FA")"
INPUT_NOEXT="${INPUT_BASE%.*}"

DO_14=0; DO_15=0; DO_16=0
case "$KMER_OPT" in
    14)      DO_14=1 ;;
    15)      DO_15=1 ;;
    16)      DO_16=1 ;;
    all)     DO_14=1; DO_15=1; DO_16=1 ;;
    "14,15") DO_14=1; DO_15=1 ;;
    "14,16") DO_14=1; DO_16=1 ;;
    "15,16") DO_15=1; DO_16=1 ;;
    *) die "Invalid -j value '$KMER_OPT'. Use: 14  15  16  all  14,15  14,16  15,16" ;;
esac

[[ -z "$OUT_DIR" ]] && OUT_DIR="$SCRIPT_DIR/index/${INPUT_NOEXT}_compact"
mkdir -p "$OUT_DIR" "$SRC_DIR" "$TOOLS_DIR"

CLEAN_FILE="$OUT_DIR/${INPUT_NOEXT}_clean.txt"
TEXT_FILE="$OUT_DIR/${INPUT_NOEXT}_compact_ref.txt"
BWT_FILE="$OUT_DIR/$(basename "$TEXT_FILE").bwt"
SA_FILE="$OUT_DIR/$(basename "$TEXT_FILE").8.sa"
CP_OCC_FILE="$OUT_DIR/cp_occ_compact.bin"
REF4_FILE="$OUT_DIR/ref4_packed.bin"

# ── RSS monitor setup ─────────────────────────────────────────────────────────
PEAK_FILE="$OUT_DIR/.rss_peak.txt"
LABEL_FILE="${PEAK_FILE}.label"
MON_PID=""
MONITOR_PY="$SCRIPT_DIR/rss_monitor.py"

start_monitor() {
    rm -f "$PEAK_FILE" "$LABEL_FILE"
    if [[ -f "$MONITOR_PY" ]]; then
        python3 "$MONITOR_PY" "$$" "$PEAK_FILE" 1.0 &
        MON_PID=$!
    fi
}

set_label() {
    echo "$*" > "$LABEL_FILE" 2>/dev/null || true
}

stop_monitor() {
    [[ -n "$MON_PID" ]] && kill "$MON_PID" 2>/dev/null || true
    MON_PID=""
}

read_peak() {
    if [[ -f "$PEAK_FILE" ]]; then
        local kb label
        kb=$(cut -f1 "$PEAK_FILE")
        label=$(cut -f2 "$PEAK_FILE")
        local gb
        gb=$(python3 -c "print(f'{$kb/1048576:.1f}')")
        echo "${gb} GB  (during: ${label})"
    else
        echo "unavailable (rss_monitor.py not found)"
    fi
}

# Clean up monitor on exit
trap 'stop_monitor; rm -f "$LABEL_FILE"' EXIT

hr
echo -e "${BOLD}  RosaSeed-Compact FM-index build pipeline${RESET}"
hr
info "Input FASTA   : $INPUT_FA"
info "Output dir    : $OUT_DIR"
info "Jump tables   : k=${KMER_OPT}"
[[ $KEEP_SA -eq 1 ]] \
    && info "Keep SA       : yes" \
    || info "Keep SA       : no  (deleted after CF files, saves ~50 GB)"
echo ""

WALL_T0=$(date +%s)
start_monitor

# =============================================================================
step "[ 1/13 ] Checking dependencies"
# =============================================================================
set_label "dependency check"
for cmd in gcc make git python3; do
    command -v "$cmd" &>/dev/null \
        && ok "$cmd  →  $(command -v $cmd)" \
        || die "$cmd not found. Install: sudo apt install $cmd"
done
python3 -c "import numpy" 2>/dev/null \
    && ok "python3 numpy" \
    || die "numpy not found. Install: sudo apt install python3-numpy"

# =============================================================================
step "[ 2/13 ] Compiling tools"
# =============================================================================
set_label "compiling tools"
compile_tool() {
    local name="$1" flags="${2:--O2}"
    local src="$SRC_DIR/${name}.c" bin="$SCRIPT_DIR/$name"
    [[ -f "$src" ]] || die "$name.c not found in $SRC_DIR"
    if [[ ! -x "$bin" ]] || [[ "$src" -nt "$bin" ]]; then
        info "Compiling $name ..."
        gcc $flags -o "$bin" "$src" \
            && ok "Compiled: $bin" || die "$name compile failed."
    else
        ok "Up-to-date: $bin"
    fi
}
compile_tool preprocess_genome        "-O2"
compile_tool derive_compactref        "-O2"
compile_tool convert_sa_to_bin_CF     "-O2"
compile_tool make_cp_occ_compact      "-O2"
compile_tool make_ref4_packed_compact "-O3 -std=c11"
compile_tool make_jumptable_compact   "-O2"

# =============================================================================
step "[ 3/13 ] Locating gsufsort-64"
# =============================================================================
set_label "locating gsufsort"
if [[ -n "$GSUFSORT" ]]; then
    [[ -x "$GSUFSORT" ]] || die "Not executable: $GSUFSORT"
    ok "Using: $GSUFSORT"
else
    for candidate in \
            "$TOOLS_DIR/gsufsort/gsufsort-64" \
            "$SCRIPT_DIR/gsufsort/gsufsort-64" \
            "$(command -v gsufsort-64 2>/dev/null || true)"; do
        if [[ -n "$candidate" && -x "$candidate" ]]; then
            GSUFSORT="$candidate"; ok "Found: $GSUFSORT"; break
        fi
    done
    if [[ -z "$GSUFSORT" ]]; then
        warn "gsufsort-64 not found, cloning into tools/gsufsort ..."
        git clone --depth 1 https://github.com/felipelouza/gsufsort.git \
            "$TOOLS_DIR/gsufsort"
        (cd "$TOOLS_DIR/gsufsort" && make)
        GSUFSORT="$TOOLS_DIR/gsufsort/gsufsort-64"
        [[ -x "$GSUFSORT" ]] || die "gsufsort compile failed."
        ok "Compiled: $GSUFSORT"
    fi
fi

# =============================================================================
step "[ 4/13 ] Preprocessing FASTA → clean ACGT"
# =============================================================================
set_label "preprocess FASTA"
info "Strips headers, joins chromosomes, replaces N/ambiguity with A."
T=$(date +%s)
"$SCRIPT_DIR/preprocess_genome" "$INPUT_FA" "$CLEAN_FILE" 2>&1
ok "Done in $(( $(date +%s) - T ))s"

CLEAN_BYTES=$(stat -c%s "$CLEAN_FILE")
GENOME_LEN=$(( CLEAN_BYTES - 1 ))
[[ $(wc -l < "$CLEAN_FILE") -eq 1 ]] || die "Clean file is not single-line"
ACGT_COUNT=$(tr -cd 'ACGT' < "$CLEAN_FILE" | wc -c)
[[ $ACGT_COUNT -eq $GENOME_LEN ]] || die "Non-ACGT chars found in clean file"
ok "Genome length : $GENOME_LEN bp"

RS_MAX_GENOME_BP=4294967296          # 2^32, see file_dec.h
if [[ $GENOME_LEN -ge $RS_MAX_GENOME_BP ]]; then
    die "Genome is $GENOME_LEN bp; must be below the RosaSeed index format limit of \
$RS_MAX_GENOME_BP bp (~4.29 Gbp).
       BWT positions are stored in 33-bit fields; a larger reference would be
       silently truncated and produce incorrect alignment coordinates.
       See README, 'Genome size limits'."
fi
ok "Genome size    : within the ${RS_MAX_GENOME_BP} bp format limit ($(( 100 * GENOME_LEN / RS_MAX_GENOME_BP ))% used)"

# =============================================================================
step "[ 5/13 ] Building forward + reverse-complement text"
# =============================================================================
set_label "derive_compactref"
T=$(date +%s)
"$SCRIPT_DIR/derive_compactref" "$CLEAN_FILE" "$TEXT_FILE" 2>&1
TEXT_LEN=$(( $(stat -c%s "$TEXT_FILE") - 1 ))
[[ $TEXT_LEN -eq $(( 2 * GENOME_LEN )) ]] || die "Text length $TEXT_LEN != 2 x $GENOME_LEN"
[[ $(wc -l < "$TEXT_FILE") -eq 1 ]] || die "Text file is not single-line"
ok "Done in $(( $(date +%s) - T ))s  →  $TEXT_LEN chars (2L)"

# =============================================================================
step "[ 6/13 ] Running gsufsort-64 (SA + BWT)"
# =============================================================================
set_label "gsufsort"
info "Single-threaded, may take ~30 minutes for a 3 Gbp T2T gapless human genome."
info "Peak RAM estimate: ~$(( TEXT_LEN * 9 / 1000000000 )) GB  (~9x text size)."
T=$(date +%s)
(cd "$OUT_DIR" && "$(realpath "$GSUFSORT")" "$TEXT_FILE" --sa --bwt)
T_GSUF=$(( $(date +%s) - T ))
[[ -f "$SA_FILE"  ]] || die "SA file not found: $SA_FILE"
[[ -f "$BWT_FILE" ]] || die "BWT file not found: $BWT_FILE"
SA_N=$(( $(stat -c%s "$SA_FILE") / 8 ))
ok "SA : $SA_N entries  ($(( $(stat -c%s "$SA_FILE") / 1073741824 )) GiB)"
ok "BWT: $(stat -c%s "$BWT_FILE") bytes"
ok "gsufsort time: ${T_GSUF}s"
[[ $SA_N -eq $(( TEXT_LEN + 2 )) ]] \
    && ok "n = TEXT_LEN + 2 ✓" \
    || warn "n=$SA_N but expected $((TEXT_LEN+2))"

# =============================================================================
step "[ 7/13 ] Trimming first entry from SA and BWT (in-place)"
# =============================================================================
set_label "trim SA+BWT in-place"
info "SA: remove first 8 bytes  |  BWT: remove first 1 byte"
info "Do NOT interrupt, files will be corrupted if stopped mid-way."
python3 - "$OUT_DIR/$(basename "$TEXT_FILE")" << 'PYEOF'
import os, sys
def trim_inplace(path, skip, label):
    sz=os.path.getsize(path); chunk=64*1024*1024; new=sz-skip
    print(f'  Trimming {label}: {sz/1073741824:.2f} GiB → {new/1073741824:.2f} GiB')
    rp=skip; wp=0
    with open(path,'r+b') as f:
        while rp<sz:
            f.seek(rp); d=f.read(min(chunk,sz-rp))
            f.seek(wp); f.write(d); rp+=len(d); wp+=len(d)
            print(f'\r    {wp/1073741824:.2f}/{new/1073741824:.2f} GiB ({100*wp/new:.1f}%)',end='',flush=True)
        f.truncate(new)
    print(f'\r    done, {new} bytes                    ')
b=sys.argv[1]
trim_inplace(b+'.8.sa', 8, 'SA ')
trim_inplace(b+'.bwt',  1, 'BWT')
PYEOF
EXPECTED_TRIM=$(( TEXT_LEN + 1 ))
SA_N_TRIM=$(( $(stat -c%s "$SA_FILE") / 8 ))
BWT_N_TRIM=$(stat -c%s "$BWT_FILE")
[[ $SA_N_TRIM -eq $EXPECTED_TRIM && $BWT_N_TRIM -eq $EXPECTED_TRIM ]] \
    && ok "Trim counts match: $SA_N_TRIM entries ✓" \
    || die "Trim mismatch: SA=$SA_N_TRIM BWT=$BWT_N_TRIM expected=$EXPECTED_TRIM"

# =============================================================================
step "[ 8/13 ] Verifying BWT character counts"
# =============================================================================
set_label "BWT verification"
python3 - "$BWT_FILE" "$TEXT_LEN" << 'PYEOF'
import numpy as np, sys
bwt, tl = sys.argv[1], int(sys.argv[2])
counts = {b: 0 for b in (65, 67, 71, 84, 0, 1)}; other = 0
with open(bwt, 'rb') as f:
    while True:
        d = np.frombuffer(f.read(256 * 1024 * 1024), dtype=np.uint8)
        if not len(d): break
        for b in counts: counts[b] += int(np.sum(d == b))
        other += int(np.sum(~np.isin(d, list(counts))))
A, C, G, T = counts[65], counts[67], counts[71], counts[84]
OK = '\033[0;32m[ OK ]\033[0m'; FAIL = '\033[0;31m[FAIL]\033[0m'; errs = 0
def chk(l, c):
    global errs
    print(f'  {OK if c else FAIL}  {l}'); errs += (not c)
chk(f'ACGT chars = {tl}', A + C + G + T == tl)
chk('A == T and C == G (forward + reverse complement)', A == T and C == G)
chk('null(0x00) = 1', counts[0] == 1)
chk('sep(0x01)  = 0', counts[1] == 0)
chk('other      = 0', other == 0)
if errs: print(f'\033[0;31m{errs} FAILED\033[0m'); sys.exit(1)
print('\033[0;32m  All BWT checks passed.\033[0m')
PYEOF

# =============================================================================
step "[ 9/13 ] Building compressed SA files (CF = 1, 2, 4, 8)"
# =============================================================================
set_label "compress SA (CF=1,2,4,8)"
T=$(date +%s)
"$SCRIPT_DIR/convert_sa_to_bin_CF" "$SA_FILE" "$OUT_DIR"
ok "Done in $(( $(date +%s) - T ))s"
for cf in 1 2 4 8; do
    ok "  sa_ls_word_cf${cf}.bin  $(du -h "$OUT_DIR/sa_ls_word_cf${cf}.bin" | cut -f1)"
    ok "  sa_ms_byte_cf${cf}.bin  $(du -h "$OUT_DIR/sa_ms_byte_cf${cf}.bin" | cut -f1)"
done

# =============================================================================
step "[ 10/13 ] Cleaning up original SA file"
# =============================================================================
set_label "delete original SA"
if [[ $KEEP_SA -eq 0 ]]; then
    SA_SZ=$(du -h "$SA_FILE" | cut -f1)
    info "Deleting original SA ($SA_SZ), CF files contain all needed data."
    rm -f "$SA_FILE"
    ok "Deleted: $SA_FILE"
else
    info "Keeping SA file (-k): $SA_FILE"
fi

# =============================================================================
step "[ 11/13 ] Building cp_occ_compact.bin"
# =============================================================================
set_label "make_cp_occ_compact"
T=$(date +%s)
"$SCRIPT_DIR/make_cp_occ_compact" "$BWT_FILE" "$CP_OCC_FILE" 2>&1
ok "Done in $(( $(date +%s) - T ))s  →  cp_occ_compact.bin $(du -h "$CP_OCC_FILE" | cut -f1)"

# =============================================================================
step "[ 12/13 ] Building ref4_packed.bin"
# =============================================================================
set_label "make_ref4_packed_compact"
T=$(date +%s)
"$SCRIPT_DIR/make_ref4_packed_compact" "$TEXT_FILE" "$REF4_FILE" "$TEXT_LEN" 2>&1
ok "Done in $(( $(date +%s) - T ))s  →  ref4_packed.bin $(du -h "$REF4_FILE" | cut -f1)"

# =============================================================================
step "[ 13/13 ] Building jump table(s)"
# =============================================================================
SA_LS_CF1="$OUT_DIR/sa_ls_word_cf1.bin"
SA_MS_CF1="$OUT_DIR/sa_ms_byte_cf1.bin"
[[ -f "$SA_LS_CF1" ]] || die "sa_ls_word_cf1.bin not found"
[[ -f "$SA_MS_CF1" ]] || die "sa_ms_byte_cf1.bin not found"
info "Memory: cp_occ + SA (CF=1), about 37 GB for a human genome (fm_pipeline estimate)"

if [[ $DO_14 -eq 1 ]]; then
    set_label "jumptable_14nt"
    info "Generating jumptable_14nt.bin  (4^14=268M entries, 2 GiB) ..."
    T=$(date +%s)
    "$SCRIPT_DIR/make_jumptable_compact" \
        "$CP_OCC_FILE" "$SA_LS_CF1" "$SA_MS_CF1" "$OUT_DIR" 14 2>&1
    ok "jumptable_14nt.bin: $(du -h "$OUT_DIR/jumptable_14nt.bin" | cut -f1)  in $(( $(date +%s)-T ))s"
fi

if [[ $DO_15 -eq 1 ]]; then
    set_label "jumptable_15nt"
    info "Generating jumptable_15nt.bin  (4^15=1.07B entries, 8 GiB) ..."
    T=$(date +%s)
    "$SCRIPT_DIR/make_jumptable_compact" \
        "$CP_OCC_FILE" "$SA_LS_CF1" "$SA_MS_CF1" "$OUT_DIR" 15 2>&1
    ok "jumptable_15nt.bin: $(du -h "$OUT_DIR/jumptable_15nt.bin" | cut -f1)  in $(( $(date +%s)-T ))s"
fi

if [[ $DO_16 -eq 1 ]]; then
    set_label "jumptable_16nt"
    info "Generating jumptable_16nt.bin  (4^16=4.3B entries, 32 GiB) ..."
    T=$(date +%s)
    "$SCRIPT_DIR/make_jumptable_compact" \
        "$CP_OCC_FILE" "$SA_LS_CF1" "$SA_MS_CF1" "$OUT_DIR" 16 2>&1
    ok "jumptable_16nt.bin: $(du -h "$OUT_DIR/jumptable_16nt.bin" | cut -f1)  in $(( $(date +%s)-T ))s"
fi

stop_monitor
set_label ""

# =============================================================================
# ── Final summary ─────────────────────────────────────────────────────────────
# =============================================================================
WALL_ELAPSED=$(( $(date +%s) - WALL_T0 ))
hr
echo -e "${BOLD}  Index build complete${RESET}"
hr
echo ""
echo -e "${BOLD}  Output: $OUT_DIR${RESET}"
echo ""

print_file() {
    local f="$OUT_DIR/$1"
    [[ -f "$f" ]] \
        && printf "  ✓  %-40s %s\n" "$1" "$(du -h "$f" | cut -f1)" \
        || printf "  –  %-40s %s\n" "$1" "(not generated)"
}

print_file "${INPUT_NOEXT}_clean.txt"
print_file "${INPUT_NOEXT}_compact_ref.txt"
print_file "$(basename "$TEXT_FILE").bwt"
echo ""
print_file "cp_occ_compact.bin"
print_file "ref4_packed.bin"
echo ""
for cf in 1 2 4 8; do
    print_file "sa_ls_word_cf${cf}.bin"
    print_file "sa_ms_byte_cf${cf}.bin"
done
echo ""
for k in 14 15 16; do
    print_file "jumptable_${k}nt.bin"
done

echo ""
hr
echo -e "${BOLD}  Resource summary${RESET}"
hr
echo -e "  Wall time : $(( WALL_ELAPSED/3600 ))h $(( (WALL_ELAPSED%3600)/60 ))m $(( WALL_ELAPSED%60 ))s"
echo -e "  Peak RSS  : $(read_peak)"
echo ""
info "Index ready for RosaSeed-Compact (build with ROSASEED_1STEP=1)."
info "Pass to the aligner with: --rs-index $OUT_DIR"
echo ""
rm -f "$PEAK_FILE" "$LABEL_FILE"
