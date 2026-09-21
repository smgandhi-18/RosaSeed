#!/usr/bin/env bash
# =============================================================================
# build_2step_pipeline.sh
#
# Complete FM-index build pipeline for the 2-step RosaSeed aligner.
#
# Steps:
#   1.  Check dependencies
#   2.  Compile all tools  (from src/)
#   3.  Locate or clone gsufsort-64  (into tools/)
#   4.  Preprocess FASTA  → clean single-line ACGT .txt
#   5.  Build 2-step base-16 reference  (FE+FO+RCE+RCO)
#   6.  Validate 2-step reference
#   7.  Convert to single-line for gsufsort
#   8.  Run gsufsort-64  → .8.sa  .bwt
#   9.  Trim first entry from SA and BWT  (in-place)
#  10.  Verify BWT character counts
#  11.  Build compressed SA files  (CF = 1, 2, 4, 8)
#  12.  Delete original .8.sa  (saves ~50 GB, -k to keep)
#  13.  Build cp_occ_full.bin + c_vector.txt
#  14.  Build ref16_packed.bin
#  15.  Build jump table(s)
#
# Usage:
#   ./build_2step_pipeline.sh [OPTIONS] <genome.fa>
#
# Options:
#   -o <dir>   Output directory  (default: ./index/<genome_basename>/)
#   -g <path>  Path to gsufsort-64 binary  (default: auto-install to tools/)
#   -k         Keep original .8.sa file  (default: delete after CF files)
#   -j <k>     Jump table k-mer size:
#                14      14-mer only  ( 2 GiB, ~10 min)
#                15      15-mer only  ( 8 GiB, ~25 min)
#                16      16-mer only  (32 GiB, ~20 min)
#                all     14 + 15 + 16
#                14,15   default — recommended for most read lengths
#   -h         Help
#
# Examples:
#   ./build_2step_pipeline.sh genome.fna
#   ./build_2step_pipeline.sh -o /data/index -j all genome.fna
#   ./build_2step_pipeline.sh -k -j 15 -o /scratch/index hg38.fa
#
# Requirements: gcc  make  git  python3  numpy
# RAM:  ~57 GB peak for T2T Human Genome (gsufsort)  ~54 GB (jump table generation)
# Disk: ~200 GB during build, ~120 GB final index (includes all SA compression factors files)
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

[[ -z "$OUT_DIR" ]] && OUT_DIR="$SCRIPT_DIR/index/$INPUT_NOEXT"
mkdir -p "$OUT_DIR" "$SRC_DIR" "$TOOLS_DIR"

CLEAN_FILE="$OUT_DIR/${INPUT_NOEXT}_clean.txt"
REF2STEP_FILE="$OUT_DIR/${INPUT_NOEXT}_2step_ref.txt"
SL_FILE="$OUT_DIR/${INPUT_NOEXT}_2step_ref_singleline.txt"
BWT_FILE="$OUT_DIR/$(basename "$SL_FILE").bwt"
SA_FILE="$OUT_DIR/$(basename "$SL_FILE").8.sa"
CP_OCC_FILE="$OUT_DIR/cp_occ_full.bin"
CVEC_FILE="$OUT_DIR/c_vector.txt"
REF16_FILE="$OUT_DIR/ref16_packed.bin"

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

# ── banner ────────────────────────────────────────────────────────────────────
hr
echo -e "${BOLD}  2-step RosaSeed FM-index build pipeline${RESET}"
hr
info "Input FASTA   : $INPUT_FA"
info "Output dir    : $OUT_DIR"
info "Jump tables   : k=${KMER_OPT}"
[[ $KEEP_SA -eq 1 ]] \
    && info "Keep SA       : yes" \
    || info "Keep SA       : no  (deleted after CF files — saves ~50 GB)"
echo ""

WALL_T0=$(date +%s)
start_monitor

# =============================================================================
step "[ 1/15 ] Checking dependencies"
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
step "[ 2/15 ] Compiling tools"
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
compile_tool preprocess_genome    "-O2"
compile_tool derive_2stepref      "-O2"
compile_tool convert_sa_to_bin_CF "-O2"
compile_tool make_cp_occ_2step    "-O2"
compile_tool make_ref16_packed    "-O3 -std=c11"
compile_tool make_jumptable_2step "-O2"

# =============================================================================
step "[ 3/15 ] Locating gsufsort-64"
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
        warn "gsufsort-64 not found — cloning into tools/gsufsort ..."
        git clone --depth 1 https://github.com/felipelouza/gsufsort.git \
            "$TOOLS_DIR/gsufsort"
        (cd "$TOOLS_DIR/gsufsort" && make)
        GSUFSORT="$TOOLS_DIR/gsufsort/gsufsort-64"
        [[ -x "$GSUFSORT" ]] || die "gsufsort compile failed."
        ok "Compiled: $GSUFSORT"
    fi
fi

# =============================================================================
step "[ 4/15 ] Preprocessing FASTA → clean ACGT"
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

RS_MAX_GENOME_BP=4294967296          # 2^32 — see file_dec.h
if [[ $GENOME_LEN -gt $RS_MAX_GENOME_BP ]]; then
    die "Genome is $GENOME_LEN bp, exceeding the RosaSeed index format limit of \
$RS_MAX_GENOME_BP bp (~4.29 Gbp).
       BWT positions are stored in 33-bit fields; a larger reference would be
       silently truncated and produce incorrect alignment coordinates.
       See README, 'Genome size limits'."
fi
ok "Genome size    : within the ${RS_MAX_GENOME_BP} bp format limit ($(( 100 * GENOME_LEN / RS_MAX_GENOME_BP ))% used)"

# =============================================================================
step "[ 5/15 ] Building 2-step base-16 reference (FE+FO+RCE+RCO)"
# =============================================================================
set_label "derive_2stepref (~7 GB RAM)"
info "Encodes ACGT pairs → base-16 symbols, builds FE+FO+RCE+RCO"
DERIVE_TMP="$OUT_DIR/GaplessGenomeT2T_2step_ref_verify.fasta"
T=$(date +%s)
(cd "$OUT_DIR" && "$(realpath "$SCRIPT_DIR/derive_2stepref")" "$CLEAN_FILE") 2>&1
[[ -f "$DERIVE_TMP" ]] || die "derive_2stepref output not found"
tail -n +2 "$DERIVE_TMP" > "$REF2STEP_FILE"
rm -f "$DERIVE_TMP"
ok "Done in $(( $(date +%s) - T ))s  →  $(du -h "$REF2STEP_FILE" | cut -f1)"

# =============================================================================
step "[ 6/15 ] Validating 2-step reference"
# =============================================================================
set_label "validating reference"
HEX_COUNT=$(tr -cd '0-9A-Fa-f' < "$REF2STEP_FILE" | wc -c)
TOTAL_CHARS=$(tr -d '\n' < "$REF2STEP_FILE" | wc -c)
python3 - "$TOTAL_CHARS" "$GENOME_LEN" "$HEX_COUNT" << 'PYEOF'
import sys
total, L, hexc = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
fe = L//2; fo = (L-1)//2; expected = 2*fe+2*fo
OK='\033[0;32m[ OK ]\033[0m'; FAIL='\033[0;31m[FAIL]\033[0m'
errs=0
def chk(l,c):
    global errs
    print(f'  {OK if c else FAIL}  {l}')
    if not c: errs+=1
chk(f'all chars valid base-16 ({hexc})', hexc==total)
chk(f'total = {expected}  (2xFE + 2xFO)', total==expected)
if errs: print(f'\033[0;31m  {errs} check(s) FAILED\033[0m'); exit(1)
print(f'\033[0;32m  Checks passed.\033[0m  L={L} ({"odd" if L%2 else "even"})  FE={fe}  FO={fo}  total={total}')
PYEOF

# =============================================================================
step "[ 7/15 ] Converting to single-line for gsufsort"
# =============================================================================
set_label "single-line conversion"
LINE_COUNT=$(wc -l < "$REF2STEP_FILE")
if [[ $LINE_COUNT -eq 1 ]]; then
    ok "Already single-line — using directly."
    SL_FILE="$REF2STEP_FILE"
else
    warn "Input has $LINE_COUNT lines — converting to single line."
    tr -d '\n' < "$REF2STEP_FILE" > "$SL_FILE"
    printf '\n' >> "$SL_FILE"
    SL_HEX=$(tr -cd '0-9A-Fa-f' < "$SL_FILE" | wc -c)
    [[ $SL_HEX -eq $HEX_COUNT ]] \
        || die "Char count mismatch after conversion: $SL_HEX vs $HEX_COUNT"
    ok "Converted: $SL_FILE  ($SL_HEX chars ✓)"
fi
SL_BASE="$(basename "$SL_FILE")"
BWT_FILE="$OUT_DIR/${SL_BASE}.bwt"
SA_FILE="$OUT_DIR/${SL_BASE}.8.sa"

# =============================================================================
step "[ 8/15 ] Running gsufsort-64 (SA + BWT)"
# =============================================================================
set_label "gsufsort (~57 GB RAM peak for T2T Gapless Human Genome)"
info "Single-threaded — may take ~30 minutes for a 3 Gbp T2T gapless human genome."
info "Peak RAM Estimation: ~$(( HEX_COUNT * 9 / 1000000000 )) GB  (~9x reference file size)."
T=$(date +%s)
(cd "$OUT_DIR" && "$(realpath "$GSUFSORT")" "$SL_FILE" --sa --bwt)
T_GSUF=$(( $(date +%s) - T ))
[[ -f "$SA_FILE"  ]] || die "SA file not found: $SA_FILE"
[[ -f "$BWT_FILE" ]] || die "BWT file not found: $BWT_FILE"
SA_N=$(( $(stat -c%s "$SA_FILE") / 8 ))
ok "SA : $SA_N entries  ($(( $(stat -c%s "$SA_FILE") / 1073741824 )) GiB)"
ok "BWT: $(stat -c%s "$BWT_FILE") bytes"
ok "gsufsort time: ${T_GSUF}s"
[[ $SA_N -eq $(( HEX_COUNT + 2 )) ]] \
    && ok "n = HEX_COUNT + 2 ✓" \
    || warn "n=$SA_N but expected $((HEX_COUNT+2))"

# =============================================================================
step "[ 9/15 ] Trimming first entry from SA and BWT (in-place)"
# =============================================================================
set_label "trim SA+BWT in-place"
info "SA: remove first 8 bytes  |  BWT: remove first 1 byte"
info "Do NOT interrupt — files will be corrupted if stopped mid-way."
python3 - "$OUT_DIR/$SL_BASE" << 'PYEOF'
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
    print(f'\r    done — {new} bytes                    ')
b=sys.argv[1]
trim_inplace(b+'.8.sa', 8, 'SA ')
trim_inplace(b+'.bwt',  1, 'BWT')
PYEOF
EXPECTED_TRIM=$(( HEX_COUNT + 1 ))
SA_N_TRIM=$(( $(stat -c%s "$SA_FILE") / 8 ))
BWT_N_TRIM=$(stat -c%s "$BWT_FILE")
[[ $SA_N_TRIM -eq $EXPECTED_TRIM && $BWT_N_TRIM -eq $EXPECTED_TRIM ]] \
    && ok "Trim counts match: $SA_N_TRIM entries ✓" \
    || die "Trim mismatch: SA=$SA_N_TRIM BWT=$BWT_N_TRIM expected=$EXPECTED_TRIM"

# =============================================================================
step "[ 10/15 ] Verifying BWT character counts"
# =============================================================================
set_label "BWT verification"
python3 - "$BWT_FILE" "$HEX_COUNT" << 'PYEOF'
import numpy as np, sys
bwt,hc=sys.argv[1],int(sys.argv[2])
valid=list(range(48,58))+list(range(65,71))
counts={b:0 for b in valid+[0,1]}; counts['other']=0
with open(bwt,'rb') as f:
    while True:
        d=np.frombuffer(f.read(256*1024*1024),dtype=np.uint8)
        if not len(d): break
        for b in valid+[0,1]: counts[b]+=int(np.sum(d==b))
        counts['other']+=int(np.sum(~np.isin(d,valid+[0,1])))
total=sum(counts.values()); hex_total=sum(counts[b] for b in valid)
OK='\033[0;32m[ OK ]\033[0m'; FAIL='\033[0;31m[FAIL]\033[0m'
errs=0
def chk(l,c):
    global errs
    print(f'  {OK if c else FAIL}  {l}')
    if not c: errs+=1
chk(f'base-16 chars = {hc}', hex_total==hc)
chk('null(0x00) = 1', counts[0]==1)
chk('sep(0x01)  = 0', counts[1]==0)
chk('other      = 0', counts['other']==0)
chk(f'total      = {hc+1}', total==hc+1)
if errs: print(f'\033[0;31m{errs} FAILED\033[0m'); exit(1)
print('\033[0;32m  All BWT checks passed.\033[0m')
PYEOF

# =============================================================================
step "[ 11/15 ] Building compressed SA files (CF = 1, 2, 4, 8)"
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
step "[ 12/15 ] Cleaning up original SA file"
# =============================================================================
set_label "delete original SA"
if [[ $KEEP_SA -eq 0 ]]; then
    SA_SZ=$(du -h "$SA_FILE" | cut -f1)
    info "Deleting original SA ($SA_SZ) — CF files contain all needed data."
    rm -f "$SA_FILE"
    ok "Deleted: $SA_FILE"
else
    info "Keeping SA file (-k): $SA_FILE"
fi

# =============================================================================
step "[ 13/15 ] Building cp_occ_full.bin + c_vector.txt"
# =============================================================================
set_label "make_cp_occ_2step"
info "Single-pass over binary BWT — ~1 min.  Output: ~25 GB"
T=$(date +%s)
"$SCRIPT_DIR/make_cp_occ_2step" "$BWT_FILE" "$CP_OCC_FILE" "$CVEC_FILE" 2>&1
ok "Done in $(( $(date +%s) - T ))s"
ok "cp_occ_full.bin : $(du -h "$CP_OCC_FILE" | cut -f1)"
ok "c_vector.txt    : $(cat "$CVEC_FILE")"

# =============================================================================
step "[ 14/15 ] Building ref16_packed.bin"
# =============================================================================
set_label "make_ref16_packed (~3 GB RAM)"
info "Two-pass over base-16 reference.  Output: ~3 GB"
T=$(date +%s)
"$SCRIPT_DIR/make_ref16_packed" "$REF2STEP_FILE" "$REF16_FILE" 2>&1
ok "Done in $(( $(date +%s) - T ))s"
ok "ref16_packed.bin : $(du -h "$REF16_FILE" | cut -f1)"

# =============================================================================
step "[ 15/15 ] Building jump table(s)"
# =============================================================================
SA_LS_CF1="$OUT_DIR/sa_ls_word_cf1.bin"
SA_MS_CF1="$OUT_DIR/sa_ms_byte_cf1.bin"
[[ -f "$SA_LS_CF1" ]] || die "sa_ls_word_cf1.bin not found"
[[ -f "$SA_MS_CF1" ]] || die "sa_ms_byte_cf1.bin not found"
info "Memory: ~54 GB  (cp_occ 25 + sa_ls 24 + sa_ms 6)"

if [[ $DO_14 -eq 1 ]]; then
    set_label "jumptable_14nt (~54 GB RAM)"
    info "Generating jumptable_14nt.bin  (16^7=268M entries, 2 GiB) ..."
    T=$(date +%s)
    "$SCRIPT_DIR/make_jumptable_2step" \
        "$CP_OCC_FILE" "$SA_LS_CF1" "$SA_MS_CF1" "$CVEC_FILE" "$OUT_DIR" 14 2>&1
    ok "jumptable_14nt.bin: $(du -h "$OUT_DIR/jumptable_14nt.bin" | cut -f1)  in $(( $(date +%s)-T ))s"
fi

if [[ $DO_15 -eq 1 ]]; then
    set_label "jumptable_15nt (~54 GB RAM)"
    info "Generating jumptable_15nt.bin  (4^15=1.07B entries, 8 GiB) ..."
    T=$(date +%s)
    "$SCRIPT_DIR/make_jumptable_2step" \
        "$CP_OCC_FILE" "$SA_LS_CF1" "$SA_MS_CF1" "$CVEC_FILE" "$OUT_DIR" 15 2>&1
    ok "jumptable_15nt.bin: $(du -h "$OUT_DIR/jumptable_15nt.bin" | cut -f1)  in $(( $(date +%s)-T ))s"
fi

if [[ $DO_16 -eq 1 ]]; then
    set_label "jumptable_16nt (~54 GB RAM)"
    info "Generating jumptable_16nt.bin  (4^16=4.3B entries, 32 GiB) ..."
    T=$(date +%s)
    "$SCRIPT_DIR/make_jumptable_2step" \
        "$CP_OCC_FILE" "$SA_LS_CF1" "$SA_MS_CF1" "$CVEC_FILE" "$OUT_DIR" 16 2>&1
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
print_file "${INPUT_NOEXT}_2step_ref.txt"
print_file "${SL_BASE}.bwt"
echo ""
print_file "cp_occ_full.bin"
print_file "c_vector.txt"
print_file "ref16_packed.bin"
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
info "Index ready for 2-step RosaSeed."
info "Point your aligner config to: $OUT_DIR"
echo ""
rm -f "$PEAK_FILE" "$LABEL_FILE"
