#!/usr/bin/env bash
set -euo pipefail

# STAR alignment for 10x-like scRNA FASTQ files (WSL/Linux).
# Usage:
#   ./run_star_wsl.sh SRR14762239 SRR14762240 ...
# If no args are given, defaults to all `*_3.fastq` under FASTQ_DIR.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

STAR_BIN="${STAR_BIN:-STAR}"
GENOME_DIR="${STAR_INDEX_DIR:-$ROOT/ref/star_index}"
FASTQ_DIR="${FASTQ_DIR:-$ROOT/data/fastq}"
OUT_ROOT="${OUT_ROOT:-$ROOT}"
THREADS="${THREADS:-8}"
BAM_DIR="${BAM_DIR:-$OUT_ROOT/data/raw}"
LOG_DIR="${LOG_DIR:-$OUT_ROOT/logs}"

mkdir -p "$BAM_DIR" "$LOG_DIR"

if ! command -v "$STAR_BIN" >/dev/null 2>&1; then
  echo "STAR not found. Set STAR_BIN or add STAR to PATH." >&2
  exit 1
fi
if [[ ! -d "$GENOME_DIR" ]]; then
  echo "Missing STAR index dir: $GENOME_DIR" >&2
  echo "Build it with: bash scripts/build_star_index.sh" >&2
  exit 1
fi

if [[ $# -eq 0 ]]; then
  mapfile -t SAMPLES < <(cd "$FASTQ_DIR" && ls *_3.fastq 2>/dev/null | sed 's/_3.fastq//')
else
  SAMPLES=("$@")
fi

for s in "${SAMPLES[@]}"; do
  fq1="${FASTQ_DIR}/${s}_1.fastq"
  fq2="${FASTQ_DIR}/${s}_2.fastq"
  fq3="${FASTQ_DIR}/${s}_3.fastq"
  if [[ ! -f "$fq2" || ! -f "$fq3" ]]; then
    echo "Skipping ${s}: FASTQ missing" >&2
    continue
  fi

  prefix="${OUT_ROOT}/starsolo/${s}/"
  bam_out="${BAM_DIR}/${s}.bam"
  log_out="${LOG_DIR}/star_${s}.log"
  mkdir -p "$prefix"

  echo "[STAR] ${s} -> ${bam_out}"
  "$STAR_BIN" \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$fq3" "$fq2" \
    --soloType CB_UMI_Simple \
    --soloCBlen 16 \
    --soloUMIlen 12 \
    --soloFeatures Gene \
    --soloCBwhitelist None \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes CB UB GX GN NM NH HI AS nM CR UR \
    --limitBAMsortRAM 32000000000 \
    --outFileNamePrefix "$prefix" \
    >"$log_out" 2>&1

  star_bam="${prefix}Aligned.sortedByCoord.out.bam"
  if [[ -f "$star_bam" ]]; then
    mv "$star_bam" "$bam_out"
    samtools index -@ "$THREADS" "$bam_out"
    echo "[OK] ${s} BAM + BAI at ${bam_out}"
  else
    echo "[WARN] ${s} STAR BAM missing at ${star_bam}" >&2
  fi
done
