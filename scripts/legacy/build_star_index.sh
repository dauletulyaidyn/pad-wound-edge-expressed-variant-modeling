#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

STAR_BIN="${STAR_BIN:-STAR}"
REF_FASTA="${REF_FASTA:-$ROOT/ref/GRCh38.primary_assembly.genome.fa}"
GTF_GZ="${GTF_GZ:-$ROOT/ref/gencode.v43.annotation.gtf.gz}"
STAR_INDEX_DIR="${STAR_INDEX_DIR:-$ROOT/ref/star_index}"
THREADS="${THREADS:-8}"

if ! command -v "$STAR_BIN" >/dev/null 2>&1; then
  echo "STAR not found. Set STAR_BIN or add STAR to PATH." >&2
  exit 1
fi
if [[ ! -s "$REF_FASTA" ]]; then
  echo "Missing REF_FASTA: $REF_FASTA" >&2
  exit 1
fi
if [[ ! -s "$GTF_GZ" ]]; then
  echo "Missing GTF_GZ: $GTF_GZ" >&2
  exit 1
fi

mkdir -p "$STAR_INDEX_DIR"

GTF="$GTF_GZ"
if [[ "$GTF_GZ" == *.gz ]]; then
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' EXIT
  GTF="$tmp/annotation.gtf"
  zcat "$GTF_GZ" > "$GTF"
fi

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX_DIR" \
  --genomeFastaFiles "$REF_FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang 100

echo "Done. STAR index: $STAR_INDEX_DIR"

