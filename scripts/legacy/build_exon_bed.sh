#!/usr/bin/env bash
set -euo pipefail

# Build an exon BED from GENCODE GTF for gene-overlap annotation.
#
# Output columns:
#   chrom  start0  end  gene_name  gene_id  transcript_id
#
# Usage (WSL):
#   bash scripts/build_exon_bed.sh

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="${REF_DIR:-$ROOT/ref}"
GTF_GZ="${GTF_GZ:-$REF_DIR/gencode.v43.annotation.gtf.gz}"
OUT="${OUT_BED_GZ:-$REF_DIR/gencode_v43_exons.bed.gz}"

if [[ ! -f "$GTF_GZ" ]]; then
  echo "Missing: $GTF_GZ" >&2
  exit 1
fi

tmp="$(mktemp)"

zcat "$GTF_GZ" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="exon" {print $1, $4-1, $5, $9}' \
  | awk -F'\t' 'BEGIN{OFS="\t"}
      {
        chrom=$1; start=$2; end=$3; attr=$4;
        gene_name=""; gene_id=""; transcript_id="";
        if (match(attr, /gene_name \"([^\\\"]+)\"/, m)) gene_name=m[1];
        if (match(attr, /gene_id \"([^\\\"]+)\"/, m)) gene_id=m[1];
        if (match(attr, /transcript_id \"([^\\\"]+)\"/, m)) transcript_id=m[1];
        if (gene_name!="") print chrom, start, end, gene_name, gene_id, transcript_id;
      }' \
  | sort -k1,1 -k2,2n > "$tmp"

gzip -c "$tmp" > "$OUT"
rm -f "$tmp"

echo "Wrote: $OUT"
