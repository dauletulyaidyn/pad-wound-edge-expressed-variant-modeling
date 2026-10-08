#!/usr/bin/env bash
set -euo pipefail

# Build:
# - results/tables/variants_long.tsv.gz (per-variant per-sample metrics)
# - results/tables/variant_gene_long.tsv.gz (variant rows joined to overlapping exon gene_name via bedtools)
#
# Requires (WSL): bcftools, bedtools, gzip
#
# Usage:
#   bash scripts/build_variant_gene_tables.sh
#   bash scripts/build_variant_gene_tables.sh --vcf-glob 'results/SRR*.filtered.novel.cohort.vcf' --tag cohort

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RES="${RES_DIR:-$ROOT/results}"
TABLES="${TABLES_DIR:-$RES/tables}"
EXON_BED_GZ="${EXON_BED_GZ:-$ROOT/ref/gencode_v43_exons.bed.gz}"

VCF_GLOB="${RES}/*.filtered.vcf"
VCFS=()
TAG=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf-glob) VCF_GLOB="$2"; shift 2 ;;
    --vcf) VCFS+=("$2"); shift 2 ;;
    --tag) TAG="$2"; shift 2 ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# //'
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$TABLES"

if [[ ! -f "$EXON_BED_GZ" ]]; then
  echo "Missing exon bed: $EXON_BED_GZ" >&2
  echo "Run: bash ${BASE}/scripts/build_exon_bed.sh" >&2
  exit 1
fi

tmp_variants="$(mktemp)"
tmp_bed="$(mktemp)"
tmp_intersect="$(mktemp)"

out_variants="${TABLES}/variants_long.tsv.gz"
out_variant_gene="${TABLES}/variant_gene_long.tsv.gz"
if [[ -n "$TAG" ]]; then
  out_variants="${TABLES}/variants_long.${TAG}.tsv.gz"
  out_variant_gene="${TABLES}/variant_gene_long.${TAG}.tsv.gz"
fi

echo -e "sample\tchrom\tpos\tref\talt\tqual\tdp\tad_ref\tad_alt\tvaf" > "$tmp_variants"

if [[ ${#VCFS[@]} -eq 0 ]]; then
  shopt -s nullglob
  vcfs=( $VCF_GLOB )
  shopt -u nullglob
else
  vcfs=( "${VCFS[@]}" )
fi

for vcf in "${vcfs[@]}"; do
  [[ -f "$vcf" ]] || continue
  base="$(basename "$vcf")"
  [[ "$base" == *"_chr1."* ]] && continue
  sample="${base%%.*}"
  # bcftools query: DP and AD in FORMAT; AD is ref,alt for biallelic.
  bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%AD{0}]\t[%AD{1}]\n" "$vcf" \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        {
          dp=$7+0; ar=$8+0; aa=$9+0;
          vaf=(dp>0 ? aa/dp : 0);
          printf("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.6f\n",$1,$2,$3,$4,$5,$6,dp,ar,aa,vaf);
        }' >> "$tmp_variants"
done

gzip -c "$tmp_variants" > "$out_variants"
echo "Wrote: $out_variants"

# Build BED (0-based, half-open) for intersections
zcat "$out_variants" \
  | awk -F'\t' 'NR>1 {OFS="\t"; print $2, $3-1, $3, $1":"$2":"$3":"$4":"$5, $1, $6, $7, $10}' > "$tmp_bed"

bedtools intersect -a "$tmp_bed" -b <(zcat "$EXON_BED_GZ") -wa -wb > "$tmp_intersect"

# Emit joined table:
# a: chrom start end var_id sample qual dp vaf
# b: chrom start end gene_name gene_id transcript_id
{
  echo -e "chrom\tpos\tref\talt\tsample\tqual\tdp\tvaf\tgene_name\tgene_id\ttranscript_id"
  awk -F'\t' 'BEGIN{OFS="\t"}
    {
      # var_id: sample:chrom:pos:ref:alt
      split($4, a, ":");
      sample=a[1]; chrom=a[2]; pos=a[3]; ref=a[4]; alt=a[5];
      qual=$6; dp=$7; vaf=$8;
      gene_name=$12; gene_id=$13; transcript_id=$14;
      print chrom, pos, ref, alt, sample, qual, dp, vaf, gene_name, gene_id, transcript_id;
    }' "$tmp_intersect" \
    | sort -k5,5 -k1,1 -k2,2n
} | gzip -c > "$out_variant_gene"

echo "Wrote: $out_variant_gene"

rm -f "$tmp_variants" "$tmp_bed" "$tmp_intersect"
