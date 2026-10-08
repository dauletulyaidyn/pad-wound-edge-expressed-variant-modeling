#!/usr/bin/env bash
set -euo pipefail

VCF_GLOB=""
VCFS=()
MIN_SAMPLES=${MIN_SAMPLES:-3}
OUT_REGIONS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf-glob) VCF_GLOB="$2"; shift 2 ;;
    --vcf) VCFS+=("$2"); shift 2 ;;
    --min-samples) MIN_SAMPLES="$2"; shift 2 ;;
    --out-regions) OUT_REGIONS="$2"; shift 2 ;;
    -h|--help)
      echo "Build cohort-common CHROM\\tFROM\\tTO regions from per-sample VCFs."
      echo "Usage: $0 --vcf-glob 'results/*.filtered.vcf' --min-samples 3 --out-regions results/tables/cohort_common_regions.tsv"
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$OUT_REGIONS" ]]; then
  echo "Missing required arg: --out-regions. See --help." >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_REGIONS")"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# Resolve VCF list
if [[ ${#VCFS[@]} -eq 0 ]]; then
  if [[ -z "$VCF_GLOB" ]]; then
    echo "Provide --vcf-glob or one/more --vcf paths." >&2
    exit 1
  fi
  # Expand glob safely inside this script.
  mapfile -t VCFS < <(compgen -G "$VCF_GLOB" || true)
fi

if [[ ${#VCFS[@]} -eq 0 ]]; then
  echo "No VCFs provided/matched. Glob: $VCF_GLOB" >&2
  exit 1
fi

# 1) For each VCF: collect unique CHROM:POS positions.
IDX=0

for vcf in "${VCFS[@]}"; do
  IDX=$((IDX+1))
  bcftools query -f '%CHROM\t%POS\n' "$vcf" 2>/dev/null | sort -u > "${TMP_DIR}/pos_${IDX}.tsv"
done

# 2) Count in how many samples each position appears.
cat "${TMP_DIR}"/pos_*.tsv | sort | uniq -c | awk -v min="$MIN_SAMPLES" '
  BEGIN{OFS="\t"}
  $1>=min{
    chrom=$2; pos=$3;
    print chrom, pos, pos
  }' > "$OUT_REGIONS"

echo "Wrote: $OUT_REGIONS"
echo "Min samples: $MIN_SAMPLES"
echo "VCFs: ${#VCFS[@]}"
