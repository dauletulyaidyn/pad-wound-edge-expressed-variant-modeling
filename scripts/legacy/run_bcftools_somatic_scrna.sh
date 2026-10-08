#!/usr/bin/env bash
# Minimal bcftools pipeline for exploratory somatic calling on scRNA BAMs.
# Usage:
#   ./run_bcftools_somatic_scrna.sh \
#     --bam data/raw/SAMPLE.bam \
#     --ref ref/GRCh38.fa \
#     --sample SAMPLE_ID \
#     --out-prefix results/SAMPLE
#
# Notes:
# - Intended for expressed variants only; expect sparse coverage.
# - Adjust filters below to match desired depth/VAF/quality thresholds.
# - Optional: add --targets <BED> to limit to exonic/transcript regions.

set -euo pipefail

BAM=""
REF=""
SAMPLE=""
OUT_PREFIX=""
THREADS=${THREADS:-4}
TARGETS=""
REGIONS_FILE=""
KNOWN_SITES_VCF=""
COHORT_COMMON_REGIONS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --targets) TARGETS="--targets-file $2"; shift 2 ;;
    --regions-file) REGIONS_FILE="$2"; shift 2 ;;
    --known-sites-vcf) KNOWN_SITES_VCF="$2"; shift 2 ;;
    --cohort-common-regions) COHORT_COMMON_REGIONS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help)
      grep '^#' "$0" | sed 's/^# //'
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$BAM" || -z "$REF" || -z "$SAMPLE" || -z "$OUT_PREFIX" ]]; then
  echo "Missing required arguments. See --help." >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")"

LOG="${OUT_PREFIX}.bcftools.log"

# If BAM contains contigs not present in REF, bcftools can spam stderr.
# Default: restrict mpileup to contigs present in REF (from REF.fai).
if [[ -z "$REGIONS_FILE" ]]; then
  if [[ ! -f "${REF}.fai" ]]; then
    if command -v samtools >/dev/null 2>&1; then
      echo "[prep] indexing reference (samtools faidx) -> ${REF}.fai"
      samtools faidx "$REF" 2>>"$LOG"
    else
      echo "Missing ${REF}.fai and samtools not found to create it." >&2
      exit 1
    fi
  fi
  if command -v samtools >/dev/null 2>&1; then
    REGIONS_FILE="${OUT_PREFIX}.ref_bam_regions.tsv"
    # bcftools -R: regions file as CHROM<TAB>FROM<TAB>TO (1-based, inclusive).
    # Build intersection of contigs present in BAM and REF to avoid:
    # - REF missing contigs present in BAM (faidx errors)
    # - BAM missing contigs present in REF (mpileup errors)
    BAM_SNS="$(mktemp)"
    samtools view -H "$BAM" 2>>"$LOG" | awk -F'\t' '
      $1=="@SQ"{
        for(i=2;i<=NF;i++){
          if($i ~ /^SN:/){
            sn=$i; sub(/^SN:/,"",sn); print sn
          }
        }
      }' > "$BAM_SNS"
    awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=1;next} ($1 in a){print $1, 1, $2}' "$BAM_SNS" "${REF}.fai" > "$REGIONS_FILE"
    rm -f "$BAM_SNS"
  else
    REGIONS_FILE="${OUT_PREFIX}.ref_regions.tsv"
    awk 'BEGIN{OFS="\t"} {print $1, 1, $2}' "${REF}.fai" > "$REGIONS_FILE"
  fi
fi

if [[ -s "${OUT_PREFIX}.vcf" ]]; then
  echo "[1/3] raw VCF exists, skipping mpileup+call -> ${OUT_PREFIX}.vcf"
else
  echo "[1/3] mpileup + call -> ${OUT_PREFIX}.vcf"
  bcftools mpileup \
    -f "$REF" \
    -Ou \
    -a "FORMAT/DP,FORMAT/AD" \
    -R "$REGIONS_FILE" \
    $TARGETS \
    "$BAM" 2>>"$LOG" |
    bcftools call -mv -Ov -o "${OUT_PREFIX}.vcf" 2>>"$LOG"
fi

if [[ -s "${OUT_PREFIX}.filtered.vcf" ]]; then
  echo "[2/3] filtered VCF exists, skipping -> ${OUT_PREFIX}.filtered.vcf"
else
  echo "[2/3] filtering -> ${OUT_PREFIX}.filtered.vcf"
  # Basic scRNA-friendly filters: depth, qual, alt reads, VAF
  bcftools filter \
    -e 'FMT/DP<10 || QUAL<30 || (FMT/AD[0:1]/FMT/DP)<0.1 || FMT/AD[0:1]<3' \
    "${OUT_PREFIX}.vcf" \
    -Ov -o "${OUT_PREFIX}.filtered.vcf" 2>>"$LOG"
fi

FINAL_VCF="${OUT_PREFIX}.filtered.vcf"

# Optional: exclude known polymorphisms (dbSNP/gnomAD/etc).
# This does NOT prove somatic, but reduces obvious germline.
if [[ -n "$KNOWN_SITES_VCF" ]]; then
  if [[ -s "${OUT_PREFIX}.filtered.novel.vcf" ]]; then
    echo "[2b] novel VCF exists, skipping -> ${OUT_PREFIX}.filtered.novel.vcf"
  else
    echo "[2b] excluding known sites -> ${OUT_PREFIX}.filtered.novel.vcf"
    bcftools view -T "^${KNOWN_SITES_VCF}" "${OUT_PREFIX}.filtered.vcf" -Ov -o "${OUT_PREFIX}.filtered.novel.vcf" 2>>"$LOG"
  fi
  FINAL_VCF="${OUT_PREFIX}.filtered.novel.vcf"
fi

# Optional: exclude cohort-common positions (cross-sample "likely germline/artefact").
if [[ -n "$COHORT_COMMON_REGIONS" ]]; then
  if [[ -s "${OUT_PREFIX}.filtered.novel.cohort.vcf" ]]; then
    echo "[2c] cohort-filtered VCF exists, skipping -> ${OUT_PREFIX}.filtered.novel.cohort.vcf"
  else
    echo "[2c] excluding cohort-common regions -> ${OUT_PREFIX}.filtered.novel.cohort.vcf"
    bcftools view -T "^${COHORT_COMMON_REGIONS}" "$FINAL_VCF" -Ov -o "${OUT_PREFIX}.filtered.novel.cohort.vcf" 2>>"$LOG"
  fi
  FINAL_VCF="${OUT_PREFIX}.filtered.novel.cohort.vcf"
fi

if [[ -s "${FINAL_VCF}.stats.txt" ]]; then
  echo "[3/3] stats exists, skipping -> ${FINAL_VCF}.stats.txt"
else
  echo "[3/3] stats -> ${FINAL_VCF}.stats.txt"
  bcftools stats "$FINAL_VCF" > "${FINAL_VCF}.stats.txt" 2>>"$LOG"
fi

# Standard summary across samples
grep '^SN' "${FINAL_VCF}.stats.txt" > "${FINAL_VCF}.stats.summary.txt" || true

echo "Done. VCF: ${FINAL_VCF}"
