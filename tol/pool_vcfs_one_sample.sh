#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
  echo "Usage: $0 <vcf_dir> <extra_vcf> [output_prefix] [sample_name]" >&2
  exit 1
fi

VCF_DIR="$1"
EXTRA_VCF="$2"
OUT_PREFIX="${3:-batch3_batch4_pool}"
SAMPLE_NAME="${4:-batch3_batch4_pool}"
TMPDIR="${OUT_PREFIX}.tmp"

mkdir -p "${TMPDIR}"

collect_bgzip_index() {
  local in="$1"
  local base
  base="$(basename "$in")"
  local staged="${TMPDIR}/${base}"
  local out

  if [[ "$in" == *.vcf.gz ]]; then
    cp "$in" "$staged"
  elif [[ "$in" == *.vcf ]]; then
    staged="${staged}.gz"
    bgzip -c "$in" > "$staged"
  else
    return 0
  fi

  out="${TMPDIR}/sorted_${base%.vcf.gz}.vcf.gz"
  out="${out%.vcf.vcf.gz}.vcf.gz"
  bcftools sort -Oz -o "$out" "$staged"
  bcftools index -t "$out"

  echo "$out"
}

INPUT_LIST="${TMPDIR}/inputs.txt"
: > "$INPUT_LIST"

while IFS= read -r f; do
  collect_bgzip_index "$f" >> "$INPUT_LIST"
done < <(find "$VCF_DIR" -type f \( -name '*.vcf' -o -name '*.vcf.gz' \) | sort)

collect_bgzip_index "$EXTRA_VCF" >> "$INPUT_LIST"

SITES_LIST="${TMPDIR}/sites_only.txt"
: > "$SITES_LIST"

i=0
while IFS= read -r vcf; do
  i=$((i + 1))
  out="${TMPDIR}/site_${i}.vcf.gz"
  bcftools view -G -v snps -Oz -o "$out" "$vcf"
  bcftools index -t "$out"
  echo "$out" >> "$SITES_LIST"
done < "$INPUT_LIST"

bcftools concat -a -Oz -o "${TMPDIR}/concat.vcf.gz" -f "$SITES_LIST"
bcftools sort -Oz -o "${TMPDIR}/sorted.vcf.gz" "${TMPDIR}/concat.vcf.gz"
bcftools index -t "${TMPDIR}/sorted.vcf.gz"

bcftools norm -d exact -Oz -o "${TMPDIR}/dedup.vcf.gz" "${TMPDIR}/sorted.vcf.gz"
bcftools index -t "${TMPDIR}/dedup.vcf.gz"

bcftools view -h "${TMPDIR}/dedup.vcf.gz" > "${TMPDIR}/header.vcf"
grep -v '^##FORMAT=<ID=GT,' "${TMPDIR}/header.vcf" | grep -v '^#CHROM' > "${TMPDIR}/header.nochrom.vcf"
printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' >> "${TMPDIR}/header.nochrom.vcf"
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$SAMPLE_NAME" >> "${TMPDIR}/header.nochrom.vcf"

bcftools view -H "${TMPDIR}/dedup.vcf.gz" \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,"GT","0/1"}' \
> "${TMPDIR}/body.vcf"

cat "${TMPDIR}/header.nochrom.vcf" "${TMPDIR}/body.vcf" \
| bgzip -c > "${OUT_PREFIX}.vcf.gz"

bcftools index -t "${OUT_PREFIX}.vcf.gz"

echo "Wrote ${OUT_PREFIX}.vcf.gz"
