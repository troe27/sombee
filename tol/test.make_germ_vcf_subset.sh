#!/usr/bin/env bash

#SBATCH --job-name=germ.ss
#SBATCH --account=UPPMAX2025-2-481
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=logs/germ.ss_%A_%a.out
#SBATCH --error=logs/germ.ss_%A_%a.err


set -euo pipefail
module load BCFtools

IN_VCF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.vcf.gz"
SAMPLE_LIST="./germ.sample_names.txt"
N_SAMPLES=5
OUT_PREFIX="test5"

TEST_SAMPLES="${OUT_PREFIX}.samples.txt"
OUT_VCF="${OUT_PREFIX}.vcf.gz"

head -n "${N_SAMPLES}" "${SAMPLE_LIST}" > "${TEST_SAMPLES}"

bcftools view \
  --samples-file "${TEST_SAMPLES}" \
  --output-type z \
  --output-file "${OUT_VCF}" \
  "${IN_VCF}"

bcftools index -t "${OUT_VCF}"

echo "Wrote sample list: ${TEST_SAMPLES}"
echo "Wrote subset VCF: ${OUT_VCF}"
echo "Samples included:"
cat "${TEST_SAMPLES}"
