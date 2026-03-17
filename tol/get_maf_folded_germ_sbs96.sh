#!/usr/bin/env bash
#SBATCH --job-name=maf_folded_germ
#SBATCH --account=UPPMAX2025-2-481
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --output=logs/maf_folded_germ_%A_%a.out
#SBATCH --error=logs/maf_folded_germ_%A_%a.err


set -euo pipefail
module load BCFtools
bcftools index -t /home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.major_allele_relative.vcf.gz

VCF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.major_allele_relative.vcf.gz"
REF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/GCF_003254395.2_Amel_HAv3.1_genomic.BetterB_majorAllele.fna"
OUT="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.major_allele_relative.sbs96.tsv}"

mamba run -n tol3.9 python "get_sbs96_counts.py" -i "${VCF}" --ref-fasta "${REF}" -o "${OUT}"
