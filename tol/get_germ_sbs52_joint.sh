#!/usr/bin/env bash
#SBATCH --job-name=germ_sbs52
#SBATCH --account=UPPMAX2025-2-481
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --output=logs/germ_sbs52_%A_%a.out
#SBATCH --error=logs/germ_sbs52_%A_%a.err

VCF_IN="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.vcf.gz"
REF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic_aux_data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
OUT_TSV="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/germ_sbs52_counts.joint.tsv"

mamba run -n tol3.9 python scripts/count_germline_sbs52.py -i ${VCF_IN} --ref-fasta ${REF}  -o ${OUT_TSV}
