#!/usr/bin/env bash
#SBATCH --job-name=fold_germ_maf
#SBATCH --account=UPPMAX2025-2-481
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --output=logs/fold_germ_maf_%A_%a.out
#SBATCH --error=logs/fold_germ_maf_%A_%a.err



IN_VCF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.vcf.gz"
TMP_VCF='./tmp.major_alt.vcf.gz'
IN_FASTA="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic_aux_data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
OUT_VCF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.major_allele_relative.vcf.gz"
OUT_FASTA="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/GCF_003254395.2_Amel_HAv3.1_genomic.BetterB_majorAllele.fna"

mkdir -p logs

module load BCFtools


bcftools +fill-tags ${IN_VCF} -- -t AF | bcftools view -i 'AF>0.5 && N_ALT=1 && TYPE="snp"' -Oz -o ${TMP_VCF}
tabix ${TMP_VCF}

bcftools consensus -f ${IN_FASTA} ${TMP_VCF} > ${OUT_FASTA}

bcftools norm -f ${OUT_FASTA} -c s ${IN_VCF} -Oz -o ${OUT_VCF} 

#rm ${TMP_VCF}

