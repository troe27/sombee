#!/usr/bin/env bash
set -euo pipefail

VCF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic.pool.batch_3_batch4.vcf.gz"
REF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic_aux_data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
OUT="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic.pool.batch3_batch4.sbs96.tsv"

mamba run -n tol3.9 python get_sbs96_counts.py -i "${VCF}" --ref-fasta "${REF}" -o "${OUT}"
