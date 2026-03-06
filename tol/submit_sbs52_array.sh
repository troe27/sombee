#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="germ.sample_names.txt"
MAX_CONCURRENT=50                 # throttle to be nice to the filesystem
VCF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/BetterB_Batch1_2_3_full_geno_HQ.vcf.gz"
REF="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic_aux_data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
OUTDIR="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/germ_sbs52_counts"
SCRIPT="./get_germ_sbs52_counts_per_sample.py"

mkdir -p "$OUTDIR"

NSAMPLES=$(wc -l < "$SAMPLES_FILE" | tr -d ' ')
if [[ "$NSAMPLES" -lt 1 ]]; then
  echo "ERROR: $SAMPLES_FILE is empty" >&2
  exit 1
fi

sbatch \
  --export=ALL,SAMPLES_FILE="$SAMPLES_FILE",VCF="$VCF",REF="$REF",OUTDIR="$OUTDIR",SCRIPT="$SCRIPT" \
  --array="1-${NSAMPLES}%${MAX_CONCURRENT}" \
  sbs52_array_task.sbatch
