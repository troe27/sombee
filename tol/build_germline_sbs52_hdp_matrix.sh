#!/usr/bin/env bash
#SBATCH --job-name=hdp_test
#SBATCH --account=UPPMAX2025-2-481
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=logs/hdp_%A_%a.out
#SBATCH --error=logs/hdp_%A_%a.err

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mamba_env="tol3.9"
infile="${ROOT_DIR}/germ.test.sbs52.files.txt"
samplefile="${ROOT_DIR}/test20.germ.samples.txt"
#outdir="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/hdp"
outdir="${ROOT_DIR}"
outfile=${outdir}"/test.germ.hdp.csv"

mkdir -p ${outdir}
mkdir -p "${ROOT_DIR}/logs"
mamba run -n $mamba_env \
      python "${ROOT_DIR}/scripts/build_germline_sbs52_hdp_matrix.py" \
             --input $infile \
             --samples $samplefile \
             --output  $outfile
