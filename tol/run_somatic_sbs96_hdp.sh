#!/usr/bin/env bash
set -euo pipefail

module load R/4.4.2-gfbf-2024a
export R_LIBS_USER="${HOME}/Rlibs"


INPUT_MATRIX="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/tab10_raw_somatic_mutcount_tol.plus_batch3_batch4.hdp.tsv"
OUTPUT_DIR="/home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic_sbs96_hdp"
PREFIX="batch3_batch4.somatic.sbs96"
CHAINS=10


HDP_CHAIN_SCRIPT="./hdp_noprior_SBS96.R"
HDP_EXTRACT_SCRIPT="./hdp_extraction_SBS96.R"

mkdir -p "${OUTPUT_DIR}/chains"
mkdir -p "${OUTPUT_DIR}/results"



CHAIN_PREFIX="${OUTPUT_DIR}/chains/${PREFIX}_chain_"

for i in $(seq 1 "${CHAINS}"); do
  echo "Running SBS96 HDP chain ${i}/${CHAINS}"
  mamba run -n tol3.9 Rscript "${HDP_CHAIN_SCRIPT}" \
    "${INPUT_MATRIX}" \
    "${i}" \
    "${CHAIN_PREFIX}"
done

echo "Running SBS96 HDP extraction"

mamba run -n tol3.9 Rscript "${HDP_EXTRACT_SCRIPT}" \
  "${CHAIN_PREFIX}" \
  "${INPUT_MATRIX}" \
  "${OUTPUT_DIR}/results" \
  "${PREFIX}"

echo "Done"
