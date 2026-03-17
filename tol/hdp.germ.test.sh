#!/usr/bin/env bash
#SBATCH --job-name=hdp_chains
#SBATCH --account=UPPMAX2025-2-481
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --output=logs/hdp_chain_%A_%a.out
#SBATCH --error=logs/hdp_chain_%A_%a.err


set -euo pipefail
module load R/4.4.2-gfbf-2024a
export R_LIBS_USER=$HOME/Rlibs

# -------- SETTINGS --------
INPUT_MATRIX=./test.germ.hdp.csv
OUTPUT_DIR=/home/tilman/bees1/private/tilman/nanoseq_batch4/data/hdp_chains
script_loc=.
PREFIX=test20.germ.hdp_

CHAINS=10

mkdir -p "$OUTPUT_DIR/chains"
mkdir -p "$OUTPUT_DIR/results"

CHAIN_PREFIX="$OUTPUT_DIR/chains/${PREFIX}_chain_"

echo "Running $CHAINS HDP chains..."

for i in $(seq 1 $CHAINS); do
    echo "Running chain $i"
    
    mamba run -n tol3.9 Rscript ${script_loc}/hdp_noprior_SBS52.R \
        "$INPUT_MATRIX" \
        "$i" \
        "$CHAIN_PREFIX"
done

echo "Running extraction step..."

mamba run -n tol3.9 Rscript ${script_loc}/hdp_extraction_SBS52.R \
    "$CHAIN_PREFIX" \
    "$INPUT_MATRIX" \
    "$OUTPUT_DIR/results" \
    "$PREFIX"

echo "Done."
