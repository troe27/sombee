#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
from typing import Dict, Any

import pandas as pd
from Bio import SeqIO


def load_matched_normal_dict(sample_id_csv: str) -> Dict[str, Dict[str, Any]]:
    """
    Build a dict:
      outer-key: 'Sample' (logical sample name)
      value: {'sample': NGI Sample ID (diluted),
              'matched_normal': NGI Sample ID (undiluted),
              'caste': caste}
    """
    df = pd.read_csv(sample_id_csv)
    mn_dict: Dict[str, Dict[str, Any]] = {}

    for i in df["Sample"].drop_duplicates():
        diluted = df[(df["Sample"] == i) & (df["Dilution"] == "diluted")]
        undiluted = df[(df["Sample"] == i) & (df["Dilution"] == "undiluted")]

        if diluted.empty or undiluted.empty:
            print(f"[WARN] Sample {i}: missing diluted or undiluted entry, skipping")
            continue

        s = diluted["NGI Sample ID"].iloc[0]
        caste = diluted["Caste"].iloc[0]
        mn = undiluted["NGI Sample ID"].iloc[0]

        mn_dict[i] = {"sample": s, "matched_normal": mn, "caste": caste}

    return mn_dict


def generate_panel_if_missing(ref_fasta: str, panel_bed: str) -> None:
    """
    If panel_bed does not exist, generate a BED that covers all chromosomes
    in ref_fasta: one row per record: <chrom>  0  len(seq)
    """
    panel_path = Path(panel_bed)
    if panel_path.exists():
        print(f"[INFO] Panel BED exists: {panel_path}")
        return

    print(f"[INFO] Generating panel BED at {panel_path}")
    panel_path.parent.mkdir(parents=True, exist_ok=True)

    with panel_path.open("w") as handle:
        for record in SeqIO.parse(ref_fasta, "fasta"):
            chrom_id = record.id.split(" ")[0]
            stop = len(record.seq)
            handle.write("\t".join([chrom_id, "0", str(stop)]) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare DupCaller tumor/normal pairs and submit Slurm array for 'call' step."
    )
    parser.add_argument("--sample-id-csv", required=True,
                        help="Beetle_IDs.csv with Sample / Dilution / NGI Sample ID / Caste.")
    parser.add_argument("--bam-folder", required=True,
                        help="Folder with *.mkdped.bam files from per-sample pipeline.")
    parser.add_argument("--dupcaller-path", required=True,
                        help="Path to DupCaller.py")
    parser.add_argument("--ref-fasta", required=True,
                        help="Reference FASTA used for alignment.")
    parser.add_argument("--germline-vcf", required=True,
                        help="Germline VCF for DupCaller.")
    parser.add_argument("--panel-bed", required=True,
                        help="BED file with panel regions; will be created if missing.")
    parser.add_argument("--output-folder", required=True,
                        help="Folder where DupCaller call outputs should go.")

    parser.add_argument("--pairs-tsv", default="pairs.tsv",
                        help="Output TSV listing logical_sample, tumor_bam, normal_bam.")
    parser.add_argument("--job-name", default="dupcaller_call",
                        help="Slurm job name.")
    parser.add_argument("--max-parallel", type=int, default=6,
                        help="Max concurrent array tasks (% limit).")
    parser.add_argument("--slurm-time", default="24:00:00")
    parser.add_argument("--slurm-cpus-per-task", type=int, default=16)
    parser.add_argument(
        "--project-id",
        default=None,
        help="Slurm project/account ID to charge (adds #SBATCH --account=ID).",
    )

    args = parser.parse_args()

    bam_folder = Path(args.bam_folder)
    bam_folder.mkdir(parents=True, exist_ok=True)
    output_folder = Path(args.output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    # 1) load tumor/normal mapping
    mn_dict = load_matched_normal_dict(args.sample_id_csv)
    print(f"[INFO] Loaded {len(mn_dict)} logical samples with matched normals")

    # 2) generate panel bed if needed
    generate_panel_if_missing(args.ref_fasta, args.panel_bed)

    # 3) build pairs.tsv
    pairs_path = Path(args.pairs_tsv)
    with pairs_path.open("w") as handle:
        handle.write("logical_sample\ttumor_bam\tnormal_bam\n")
        n_pairs = 0

        for logical_sample, item in mn_dict.items():
            tumor_bam = bam_folder / f"{item['sample']}.mkdped.bam"
            normal_bam = bam_folder / f"{item['matched_normal']}.mkdped.bam"

            if not tumor_bam.exists():
                print(f"[WARN] Tumor BAM missing for {logical_sample}: {tumor_bam}, skipping")
                continue
            if not normal_bam.exists():
                print(f"[WARN] Normal BAM missing for {logical_sample}: {normal_bam}, skipping")
                continue

            handle.write(f"{logical_sample}\t{tumor_bam}\t{normal_bam}\n")
            n_pairs += 1

    if n_pairs == 0:
        print("[ERROR] No valid tumor/normal pairs found; aborting.")
        return

    print(f"[INFO] Wrote {n_pairs} pairs to {pairs_path}")

    array_spec = f"0-{n_pairs - 1}%{args.max_parallel}"
    account_line = f"#SBATCH --account={args.project_id}\n" if args.project_id else ""

    # 4) build sbatch script
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={args.job_name}
{account_line}#SBATCH --cpus-per-task={args.slurm_cpus_per_task}
#SBATCH --time={args.slurm_time}
#SBATCH --array={array_spec}
#SBATCH --output=logs/{args.job_name}_%A_%a.out
#SBATCH --error=logs/{args.job_name}_%A_%a.err

set -euo pipefail

module load GATK
module load BWA-mem2

TASK_ID=${{SLURM_ARRAY_TASK_ID}}
PAIRS_TSV="{pairs_path.resolve()}"
DUPCALLER_PATH="{Path(args.dupcaller_path).resolve()}"
REF_FASTA="{Path(args.ref_fasta).resolve()}"
GERMLINE_VCF="{Path(args.germline_vcf).resolve()}"
PANEL_BED="{Path(args.panel_bed).resolve()}"
OUTPUT_FOLDER="{output_folder.resolve()}"
THREADS={args.slurm_cpus_per_task}

# get the line for this task (skip header)
LINE=$(awk -v line=$((TASK_ID+2)) 'NR==line' "$PAIRS_TSV")
if [ -z "$LINE" ]; then
    echo "No line for TASK_ID=$TASK_ID" >&2
    exit 1
fi

logical_sample=$(echo "$LINE" | cut -f1)
tumor_bam=$(echo "$LINE" | cut -f2)
normal_bam=$(echo "$LINE" | cut -f3)

echo "Running DupCaller call for logical sample: $logical_sample"
echo " Tumor:  $tumor_bam"
echo " Normal: $normal_bam"

# regions: all chrom IDs from panel BED (col 1)
regions=$(cut -f1 "$PANEL_BED" | tr '\\n' ' ')

out_prefix="$OUTPUT_FOLDER/${{logical_sample}}.dupcaller"

# activate env here if needed
# module load ...
# conda activate ...

python "$DUPCALLER_PATH" call \\
  -b "$tumor_bam" \\
  -f "$REF_FASTA" \\
  -r $regions \\
  -o "$out_prefix" \\
  -p "$THREADS" \\
  -n "$normal_bam" \\
  -g "$GERMLINE_VCF"
"""

    Path("logs").mkdir(exist_ok=True)
    sbatch_file = Path("run_dupcaller_array.sbatch")
    sbatch_file.write_text(sbatch_script)
    print(f"[INFO] Wrote {sbatch_file}")

    # 5) submit
    cmd = ["sbatch", str(sbatch_file)]
    print("[INFO] Submitting:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
