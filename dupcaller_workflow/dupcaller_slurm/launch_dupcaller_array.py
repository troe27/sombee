#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Tuple

import pandas as pd
from Bio import SeqIO


def load_sample_normal_dict(sample_id_tsv: str) -> Dict[str, Dict[str, str]]:
    """
    Input TSV must contain:
      - NGI_ID (used to locate BAM: <NGI_ID>.mkdped.bam)
      - sample_type ('sample' or 'matched_normal')
      - sample_ID (logical pairing id)

    Returns:
      { sample_ID: {'sample': <NGI_ID>, 'normal': <NGI_ID>} }
    """
    df = pd.read_csv(sample_id_tsv, sep="\t", dtype=str)

    required = {"NGI_ID", "sample_type", "sample_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {sample_id_tsv}: {sorted(missing)}")

    df = df.copy()
    df["NGI_ID"] = df["NGI_ID"].str.strip()
    df["sample_type"] = df["sample_type"].str.strip()
    df["sample_ID"] = df["sample_ID"].str.strip()

    sn: Dict[str, Dict[str, str]] = {}

    for sample_id, sub in df.groupby("sample_ID", sort=False):
        sample_rows = sub[sub["sample_type"] == "sample"]
        normal_rows = sub[sub["sample_type"] == "matched_normal"]

        if sample_rows.empty or normal_rows.empty:
            print(f"[WARN] sample_ID {sample_id}: missing sample or matched_normal entry, skipping")
            continue
        if len(sample_rows) > 1:
            print(f"[WARN] sample_ID {sample_id}: multiple 'sample' rows; using first")
        if len(normal_rows) > 1:
            print(f"[WARN] sample_ID {sample_id}: multiple 'matched_normal' rows; using first")

        sn[sample_id] = {
            "sample": sample_rows["NGI_ID"].iloc[0],
            "normal": normal_rows["NGI_ID"].iloc[0],
        }

    return sn


def ensure_panel_bed(ref_fasta: str, panel_bed: str) -> None:
    """
    If panel_bed doesn't exist, create a BED that spans each contig in ref_fasta:
      <contig>  0  <len>
    """
    panel_path = Path(panel_bed)
    if panel_path.exists():
        print(f"[INFO] Panel BED exists: {panel_path}")
        return

    print(f"[INFO] Creating panel BED: {panel_path}")
    panel_path.parent.mkdir(parents=True, exist_ok=True)

    with panel_path.open("w") as out:
        for rec in SeqIO.parse(ref_fasta, "fasta"):
            contig = rec.id.split(" ")[0]
            out.write(f"{contig}\t0\t{len(rec.seq)}\n")


def build_pairs_tsv(
    sn_dict: Dict[str, Dict[str, str]],
    bam_folder: Path,
    pairs_path: Path,
    ) -> int:
    """
    Write pairs TSV: sample_id, sample_bam, normal_bam.
    Returns number of valid pairs written.
    """
    n_pairs = 0
    with pairs_path.open("w") as out:
        out.write("sample_id\tsample_bam\tnormal_bam\n")

        for sample_id, ids in sn_dict.items():
                sample_matches = list(bam_folder.glob(f"{ids['sample']}_*.mkdped.bam"))
                normal_matches = list(bam_folder.glob(f"{ids['normal']}_*.mkdped.bam"))

                if len(sample_matches) == 0:
                        print(f"[WARN] No sample BAM found for {sample_id} (pattern {ids['sample']}_*.mkped.bam), skipping")
                        continue

                if len(normal_matches) == 0:
                        print(f"[WARN] No normal BAM found for {sample_id} (pattern {ids['normal']}_*.mkped.bam), skipping")
                        continue

                if len(sample_matches) > 1:
                       print(f"[WARN] Multiple sample BAMs found for {sample_id}, using first: {sample_matches[0]}")

                if len(normal_matches) > 1:
                      print(f"[WARN] Multiple normal BAMs found for {sample_id}, using first: {normal_matches[0]}")

                sample_bam = sample_matches[0]
                normal_bam = normal_matches[0]

                out.write(f"{sample_id}\t{sample_bam}\t{normal_bam}\n")
                n_pairs += 1

    return n_pairs


def write_sbatch_script(
    sbatch_file: Path,
    *,
    job_name: str,
    account_line: str,
    cpus_per_task: int,
    slurm_time: str,
    array_spec: str,
    logdir: str,
    pairs_path: Path,
    dupcaller_path: Path,
    ref_fasta: Path,
    germline_vcf: Path,
    panel_bed: Path,
    output_folder: Path,
    env: str
) -> None:
    """
    Create a Slurm array script that reads pairs.tsv and runs DupCaller call.
    """
    if env is None:
        python_cmd = 'python'
    else:
        python_cmd = f'mamba run -n {env} python'
    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
{account_line}#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --time={slurm_time}
#SBATCH --array={array_spec}
#SBATCH --output={logdir}/{job_name}_%A_%a.out
#SBATCH --error={logdir}/{job_name}_%A_%a.err

set -euo pipefail


TASK_ID=${{SLURM_ARRAY_TASK_ID}}
PAIRS_TSV="{pairs_path}"
DUPCALLER_PATH="{dupcaller_path}"
REF_FASTA="{ref_fasta}"
GERMLINE_VCF="{germline_vcf}"
PANEL_BED="{panel_bed}"
OUTPUT_FOLDER="{output_folder}"
THREADS={cpus_per_task}

# line for this task (skip header)
LINE=$(awk -v line=$((TASK_ID+2)) 'NR==line' "$PAIRS_TSV")
if [ -z "$LINE" ]; then
    echo "No line for TASK_ID=$TASK_ID" >&2
    exit 1
fi

sample_id=$(echo "$LINE" | cut -f1)
sample_bam=$(echo "$LINE" | cut -f2)
normal_bam=$(echo "$LINE" | cut -f3)

echo "Running DupCaller call for sample: $sample_id"
echo " Sample BAM: $sample_bam"
echo " Normal BAM: $normal_bam"

# regions: contig IDs from panel BED (col 1)
regions=$(cut -f1 "$PANEL_BED" | tr '\\n' ' ')

out_prefix="$OUTPUT_FOLDER/${{sample_id}}.dupcaller"

{python_cmd} "$DUPCALLER_PATH" call \\
  -b "$sample_bam" \\
  -f "$REF_FASTA" \\
  -r $regions \\
  -o "$out_prefix" \\
  -p "$THREADS" \\
  -n "$normal_bam" \\
  -g "$GERMLINE_VCF"
"""
    sbatch_file.write_text(script)
    print(f"[INFO] Wrote {sbatch_file}")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Prepare DupCaller sample/normal pairs and submit Slurm array for 'call' step."
    )
    p.add_argument(
        "--sample-id-csv",
        required=True,
        help="TSV with NGI_ID, sample_type (sample/matched_normal), sample_ID (logical id).",
    )
    p.add_argument(
        "--bam-folder",
        required=True,
        help="Folder with *.mkdped.bam named by NGI_ID (e.g. P38907_1008.mkdped.bam).",
    )
    p.add_argument("--dupcaller-path", required=True, help="Path to DupCaller.py")
    p.add_argument("--ref-fasta", required=True, help="Reference FASTA used for alignment.")
    p.add_argument("--germline-vcf", required=True, help="Germline VCF for DupCaller.")
    p.add_argument("--panel-bed", required=True, help="BED file with panel regions; created if missing.")
    p.add_argument("--output-folder", required=True, help="Folder for DupCaller call outputs.")

    p.add_argument("--pairs-tsv", default="pairs.tsv", help="Output TSV listing sample_id, sample_bam, normal_bam.")
    p.add_argument("--job-name", default="dupcaller_call", help="Slurm job name.")
    p.add_argument("--max-parallel", type=int, default=6, help="Max concurrent array tasks (% limit).")
    p.add_argument("--slurm-time", default="24:00:00")
    p.add_argument("--slurm-cpus-per-task", type=int, default=16)
    p.add_argument("--project-id", default=None, help="Slurm project/account ID (#SBATCH --account=ID).")
    p.add_argument("--logdir", default="logs", help="Directory for Slurm stdout/stderr logs.")
    p.add_argument('--env', default=None, help='name of mamba environment to run command in')
    args = p.parse_args()

    bam_folder = Path(args.bam_folder)
    output_folder = Path(args.output_folder)
    pairs_path = Path(args.pairs_tsv)
    sbatch_file = Path("run_dupcaller_array.sbatch")

    bam_folder.mkdir(parents=True, exist_ok=True)
    output_folder.mkdir(parents=True, exist_ok=True)
    Path(args.logdir).mkdir(parents=True, exist_ok=True)

    # 1) mapping
    sn_dict = load_sample_normal_dict(args.sample_id_csv)
    print(f"[INFO] Loaded {len(sn_dict)} logical samples with matched normals")

    # 2) panel BED
    ensure_panel_bed(args.ref_fasta, args.panel_bed)

    # 3) pairs.tsv
    n_pairs = build_pairs_tsv(sn_dict, bam_folder, pairs_path)
    if n_pairs == 0:
        print("[ERROR] No valid sample/normal pairs found; aborting.")
        return
    print(f"[INFO] Wrote {n_pairs} pairs to {pairs_path.resolve()}")

    # 4) sbatch script
    array_spec = f"0-{n_pairs - 1}%{args.max_parallel}"
    account_line = f"#SBATCH --account={args.project_id}\n" if args.project_id else ""

    write_sbatch_script(
        sbatch_file,
        job_name=args.job_name,
        account_line=account_line,
        cpus_per_task=args.slurm_cpus_per_task,
        slurm_time=args.slurm_time,
        array_spec=array_spec,
        logdir=args.logdir,
        pairs_path=pairs_path.resolve(),
        dupcaller_path=Path(args.dupcaller_path).resolve(),
        ref_fasta=Path(args.ref_fasta).resolve(),
        germline_vcf=Path(args.germline_vcf).resolve(),
        panel_bed=Path(args.panel_bed).resolve(),
        output_folder=output_folder.resolve(),
        env = args.env,
    )

    # 5) submit
    cmd = ["sbatch", str(sbatch_file)]
    print("[INFO] Submitting:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
