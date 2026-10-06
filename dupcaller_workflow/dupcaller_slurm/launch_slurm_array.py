#!/usr/bin/env python3
import subprocess
from pathlib import Path
from typing import Dict, List
import argparse
import re


def collect_fastq_pairs_from_tree(root_dir: str, pattern: str = None) -> Dict[str, Dict[str, List[str]]]:
    from pathlib import Path
    if pattern is None:
        pattern = r'^(?P<sample>.+?)[._-](?P<read>R[12])(?:[_.\-]\d+)?\.f(ast)?q(\.gz)?$'
    regex = re.compile(pattern)
    root = Path(root_dir)
    samples: Dict[str, Dict[str, List[str]]] = {}
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if not any(str(path).endswith(ext) for ext in (".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            continue
        m = regex.match(path.name)
        if not m:
            continue
        sample = m.group("sample")
        read = m.group("read")
        samples.setdefault(sample, {}).setdefault(read, []).append(str(path.resolve()))
    return samples


def main():
    parser = argparse.ArgumentParser(description="Discover samples and submit Slurm array for per-sample pipeline.")
    parser.add_argument("--raw-data-path", required=True)
    parser.add_argument("--detag-path", required=True)
    parser.add_argument("--bam-folder", required=True)
    parser.add_argument("--dupcaller-path", required=True)
    parser.add_argument("--ref-fasta", required=True)
    parser.add_argument("--per-sample-script", default="per_sample_pipeline.py")
    parser.add_argument("--max-parallel", type=int, default=10, help="Max concurrent array tasks (% in Slurm)")
    parser.add_argument("--slurm-time", default="24:00:00")
    parser.add_argument("--slurm-cpus-per-task", type=int, default=16)
    parser.add_argument("--job-name", default="nanoseq")
    parser.add_argument("--samples-tsv", default="samples.tsv")
    parser.add_argument("--project-id",default=None,help="Slurm project/account ID to charge (adds #SBATCH --account=ID).",)
    parser.add_argument('--env', default=None, help='specify mamba/conda environment to use')
    parser.add_argument('--logdir', default='logs', help='directory to deposit the logfiles from slurm')
    args = parser.parse_args()

    # 1) discover samples
    samples = collect_fastq_pairs_from_tree(args.raw_data_path)
    print(f"Discovered {len(samples)} samples")

    # 2) write samples.tsv
    tsv_path = Path(args.samples_tsv)
    with tsv_path.open("w") as handle:
        handle.write("sample_id\tr1\tr2\n")
        for sample_id, reads in samples.items():
            r1 = reads.get("R1", [None])[0]
            r2 = reads.get("R2", [None])[0]
            if not (r1 and r2):
                print(f"Skipping {sample_id}: missing R1/R2")
                continue
            handle.write(f"{sample_id}\t{r1}\t{r2}\n")

    # count lines (minus header) to set array range
    n_samples = sum(1 for _ in tsv_path.open()) - 1
    if n_samples <= 0:
        print("No valid samples in samples.tsv; aborting")
        return

    array_spec = f"0-{n_samples - 1}%{args.max_parallel}"

    account_line = f"#SBATCH --account={args.project_id}\n" if args.project_id else ""
    if args.env:
       python_cmd = f"mamba run -n {args.env} python"
    else:
       python_cmd = "python"

    # 3) build sbatch script content
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={args.job_name}
{account_line}#SBATCH --cpus-per-task={args.slurm_cpus_per_task}
#SBATCH --time={args.slurm_time}
#SBATCH --array={array_spec}
#SBATCH --output={args.logdir}/{args.job_name}_%A_%a.out
#SBATCH --error={args.logdir}/{args.job_name}_%A_%a.err

set -euo pipefail

TASK_ID=${{SLURM_ARRAY_TASK_ID}}

SAMPLES_TSV="{tsv_path.resolve()}"

# skip header, get line TASK_ID+2 (because NR starts at 1 and line 1 is header)
LINE=$(awk -v line=$((TASK_ID+2)) 'NR==line' "$SAMPLES_TSV")
if [ -z "$LINE" ]; then
    echo "No line for TASK_ID=$TASK_ID" >&2
    exit 1
fi

sample_id=$(echo "$LINE" | cut -f1)
r1=$(echo "$LINE" | cut -f2)
r2=$(echo "$LINE" | cut -f3)

echo "Running sample: $sample_id"
echo " R1: $r1"
echo " R2: $r2"

# activate env here if needed:
module load GATK
module load bwa-mem2 
module load SAMtools

{python_cmd} {Path(args.per_sample_script).resolve()} \\
  --sample-id "$sample_id" \\
  --r1 "$r1" \\
  --r2 "$r2" \\
  --detag-path "{Path(args.detag_path).resolve()}" \\
  --bam-folder "{Path(args.bam_folder).resolve()}" \\
  --dupcaller-path "{Path(args.dupcaller_path).resolve()}" \\
  --ref-fasta "{Path(args.ref_fasta).resolve()}"
"""

    # 4) write sbatch script to disk
    Path(args.logdir).mkdir(exist_ok=True)
    sbatch_file = Path("run_array.sbatch")
    sbatch_file.write_text(sbatch_script)
    print(f"Wrote {sbatch_file}")

    # 5) submit it
    cmd = ["sbatch", str(sbatch_file)]
    print("Submitting:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
