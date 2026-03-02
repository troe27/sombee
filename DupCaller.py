import os
import re
import time
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

import pandas as pd
from Bio import SeqIO


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

@dataclass
class PipelineConfig:
    # paths - software
    dupcaller_path: str
    ref_fasta: str

    # input
    raw_data_path: str
    sample_id_csv: str
    germline_vcf: str

    # output
    detag_path: str
    bam_folder: str
    output_folder: str
    panel_bed: str

    # UMI / DupCaller trim
    umi_pattern: str = "NNNXXXX"

    # pigz
    pigz_threads: int = 16

    # bwa / samtools
    bwa_threads_per_sample: int = 12
    max_parallel_samples: int = 6   # how many samples to process concurrently

    # gatk MarkDuplicates
    gatk_max_jobs: int = 4
    gatk_java_mem_gb: int = 24

    # samtools index for mkdp bams
    samtools_threads_per_index: int = 8
    samtools_index_max_jobs: int = 3

    # dupcaller call
    dupcaller_threads_per_job: int = 16
    dupcaller_max_jobs: int = 6


# ----------------------------------------------------------------------
# Utility: run commands with logging
# ----------------------------------------------------------------------

def run_cmd(cmd, shell: bool = False, cwd: Optional[str] = None) -> None:
    """
    Run a command via subprocess.run with basic logging and error checking.
    """
    if isinstance(cmd, list):
        printable = " ".join(str(c) for c in cmd)
    else:
        printable = cmd
    print(f"[CMD] {printable}")

    subprocess.run(cmd, shell=shell, cwd=cwd, check=True)


# ----------------------------------------------------------------------
# FASTQ discovery
# ----------------------------------------------------------------------

def collect_fastq_pairs_from_tree(
    root_dir: str,
    pattern: Optional[str] = None,
) -> Dict[str, Dict[str, List[str]]]:
    """
    Recursively walk root_dir, find FASTQ files, and group them into R1/R2
    per sample based ONLY on the file name.

    Returns:
        {
          "SAMPLE1": {
              "R1": ["/path/to/SAMPLE1_R1_001.fastq.gz", ...],
              "R2": ["/path/to/SAMPLE1_R2_001.fastq.gz", ...],
          },
          ...
        }
    """
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
        read = m.group("read")   # 'R1' or 'R2'

        samples.setdefault(sample, {}).setdefault(read, []).append(str(path.resolve()))

    return samples


# ----------------------------------------------------------------------
# Sample ID mapping (tumor/matched normal)
# ----------------------------------------------------------------------

def load_matched_normal_dict(sample_id_csv: str) -> Dict[str, Dict[str, str]]:
    """
    Build a dict:
      outer-key: 'Sample' (logical sample name)
      value: {'sample': NGI Sample ID (diluted),
              'matched_normal': NGI Sample ID (undiluted),
              'caste': caste}
    """
    df = pd.read_csv(sample_id_csv)
    mn_dict = {}

    for i in df["Sample"].drop_duplicates():
        diluted = df[(df["Sample"] == i) & (df["Dilution"] == "diluted")]
        undiluted = df[(df["Sample"] == i) & (df["Dilution"] == "undiluted")]

        s = diluted["NGI Sample ID"].iloc[0]
        caste = diluted["Caste"].iloc[0]
        mn = undiluted["NGI Sample ID"].iloc[0]

        mn_dict[i] = {"sample": s, "matched_normal": mn, "caste": caste}

    return mn_dict


# ----------------------------------------------------------------------
# Reference indexing & panel generation
# ----------------------------------------------------------------------

def index_reference(cfg: PipelineConfig) -> None:
    """
    Run DupCaller index on the reference once.
    """
    cmd = ["python", cfg.dupcaller_path, "index", "-f", cfg.ref_fasta]
    run_cmd(cmd)


def generate_panel_and_regions(cfg: PipelineConfig) -> List[str]:
    """
    Generate a BED panel that covers all chromosomes and return a list
    of region IDs (chromosome names).
    """
    regions_list: List[str] = []
    ref = cfg.ref_fasta
    panel = cfg.panel_bed

    Path(panel).parent.mkdir(parents=True, exist_ok=True)

    with open(panel, "wt") as handle:
        for record in SeqIO.parse(ref, "fasta"):
            chrom_id = record.id.split(" ")[0]
            stop = len(record.seq)
            regions_list.append(chrom_id)
            handle.write("\t".join([chrom_id, "0", str(stop)]) + "\n")

    return regions_list


# ----------------------------------------------------------------------
# Per-sample steps
# ----------------------------------------------------------------------

def dupcaller_trim_sample(
    sample_id: str,
    r1: str,
    r2: str,
    cfg: PipelineConfig,
) -> None:
    """
    Run DupCaller trim for one sample.
    """
    detag_prefix = str(Path(cfg.detag_path) / sample_id)
    cmd = [
        "python", cfg.dupcaller_path, "trim",
        "-i", r1,
        "-i2", r2,
        "-p", cfg.umi_pattern,
        "-o", detag_prefix,
    ]
    run_cmd(cmd)


def gzip_detagged_fastqs(sample_id: str, cfg: PipelineConfig) -> None:
    """
    Gzip the DupCaller-trimmed FASTQs for a sample using pigz.
    Assumes DupCaller produced <prefix>_1.fastq and <prefix>_2.fastq.
    """
    prefix = Path(cfg.detag_path) / sample_id
    for suffix in ("_1.fastq", "_2.fastq"):
        fq = prefix.with_name(prefix.name + suffix)
        if not fq.exists():
            raise FileNotFoundError(f"Expected detagged FASTQ not found: {fq}")
        cmd = ["pigz", "-p", str(cfg.pigz_threads), str(fq)]
        run_cmd(cmd)


def align_and_index_sample(
    sample_id: str,
    cfg: PipelineConfig,
    bwa_semaphore: threading.Semaphore,
) -> Path:
    """
    Run bwa mem -> samtools sort -> samtools index for one sample.
    Returns path to plain BAM.
    """
    detag_r1 = Path(cfg.detag_path) / f"{sample_id}_1.fastq.gz"
    detag_r2 = Path(cfg.detag_path) / f"{sample_id}_2.fastq.gz"

    if not detag_r1.exists() or not detag_r2.exists():
        raise FileNotFoundError(f"Detagged gz FASTQs not found for sample {sample_id}")

    bam_path = Path(cfg.bam_folder) / f"{sample_id}.bam"
    bam_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"bwa mem -C -t {cfg.bwa_threads_per_sample} "
        f"-R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
        f"{cfg.ref_fasta} {detag_r1} {detag_r2} "
        f"| samtools sort -@ {cfg.bwa_threads_per_sample} -o {bam_path} - "
        f"&& samtools index -@ {cfg.bwa_threads_per_sample} {bam_path}"
    )

    # limit number of bwa jobs at once
    with bwa_semaphore:
        run_cmd(cmd, shell=True)

    return bam_path


def mark_duplicates_sample(
    sample_id: str,
    cfg: PipelineConfig,
    gatk_semaphore: threading.Semaphore,
    index_semaphore: threading.Semaphore,
) -> Path:
    """
    Run GATK MarkDuplicates and index the mkdped BAM for one sample.
    Returns path to mkdped.bam.
    """
    bam_in = Path(cfg.bam_folder) / f"{sample_id}.bam"
    bam_out = Path(cfg.bam_folder) / f"{sample_id}.mkdped.bam"
    metrics = Path(cfg.bam_folder) / f"{sample_id}.mkdp_metrics.txt"

    if not bam_in.exists():
        raise FileNotFoundError(f"Input BAM for MarkDuplicates not found: {bam_in}")

    cmd_gatk = [
        "gatk",
        "--java-options", f"-Xmx{cfg.gatk_java_mem_gb}g",
        "MarkDuplicates",
        "-I", str(bam_in),
        "-O", str(bam_out),
        "-M", str(metrics),
        "--READ_NAME_REGEX",
        "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$",
        "--DUPLEX_UMI",
        "--TAGGING_POLICY", "OpticalOnly",
        "--BARCODE_TAG", "DB",
    ]

    with gatk_semaphore:
        run_cmd(cmd_gatk)

    # index mkdped bam with samtools
    cmd_index = [
        "samtools", "index",
        "-@", str(cfg.samtools_threads_per_index),
        str(bam_out),
    ]
    with index_semaphore:
        run_cmd(cmd_index)

    return bam_out


def process_single_sample(
    sample_id: str,
    reads: Dict[str, List[str]],
    cfg: PipelineConfig,
    bwa_semaphore: threading.Semaphore,
    gatk_semaphore: threading.Semaphore,
    index_semaphore: threading.Semaphore,
) -> None:
    """
    Full per-sample pipeline:
        DupCaller trim -> pigz -> bwa+samtools -> MarkDuplicates + index
    """
    print(f"[SAMPLE] Starting {sample_id}")

    # For simplicity we assume one R1/R2 file per sample (no lane splits).
    try:
        r1 = reads["R1"][0]
        r2 = reads["R2"][0]
    except (KeyError, IndexError):
        raise ValueError(f"Sample {sample_id} does not have both R1 and R2 entries")

    dupcaller_trim_sample(sample_id, r1, r2, cfg)
    gzip_detagged_fastqs(sample_id, cfg)
    align_and_index_sample(sample_id, cfg, bwa_semaphore)
    mark_duplicates_sample(sample_id, cfg, gatk_semaphore, index_semaphore)

    print(f"[SAMPLE] Finished {sample_id}")


# ----------------------------------------------------------------------
# DupCaller call (tumor/normal pairs)
# ----------------------------------------------------------------------

def run_dupcaller_calls(
    mn_dict: Dict[str, Dict[str, str]],
    regions_list: List[str],
    cfg: PipelineConfig,
) -> None:
    """
    Run DupCaller 'call' for each logical sample (tumor + matched normal),
    in parallel with a simple max-jobs throttling.
    """
    procs: List[subprocess.Popen] = []

    for logical_sample, item in mn_dict.items():
        sample_bam = Path(cfg.bam_folder) / f"{item['sample']}.mkdped.bam"
        normal_bam = Path(cfg.bam_folder) / f"{item['matched_normal']}.mkdped.bam"
        out_prefix = Path(cfg.output_folder) / f"{logical_sample}.dupcaller"

        cmd = [
            "python", cfg.dupcaller_path, "call",
            "-b", str(sample_bam),
            "-f", cfg.ref_fasta,
            "-r", *regions_list,
            "-o", str(out_prefix),
            "-p", str(cfg.dupcaller_threads_per_job),
            "-n", str(normal_bam),
            "-g", cfg.germline_vcf,
        ]

        print("[DupCaller call]", " ".join(cmd))
        p = subprocess.Popen(cmd)
        procs.append(p)

        # throttle dupcaller jobs
        while len(procs) >= cfg.dupcaller_max_jobs:
            procs = [pr for pr in procs if pr.poll() is None]
            if len(procs) >= cfg.dupcaller_max_jobs:
                time.sleep(5)

    # wait for remaining
    for p in procs:
        p.wait()


# ----------------------------------------------------------------------
# Top-level orchestrator
# ----------------------------------------------------------------------

def run_full_pipeline(cfg: PipelineConfig) -> None:
    """
    Orchestrate the whole pipeline:
      - create output dirs
      - find FASTQs
      - index reference + panel
      - per-sample pipeline in parallel (trim -> gzip -> map -> mkdup)
      - DupCaller call on tumor/normal pairs
    """
    # ensure output dirs exist
    Path(cfg.detag_path).mkdir(parents=True, exist_ok=True)
    Path(cfg.bam_folder).mkdir(parents=True, exist_ok=True)
    Path(cfg.output_folder).mkdir(parents=True, exist_ok=True)

    # discover FASTQs
    samples = collect_fastq_pairs_from_tree(cfg.raw_data_path)
    print(f"Discovered {len(samples)} samples")

    # load tumor/normal mapping
    mn_dict = load_matched_normal_dict(cfg.sample_id_csv)

    # index reference & build panel
    index_reference(cfg)
    regions_list = generate_panel_and_regions(cfg)

    # semaphores for step-specific concurrency
    bwa_semaphore = threading.Semaphore(cfg.max_parallel_samples)        # for bwa/samtools pipeline
    gatk_semaphore = threading.Semaphore(cfg.gatk_max_jobs)             # MarkDuplicates
    index_semaphore = threading.Semaphore(cfg.samtools_index_max_jobs)  # mkdped index

    # per-sample processing in parallel
    with ThreadPoolExecutor(max_workers=cfg.max_parallel_samples) as ex:
        futures = {
            ex.submit(
                process_single_sample,
                sample_id,
                reads,
                cfg,
                bwa_semaphore,
                gatk_semaphore,
                index_semaphore,
            ): sample_id
            for sample_id, reads in samples.items()
        }

        for fut in as_completed(futures):
            sid = futures[fut]
            try:
                fut.result()
            except Exception as e:
                print(f"[ERROR] Sample {sid} failed: {e}")

    # after all mkdped BAMs exist, run DupCaller call on pairs
    run_dupcaller_calls(mn_dict, regions_list, cfg)


# ----------------------------------------------------------------------
# Example usage (as a script)
# ----------------------------------------------------------------------

if __name__ == "__main__":
    cfg = PipelineConfig(
        dupcaller_path="/home/tilman/osmia_store/dupcaller-test/DupCaller/build/scripts-3.10/DupCaller.py",
        ref_fasta="/home/tilman/nanoseq_test/data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.fna",
        raw_data_path="/home/tilman/osmia_store/20250804_nanoseq/files/P35208/",
        sample_id_csv="/home/tilman/osmia_store/20250804_nanoseq/Beetle_IDs.csv",
        germline_vcf="/media/osmia/taliadoros/raw_data/aligned/sorted/RGed2/merged/BetterB_Batch1_2_3_full_geno_HQ.vcf.gz",
        detag_path="./test-data/detag",
        bam_folder="./test-data/bams",
        output_folder="./test-data/output",
        panel_bed="/home/tilman/nanoseq_test/data/ref/ref.bed",
        umi_pattern="NNNXXXX",
        pigz_threads=16,
        bwa_threads_per_sample=12,
        max_parallel_samples=6,
        gatk_max_jobs=4,
        gatk_java_mem_gb=24,
        samtools_threads_per_index=8,
        samtools_index_max_jobs=3,
        dupcaller_threads_per_job=16,
        dupcaller_max_jobs=6,
    )

    run_full_pipeline(cfg)