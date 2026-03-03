#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
from pathlib import Path
from dataclasses import dataclass
from typing import Optional
from Bio import SeqIO   # only needed if you want to make the panel here as well


@dataclass
class Config:
    dupcaller_path: str
    ref_fasta: str
    detag_path: str
    bam_folder: str

    umi_pattern: str = "NNNXXXX"
    pigz_threads: int = 16
    bwa_threads_per_sample: int = 12
    gatk_java_mem_gb: int = 24
    samtools_threads_per_index: int = 8


def run_cmd(cmd, shell: bool = False) -> None:
    if isinstance(cmd, list):
        printable = " ".join(str(c) for c in cmd)
    else:
        printable = cmd
    print(f"[CMD] {printable}")
    subprocess.run(cmd, shell=shell, check=True)


def dupcaller_trim(sample_id: str, r1: str, r2: str, cfg: Config) -> None:
    out_prefix = str(Path(cfg.detag_path) / sample_id)
    cmd = [
        "python", cfg.dupcaller_path, "trim",
        "-i", r1,
        "-i2", r2,
        "-p", cfg.umi_pattern,
        "-o", out_prefix,
    ]
    run_cmd(cmd)


def gzip_detagged(sample_id: str, cfg: Config) -> None:
    prefix = Path(cfg.detag_path) / sample_id
    for suffix in ("_1.fastq", "_2.fastq"):
        fq = prefix.with_name(prefix.name + suffix)
        if not fq.exists():
            raise FileNotFoundError(f"Expected detagged FASTQ not found: {fq}")
        cmd = ["pigz", "-p", str(cfg.pigz_threads), str(fq)]
        run_cmd(cmd)


def align_and_index(sample_id: str, cfg: Config) -> Path:
    r1 = Path(cfg.detag_path) / f"{sample_id}_1.fastq.gz"
    r2 = Path(cfg.detag_path) / f"{sample_id}_2.fastq.gz"
    if not r1.exists() or not r2.exists():
        raise FileNotFoundError(f"Detagged gz FASTQs not found for {sample_id}")

    bam = Path(cfg.bam_folder) / f"{sample_id}.bam"
    bam.parent.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"bwa-mem2 mem -C -t {cfg.bwa_threads_per_sample} "
        f"-R '@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA' "
        f"{cfg.ref_fasta} {r1} {r2} "
        f"| samtools sort -@ {cfg.bwa_threads_per_sample} -o {bam} - "
        f"&& samtools index -@ {cfg.bwa_threads_per_sample} {bam}"
    )
    run_cmd(cmd, shell=True)
    return bam


def mark_duplicates(sample_id: str, cfg: Config) -> Path:
    bam_in = Path(cfg.bam_folder) / f"{sample_id}.bam"
    bam_out = Path(cfg.bam_folder) / f"{sample_id}.mkdped.bam"
    metrics = Path(cfg.bam_folder) / f"{sample_id}.mkdp_metrics.txt"

    cmd = [
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
    run_cmd(cmd)

    cmd_idx = [
        "samtools", "index",
        "-@", str(cfg.samtools_threads_per_index),
        str(bam_out),
    ]
    run_cmd(cmd_idx)
    return bam_out


def main(argv: Optional[list] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Per-sample pipeline: DupCaller trim -> gzip -> bwa/samtools -> MarkDuplicates + index"
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--r1", required=True, help="R1 FASTQ (raw)")
    parser.add_argument("--r2", required=True, help="R2 FASTQ (raw)")

    # you can make these CLI args or just hardcode / env vars
    parser.add_argument("--detag-path", required=True)
    parser.add_argument("--bam-folder", required=True)
    parser.add_argument("--dupcaller-path", required=True)
    parser.add_argument("--ref-fasta", required=True)

    args = parser.parse_args(argv)

    cfg = Config(
        dupcaller_path=args.dupcaller_path,
        ref_fasta=args.ref_fasta,
        detag_path=args.detag_path,
        bam_folder=args.bam_folder,
    )

    Path(cfg.detag_path).mkdir(parents=True, exist_ok=True)
    Path(cfg.bam_folder).mkdir(parents=True, exist_ok=True)

    sample_id = args.sample_id
    print(f"[SAMPLE] Starting {sample_id}")

    dupcaller_trim(sample_id, args.r1, args.r2, cfg)
    gzip_detagged(sample_id, cfg)
    align_and_index(sample_id, cfg)
    mark_duplicates(sample_id, cfg)

    print(f"[SAMPLE] Finished {sample_id}")


if __name__ == "__main__":
    main()
