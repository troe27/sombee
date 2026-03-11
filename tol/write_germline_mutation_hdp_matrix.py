#!/usr/bin/env python

import argparse
import csv
import math
import os
from pathlib import Path
from typing import List
import sys

import natsort

SBS52_COLUMNS = [
    "C>A,A-C", "C>A,A-T", "C>A,C-A", "C>A,C-C", "C>A,C-G", "C>A,C-T", "C>A,G-C", "C>A,T-A", "C>A,T-C", "C>A,T-T",
    "C>G,A-C", "C>G,A-T", "C>G,C-A", "C>G,C-C", "C>G,C-G", "C>G,C-T", "C>G,G-C", "C>G,T-A", "C>G,T-C", "C>G,T-T",
    "C>T,A-A", "C>T,A-C", "C>T,A-G", "C>T,A-T", "C>T,C-A", "C>T,C-C", "C>T,C-G", "C>T,C-T", "C>T,G-A", "C>T,G-C", "C>T,G-G", "C>T,G-T", "C>T,T-A", "C>T,T-C", "C>T,T-G", "C>T,T-T",
    "T>A,A-C", "T>A,A-T", "T>A,C-A", "T>A,C-C", "T>A,C-G", "T>A,C-T", "T>A,G-C", "T>A,T-A", "T>A,T-C", "T>A,T-T",
    "T>G,A-C", "T>G,C-A", "T>G,C-C", "T>G,C-T", "T>G,T-C", "T>G,T-T",
]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Generate an HDP input matrix from per-sample germline mutation count files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="File containing full paths to germline mutation count files, one per line.",
    )
    parser.add_argument(
        "--samples",
        type=Path,
        required=True,
        help="File containing sample names to include in the HDP matrix, one per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output HDP matrix CSV file.",
    )
    args = args[1:]
    return parser.parse_args(args)


def get_common_suffix(strings: List[str]) -> str:
    if not strings:
        return ""
    rev_strings = [s[::-1] for s in strings]
    rev_prefix = os.path.commonprefix(rev_strings)
    return rev_prefix[::-1]


def read_nonempty_lines(path: Path) -> List[str]:
    with open(path) as infile:
        return [line.strip() for line in infile if line.strip()]


def write_hdp_matrix(input_path: Path, sample_path: Path, output_path: Path):
    # Load file paths
    file_paths = read_nonempty_lines(input_path)
    if not file_paths:
        raise ValueError(f"No input file paths found in {input_path}")

    # Work on basenames, not full paths
    file_names = [Path(p).name for p in file_paths]
    suffix = get_common_suffix(file_names)

    # Load target samples
    target_samples = set(read_nonempty_lines(sample_path))
    if not target_samples:
        raise ValueError(f"No sample names found in {sample_path}")

    # Aggregate counts
    sample_sbs52_counts = {}
    matched_samples = set()
    observed_classes = set()

    for file_path_str in file_paths:
        file_path = Path(file_path_str)
        file_name = file_path.name

        if suffix:
            query_sample = file_name.removesuffix(suffix)
        else:
            query_sample = file_name

        if query_sample not in target_samples:
            continue

        matched_samples.add(query_sample)

        with open(file_path) as infile:
            for line in infile:
                if line.startswith("SUBSTITUTION"):
                    continue

                fields = line.rstrip().split("\t")
                if len(fields) < 4:
                    continue

                sub = fields[0]
                tri = fields[1].split()[-1]   # e.g. ACC
                count = fields[3]

                if len(tri) != 3:
                    continue

                upstream, _, downstream = tri
                sbs52 = f"{sub},{upstream}-{downstream}"
                if sbs52 not in SBS52_COLUMNS:
                    raise ValueError(f"Unrecognized SBS52 channel '{sbs52}' in {file_path}")

                sample_sbs52_counts[f"{query_sample}:{sbs52}"] = count
                observed_classes.add(sbs52)

    output_samples = natsort.natsorted(list(target_samples))

    # Write CSV matrix
    with open(output_path, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["sample"] + SBS52_COLUMNS)

        for sample in output_samples:
            row_values = [sample] + [
                str(math.ceil(float(sample_sbs52_counts.get(f"{sample}:{sbs52}", 0))))
                for sbs52 in SBS52_COLUMNS
            ]
            writer.writerow(row_values)

    # Helpful stderr summary
    missing_samples = sorted(target_samples - matched_samples)
    print(f"Loaded {len(file_paths)} input files", file=sys.stderr)
    print(f"Matched {len(matched_samples)} / {len(target_samples)} requested samples", file=sys.stderr)
    print(f"Observed {len(observed_classes)} / {len(SBS52_COLUMNS)} SBS52 classes", file=sys.stderr)
    print(f"Wrote matrix to: {output_path}", file=sys.stderr)

    if missing_samples:
        print(
            "Warning: the following requested samples were not matched to any input file:",
            file=sys.stderr,
        )
        for sample in missing_samples:
            print(f"  {sample}", file=sys.stderr)


def main():
    options = parse_args(sys.argv)
    write_hdp_matrix(options.input, options.samples, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
