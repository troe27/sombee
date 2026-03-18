#!/usr/bin/env python

import argparse
import csv
from pathlib import Path
import sys


NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_COLUMNS = [f"{sub},{up}-{down}" for sub in SUBS for up in NTS for down in NTS]


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Append one SBS96 TSV sample to the ToL raw somatic count table and write an HDP-ready matrix."
    )
    parser.add_argument("--reference-csv", type=Path, required=True, help="Large ToL somatic count CSV")
    parser.add_argument("--single-sample-tsv", type=Path, required=True, help="Single-sample SBS96 TSV")
    parser.add_argument("--sample-name", type=str, required=True, help="Sample name to use for the appended row")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output HDP matrix path")
    return parser.parse_args(argv[1:])


def load_reference_rows(path: Path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        rows = [row for row in reader if row]

    data_rows = [row for row in rows if not row[0].startswith("#")]
    if not data_rows:
        raise ValueError(f"No data rows found in {path}")

    header = data_rows[0]
    if header[0] != "Sample":
        raise ValueError(f"Unexpected first header column in {path}: {header[0]}")
    if header[1:] != SBS96_COLUMNS:
        raise ValueError("Reference CSV SBS96 columns do not match the expected canonical order")

    sample_rows = data_rows[1:]
    return sample_rows


def load_single_sample_counts(path: Path):
    counts = {column: 0 for column in SBS96_COLUMNS}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = ["SUBSTITUTION", "TRINUCLEOTIDE", "SBS96", "COUNT"]
        if reader.fieldnames != expected_fields:
            raise ValueError(f"Unexpected SBS96 TSV header in {path}: {reader.fieldnames}")
        for row in reader:
            sub = row["SUBSTITUTION"]
            tri = row["TRINUCLEOTIDE"]
            if len(tri) != 3:
                raise ValueError(f"Unexpected trinucleotide '{tri}' in {path}")
            key = f"{sub},{tri[0]}-{tri[2]}"
            if key not in counts:
                raise ValueError(f"Unexpected SBS96 key '{key}' derived from {path}")
            counts[key] = int(float(row["COUNT"]))
    return counts


def write_hdp_matrix(reference_rows, sample_name: str, sample_counts, output_path: Path):
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(SBS96_COLUMNS)
        for row in reference_rows:
            writer.writerow([row[0]] + row[1:])
        writer.writerow([sample_name] + [sample_counts[column] for column in SBS96_COLUMNS])


def main():
    args = parse_args(sys.argv)
    reference_rows = load_reference_rows(args.reference_csv)
    sample_counts = load_single_sample_counts(args.single_sample_tsv)
    write_hdp_matrix(reference_rows, args.sample_name, sample_counts, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
