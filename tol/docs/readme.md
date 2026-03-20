# README

This folder contains local workflow scripts in `scripts/`, runnable wrappers in the repository root, notes in `docs/`, and the upstream reference code in `reference/`.

## Germline

### SBS52

Main local scripts:

- `scripts/count_germline_sbs52.py`
- `scripts/build_germline_sbs52_hdp_matrix.py`
- `scripts/hdp_germline_sbs52.R`
- `scripts/extract_germline_sbs52.R`

Typical order:

1. Count per-sample SBS52 from a multi-sample VCF:

```bash
python ./scripts/count_germline_sbs52.py \
  -i cohort.vcf.gz \
  --ref-fasta ref.fna \
  --sample SAMPLE_ID \
  -o SAMPLE_ID.sbs52_counts.tsv
```

2. Build the HDP matrix:

```bash
python ./scripts/build_germline_sbs52_hdp_matrix.py \
  --input germ.sbs52.files.txt \
  --samples germ.samples.txt \
  --output test.germ.hdp.csv
```

3. Run HDP:

```bash
bash ./run_germline_sbs52_hdp.sh
```

Helper wrappers:

- `submit_germline_sbs52_array.sh`
- `build_germline_sbs52_hdp_matrix.sh`
- `run_germline_sbs52_hdp.sh`

Notes:

- The SBS52 matrix is a quoted CSV with `sample` as the first column.
- The R scripts validate and reorder the 52 channels by name.

### SBS96

There is no maintained germline SBS96 HDP workflow here. Germline SBS96 generation is currently used mainly for exploratory counting and comparison.

Useful wrappers:

- `run_maf_folded_germ_sbs96.sh`
- `get_maf_folded_germ_sbs96.sh`

Main counter:

- `scripts/count_somatic_sbs96.py`

Note:

- The SBS96 counter was patched to uppercase reference context, which matters when the FASTA is lowercase-masked.

## Somatic

### SBS96

Main local scripts:

- `scripts/count_somatic_sbs96.py`
- `scripts/build_somatic_sbs96_hdp_matrix.py`
- `run_somatic_sbs96_hdp.sh`

Reference HDP scripts currently used by the wrapper:

- `reference/treeoflife/code/hdp_noprior_SBS96.R`
- `reference/treeoflife/code/hdp_extraction_SBS96.R`

Typical order:

1. Count SBS96 from a somatic VCF:

```bash
python ./scripts/count_somatic_sbs96.py \
  -i sample.vcf.gz \
  --ref-fasta ref.fna \
  -o sample.sbs96.tsv
```

2. If needed, append one sample to the large ToL somatic count matrix and write an HDP-ready table:

```bash
python ./scripts/build_somatic_sbs96_hdp_matrix.py \
  --reference-csv tab10_raw_somatic_mutcount_tol.csv \
  --single-sample-tsv sample.sbs96.tsv \
  --sample-name my_sample \
  -o somatic.plus_sample.hdp.tsv
```

3. Run HDP:

```bash
bash ./run_somatic_sbs96_hdp.sh
```

Notes:

- The SBS96 HDP matrix is tab-delimited with `Sample` as the first column.
- The large ToL reference count table contains comma thousands separators; the matrix builder strips them.
- The reference SBS96 R scripts still need local patching if their input parsing changes.

### SBS52

Somatic SBS52 is currently treated as derived output from SBS96 rather than a primary maintained workflow.

Useful helper scripts:

- `get_somatic_pool_sbs52.sh`
- `reference/treeoflife/scripts/get_sbs96_mapped_to_sbs52.py`

## Practical checks

- Per-sample count files should not all be identical.
- HDP matrices should parse cleanly with no `NA` values after numeric conversion.
- Test on a small subset before running the full cohort.
