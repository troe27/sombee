# sombee

Analysis code for measuring somatic mutation rates in honeybees (*Apis mellifera*) with NanoSeq duplex sequencing.


1. **`dupcaller_workflow/`**:

- fastq2vcf using DupCaller including slurm setup
- Importable functions for handling and formatting the output of DupCaller in dupIO.py

2. **`downstream_processing/`**: filters for somatic calls - output of the NanoSeq pipeline: germline overlap, discarded sites,  double variants, variants within repeats. Import as `downstream_processing`,
   
3. **`correction/`**: False-positive correction.
    - Estimate completeness of reference germline VCF using variant discovery saturation
    - Combine estimated completeness of reference dataset with known false positives to get adjustment factor for mutation rate.
   -  builds the annotated variant tables from NanoSeq output.
   -  positive-unlabelled classifier (Elkan–Noto)/Random Forest ensemble to rate somatic variants on likelihood to be false positive - use with threshold from variant saturation discovery as cutoff
     
4. **`tol/`**: mutational signatures (SBS52 germline, SBS96 somatic) and HDP extraction using the treeoflife code. See `tol/docs/readme.md`.

## Setup

See `install_dependencies.txt` (mamba env `sombee`).

Data is not in the repo. `data/` is a gitignored symlink to the data folder on each machine:

    ln -s /path/to/sombee_data data

Notebooks read and write via `data/` or `../data/`.
