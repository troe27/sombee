#!/bin/bash

mamba run -n tol3.9  python ../../treeoflife/scripts/get_sbs96_mapped_to_sbs52.py \
    -i /home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic.pool.batch3_batch4.sbs96.tsv  \
    --sbs96-to-sbs52 ../../treeoflife/scripts/sbs96_to_sbs52_lookup_table.tsv \
    -o /home/tilman/bees1/private/tilman/nanoseq_batch4/data/somatic.pool.batch3_batch4.sbs52.tsv
