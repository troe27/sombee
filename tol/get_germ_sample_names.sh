#!/bin/bash
module load BCFtools
bcftools query -l ../../data/BetterB_Batch1_2_3_full_geno_HQ.vcf.gz > germ.sample_names.txt
