#!/bin/bash

python3 ~/pype/pype/run_mr.py \
    --exposure_variants ~/pype/test_data/input_ex/mr_ready_SNPs.tsv \
    --exposure_phenotypes test \
    --traits bmi \
    --output ~//pype/test_data/mr_results.tsv \
    --all_studies_path ~/pype/pype/cached/all_studies.pkl \
    --mr_type ivw \
    --mr_type egger \
    --mr_type simple_median 