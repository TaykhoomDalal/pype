#!/bin/bash

python ../../../pype/run_mr.py \
    --exposure_variants ../example_data/variant_files/example_mr_SNPs_1.tsv \
    --exposure_phenotypes random \
    --traits bmi \
    --output ../example_output/mr_results.tsv \
    --all_studies_path ../../../pype/cached/all_studies.pkl \
    --run_all_mr


python ../../../pype/run_mr.py \
    --exposure_variants ../example_data/variant_files/example_mr_SNPs_1.tsv \
    --exposure_phenotypes example1 \
    --outcome_variants ../example_data/variant_files/example_mr_SNPs_2.tsv \
    --outcome_phenotype example2 \
    --output ../example_output/example_mr_results.tsv \
    --run_all_mr