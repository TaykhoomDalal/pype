#!/bin/bash

python3 ~/pype/pype/annotate_categories.py \
    --input_file ~/pype/test_data/geno_phewas_1006/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1006 \
    --singulars False \
    --input_file ~/pype/test_data/geno_phewas_1019/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1019 \
    --singulars False \
    --output_file ~/pype/test_data/input_ex/ukbb_categories_aggregated.tsv