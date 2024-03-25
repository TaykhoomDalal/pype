#!/bin/bash

python ../../../pype/annotate_ukbb_categories.py \
    --input_file ../example_output/genotype_pheWAS_1006/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1006 \
    --singulars False \
    --input_file ../example_output/genotype_phewas_1019/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1019 \
    --singulars False \
    --input_file ../example_output/genotype_phewas_17518/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 17518 \
    --singulars True \
    --output_file ../example_output/genotype_pheWAS_categories_aggregated.tsv

python ../../../pype/annotate_ukbb_categories.py \
    --input_file ../example_output/phenotype_pheWAS_1006/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1006 \
    --singulars False \
    --input_file ../example_output/phenotype_phewas_1019/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1019 \
    --singulars False \
    --input_file ../example_output/phenotype_phewas_17518/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 17518 \
    --singulars True \
    --output_file ../example_output/phenotype_pheWAS_categories_aggregated.tsv