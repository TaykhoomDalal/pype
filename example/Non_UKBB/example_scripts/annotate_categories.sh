#!/bin/bash

python ../../../pype/annotate_external_categories.py \
    --input_file ../example_output/genotype_pheWAS/pheWAS/aggregate_pheWAS_results.tsv \
    --mapping_file ../example_data/phenotype_files/example_mapper.tsv \
    --output_file ../example_output/genotype_pheWAS_categories_aggregated.tsv

python ../../../pype/annotate_external_categories.py \
    --input_file ../example_output/phenotype_pheWAS/pheWAS/aggregate_pheWAS_results.tsv \
    --mapping_file ../example_data/phenotype_files/example_mapper.tsv \
    --output_file ../example_output/phenotype_pheWAS_categories_aggregated.tsv