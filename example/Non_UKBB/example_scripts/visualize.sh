#!/bin/bash

python ../../../pype/visualize.py \
    --phewas_results ../example_output/genotype_pheWAS_categories_aggregated.tsv \
    --directory_name ../example_output/genotype_pheWAS_plots \
    --output_prefix 'example' \
    --output_format 'png' \
    --dpi 300 \
    --height 6 \
    --width 10 \
    --phewas_type 'genotype' \
    --clear_old_files \
    --single_category \
    --variant_files ../example_data/variant_files/example_pheWAS_SNPS.tsv \
    --seed 0 \
    --plot_manhattan \
    --plot_top_data_field_categories \
    --plot_volcano \
    --annotate_top_N_volcano 10 \
    --compare_original_betas \
    --plot_category_enrichment \
    --annotate

python ../../../pype/visualize.py \
    --phewas_results ../example_output/phenotype_pheWAS_categories_aggregated.tsv \
    --directory_name ../example_output/phenotype_pheWAS_plots \
    --output_prefix 'example' \
    --output_format 'png' \
    --dpi 300 \
    --height 6 \
    --width 10 \
    --phewas_type 'phenotype' \
    --clear_old_files \
    --single_category \
    --seed 0 \
    --plot_manhattan \
    --plot_top_data_field_categories \
    --plot_volcano \
    --annotate_top_N_volcano 10 \
    --plot_category_enrichment \
    --annotate