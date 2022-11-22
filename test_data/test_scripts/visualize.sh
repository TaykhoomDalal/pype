#!/bin/bash

python3 ~/pype/pype/visualize.py \
    --phewas_results ~/pype/test_data/input_ex/ukbb_categories_aggregated.tsv \
    --directory_name ~/pype/test_data/test_images \
    --output_prefix 'testVis' \
    --output_extension '.png' \
    --phewas_type 'genotype' \
    --mapping ~/pype/pype/category_mappings.txt \
    --clear_old_files \
    --single_category \
    --variant_files ~/pype/test_data/input_ex/test_Top_SNPS.tsv \
    --seed 0 \
    --plot_manhattan \
    --plot_top_data_field_categories \
    --plot_volcano \
    --annotate_top_N_volcano 10 \
    --compare_original_betas \
    --plot_category_enrichment \
    --annotate