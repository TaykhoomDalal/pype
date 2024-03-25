#!/bin/bash

python ../../../pype/run_phewas.py \
  --aggregate_phenotype_file ../example_data/phenotype_files/example.tsv \
  --independent_phenotypes_list "Overall acceleration average" "Body mass index (BMI)" \
  --dependent_phenotypes_file ../example_data/phenotype_files/example_dependent_phenotypes.txt \
  --covariates_file ../example_data/example_covariates.tsv \
  --sample_file ../example_data/example_excluded_sample_IDs.txt \
  --output_prefix example \
  --ukbiobank_url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
  --threads 4 \
  --memory 12000 \
  --directory_name ../example_output/phenotype_pheWAS \
  --save_phenotypes \
  --pickle_intermediates \
  --correction 'bonferroni' \
  --sample_threshold 100