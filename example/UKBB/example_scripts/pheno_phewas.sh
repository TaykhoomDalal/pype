#!/bin/bash

# the $1 should be the category of phenotype you want to run the pheWAS on (1006, 1019, etc)

python ../../../pype/run_phewas.py \
  --UKBB \
  --aggregate_phenotype_file ../example_data/phenotype_files/example.tab \
  --independent_phenotypes_list 90012 21001 \
  --covariates_file ../example_data/example_covariates.tsv \
  --sample_file ../example_data/example_excluded_sample_IDs.txt \
  --output_prefix example \
  --ukbiobank_url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
  --threads 4 \
  --memory 12000 \
  --directory_name ../example_output/phenotype_pheWAS_$1 \
  --save_phenotypes \
  --pickle_intermediates \
  --correction 'bonferroni' \
  --sample_threshold 100