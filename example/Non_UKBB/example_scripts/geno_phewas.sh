#!/bin/bash

python ../../../pype/run_phewas.py \
  --aggregate_phenotype_file ../example_data/phenotype_files/example.tsv \
  --bfiles_directory ../example_data/genotype_files \
  --variant_file ../example_data/variant_files/example_pheWAS_SNPS.tsv \
  --covariates_file ../example_data/example_covariates.tsv \
  --dependent_phenotypes_file ../example_data/phenotype_files/example_dependent_phenotypes.txt \
  --sample_file ../example_data/example_excluded_sample_IDs.txt \
  --output_prefix example \
  --directory_name ../example_output/genotype_pheWAS \
  --threads 4 \
  --memory 12000 \
  --plink2_installation 'local' \
  --save_phenotypes \
  --pickle_intermediates \
  --correction 'bonferroni' \
  --sample_threshold 100 \
  --gene_file ../../../pype/cached/hg19_gene_list \
  --annotate