#!/bin/bash

python3 ~/pype/pype/full_pheWAS_pipeline.py \
  --ukbiobank_url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
  --ukbiobank_phenotype_file ~/pype/test_data/input_ex/test.tab \
  --bfiles_directory ~/pype/test_data/input_ex/test_bfiles \
  --phenotypes_list 90012 21001 \
  --sample_file ~/pype/test_data/input_ex/test_in_sample_IDs.txt \
  --output_prefix testResults \
  --threads 1 \
  --memory 8000 \
  --directory_name ~/pype/test_data/pheno_phewas_$1 \
  --save_phenotypes \
  --pickle_intermediates \
  --correction 'bonferroni' \
  --sample_threshold 100 \
  --covariates ~/pype/test_data/input_ex/test_covariates.tsv \
  # --reuse_phenotypes ~/pype/test_data/[dir]/pheno \
  # --old_data_fields_dir ~/pype/test_data/[dir]/data_fields \