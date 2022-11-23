#!/bin/bash

python3 ~/pype/pype/run_phewas.py \
  --ukbiobank_url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
  --ukbiobank_phenotype_file ~/pype/test_data/input_ex/test.tab \
  --bfiles_directory ~/pype/test_data/input_ex/test_bfiles \
  --variant_file ~/pype/test_data/input_ex/test_Top_SNPS.tsv \
  --sample_file ~/pype/test_data/input_ex/test_in_sample_IDs.txt \
  --output_prefix test_results \
  --threads 1 \
  --memory 8000 \
  --directory_name ~/pype/test_data/geno_phewas_$1 \
  --save_phenotypes \
  --pickle_intermediates \
  --plink2_installation 'module' \
  --plink2_module_load 'module load plink2' \
  --correction 'bonferroni' \
  --sample_threshold 100 \
  --covariates ~/pype/test_data/input_ex/test_covariates.tsv \
  --gene_file ~/pype/pype/cached/hg19_gene_list \
  --annotate \
  # --reuse_genos ~/pype/test_data/geno_phewas_$1/geno \
  # --reuse_raw ~/pype/test_data/geno_phewas_$1/raw \
  # --reuse_phenotypes ~/pype/test_data/geno_phewas_$1/pheno \
  # --old_data_fields_dir ~/pype/test_data/geno_phewas_$1/data_fields \