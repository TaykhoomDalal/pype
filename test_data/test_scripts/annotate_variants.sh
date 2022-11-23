#!/bin/bash

python3 ~/pype/pype/annotate_variants.py \
	--significant_phewas_results ~/pype/test_data/geno_phewas_1006/pheWAS/TOP_SNPS_bonferroni__aggregate_pheWAS_results.tsv \
	--out_dir ~/pype/test_data/test_annotations