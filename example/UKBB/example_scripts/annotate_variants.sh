#!/bin/bash

python ../../../pype/annotate_variants.py \
	--significant_phewas_results ../example_output/genotype_pheWAS_$1/pheWAS/bonferroni_significant_aggregate_pheWAS_results.tsv \
	--out_dir ../example_output/gene_annotations_$1