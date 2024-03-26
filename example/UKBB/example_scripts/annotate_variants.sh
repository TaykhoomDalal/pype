#!/bin/bash

# the $1 should specify which directory to look for to annotate its significant variants
# with nearby genes and functional information (1006, 1019, etc)

python ../../../pype/annotate_variants.py \
	--significant_phewas_results ../example_output/genotype_pheWAS_$1/pheWAS/bonferroni_significant_aggregate_pheWAS_results.tsv \
	--out_dir ../example_output/gene_annotations_$1