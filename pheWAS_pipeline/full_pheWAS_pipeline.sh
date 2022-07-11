#!/bin/bash
#SBATCH -c 8                                                                                        # Request one core
#SBATCH -t 0-4:00:00                                                                               # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=128000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/abdomen_pheWAS_pipeline%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/abdomen_pheWAS_pipeline%j.err   # File to which STDERR will be written, including job ID (%j)

python3 ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/full_pheWAS_pipeline.py \
                                --phenotype ~/data_dir/pheno/pheno.tab \
                                --bfiles ~/data_dir/bfiles \
                                --covariates ~/data_dir/pheno/covariates_data_field_list.txt \
                                --url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Abdomen_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Abdomen_Top_SNPS.tsv \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Pancreas_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Pancreas_Top_SNPS.tsv \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Liver_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Liver_Top_SNPS.tsv \
                                --threads 16 \
                                --memory 128000 \
                                --location ~/data_dir \
                                --name ukbb_category_$1