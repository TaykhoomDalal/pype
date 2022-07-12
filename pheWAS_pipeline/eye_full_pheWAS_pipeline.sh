#!/bin/bash
#SBATCH -c 4                                                                                        # Request four cores
#SBATCH -t 0-36:00:00                                                                               # Runtime in D-HH:MM format
#SBATCH -p medium                                                                                    # Partition to run in
#SBATCH --mem=64000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/eyes_pheWAS_pipeline%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/eyes_pheWAS_pipeline%j.err   # File to which STDERR will be written, including job ID (%j)

python3 ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/full_pheWAS_pipeline.py \
                                --phenotype /home/tad368/orig_phenos/ukb41230.tab \
                                --bfiles ~/data_dir/bfiles \
                                --covariates ~/data_dir/data_fields/cov_fields.txt \
                                --url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Eyes_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Eyes/Eyes_Top_SNPS.tsv \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/EyesFundus_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Eyes/EyesFundus_Top_SNPS.tsv \
                                --threads 4 \
                                --memory 64000 \
                                --location ~/data_dir \
                                --name ukbb_category_$1_eyes