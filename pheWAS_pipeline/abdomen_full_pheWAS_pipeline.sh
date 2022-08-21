#!/bin/bash
#SBATCH -c 4                                                                                       # Request four cores
#SBATCH -t 1-1:00                                                                               # Runtime in D-HH:MM format
#SBATCH -p medium                                                                                    # Partition to run in
#SBATCH --mem-per-cpu=64000                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/abdomen_pheWAS_pipeline%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/abdomen_pheWAS_pipeline%j.err   # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                                                                           # Type of email notification- BEGIN,END,FAIL,ALL

python3 ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/full_pheWAS_pipeline.py \
                                --phenotype /home/tad368/orig_phenos/ukb41230.tab \
                                --bfiles ~/data_dir/bfiles \
                                --covariates /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/cov_fields.txt \
                                --url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Abdomen_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Abdomen_Top_SNPS.tsv \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Pancreas_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Pancreas_Top_SNPS.tsv \
                                --samples ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Liver_in_sample_pop_IDs.txt \
                                --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Liver_Top_SNPS.tsv \
                                --threads 4 \
                                --memory 64000 \
                                --save_phenotypes \
                                --location ~/data_dir/AC \
                                --name ukbb_category_$1_abdomen \
                                --pickle_intermediates \
                                --reuse_genos /home/tad368/data_dir/ABDOMEN_COMMON_FILES/geno \
                                --reuse_raw /home/tad368/data_dir/ABDOMEN_COMMON_FILES/raw \
                                # --reuse_phenotypes /home/tad368/data_dir/ABDOMEN_COMMON_FILES/AC/$1/pheno \
                                # --old_data_fields_dir /home/tad368/data_dir/ABDOMEN_COMMON_FILES/AC/$1/data_fields \