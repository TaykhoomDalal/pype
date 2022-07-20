#!/bin/bash
#SBATCH -c 1                                                                                        # Request one core
#SBATCH -t 0-00:10:00                                                                                # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=16000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/annotate_results_%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/annotate_results_%j.err   # File to which STDERR will be written, including job ID (%j)

python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/annotate_results.py \
    --input_file /home/tad368/data_dir/ukbb_category_1003_$1/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1003 \
    --singulars False \
    --input_file /home/tad368/data_dir/ukbb_category_1006_$1/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1006 \
    --singulars False \
    --input_file /home/tad368/data_dir/ukbb_category_1019_$1/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1019 \
    --singulars False \
    --input_file /home/tad368/data_dir/ukbb_category_1307_$1/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 1307 \
    --singulars True \
    --input_file /home/tad368/data_dir/ukbb_category_17518_$1/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 17518 \
    --singulars True \
    --input_file /home/tad368/data_dir/ukbb_category_100081_$1/pheWAS/aggregate_pheWAS_results.tsv \
    --categories 100081 \
    --singulars True \
    --output_file /home/tad368/data_dir/PheWAS_Visualization/ukbb_categories_aggregated_$1.tsv
