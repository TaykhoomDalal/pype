#!/bin/bash
#SBATCH -c 4                                                                                        # Request one core
#SBATCH -t 0-00:30:00                                                                                # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=64000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/extract_%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/extract_%j.err   # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                                                                           # Type of email notification- BEGIN,END,FAIL,ALL

python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/2_extract_variants_and_samples.py \
        --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Abdomen_Top_SNPS.tsv \
        --samples /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Abdomen_in_sample_pop_IDs.txt \
        --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Liver_Top_SNPS.tsv \
        --samples /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Liver_in_sample_pop_IDs.txt \
        --variants /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Pancreas_Top_SNPS.tsv \
        --samples /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Pancreas_in_sample_pop_IDs.txt \
        --bfiles /home/tad368/data_dir/bfiles \
        --geno /home/tad368/data_dir/geno \
        --threads 4 \
        --memory 64000