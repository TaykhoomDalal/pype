#!/bin/bash
#SBATCH -c 1                                                                                        # Request one core
#SBATCH -t 0-00:30:00                                                                                # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=64000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/save_top_variants_%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/save_top_variants_%j.err   # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                                                                           # Type of email notification- BEGIN,END,FAIL,ALL

python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/7_save_top_variants.py \
        --url https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=$1 \
        --out_dir /home/tad368/pheWAS