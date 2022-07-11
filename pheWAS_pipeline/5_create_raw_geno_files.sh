#!/bin/bash
#SBATCH -c 4                                                                                        # Request one core
#SBATCH -t 0-3:00:00                                                                                # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=64000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/create_raw_geno_files_%j.out  # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/create_raw_geno_files_%j.err  # File to which STDERR will be written, including job ID (%j)

python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/5_create_raw_geno_files.py \
        --geno_dir /home/tad368/data_dir/geno \
        --raw_dir /home/tad368/data_dir/raw \
