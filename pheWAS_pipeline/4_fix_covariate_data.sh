#!/bin/bash
#SBATCH -c 4                                                                                        # Request one core
#SBATCH -t 0-3:00:00                                                                                # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=64000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/fix_covariate_data_%j.out  # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/fix_covariate_data_%j.err  # File to which STDERR will be written, including job ID (%j)

python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/4_fix_covariate_data.py \
        --pheno_files /home/tad368/data_dir/pheno/Abdomen_subsetted_pheno.tab \
        --fields_file /home/tad368/data_dir/data_fields/ukb_data_fields.txt \
        --out_dir /home/tad368/data_dir/pheno \
        --data_fields_dir /home/tad368/data_dir/data_fields
