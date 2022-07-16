#!/bin/bash
#SBATCH -c 4                                                                                        # Request one core
#SBATCH -t 0-3:00:00                                                                                # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=64000M                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/make_pheno_%j.out  # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/make_pheno_%j.err  # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                                                                           # Type of email notification- BEGIN,END,FAIL,ALL

python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/3_make_out_of_sample_pheno_file.py \
        --sample_files /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/HeartMRI_in_sample_pop_IDs.txt \
        --pheno_file /home/tad368/orig_phenos/ukb41230.tab \
        --fields /home/tad368/data_dir/data_fields/data_fields.txt \
        --covariates_file /home/tad368/data_dir/data_fields/cov_fields.txt \
        --out_dir /home/tad368/data_dir/pheno \
        --data_fields_dir /home/tad368/data_dir/data_fields
        # --samples /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Liver_in_sample_pop_IDs.txt \
        # --samples /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/in_pop_IDs/Pancreas_in_sample_pop_IDs.txt \
# 31 is sex
# 34 is year of birth
# 52 is month of birth 
# 53 is date of attending recruitment center
# 21000 is ethnic background
# 22001 is genetic sex
# 22009 is principal component 1-40