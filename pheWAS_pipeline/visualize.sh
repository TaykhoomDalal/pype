#!/bin/bash
#SBATCH -c 1                                                                                        # Request four cores
#SBATCH -t 0-00:50                                                                               # Runtime in D-HH:MM format
#SBATCH -p short                                                                                    # Partition to run in
#SBATCH --mem=8GB                                                                                # Memory total in MiB (for all cores)
#SBATCH -o /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/visualize_%j.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/batch_logs/visualize_%j.err   # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                                                                           # Type of email notification- BEGIN,END,FAIL,ALL

if [ $1 == "-abdomen" ]; then
	python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/visualize.py \
		--input ~/data_dir/AC_included/PheWAS_Visualization/ukbb_categories_aggregated_abdomen.tab \
		--output ~/data_dir/AC_included/PheWAS_Visualization/abdomen/$2/abdomen_all.png \
		--mapping category_mappings.txt \
		--gene_file ~/data_dir/AC_included/PheWAS_Visualization/hg19_gene_list \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Abdomen_Top_SNPS.tsv \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Liver_Top_SNPS.tsv \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Abdomen/Pancreas_Top_SNPS.tsv \
		--plot_top_categories \
		--number_of_top_results 20 \
		$3 \
		$4
elif [ $1 == '-eyes' ]; then
	python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/visualize.py \
		--input ~/data_dir/AC_included/PheWAS_Visualization/ukbb_categories_aggregated_eyes.tab \
		--output ~/data_dir/AC_included/PheWAS_Visualization/eyes/$2/eyes_all.png \
		--mapping category_mappings.txt \
		--gene_file ~/data_dir/AC_included/PheWAS_Visualization/hg19_gene_list \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Eyes/EyesFundus_Top_SNPS.tsv \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Eyes/EyesOCT_Top_SNPS.tsv \
		--plot_top_categories \
		--number_of_top_results 20 \
		$3 \
		$4
elif [ $1 == '-heart' ]; then
	python3 /home/tad368/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/pheWAS_pipeline/visualize.py \
		--input ~/data_dir/AC_included/PheWAS_Visualization/ukbb_categories_aggregated_heart.tab \
		--output ~/data_dir/AC_included/PheWAS_Visualization/heart/$2/heart_all.png \
		--mapping category_mappings.txt \
		--gene_file ~/data_dir/AC_included/PheWAS_Visualization/hg19_gene_list \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Heart/Heart_Top_SNPS.tsv \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Heart/HeartECG_Top_SNPS.tsv \
		--rsid_files ~/PheWAS-and-GWAS-Diabetes-Correlation-pipeline/variants_list/Heart/HeartMRI_Top_SNPS.tsv \
		--plot_top_categories \
		--number_of_top_results 20 \
		$3 \
		$4
fi
