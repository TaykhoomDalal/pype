import pandas as pd
import argparse
import pheWAS
import os

# major ethnic groups in UKBB
ethnicity_groups = ['Ethnicity.White']#, 'Ethnicity.Asian', 'Ethnicity.Black', 'Ethnicity.Other']

def main():
	# Parse the command line arguments
	parser = argparse.ArgumentParser(description='Pass the input html file')
	parser.add_argument('-p', '--pheno_files', required=True, help='Phenotype files', action = 'append')
	parser.add_argument('-o', '--out_dir', required=True, help = 'PheWAS Results directory')
	parser.add_argument('-r', '--raw_dir', required=True, help = 'Raw data directory')
	parser.add_argument('-d', '--data_fields_dir', required=True, help = 'Data fields directory')
	parser.add_argument('-l', '--l1_regularization', help = 'Regularization', action = 'store_true')
	args = parser.parse_args()

	# parse arguments
	pheno_files = args.pheno_files
	out_dir = args.out_dir
	raw_dir = args.raw_dir
	data_fields_dir = args.data_fields_dir
	reg = args.l1_regularization

	for file in pheno_files:

		file_prefix = os.path.basename(file).split('_')[0]

		pheno = pd.read_csv(file, sep='\t', low_memory = False, header = 0)

		# find the raw genotype files with the same prefix
		raw_geno_files = [os.path.join(raw_dir, f) for f in os.listdir(raw_dir) if f.endswith('.raw') and os.path.basename(f).startswith(file_prefix)]
		
		# read covariate fixed headers file into list
		with open(data_fields_dir + '/' +  file_prefix + '_fixed_cov_fields.txt', "r") as f:
				cov_headers = f.read().splitlines()

		for geno_file in raw_geno_files:
			
			for ethnicity_group in ethnicity_groups:

				# get the prefix of the geno_file name (before the .raw)
				geno_file_prefix = geno_file.split('/')[-1].split('.')[0]
				
				results_file = out_dir + '/' + geno_file_prefix + '_'+ ethnicity_group.replace('.', '_') + '_results.tsv'

				print('Running PheWAS for ' + results_file)

				# set index column to the IID and then rearrange the rows according to the order of the phenotype file by the f.eid column
				geno_data = pd.read_csv(geno_file, sep='\t', index_col = 'IID').reindex(pheno.loc[:, 'f.eid'].tolist())

				# reset the index first to get rid of IID column, and then since after the 'PHENOTYPE' column, the variant data starts, select everything after, and drop the phenotype column
				geno_data = geno_data.reset_index().loc[:, 'PHENOTYPE':].drop('PHENOTYPE', axis = 1)

				# get list of columns starting with 'Ethnicity'
				cols_to_remove = [col for col in cov_headers if col.startswith('Ethnicity')]

				# remove all the columns starting with ethnicity since we will only be looking at once ethnicity (ethnicity covariate is uneeded)
				cov_headers_ethnicity_fixed = [cov for cov in cov_headers if cov not in cols_to_remove]

				# only keep individuals from the chosen population
				pheno_fixed = pheno[pheno[ethnicity_group] == 1]
				
				# drop the f.eid column and all the ethnicities in the phenotype data since we no longer need them
				pheno_fixed = pheno_fixed.drop(cols_to_remove + ['f.eid'], axis = 1)

				# Take the difference between the covariate headers and the phenotype headers
				diff_cov_pheno = list(set(pheno_fixed.columns) - set(cov_headers_ethnicity_fixed))

				results = pheWAS.run_phewas(phenotypes = pheno_fixed, 
											genotypes = geno_data, 
											non_cov_pheno_list= diff_cov_pheno, 
											reg = reg, 
											covariates = cov_headers_ethnicity_fixed)

				# save the results to a file
				results.to_csv(results_file, sep='\t', index=False)



if __name__ == "__main__":
	main()