import mr
import pandas as pd
import mr_utils
import argparse


def main():
	parser = argparse.ArgumentParser(description="Run Mendelian Randomization analyses with the variants used in the PheWAS")

	# Required Arguments
	parser.add_argument('--exposure_variants', help = 'File containing variants from PheWAS', required = True, type = str, action = 'append')
	parser.add_argument('--exposure_phenotypes', help = 'Name of the phenotype(s) that the exposure variants are associated with (must be in the same order as the files)', required = True, type = str, action = 'append')
	parser.add_argument('--outcome_variants', help = 'File containing variants to run MR against using exposures', required = False, default = None, type = str)
	parser.add_argument('--outcome_phenotype', help = 'Name of the phenotype that the outcome variants are associated with', required = False, default = None, type = str)
	parser.add_argument('--traits', help = 'List of traits to search for in OpenGWAS to use for MR', required = False, type = str, action ='append', nargs = '+')
	parser.add_argument('--output', help ='File to store output of MR', required = True, type = str)
	parser.add_argument('--mr_type', help = 'MR test to run', choices = ['ivw', 'egger', 'simple_median'], required = True, type = str, action='append')

	# Optional Arguments for MR
	parser.add_argument('--similarity_func', help = 'Similarity function to use for matching trait strings in OpenGWAS', choices = ['jaccard', 'levenshtein'], default = 'jaccard', required = False, type = str)
	parser.add_argument('--all_studies_path', help = 'Path to file containing all studies in OpenGWAS', default = None, required = False, type = str)
	parser.add_argument('--cache_all_studies', help = 'Cache all studies in OpenGWAS', required = False, action = 'store_true')
	parser.add_argument('--batches', help = 'GWAS summary dataset batches from Open GWAS', choices=mr_utils.open_gwas_batches.keys(), required = False, type = str, nargs = '+')
	parser.add_argument('--strip_names', help = 'If trying to match traits is not working well, try standardizing to strip all non alpha characters from the traits', required = False, action = 'store_true')

	args = parser.parse_args()

	# Required Arguments
	exposure_variants = args.exposure_variants
	exposure_phenotypes = args.exposure_phenotypes
	outcome_variants = args.outcome_variants
	outcome_phenotype = args.outcome_phenotype
	traits = args.traits
	output = args.output
	mr_type = args.mr_type

	# Optional Arguments
	similarity_func = args.similarity_func
	all_studies_path = args.all_studies_path
	cache_all_studies = args.cache_all_studies
	batches = args.batches
	strip_names = args.strip_names

	# ---------------------------------------VERIFY ARGS--------------------------------------- #
	
	if traits is None and outcome_variants is None:
		print("Error: You must provide either a list of traits or a file containing outcome variants to run the MR")
		exit()
	
	if traits is not None and outcome_variants is not None:
		print("Error: You must provide either a list of traits or a file containing outcome variants to run the MR, but not both")
		exit()
	
	if outcome_variants is not None and outcome_phenotype is None:
		print("Error: You must provide the phenotype name for the outcome variants")
		exit()
	
	if outcome_variants is None and outcome_phenotype is not None:
		print("Warning: You provided an outcome phenotype, but no outcome variants. This will not be used")

	if len(exposure_variants) != len(exposure_phenotypes):
		print("Error: The number of exposure variants files must match the number of exposure phenotypes")
		exit()

	for exp_pheno in exposure_phenotypes:
		if exp_pheno == outcome_phenotype:
			print("Error: The exposure and outcome phenotypes cannot be the same")
			exit()

	if len(mr_type) == 0:
		print("Error: You must provide at least one MR test to run")
		exit()

	# ----------------------------------------------------------------------------------------- #
	
	# --------------------------------------SETUP MR CODE-------------------------------------- #
	
	if all_studies_path is not None:
		mr_utils.setAllStudiesPath(all_studies_path)

	# if the user wants to re-cache all studies
	if cache_all_studies:
		mr_utils.cache_all_studies()
	
	# if the user wants to change the similarity function
	similarity_func = mr_utils.get_similarity_func(similarity_func)

	# get the output file base name and extension
	output_base, output_ext = output.rsplit('.', 1)

	# ----------------------------------------------------------------------------------------- #

	# --------------------------------------RUN MR CODE---------------------------------------- #

	# if the user wants to run MR against a list of traits
	if traits is not None:
		
		for index, exposure_pheno in enumerate(exposure_phenotypes):

			# read in the exposure variants
			exposure_variants_i = pd.read_csv(exposure_variants[index], sep = '\t')

			temp_data = []
			for trait_family in traits:
				for trait in trait_family:

					traits_to_study_mapping = mr_utils.find_studies_based_on_traits([[trait]], similarity = similarity_func, batch_list = batches, strip = strip_names)

					external_gwas_data = mr_utils.extract_snps_from_outcomes(exposure_variants_i['rsID'].to_list(), list(traits_to_study_mapping.values())[0], proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
					
					data = {}

					for phenotype in external_gwas_data['phenotype'].unique():
						external_gwas_data_i = external_gwas_data.loc[external_gwas_data['phenotype'] == phenotype].copy()
						data[phenotype] = mr_utils.harmonize(exposure_variants_i, external_gwas_data_i, '_' + exposure_pheno, '_' + phenotype)
						
						print('Running MR ({}) for Exposure ({}) against Outcome ({})'.format(", ".join(mr_type), exposure_pheno, phenotype))
						
						mr_res = mr.run_mr(mr_type, data[phenotype], 'BETA_' + exposure_pheno, 'BETA_' + phenotype, 'SE_' + exposure_pheno, 'SE_' + phenotype)

						for mr_type_i, results in mr_res.items():
							if results is None:
								print('{} MR failed for Exposure ({}) against Outcome ({})'.format(mr_type_i.capitalize(), exposure_pheno, phenotype))
							else:
								pval, beta, std_err, *additional = results

								temp_data.append([mr_type_i, phenotype, pval, beta, std_err, data[phenotype].shape[0]])
			
			mr_results = pd.DataFrame(temp_data, columns = ['MR_Method', 'Outcome', 'P_value', 'Effect_Size', 'Standard_Error', 'Number_SNPs'])

			# write the results to a file
			mr_results.to_csv(output_base + '_' + exposure_pheno + '.' + output_ext, sep = '\t', index = False)

	elif outcome_variants != None:
		
		for index, exposure_pheno in enumerate(exposure_phenotypes):
			
			print('Running MR ({}) for Exposure ({}) against Outcome ({})'.format(", ".join(mr_type), exposure_pheno, outcome_phenotype))

			# read in the exposure variants
			exposure_variants_i = pd.read_csv(exposure_variants[index], sep = '\t')

			# read in the outcome variants
			outcome_variants = pd.read_csv(outcome_variants, sep = '\t')

			# harmonize the exposure and outcome variants
			harmonized_data = mr_utils.harmonize(exposure_variants_i, outcome_variants, '_' + exposure_pheno, '_' + outcome_phenotype)

			# run MR
			mr_res = mr.run_mr(mr_type, harmonized_data, 'BETA_' + exposure_pheno, 'BETA_' + outcome_phenotype, 'SE_' + exposure_pheno, 'SE_' + outcome_phenotype)

			temp_data = []
			for mr_type_i, results in mr_res.items():
				if results is None:
					print('{} MR failed for Exposure ({}) against Outcome ({})'.format(mr_type_i.capitalize(), exposure_pheno, outcome_phenotype))
				else:
					pval, beta, std_err, *additional = results

					temp_data.append([mr_type_i, outcome_phenotype, pval, beta, std_err, harmonized_data.shape[0]])

			# convert the results to a dataframe
			mr_results = pd.DataFrame(data = temp_data, columns = ['Outcome', 'P_value', 'Effect_Size', 'Standard_Error', 'Number_SNPs'])

			# write the results to a file
			mr_results.to_csv(output_base + '_' + exposure_pheno + '.' + output_ext, sep = '\t', index = False)

if __name__ == '__main__':
	main()