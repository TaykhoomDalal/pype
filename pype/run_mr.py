import os
import mr
import json
import inspect
import mr_utils
import argparse
import pandas as pd

expected_ieugwas_column_names = ["rsID", "CHR", "Effect", "BETA", "P", "SE", "N"]

def save_presso_extras(presso_file_path, results, presso_args):

	with open(presso_file_path, 'w') as file:

		global_test = results["Global Test"]
		rssobs, p, p_sig = global_test["RSSobs"], global_test["P_value"], global_test["P_value_significance"]
		print("MR Presso Global Test:\n\tRSSobs: {}\n\tP_value: {}\n\tP_value Significance: {}".format(rssobs, p, p_sig))

		file.writelines("MR Presso Global Test:\n\tRSSobs: {}\n\tP_value: {}\n\tP_value Significance: {}".format(rssobs, p, p_sig))

		if len(results) > 2:
			# we have the outlier test and the distortion test or just the outlier test
			outlier_test = presso_args['outlier_test']
			distortion_test = presso_args['distortion_test']

			if outlier_test:
				file.writelines("\n\n")
				outlier_df = results['Outlier Test']
				outlier_df.index.name = 'Index'
				file.write(outlier_df.to_string(header=True))

			if distortion_test:
				distortion_results = results['Distortion Test']
				indices, d_coeff, p = ",".join([str(i) for i in distortion_results["Outliers Indices"]]), distortion_results["Distortion Coefficient"], distortion_results["P_value"]
				file.writelines("\n\nMR Presso Distortion Test:\n\tOutlier Indices: {}\n\tDistortion Coefficient: {}\n\tP_value: {}".format(indices, d_coeff, p.values[0]))


def main():
	parser = argparse.ArgumentParser(description="Run Mendelian Randomization analyses with the variants used in the PheWAS")

	# Required Arguments
	parser.add_argument('--exposure_variants', help = 'File containing variants from PheWAS', required = True, type = str, action = 'append')
	parser.add_argument('--exposure_phenotypes', help = 'Name of the phenotype(s) that the exposure variants are associated with (must be in the same order as the files)', required = True, type = str, action = 'append')
	parser.add_argument('--outcome_variants', help = 'File containing variants to run MR against using exposures', required = False, default = None, type = str)
	parser.add_argument('--outcome_phenotype', help = 'Name of the phenotype that the outcome variants are associated with', required = False, default = None, type = str)
	parser.add_argument('--traits', help = 'List of traits to search for in OpenGWAS to use for MR', required = False, type = str, action ='append', nargs = '+')
	parser.add_argument('--output', help ='File to store output of MR', required = True, type = str)
	parser.add_argument('--mr_type', help = 'MR test to run', choices = ['ivw', 'egger', 'simple_median', 
																		'weighted_median', 'penalized_weighted_median',
																		'simple_mode', 'simple_mode_nome', 'weighted_mode',
																		'penalized_mode', 'weighted_mode_nome', 'presso'], required = False, type = str, action='append')
	parser.add_argument('--run_all_mr', help = 'Whether or not to run all the MR tests.', required=False, action = 'store_true')

	# Optional if headers are standard with what the script expects
	parser.add_argument('--exp_rsID', help = 'Name of the column indicating the rsID for the exposure variants', required = False, default = 'rsID', type = str)
	parser.add_argument('--out_rsID', help = 'Name of the column indicating the rsID for the outcome variants', required = False, default = 'rsID', type = str)
	parser.add_argument('--exp_CHR', help = 'Name of the column indicating the chromosome for the exposure variants', required = False, default = 'CHR', type = str)
	parser.add_argument('--out_CHR', help = 'Name of the column indicating the chromosome for the outcome variants', required = False, default = 'CHR', type = str)
	parser.add_argument('--exp_B', help = 'Name of the column indicating the beta for the exposure variants', required = False, default = 'BETA', type = str)
	parser.add_argument('--out_B', help = 'Name of the column indicating the beta for the outcome variants', required = False, default = 'BETA', type = str)
	parser.add_argument('--exp_P', help = 'Name of the column indicating the p-value for the exposure variants', required = False, default = 'P', type = str)
	parser.add_argument('--out_P', help = 'Name of the column indicating the p-value for the outcome variants', required = False, default = 'P', type = str)
	parser.add_argument('--exp_SE', help = 'Name of the column indicating the standard error for the exposure variants', required = False, default = 'SE', type = str)
	parser.add_argument('--out_SE', help = 'Name of the column indicating the standard error for the outcome variants', required = False, default = 'SE', type = str)
	parser.add_argument('--exp_N', help = 'Name of the column indicating the number of samples for the exposure variants', required = False, default = 'N', type = str)
	parser.add_argument('--out_N', help = 'Name of the column indicating the number of samples for the outcome variants', required = False, default = 'N', type = str)
	parser.add_argument('--exp_EFFECT', help = 'Name of the column indicating the effect allele for the exposure variants', required = False, default = 'Effect_Allele', type = str)
	parser.add_argument('--out_EFFECT', help = 'Name of the column indicating the effect allele for the outcome variants', required = False, default = 'Effect_Allele', type = str)

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
	run_all_mr = args.run_all_mr

	# Headers
	exp_rsID = args.exp_rsID
	out_rsID = args.out_rsID
	exp_CHR = args.exp_CHR
	out_CHR = args.out_CHR	
	exp_B = args.exp_B
	out_B = args.out_B
	exp_P = args.exp_P
	out_P = args.out_P
	exp_SE = args.exp_SE
	out_SE = args.out_SE
	exp_N = args.exp_N
	out_N = args.out_N
	exp_EFFECT = args.exp_EFFECT
	out_EFFECT = args.out_EFFECT

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

	if mr_type is None and run_all_mr is False:
		print("Error: You must provide at least one MR test to run or specify to run all MR tests.")
		exit()

	if mr_type is not None and run_all_mr is True:
		print("Warning: You specified to run all MR tests but also specified a specific MR test to run. All MR tests will be run.")
		mr_type = None

	if cache_all_studies and all_studies_path is None:
		print("Error: You must provide a path where the studies from OpenGWAS will be cached")
		exit()

	try:
		filename = inspect.getframeinfo(inspect.currentframe()).filename
		path = os.path.dirname(os.path.abspath(filename))
		# Attempt to open the JSON file (is in current directory)
		print("Attempting to open the JSON file containing argument values for MR (located at {})".format(path + '/mr_args.json'))
		with open(path + "/mr_args.json", 'r') as file:
			# Load its content as a Python dictionary
			mr_args = json.load(file)
			print("JSON data loaded successfully.")
	except FileNotFoundError:
		print("mr_args.json file not found. Please copy the file from the github")
	except json.JSONDecodeError:
		print("Error decoding JSON from mr_args.json file. Please fix the errors in the JSON file and rerun this script.")

	# ----------------------------------------------------------------------------------------- #
	
	# --------------------------------------SETUP MR CODE-------------------------------------- #
	
	if all_studies_path is not None:
		mr_utils.setAllStudiesPath(all_studies_path)

	# if the user wants to re-cache all studies
	if cache_all_studies:
		print("Caching all studies in OpenGWAS")
		mr_utils.cache_all_studies()
		print("Finished caching all studies in OpenGWAS")
	
	# if the user wants to change the similarity function
	similarity_func = mr_utils.get_similarity_func(similarity_func)

	# get the output file base name and extension
	output_base, output_ext = output.rsplit('.', 1)

	# ----------------------------------------------------------------------------------------- #

	# --------------------------------------RUN MR CODE---------------------------------------- #

	# if the user wants to run MR against a list of traits
	if traits is not None:
		
		if run_all_mr == True:
			mr_string = "All Tests"
		else:
			mr_string = ", ".join(mr_type)

		for index, exposure_pheno in enumerate(exposure_phenotypes):

			# read in the exposure variants
			exposure_variants_i = pd.read_csv(exposure_variants[index], sep = None, engine='python')

			if exposure_variants_i.empty:
				print("File {} associated with exposure phenotype {} is empty. Skipping.".format(exposure_variants[index], exposure_pheno))
				continue
			else:
				if exp_N in exposure_variants_i.columns:
					exposure_variants_i = exposure_variants_i[[exp_rsID, exp_CHR, exp_EFFECT, exp_B, exp_P, exp_SE, exp_N]]
				else:
					exposure_variants_i = exposure_variants_i[[exp_rsID, exp_CHR, exp_EFFECT, exp_B, exp_P, exp_SE]]

			temp_data = []
			for trait_family in traits:
				for trait in trait_family:

					traits_to_study_mapping = mr_utils.find_studies_based_on_traits([[trait]], similarity = similarity_func, batch_list = batches, strip = strip_names)

					external_gwas_data = mr_utils.extract_snps_from_outcomes(exposure_variants_i['rsID'].to_list(), list(traits_to_study_mapping.values())[0], proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
					
					data = {}

					for phenotype in external_gwas_data['phenotype'].unique():
						external_gwas_data_i = external_gwas_data.loc[external_gwas_data['phenotype'] == phenotype].copy()[expected_ieugwas_column_names]

						data[phenotype] = mr_utils.harmonize(exposure_variants_i, external_gwas_data_i, '_' + exposure_pheno, '_' + phenotype)

						print('Running MR ({}) for Exposure ({}) against Outcome ({})'.format(mr_string, exposure_pheno, phenotype))
						
						mr_res = mr.run_mr(mr_type, data[phenotype], 'BETA_' + exposure_pheno, 'BETA_' + phenotype, 'SE_' + exposure_pheno, 'SE_' + phenotype, mr_args, run_all_mr)

						for mr_type_i, results in mr_res.items():
							if results is None:
								print('MR {} failed for Exposure ({}) against Outcome ({})'.format(mr_type_i.capitalize(), exposure_pheno, phenotype))
							else:

								if mr_type_i == 'PRESSO':

									presso_extras = output_base + '_' + exposure_pheno + '_presso_extras.txt'
									save_presso_extras(presso_file_path = presso_extras, results = results, presso_args = mr_args['PRESSO'])

									results = results['MR_RESULTS']

									for i in range(len(results)):
										temp_data.append([results[i][0], outcome_phenotype, results[i][1][0], results[i][2][0], results[i][3][0], data[phenotype].shape[0]])
								
								else:
									pval, beta, std_err, *additional = results
									temp_data.append([mr_type_i, phenotype, pval, beta, std_err, data[phenotype].shape[0]])
			
			mr_results = pd.DataFrame(temp_data, columns = ['MR_Method', 'Outcome', 'P_value', 'Effect_Size', 'Standard_Error', 'Number_SNPs'])

			# write the results to a file
			mr_results.to_csv(output_base + '_' + exposure_pheno + '.' + output_ext, sep = None, index = False)

	elif outcome_variants != None:
		
		# read in the outcome variants
		outcome_variants_d = pd.read_csv(outcome_variants, sep = None, engine='python')

		if outcome_variants_d.empty:
			print("File {} associated with outcome phenotype {} is empty. Exiting.".format(outcome_variants, outcome_phenotype))
			exit()
		else:
			if out_N in outcome_variants_d.columns:
				outcome_variants_d = outcome_variants_d[[out_rsID, out_CHR, out_EFFECT, out_B, out_P, out_SE, out_N]]
			else:
				outcome_variants_d = outcome_variants_d[[out_rsID, out_CHR, out_EFFECT, out_B, out_P, out_SE]]
	
		if run_all_mr == True:
			mr_string = "All Tests"
		else:
			mr_string = ", ".join(mr_type)

		for index, exposure_pheno in enumerate(exposure_phenotypes):
			
			print('Running MR ({}) for Exposure ({}) against Outcome ({})'.format(mr_string, exposure_pheno, outcome_phenotype))

			# read in the exposure variants
			exposure_variants_i = pd.read_csv(exposure_variants[index], sep = None, engine='python')

			if exposure_variants_i.empty:
				print("File {} associated with exposure phenotype {} is empty. Skipping.".format(exposure_variants[index], exposure_pheno))
				continue
			else:
				if exp_N in exposure_variants_i.columns:
					exposure_variants_i = exposure_variants_i[[exp_rsID, exp_CHR, exp_EFFECT, exp_B, exp_P, exp_SE, exp_N]]
				else:
					exposure_variants_i = exposure_variants_i[[exp_rsID, exp_CHR, exp_EFFECT, exp_B, exp_P, exp_SE]]

			# harmonize the exposure and outcome variants
			harmonized_data = mr_utils.harmonize(exposure_variants_i, outcome_variants_d, '_' + exposure_pheno, '_' + outcome_phenotype)

			# run MR
			mr_res = mr.run_mr(mr_type, harmonized_data, 'BETA_' + exposure_pheno, 'BETA_' + outcome_phenotype, 'SE_' + exposure_pheno, 'SE_' + outcome_phenotype, mr_args, run_all_mr)

			temp_data = []
			for mr_type_i, results in mr_res.items():
				if results is None:
					print('MR {} failed for Exposure ({}) against Outcome ({})'.format(mr_type_i.capitalize(), exposure_pheno, outcome_phenotype))
				else:

					if mr_type_i == 'PRESSO':
			
						presso_extras = output_base + '_' + exposure_pheno + '_presso_extras.txt'
						save_presso_extras(presso_file_path = presso_extras, results = results, presso_args = mr_args['PRESSO'])

						results = results['MR_RESULTS']

						for i in range(len(results)):
							temp_data.append([results[i][0], outcome_phenotype, results[i][1][0], results[i][2][0], results[i][3][0], harmonized_data.shape[0]])
					else:
						pval, beta, std_err, *additional = results
						temp_data.append([mr_type_i, outcome_phenotype, pval, beta, std_err, harmonized_data.shape[0]])

			# convert the results to a dataframe
			mr_results = pd.DataFrame(data = temp_data, columns = ['MR_Method', 'Outcome', 'P_value', 'Effect_Size', 'Standard_Error', 'Number_SNPs'])

			# write the results to a file
			mr_results.to_csv(output_base + '_' + exposure_pheno + '.' + output_ext, sep = '\t', index = False)

if __name__ == '__main__':
	main()