import argparse
import os
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup
import re
import requests
from os import listdir
from os.path import isfile, join
import pheWAS
import shutil
import math
from datetime import datetime

# adapted code from Alan's MI_Classes.py file
# add more mappings here if you need to fix other covariates
dict_UKB_fields_to_names =  {'f.31.0.0': 'Sex', 
							'f.34.0.0': 'Year_of_birth', 
							'f.52.0.0': 'Month_of_birth',
							'f.53.0.0': 'Date_attended_center_0', 
							'f.53.1.0': 'Date_attended_center_1',
							'f.53.2.0': 'Date_attended_center_2', 
							'f.53.3.0': 'Date_attended_center_3',
							'f.21000.0.0': 'Ethnicity', 
							'f.21000.1.0': 'Ethnicity_1', 
							'f.21000.2.0': 'Ethnicity_2',
							'f.22001.0.0': 'Sex_genetic'}

# major ethnic groups in UKBB
ethnicity_groups = ['Ethnicity.White', 'Ethnicity.Asian', 'Ethnicity.Black', 'Ethnicity.Other']

def grab_data_fields(url, out_dir):
	
	# Get the html file
	html = requests.get(url).text

	# Create beautiful soup object which can be used to parse the html file
	soup = BeautifulSoup(html, 'html.parser')
	
		# Get all tags which correspond to the data fields for ukbb
	desc_and_field_tags = soup.find_all("a", class_="subtle", href=re.compile("field.cgi"))

	# parse and grab the actual data-fields
	field_descriptions = [tag.string for tag in desc_and_field_tags]
	fields = [tag.get('href').split("=")[1] for tag in desc_and_field_tags]

	# write the data-fields to the output file
	with open(out_dir + '/data_fields.txt', 'w') as f:
		f.write("\n".join(fields))

	# create a dictionary mapping from field to description (and if the field is a date, remove it)
	field_dict = {field: description for field, description in zip(fields, field_descriptions) if 'date' not in description}

	return field_dict

def extract_variants_and_samples(variant_list, sample_IDs, bfiles_dir, geno_dir, threads, memory, keep):

	operation_on_samples = ''
	if keep:
		operation_on_samples = '--keep'
	else:
		operation_on_samples = '--remove'

	os.system('module load plink2')

	for input_file, sample_file in zip(variant_list, sample_IDs):

		# Determine the age predictor being used
		age_Predictor = input_file.split('/')[-1].split('.')[0].split('_')[0]

		# read the file of variants
		variants = pd.read_csv(input_file, sep='\t', header=0)
		chr = set(variants['CHR'].values)
		rsIDs = variants['rsID'].values

		temp_file = geno_dir + '/' + age_Predictor + '_temp_rsIDs_file.txt'

		with open(temp_file, 'w+') as f:
			f.write('\n'.join(rsIDs))

		for i in chr:

			# need to make sure that the directory used only contains bfiles and no other file, and that they all have the same prefix (before extension)
			file = ""
			for filename in os.listdir(bfiles_dir):
				if filename.startswith('chr' + str(i) + '_'): # every file in the directory starts with 'chr#_', where # is the chromosome
					file = filename.split('.')[0] # get the file name without the extension
			
			command = " ".join(["plink2", "--bfile", bfiles_dir + '/' + file, 
							"--extract", temp_file, 
							operation_on_samples, sample_file,
							"--threads", str(threads),
							"--memory", str(memory),
							"--make-bed", 
							"--out", geno_dir + '/' + age_Predictor + '_' + file])

			os.system('module load plink2 && ' + command)
			print('chr' + str(i) + ' done')

		os.remove(temp_file)

def extract_ukb_data_fields(col_names, fields_to_keep, file, write_file = False):
	"""
	This function takes a list of headers from a phenotype file and a list of fields to keep and returns a list with the UKBB data fields.
	It also writes the list to a file for future reference
	"""
	# get the column names that are in the list of columns to keep
	
	exact_col_names = []
	for field in fields_to_keep:

		# data fields are in the form of "f.field.[etc]"
		field = "f." + str(field) + "." 
		for col_name in col_names:
			if field in col_name:
				exact_col_names.append(col_name)
	
	# ensure columns are unique (use dict.fromkeys to remove duplicates - keeps order as compared to set())
	exact_col_names = list(dict.fromkeys(exact_col_names))
	
	if write_file:
		# write the list of columns to a file
		with open(file, "w") as data_file:
			data_file.write("\n".join(exact_col_names))

	# return the list with the fields to keep
	return exact_col_names

def make_out_of_sample_pheno_file(sample_files, pheno_file, data_fields, cov, keep, out_dir, data_fields_dir):
	
	data_fields = sorted(data_fields)
	data_fields_with_cov = sorted(data_fields + cov)
	exact_col_names = []
	
	# get the column names
	col_names = pd.read_csv(pheno_file, sep ='\t', nrows = 0).columns.tolist()

	# we only need to keep the phenotypes we want to run the PheWAS on
	exact_col_names = ['f.eid'] + extract_ukb_data_fields(col_names, data_fields_with_cov, "", False)

	# create file with covariate data field names
	extract_ukb_data_fields(col_names, cov, data_fields_dir + "/ukb_cov_fields.txt", True)

	# create file with data field names
	pheno_fields = extract_ukb_data_fields(col_names, data_fields, data_fields_dir + "/ukb_data_fields.txt", True)

	# full_data is the dataframe with all the data
	index = 0
	predictor_dict = {}
	sample_file_data_list = {}

	for file in sample_files:
		predictor_dict[file] = pd.DataFrame(columns = exact_col_names)
		sample_file_data_list[file] = pd.read_csv(file, sep='\t', header=None).iloc[:,0].tolist()

	# can't read entire .tab file at once, very slow, so read it in chunk by chunk
	chunksize = 1000 # number or rows we see at a time
	with pd.read_csv(pheno_file, sep = '\t', chunksize=chunksize, low_memory = False, usecols = exact_col_names) as reader:
		
		# for each chunk of data
		for chunk in reader:
			
			for file in sample_files:
				
				if keep:
					data = chunk[chunk['f.eid'].isin(sample_file_data_list[file])] # get the data for the chunk that is in the sample_IDs
				else:
					data = chunk[~chunk['f.eid'].isin(sample_file_data_list[file])] # get the data for the chunk that is not in sample_IDs  

				prev_data = predictor_dict[file]

				prev_data = pd.concat([prev_data, data]) # append the data to the full data frame

				predictor_dict[file] = prev_data

				print('chunk ' + str(index) +' for file '+ os.path.basename(file) +' done')
			index += 1

	for file in sample_files:
		# write the subsetted phenotype files to a file
		predictor_dict[file].to_csv(out_dir + '/' + os.path.basename(file).split('_')[0] + '_subsetted_pheno.tab', sep='\t', header=True, index=False)

	return predictor_dict, pheno_fields

def compute_sex(pheno):
	# Use genetic sex when available
	pheno['Sex_genetic'].fillna(pheno['Sex'], inplace = True)

	pheno.drop(['Sex'], axis=1, inplace=True)
	pheno.rename(columns={'Sex_genetic': 'Sex'}, inplace=True)
	pheno.dropna(subset=['Sex'], inplace=True)

	return pheno

def compute_age(pheno):
	# Recompute age with greater precision by leveraging the month of birth
	pheno['Year_of_birth'] = pheno['Year_of_birth'].astype(int)
	pheno['Month_of_birth'] = pheno['Month_of_birth'].astype(int)
	pheno['Date_of_birth'] = pheno.apply(
		lambda row: datetime(row.Year_of_birth, row.Month_of_birth, 15), axis=1)
	for i in range(4):
		i = str(i)
		pheno['Date_attended_center_' + i] = \
			pheno['Date_attended_center_' + i].apply(
				lambda x: pd.NaT if pd.isna(x) else datetime.strptime(x, '%Y-%m-%d'))
		pheno['Age_' + i] = pheno['Date_attended_center_' + i] - pheno['Date_of_birth']
		pheno['Age_' + i] = pheno['Age_' + i].dt.days / 365.25
		pheno.drop(['Date_attended_center_' + i], axis=1, inplace=True)
		
		# Save age as a float32 instead of float64
		pheno['Age_' + i] = np.float32(pheno['Age_' + i])

	pheno.drop(['Year_of_birth', 'Month_of_birth', 'Date_of_birth'], axis=1, inplace=True)
	pheno.dropna(how='all', subset=['Age_0', 'Age_1','Age_2', 'Age_3'], inplace=True)

	return pheno

def encode_ethnicity(pheno):
	# Fill NAs for ethnicity on instance 0 if available in other instances
	eids_missing_ethnicity = pheno['f.eid'][pheno['Ethnicity'].isna()]
	for eid in eids_missing_ethnicity:
		sample = pheno.loc[eid, :]
		if not math.isnan(sample['Ethnicity_1']):
			pheno.loc[eid, 'Ethnicity'] = pheno.loc[eid, 'Ethnicity_1']
		elif not math.isnan(sample['Ethnicity_2']):
			pheno.loc[eid, 'Ethnicity'] = pheno.loc[eid, 'Ethnicity_2']
	pheno.drop(['Ethnicity_1', 'Ethnicity_2'], axis=1, inplace=True)
	
	# One hot encode ethnicity
	dict_ethnicity_codes = {'1': 'Ethnicity.White', '1001': 'Ethnicity.British', '1002': 'Ethnicity.Irish',
							'1003': 'Ethnicity.White_Other',
							'2': 'Ethnicity.Mixed', '2001': 'Ethnicity.White_and_Black_Caribbean',
							'2002': 'Ethnicity.White_and_Black_African',
							'2003': 'Ethnicity.White_and_Asian', '2004': 'Ethnicity.Mixed_Other',
							'3': 'Ethnicity.Asian', '3001': 'Ethnicity.Indian', '3002': 'Ethnicity.Pakistani',
							'3003': 'Ethnicity.Bangladeshi', '3004': 'Ethnicity.Asian_Other',
							'4': 'Ethnicity.Black', '4001': 'Ethnicity.Caribbean', '4002': 'Ethnicity.African',
							'4003': 'Ethnicity.Black_Other',
							'5': 'Ethnicity.Chinese',
							'6': 'Ethnicity.Other_ethnicity',
							'-1': 'Ethnicity.Do_not_know',
							'-3': 'Ethnicity.Prefer_not_to_answer',
							'-5': 'Ethnicity.NA'}
	pheno['Ethnicity'] = pheno['Ethnicity'].fillna(-5).astype(int).astype(str)
	ethnicities = pd.get_dummies(pheno['Ethnicity'])
	pheno.drop(['Ethnicity'], axis=1, inplace=True)
	ethnicities.rename(columns=dict_ethnicity_codes, inplace=True)
	ethnicities['Ethnicity.White'] = ethnicities['Ethnicity.White'] + ethnicities['Ethnicity.British'] + \
									ethnicities['Ethnicity.Irish'] + ethnicities['Ethnicity.White_Other']
	ethnicities['Ethnicity.Mixed'] = ethnicities['Ethnicity.Mixed'] + \
									ethnicities['Ethnicity.White_and_Black_Caribbean'] + \
									ethnicities['Ethnicity.White_and_Black_African'] + \
									ethnicities['Ethnicity.White_and_Asian'] + \
									ethnicities['Ethnicity.Mixed_Other']
	ethnicities['Ethnicity.Asian'] = ethnicities['Ethnicity.Asian'] + ethnicities['Ethnicity.Indian'] + \
									ethnicities['Ethnicity.Pakistani'] + ethnicities['Ethnicity.Bangladeshi'] + \
									ethnicities['Ethnicity.Asian_Other']
	ethnicities['Ethnicity.Black'] = ethnicities['Ethnicity.Black'] + ethnicities['Ethnicity.Caribbean'] + \
									ethnicities['Ethnicity.African'] + ethnicities['Ethnicity.Black_Other']
	ethnicities['Ethnicity.Other'] = ethnicities['Ethnicity.Other_ethnicity'] + \
									ethnicities['Ethnicity.Do_not_know'] + \
									ethnicities['Ethnicity.Prefer_not_to_answer'] + \
									ethnicities['Ethnicity.NA']
	pheno = pheno.join(ethnicities)
	return pheno

def fix_covariate_data(pheno, phenotype_fields, sample_list, out_dir, data_fields_dir):

	for file in sample_list:
		file_prefix = os.path.basename(file).split('_')[0]

		# Formatting
		pheno[file].rename(columns=dict_UKB_fields_to_names, inplace=True)
		pheno[file].set_index('f.eid', drop=False, inplace=True)
		pheno[file].index.name = 'column_names'

		pheno[file] = compute_sex(pheno[file])
		pheno[file] = compute_age(pheno[file])
		pheno[file] = encode_ethnicity(pheno[file])

		# set of all covaraites, including the new covariate related columns added
		all_covariates = [col for col in pheno[file].columns if col not in phenotype_fields + ['f.eid']]

		# keep all columns with more than 80% non-nan values
		remove_columns = pheno[file].loc[:, all_covariates].dropna(thresh = len(pheno[file][all_covariates])*0.2, axis = 1).columns
		dropped_cols = [col for col in all_covariates if col not in remove_columns]
		cols_to_keep = [col for col in pheno[file] if col not in dropped_cols]

		print(file_prefix, dropped_cols)

		remaining_covariates = [col for col in all_covariates if col not in dropped_cols] 
		# fixed covariate file path
		covariate_file = os.path.join(data_fields_dir, file_prefix + '_fixed_cov_fields.txt')
		with open(covariate_file, 'w') as f:
			f.write('\n'.join(remaining_covariates))

		# reset the index since we changed it, and only keep a select few cols
		pheno[file] = pheno[file][cols_to_keep].reset_index().drop('column_names', axis = 1)

		pheno[file].to_csv(out_dir + '/' + file_prefix + '_subsetted_pheno_covariates_fixed.tab', index = False, sep = '\t')
	
	return pheno

def create_raw_geno_files(geno_dir, raw_dir):

	# get a list of the files in the geno_dir
	geno_files = [os.path.join(geno_dir, f.split('.bed')[0]) for f in os.listdir(geno_dir) if f.endswith('.bed')]

	for geno_file in geno_files:
		
		geno_file_prefix = geno_file.split('/')[-1].split('.')[0]
		
		# create a raw file for the genotype data
		command = " ".join(["plink2", "--export", 'A', 
							"--bfile", geno_file, 
							"--out", raw_dir + '/' + geno_file_prefix])

		os.system('module load plink2 && ' + command)

def run_pheWAS(sample_list, pheno_dict, out_dir, raw_dir, data_fields_dir, reg):

	for file in sample_list:

		file_prefix = os.path.basename(file).split('_')[0]

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
				geno_data = pd.read_csv(geno_file, sep='\t', index_col = 'IID').reindex(pheno_dict[file].loc[:, 'f.eid'].tolist())

				# reset the index first to get rid of IID column, and then since after the 'PHENOTYPE' column, the variant data starts, select everything after, and drop the phenotype column
				geno_data = geno_data.reset_index().loc[:, 'PHENOTYPE':].drop('PHENOTYPE', axis = 1)

				# get list of columns starting with 'Ethnicity'
				cols_to_remove = [col for col in cov_headers if col.startswith('Ethnicity')]

				# remove all the columns starting with ethnicity since we will only be looking at once ethnicity (ethnicity covariate is uneeded)
				cov_headers_ethnicity_fixed = [cov for cov in cov_headers if cov not in cols_to_remove]

				# only keep individuals from the chosen population
				pheno = pheno_dict[file][pheno_dict[file][ethnicity_group] == 1]
				
				# drop the f.eid column and all the ethnicities in the phenotype data since we no longer need them
				pheno = pheno.drop(cols_to_remove + ['f.eid'], axis = 1)

				# Take the difference between the covariate headers and the phenotype headers
				diff_cov_pheno = list(set(pheno.columns) - set(cov_headers_ethnicity_fixed))

				results = pheWAS.run_phewas(phenotypes = pheno, 
											genotypes = geno_data, 
											non_cov_pheno_list= diff_cov_pheno, 
											reg = reg, 
											covariates = cov_headers_ethnicity_fixed)

				# save the results to a file
				results.to_csv(results_file, sep='\t', index=False)

def save_top_variants(out_dir, field_dict):
	output = out_dir + '/aggregate_pheWAS_results.tsv'

	# add prefix to the output file for the top rated SNPS
	out_dir = os.path.dirname(output)
	out_file = os.path.basename(output)
	top_SNPS_file = out_dir + '/TOP_SNPS_' + out_file

	# create a list of the files in the out_dir (these are our pheWAS results)
	results_files = [f for f in listdir(out_dir) if isfile(join(out_dir, f))]

	# Create a dataframe to store the results
	# get the headers from the first file in the list 
	cols = ['Predictor','Ethnicity', 'Description'] + pd.read_csv(out_dir + '/' + results_files[0], sep='\t', nrows = 0).columns.tolist()
	results = pd.DataFrame(columns = cols)
	top_SNPS = pd.DataFrame(columns = cols)

	# Loop through all the results files and add them to the dataframe
	for file in results_files:
		# Read the file
		df = pd.read_csv(out_dir + '/' + file, sep = '\t', header = 0)
		
		data_fields = df['Data_Field'].tolist()

		# get the ethnic group that was included in the covariates 
		# (naming is like: Organ_chromosome_version_bfile_Ethnicity_Group_results.tsv)
		# first split gets us: 'Organ_chromosome_version_bfile' 'Group_results.tsv']
		# second split gets us: 'Group', ''
		group = file.split('_Ethnicity_')[1].split('_results.tsv')[0]

		# get age predictor prefix from file name # I.E. Pancreas_chr10_v3_bfile_Ethnicity_Asian_results.tsv --> Pancreas
		age_predictor_prefix = file.split('_')[0].split('Age')[0]

		for field in data_fields:
			# split the field at the period twice, once to get ['f', 'X.X.X'], and 
			# the second time to get ['f', 'X', 'X.X'], the first X represents the field
			field_without_periods = field.split(".", 2)[1]
			
			# Get the description, category, and category description for the field
			description = field_dict[field_without_periods]

			# Add the description, predictor, ethnicity to the dataframe
			df.loc[df['Data_Field'] == field, 'Description'] = description
			df.loc[df['Data_Field'] == field, 'Predictor'] = age_predictor_prefix
			df.loc[df['Data_Field'] == field, 'Ethnicity'] = group

		results = results.append(df)

		# Get the top SNPS based on bonferroni corrected significance level 
		# (only regressing on one SNP at a time so the number of tests is the same as the number of p vals)
		try:
			alpha = 0.05
			bonferroni_correction = alpha / sum(np.isfinite(df['p-val']))
			df = df[df['p-val'] <= bonferroni_correction].sort_values(by='p-val')

			top_SNPS = top_SNPS.append(df)

			# add the significance threshold used
			results['Bonferroni_correction'] = bonferroni_correction
			top_SNPS['Bonferroni_correction'] = bonferroni_correction

		except:
			print('Error: No p-vals found for ' + file)
			print(df['p-val'])
			continue

	# Save the results to a file
	results.to_csv(output, sep = '\t', index=False)

	# Save the top SNPs (based on pval) to a file
	top_SNPS.to_csv(top_SNPS_file, sep = '\t', index=False)

def main():
	# Parse the command line arguments
	parser = argparse.ArgumentParser(description='Pass the input html file')
	parser.add_argument('-b', '--bfiles', help='Directory containing bfiles', required=True, type = str)
	parser.add_argument('-p', '--phenotype', help='Tab file with all phenotypes', required=True, type = str)
	parser.add_argument('-u', '--url', help='Url to grab description of fields', required=True, type = str)
	parser.add_argument('-s', '--samples', action='append', help='File of sample IDs', required=True)
	parser.add_argument('-k', '--keep', help='Keep the samples (if ommitted, remove the samples)', action = 'store_true')
	parser.add_argument('-c', '--covariates', help='File of covariates data fields', required=True, type = str)
	parser.add_argument('-l', '--location', help='Location (directory) to store intermediate files', required=True, type = str)
	parser.add_argument('-v', '--variants', action='append', help='List of variant files to use', required=True)
	parser.add_argument('-t', '--threads', help='Number of threads to use', required = True, type=int)
	parser.add_argument('-m', '--memory', help='Memory to use', required = True, type=int)
	parser.add_argument('-r', '--regularization', action='store_true', help='Regularization parameter')
	parser.add_argument('-n', '--name', help='Name of the output directory', required = True, type=str)

	args = parser.parse_args()
	bfiles_dir = args.bfiles
	phenotype_file = args.phenotype
	url = args.url
	sample_list = args.samples # can be a list of files, or just one file
	keep = args.keep
	covariates_file = args.covariates
	directory = args.location
	variant_list = args.variants # can be a list of files, or just one file
	threads = args.threads
	memory = args.memory
	regularization = args.regularization
	name = args.name
	
	# make new directory for current run
	if directory == '':
		directory = '.'

	dir_path = directory + '/' + name #time.strftime("%Y_%m_%d_%H_%M_%S")
	shutil.rmtree(dir_path, ignore_errors=True)
	os.makedirs(dir_path)

	# make directories for intermediate files
	pheno_dir = dir_path + '/pheno'
	raw_dir = dir_path + '/raw'
	geno_dir = dir_path + '/geno'
	pheWAS_dir = dir_path + '/pheWAS'
	data_fields_dir = dir_path + '/data_fields'

	# try deleting the directory, if it exists, then directory will be deleted, if not, the error is ignored (EAFP)
	shutil.rmtree(pheno_dir, ignore_errors=True)
	shutil.rmtree(raw_dir, ignore_errors=True)
	shutil.rmtree(geno_dir, ignore_errors=True)
	shutil.rmtree(pheWAS_dir, ignore_errors=True)
	shutil.rmtree(data_fields_dir, ignore_errors=True)

	os.makedirs(pheno_dir)
	os.makedirs(raw_dir)
	os.makedirs(geno_dir)
	os.makedirs(pheWAS_dir)
	os.makedirs(data_fields_dir)
	
	field_dict = grab_data_fields(url = url, out_dir = data_fields_dir)

	# copy covariates file to the new directory to keep all data field files together
	shutil.copy(covariates_file, data_fields_dir + '/cov_fields.txt')

	# read covariates into a list
	with open(covariates_file, "r") as f:
		covariates = f.read().splitlines()

	# sort the two lists such that both lists have the same order (requires that the prefix of the files match between the lists, and suffix matches within the list)
	sample_list, variant_list = (list(t) for t in zip(*sorted(zip(sample_list, variant_list))))
	
	# create bfiles with associated samples removed or kept and target variants extracted
	extract_variants_and_samples(variant_list = variant_list, 
								sample_IDs = sample_list, 
								bfiles_dir = bfiles_dir, 
								geno_dir = geno_dir, 
								threads = threads, 
								memory = memory, 
								keep = keep)

	# get subsetted phenotype data + covariate data field names
	OOS_data_dict, phenof = make_out_of_sample_pheno_file(sample_files = sample_list, 
											pheno_file = phenotype_file, 
											data_fields = list(field_dict.keys()), 
											cov = covariates, 
											keep = keep,
											out_dir = pheno_dir,
											data_fields_dir = data_fields_dir)

	# fix the covariate data (one hot encoding, calculating age, etc)
	pheWAS_ready_phenos_dict = fix_covariate_data(pheno = OOS_data_dict, 
											phenotype_fields = phenof, 
											sample_list = sample_list,
											out_dir = pheno_dir, 
											data_fields_dir = data_fields_dir)

	# create raw genotype files for the pheWAS
	create_raw_geno_files(geno_dir = geno_dir, 
							raw_dir = raw_dir)

	# run the pheWAS across all the variants and predictors
	run_pheWAS(sample_list = sample_list, 
				pheno_dict = pheWAS_ready_phenos_dict, 
				out_dir = pheWAS_dir, 
				raw_dir = raw_dir, 
				data_fields_dir = data_fields_dir,
				reg = regularization)
	
	# create the final output file
	save_top_variants(out_dir = pheWAS_dir, 
							field_dict = field_dict)

if __name__ == '__main__':
	main()