import argparse
from distutils import extension
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
from glob import glob
from pathlib import Path
from multiprocessing import Process, Manager
from utility_funcs import multiple_testing_correction, annotate_genes
from datetime import datetime
import time

# sentinel lets us know if the process has finished
SENTINEL = None

# url for searching ukbb data fields (if updated in future, supply new URL with all fields selected)
# I.E. select all fields in Stability, Strata, Item Type, Value Type, & leave 4 XXXX in srch= part of url
current_year = str(datetime.now().year)
DEFAULT_SHOWCASE_URL = 'https://biobank.ndph.ox.ac.uk/showcase/search.cgi?wot=0&srch=XXXX&sta0=on&sta1=on&sta2=on&sta3=on&sta4=on&str0=on&str3=on&str1=on&str2=on&fit0=on&fit10=on&fit20=on&fit30=on&fvt11=on&fvt21=on&fvt22=on&fvt31=on&fvt41=on&fvt51=on&fvt61=on&fvt101=on&yfirst=2000&ylast=' + current_year

# adapted code from Alan's MI_Classes.py file
# add more mappings here if you need to fix other covariates
dict_UKB_fields_to_names =  {'f.31.0.0': 'Sex', 
							'f.34.0.0': 'Year_of_birth', 
							'f.52.0.0': 'Month_of_birth',
							'f.53.0.0': 'Date_attended_center_0', 
							'f.53.1.0': 'Date_attended_center_1',
							'f.53.2.0': 'Date_attended_center_2', 
							'f.53.3.0': 'Date_attended_center_3',
							'f.54.0.0': 'Center_of_attendance_0',
							'f.54.1.0': 'Center_of_attendance_1',
							'f.54.2.0': 'Center_of_attendance_2',
							'f.54.3.0': 'Center_of_attendance_3',
							'f.21000.0.0': 'Ethnicity', 
							'f.21000.1.0': 'Ethnicity_1', 
							'f.21000.2.0': 'Ethnicity_2',
							'f.22001.0.0': 'Sex_genetic'}

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

# mapping of ethnicity category to smaller categories
dict_ethnicity_mapping = {'Ethnicity.White' : ['Ethnicity.White', 'Ethnicity.British', 'Ethnicity.Irish', 'Ethnicity.White_Other'],
						'Ethnicity.Mixed' : ['Ethnicity.Mixed','Ethnicity.White_and_Black_Caribbean','Ethnicity.White_and_Black_African', 'Ethnicity.White_and_Asian', 'Ethnicity.Mixed_Other'],
						'Ethnicity.Asian': ['Ethnicity.Asian', 'Ethnicity.Indian', 'Ethnicity.Pakistani', 'Ethnicity.Bangladeshi','Ethnicity.Asian_Other'],
						'Ethnicity.Black': ['Ethnicity.Black', 'Ethnicity.Caribbean', 'Ethnicity.African','Ethnicity.Black_Other'],
						'Ethnicity.Other': ['Ethnicity.Other_ethnicity', 'Ethnicity.Do_not_know','Ethnicity.Prefer_not_to_answer','Ethnicity.NA']}

def query_ukbb(url, get_categories = False):
	# Get the html file
	html = requests.get(url).text

	# Create beautiful soup object which can be used to parse the html file
	soup = BeautifulSoup(html, 'html.parser')

	# Get all tags which correspond to the data fields and their descriptions
	desc_and_field_tags = soup.find_all("a", class_="subtle", href=re.compile("field.cgi"))
	
	# parse and grab the actual data-fields
	fields = [tag.get('href').split("=")[1] for tag in desc_and_field_tags]

	if get_categories:
		
		# get the descriptions of all the fields
		field_descriptions = [tag.string for tag in desc_and_field_tags]
		
		# Get all tags which correspond to the data field categories
		category_html_tags = soup.find_all("a", class_="subtle", href=re.compile("label.cgi"))

		# parse and grab the actual data-fields categories
		data_cat = [tag.string for tag in category_html_tags]
		
		return fields, field_descriptions, data_cat

	else:
		return fields

def get_descriptions_and_categories(fields):
	
	field_dict = {}

	# loop over the list of fields and get the descriptions and categories for each field
	for field in fields:
		
		url = DEFAULT_SHOWCASE_URL.replace('srch=XXXX', 'srch=' + field)

		# can't make too many requests in a row
		time.sleep(0.01)

		fields, field_descriptions, field_categories = query_ukbb(url, True)

		field_dict[fields[0]] = field_descriptions[0], field_categories[0]
	
	return field_dict

def grab_data_fields(url, out_dir):
	
	fields = query_ukbb(url)

	# write the data-fields to the output file
	with open(out_dir + '/data_fields.txt', 'w') as f:
		f.write("\n".join(fields))

	field_dict = get_descriptions_and_categories(fields)

	# remove fields if it is a date/Age
	for field, (description, _) in list(field_dict.items()):
		if "Date" in description.split() or "Age" in description.split():
			print("Removing field %s (Description: %s) from the list of fields to keep" % (field, description))
			del field_dict[field]

	return field_dict

def extract_variants_and_samples(sample_2_independent_dict, bfiles_dir, geno_dir, threads, memory, keep, delete_non_compressed):

	operation_on_samples = ''
	if keep:
		operation_on_samples = '--keep'
	else:
		operation_on_samples = '--remove'

	for sample_file, (variant_file, output_prefix) in sample_2_independent_dict.items():

		# read the file of variants
		variants = pd.read_csv(variant_file, sep='\t', header=0)
		chr = set(variants['CHR'].values)
		rsIDs = variants['rsID'].values

		temp_file = geno_dir + '/' + output_prefix + '_temp_rsIDs_file.txt'
		# temp_files.append(temp_file)

		with open(temp_file, 'w+') as f:
			f.write('\n'.join(rsIDs))

		for i in chr:

			# every file in the directory starts with 'chr#[some non numeric character]', where # is the chromosome
			# another way of doing this --> re.compile(r'chr1[^0-9]')
			regex = re.compile(r'chr' + str(i) + '\D') # \d represents digit, \D represents non-digit

			# need to make sure that the directory used only contains bfiles and no other file, and that they all have the same prefix (before extension)
			file = ""
			extension = ""
			for filename in os.listdir(bfiles_dir):
				if regex.match(filename): 
					if filename.endswith('bed') or filename.endswith('zst'):
						file = filename.split('.')[0] # get the file name without the extension
						extension = filename.split('.')[-1] # get the extension

			# if the extension is zst, we need to first extract the file
			if extension == 'zst':
				
				# check to see if only the .zst file is in the directory (if so, we need to extract it)
				if not os.path.exists(bfiles_dir + '/' + file + '.bed'):
					decompress_command = " ".join(['plink2', '--zd', bfiles_dir + '/' + file + '.bed.zst', bfiles_dir + '/' + file + '.bed'])
					os.system('module load plink2 && ' + decompress_command)
					# os.system('zstd -d ' + bfiles_dir + '/' + file + '.bed.zst')

			extraction_command = " ".join(["plink2", "--bed", bfiles_dir + '/' + file + '.bed',
							"--bim", bfiles_dir + '/' + file + '.bim',
							"--fam", bfiles_dir + '/' + file + '.fam',
							"--extract", temp_file, 
							operation_on_samples, sample_file,
							"--threads", str(threads),
							"--memory", str(memory),
							"--make-bed", 
							"--out", geno_dir + '/' + output_prefix + '_' + file])
			
			os.system('module load plink2 && ' + extraction_command)
			print('chr' + str(i) + ' done')

		os.remove(temp_file)

	if delete_non_compressed:
		for file in glob(bfiles_dir + '/*.bed'):
			os.remove(file)

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

def make_out_of_sample_pheno_file(sample_2_independent_dict, pheno_file, data_fields, cov, keep, out_dir, data_fields_dir, save_phenos, pickle_intermediates, pickle_protocol, phenotype_phewas):
	
	# get the column names from the phenotype file
	col_names = pd.read_csv(pheno_file, sep ='\t', nrows = 0).columns.tolist()

	# sort the phenotype fields + extract all the UKBB fields corresponding to the phenotype fields
	data_fields = sorted(data_fields)
	pheno_fields = extract_ukb_data_fields(col_names, data_fields, data_fields_dir + "/ukb_data_fields.txt", True)

	if cov is not None:
		# sort the data fields + the covariates
		data_fields = sorted(data_fields + cov)

		# create file with covariate data field names
		extract_ukb_data_fields(col_names, cov, data_fields_dir + "/ukb_cov_fields.txt", True)
	
	# we only need to keep the phenotypes we want to run the PheWAS on + if the covariates exist, they are added to the list of data fields
	exact_col_names = ['f.eid'] + extract_ukb_data_fields(col_names, data_fields, "", False)

	# full_data is the dataframe with all the data
	index = 0
	predictor_dict = {}
	sample_file_data_list = {}
	
	sample_2_phenos_dict = {}

	for file, (phenotypes, output_prefix) in sample_2_independent_dict.items():
		predictor_dict[file] = pd.DataFrame(columns = exact_col_names)
		sample_file_data_list[file] = pd.read_csv(file, sep='\t', header=None).iloc[:,0].tolist()

		sample_2_phenos_dict[file] = extract_ukb_data_fields(col_names, phenotypes, data_fields_dir + '/' + output_prefix + '_ukb_pheno_fields.txt', True)

	# can't read entire .tab file at once, very slow, so read it in chunk by chunk
	chunksize = 1000 # number or rows we see at a time
	with pd.read_csv(pheno_file, sep = '\t', chunksize=chunksize, low_memory = False, usecols = exact_col_names) as reader:
		
		# for each chunk of data
		for chunk in reader:
			
			for file in sample_2_phenos_dict.keys():
				
				if phenotype_phewas:
					phenotypes_to_remove = []
					for sfile, pheno_list in sample_2_phenos_dict.items():
						if sfile != file:
							phenotypes_to_remove += pheno_list

					# drop the phenotypes that are not associated with this file from the columns
					chunk = chunk.drop(phenotypes_to_remove, axis=1)

				if keep:
					data = chunk[chunk['f.eid'].isin(sample_file_data_list[file])] # get the data for the chunk that is in the sample_IDs
				else:
					data = chunk[~chunk['f.eid'].isin(sample_file_data_list[file])] # get the data for the chunk that is not in sample_IDs  

				prev_data = predictor_dict[file]

				prev_data = pd.concat([prev_data, data]) # append the data to the full data frame

				predictor_dict[file] = prev_data

				print('chunk ' + str(index) +' for file '+ os.path.basename(file) +' done')
			index += 1

	if save_phenos:
		for file, (_, output_prefix) in sample_2_independent_dict.items():

			if pickle_intermediates:
				predictor_dict[file].to_pickle(out_dir + '/' + output_prefix + '_subsetted_pheno.pkl', protocol = pickle_protocol)
			else:
				# write the subsetted phenotype files to a file
				predictor_dict[file].to_csv(out_dir + '/' + output_prefix + '_subsetted_pheno.tab', sep='\t', header=True, index=False)

	return predictor_dict, pheno_fields

def compute_sex(pheno):
	# Use genetic sex when available
	pheno['Sex_genetic'].fillna(pheno['Sex'], inplace = True)

	pheno.drop(['Sex'], axis=1, inplace=True)
	pheno.rename(columns={'Sex_genetic': 'Sex'}, inplace=True)
	pheno.dropna(subset=['Sex'], inplace=True)

	return pheno

def compute_attendance_center(pheno):
	
	# one hot encode the attendance centers
	for i in range(4):
		pheno = pd.get_dummies(pheno, prefix = 'center_' + str(i), columns = ['Center_of_attendance_' + str(i)])

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
	

	pheno['Ethnicity'] = pheno['Ethnicity'].fillna(-5).astype(int).astype(str)
	eth = pd.get_dummies(pheno['Ethnicity'])
	pheno.drop(['Ethnicity'], axis=1, inplace=True)
	eth.rename(columns=dict_ethnicity_codes, inplace=True)

	for key, value in dict_ethnicity_mapping.items():
		for eth_cat in value:
			if eth_cat not in eth:
				continue
			else:
				if key not in eth:
					eth[key] = eth[eth_cat]
				else:
					eth[key] += eth[eth_cat]

	pheno = pheno.join(eth)
	return pheno

def fix_covariate_data(pheno_dict, phenotype_fields, sample_2_independent_dict, out_dir, data_fields_dir, pickle_intermediates, pickle_protocol, covariate_name_map):

	for file, (_, output_prefix) in sample_2_independent_dict.items():

		# Formatting
		pheno_dict[file].rename(columns=covariate_name_map, inplace=True)
		pheno_dict[file].set_index('f.eid', drop=False, inplace=True)
		pheno_dict[file].index.name = 'column_names'
		
		if 'Sex' in pheno_dict[file].columns and 'Sex_genetic' in pheno_dict[file].columns:
			pheno_dict[file] = compute_sex(pheno_dict[file])

		if ('Year_of_birth' in pheno_dict[file].columns and 
			'Month_of_birth' in pheno_dict[file].columns and 
			'Date_attended_center_0' in pheno_dict[file].columns):
			pheno_dict[file] = compute_age(pheno_dict[file])
		
		if 'Ethnicity' in pheno_dict[file].columns:
			pheno_dict[file] = encode_ethnicity(pheno_dict[file])

		# set of all covariates, including the new covariate related columns added
		all_covariates = [col for col in pheno_dict[file].columns if col not in phenotype_fields + ['f.eid']]

		# keep all columns with more than 80% non-nan values
		remove_columns = pheno_dict[file].loc[:, all_covariates].dropna(thresh = len(pheno_dict[file][all_covariates])*0.2, axis = 1).columns
		dropped_cols = [col for col in all_covariates if col not in remove_columns]
		cols_to_keep = [col for col in pheno_dict[file] if col not in dropped_cols]

		print(output_prefix, dropped_cols)

		remaining_covariates = [col for col in all_covariates if col not in dropped_cols] 
		# fixed covariate file path
		covariate_file = os.path.join(data_fields_dir, output_prefix + '_fixed_cov_fields.txt')
		with open(covariate_file, 'w') as f:
			f.write('\n'.join(remaining_covariates))

		# reset the index since we changed it, and only keep a select few cols
		pheno_dict[file] = pheno_dict[file][cols_to_keep].reset_index().drop('column_names', axis = 1)

		if 'Center_of_attendance' in pheno_dict[file].columns:
			pheno_dict[file] = compute_attendance_center(pheno_dict[file])

		if pickle_intermediates:
			pheno_dict[file].to_pickle(out_dir + '/' + output_prefix + '_subsetted_pheno_covariates_fixed.pkl', protocol= pickle_protocol)
		else:
			pheno_dict[file].to_csv(out_dir + '/' + output_prefix + '_subsetted_pheno_covariates_fixed.tab', index = False, sep = '\t')
	return pheno_dict

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

def phewas_target(thresh, reg, target_args):
	
	while True:
		args_for_phewas = target_args.get()

		if args_for_phewas == SENTINEL:
			print('phewas_target: sentinel received')
			return
		else:
			pheno, geno, non_cov, cov, results_file = args_for_phewas

			results = pheWAS.run_phewas(phenotypes = pheno, 
										independent_variables = geno, 
										non_cov_pheno_list= non_cov, 
										reg = reg, 
										covariates = cov,
										thresh = thresh)

			if results is not None:
				results.to_csv(results_file, index = False, sep = '\t')
			else:
				print('No pheWAS results for ', results_file)
				return

def run_genotype_pheWAS(sample_2_independent_dict, pheno_dict, out_dir, raw_dir, data_fields_dir, reg, thresh, threads, add_covariates, ethnicity_groups):

	manager = Manager()
	process_queue = manager.Queue(threads)

	# start workers, pointing towards the function to run (phewas_target), and the arguments supplied (placed on queue)       
	pool = []
	for _ in range(threads):
		p = Process(target=phewas_target, args=(thresh, reg, process_queue))
		p.start()
		pool.append(p)

	for file, (_, output_prefix) in sample_2_independent_dict.items():

		# find the raw genotype files with the same prefix
		raw_geno_files = [os.path.join(raw_dir, f) for f in os.listdir(raw_dir) if f.endswith('.raw') and os.path.basename(f).startswith(output_prefix + '_')]
		
		cov_headers = []
		if add_covariates is not None:
			# read covariate fixed headers file into list
			with open(data_fields_dir + '/' +  output_prefix + '_fixed_cov_fields.txt', "r") as f:
				cov_headers = f.read().splitlines()

		for geno_file in raw_geno_files:
			
			# get the prefix of the geno_file name (before the .raw)
			geno_file_prefix = geno_file.split('/')[-1].split('.')[0]

			# set index column to the IID and then rearrange the rows according to the order of the phenotype file by the f.eid column
			geno_data = pd.read_csv(geno_file, sep='\t', index_col = 'IID').reindex(pheno_dict[file].loc[:, 'f.eid'].tolist())

			# reset the index first to get rid of IID column, and then since after the 'PHENOTYPE' column, the variant data starts, select everything after, and drop the phenotype column
			geno_data = geno_data.reset_index().loc[:, 'PHENOTYPE':].drop('PHENOTYPE', axis = 1)

			if add_covariates is not None and 'Ethnicity' in " ".join(cov_headers):
				for ethnicity_group in ethnicity_groups:
					
					results_file = out_dir + '/' + geno_file_prefix + '_'+ ethnicity_group.replace('.', '_') + '_results.tsv'

					print('Running PheWAS for ' + results_file)

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

					

					# add each the phewas files needed to run each phewas to the queue for the workers to use
					process_queue.put((pheno, geno_data, diff_cov_pheno, cov_headers_ethnicity_fixed, results_file))
			else:

				results_file = out_dir + '/' + geno_file_prefix + '_results.tsv'

				print('Running PheWAS for ' + results_file)

				# drop the f.eid column since we no longer need it
				pheno = pheno_dict[file].drop(['f.eid'], axis = 1)

				# add each the phewas files needed to run each phewas to the queue for the workers to use
				process_queue.put((pheno, geno_data, pheno.columns, cov_headers, results_file))
	
	for _ in range(threads):
		process_queue.put(SENTINEL)

	# wait for all the workers to finish before continuing with the main code
	for p in pool:
		p.join()

def run_phenotype_pheWAS(sample_2_independent_dict, pheno_dict, out_dir, data_fields_dir, reg, thresh, threads, add_covariates, ethnicity_groups):

	manager = Manager()
	process_queue = manager.Queue(threads)

	# start workers, pointing towards the function to run (phewas_target), and the arguments supplied (placed on queue)       
	pool = []
	for _ in range(threads):
		p = Process(target=phewas_target, args=(thresh, reg, process_queue))
		p.start()
		pool.append(p)

	for file, (_, output_prefix) in sample_2_independent_dict.items():

		cov_headers = []
		if add_covariates is not None:
			# read covariate fixed headers file into list
			with open(data_fields_dir + '/' +  output_prefix + '_fixed_cov_fields.txt', "r") as f:
				cov_headers = f.read().splitlines()
		
		independent_variable_phenotypes = sample_2_independent_dict[file][0]

		if add_covariates is not None and 'Ethnicity' in " ".join(cov_headers):
			for ethnicity_group in ethnicity_groups:
				
				results_file = out_dir + '/' + output_prefix + '_'+ ethnicity_group.replace('.', '_') + '_results.tsv'

				print('Running PheWAS for ' + results_file)

				# get list of columns starting with 'Ethnicity'
				cols_to_remove = [col for col in cov_headers if col.startswith('Ethnicity')]

				# remove all the columns starting with ethnicity since we will only be looking at once ethnicity (ethnicity covariate is uneeded)
				cov_headers_ethnicity_fixed = [cov for cov in cov_headers if cov not in cols_to_remove]

				# only keep individuals from the chosen population
				pheno = pheno_dict[file][pheno_dict[file][ethnicity_group] == 1]

				# get the independent phenotypes and then drop them from the phenotype data along with the ethnicities and the f.eid
				independent_pheno =  pheno[independent_variable_phenotypes]
				pheno = pheno.drop(cols_to_remove + ['f.eid'] + independent_variable_phenotypes, axis = 1)

				# Take the difference between the covariate headers and the phenotype headers
				diff_cov_pheno = list(set(pheno.columns) - set(cov_headers_ethnicity_fixed))

				# add each the phewas files needed to run each phewas to the queue for the workers to use
				process_queue.put((pheno, independent_pheno, diff_cov_pheno, cov_headers_ethnicity_fixed, results_file))
		else:

			results_file = out_dir + '/' + output_prefix + '_results.tsv'

			print('Running PheWAS for ' + results_file)

			# get the independent phenotypes and then drop them from the phenotype data along with the f.eid
			independent_pheno = pheno[independent_variable_phenotypes]
			pheno = pheno_dict[file].drop(['f.eid'] + independent_variable_phenotypes, axis = 1)

			# add each the phewas files needed to run each phewas to the queue for the workers to use
			process_queue.put((pheno, pheno[independent_variable_phenotypes], pheno.columns, cov_headers, results_file))
	
	for _ in range(threads):
		process_queue.put(SENTINEL)

	# wait for all the workers to finish before continuing with the main code
	for p in pool:
		p.join()

def save_top_variants(out_dir, field_dict, correction, alpha, gene_file = None, sample_2_independent_dict = None, downstream = 40, upstream = 40):
	output = out_dir + '/aggregate_pheWAS_results.tsv'
	
	# name significant results file based on correction method
	correction_method = correction + '_'
	top_SNPS_file = out_dir + '/TOP_SNPS_' + correction_method + '_aggregate_pheWAS_results.tsv'

	# create a list of the files in the out_dir (these are our pheWAS results)
	results_files = [f for f in listdir(out_dir) if isfile(join(out_dir, f))]

	if len(results_files) == 0:
		print('No pheWAS results files found in ' + out_dir)
		return

	# get the headers from the first file in the list and add extra headers depending on whether phewas were ran separately for ethniciities
	cols = []
	if 'Ethnicity' in results_files[0]:
		cols = ['PheWAS_Category','Ethnicity', 'Description']
	else:
		cols = ['PheWAS_Category', 'Description']
	
	cols += pd.read_csv(out_dir + '/' + results_files[0], sep='\t', nrows = 0).columns.tolist() + ['Category']

	# Create a dataframe to store the results
	results = pd.DataFrame(columns = cols)

	# description is 1st value
	# category is the 2nd value 
	DESC = 0
	CAT = 1

	# Loop through all the results files and add them to the dataframe
	for file in results_files:
		# Read the file
		df = pd.read_csv(out_dir + '/' + file, sep = '\t', header = 0)
		
		data_fields = df['Data_Field'].tolist()

		group = ''
		if 'Ethnicity' in file:
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
			description = field_dict[field_without_periods][DESC]
			category = field_dict[field_without_periods][CAT]

			# Add the description, predictor, ethnicity to the dataframe
			df.loc[df['Data_Field'] == field, 'Description'] = description
			df.loc[df['Data_Field'] == field, 'PheWAS_Category'] = age_predictor_prefix
			df.loc[df['Data_Field'] == field, 'Category'] = category

			if 'Ethnicity' in file:
				df.loc[df['Data_Field'] == field, 'Ethnicity'] = group

		results = results.append(df)

	if gene_file is not None:
		# create a aggregate dataframe of all the snps
		rsid_df = pd.DataFrame()

		# rename Independent_Var to rsID
		results = results.rename(columns = {'Independent_Var': 'rsID'})

		for _, (rsid_f, _) in sample_2_independent_dict.items():
			# read f and then concat to rsid_df
			rsid_df = pd.concat([rsid_df, pd.read_csv(rsid_f, sep = '\t')])
			rsid_df = rsid_df.drop_duplicates()
		
		rsid_df['rsID'] = rsid_df['rsID'].apply(lambda x: x.split('_')[0])
		old_results_rsIDs = results['rsID'].tolist()
		results['rsID'] = results['rsID'].apply(lambda x: x.split('_')[0])

		for _, (_, output_prefix) in sample_2_independent_dict.items():
			results, _ = annotate_genes(gene_file = gene_file, 
												rsid_df = rsid_df, 
												down = downstream, 
												up = upstream, 
												regressions = results, 
												output_dir = out_dir, 
												pheno_name = output_prefix)
		
		results = results.rename(columns = {'rsID': 'Independent_Var'})
		results['Independent_Var'] = old_results_rsIDs
		
	# Save the results to a file
	results.to_csv(output, sep = '\t', index=False)

	try:
		significance_threshold = multiple_testing_correction(pvalues = results['p-val'], alpha = alpha, method = correction)
		results = results[results['p-val'] <= significance_threshold].sort_values(by='p-val')
	except:
		print('Error: No p-vals found for all files')

	# Save the top SNPs (based on pval) to a file
	results.to_csv(top_SNPS_file, sep = '\t', index=False)

	return results

def annotate_function(top_results, out_dir):
	pass

def main():
	# Parse the command line arguments
	parser = argparse.ArgumentParser(description='Run PheWAS using either genotype or phenotype data from the UKBB.')

	# Arguments for PheWAS (either genotype based or phenotype)
	parser.add_argument('--bfiles_directory', help='Directory containing bfiles', required=False, type = str)
	parser.add_argument('--variant_file', help='List of files with genetic variants to use as independent variables', required = False, action='append')
	parser.add_argument('--phenotypes_file', help = 'List of files with phenotypes to use as independent variables', required = False, action='append')
	parser.add_argument('--phenotypes_list', help = 'List of phenotypes [data field #] to use as independent variables', required = False, nargs='+', action='append')
	parser.add_argument('--ukbiobank_phenotype_file', help='UKBB phenotype tab file with all phenotypes', required=True, type = str)
	parser.add_argument('--ukbiobank_url', help='Url to grab description of fields from UKBB', required=False, default = None, type = str)
	parser.add_argument('--dependent_phenotypes_file', help = 'File with list of phenotypes [data field #] to use as dependent variables', required = False, default = None, type = str)
	parser.add_argument('--dependent_phenotypes_list', help = 'List of phenotypes [data field #] to use as dependent variables', required = False, nargs='+',)
	parser.add_argument('--sample_file', help='File (or list of files) of sample IDs', required=True, action='append')
	parser.add_argument('--keep', help='Keep the samples (if ommitted, remove the samples)', required=False, default = False, action = 'store_true')
	parser.add_argument('--alpha', help = 'Alpha value to use for multiple testing', required = False, default = 0.05, type = float)
	parser.add_argument('--correction', help = 'Correction method to use', required = True, default = 'bonferroni', choices = ['bonferroni', 'sidak', 'fdr_bh', 'no_correction'])
	parser.add_argument('--covariates_file', help='File of covariate data fields', required=False, default = None, type = str)
	parser.add_argument('--ethnicity_groups', help = 'Ethnicity groups to use for pheWAS', required = False, default = ['Ethnicity.White'], choices = ['Ethnicity.White', 'Ethnicity.Asian', 'Ethnicity.Black', 'Ethnicity.Other'],  nargs = '+')
	parser.add_argument('--output_prefix', help = 'Prefix for output files (specify one for each sample file)', required = True, action = 'append')

	# Data logging / location arguments
	parser.add_argument('--directory_name', help='Directory to store intermediate files', required=True, type = str)
	parser.add_argument('--save_phenotypes', help = 'If you want to save the phenotypes, pass this flag', required = False, default = False, action = 'store_true')
	parser.add_argument('--pickle_intermediates', help = 'If you want to pickle the intermediate files, pass this flag', required = False, default = False, action = 'store_true')
	parser.add_argument('--pickle_protocol', help = 'Sets the protocol to use when pickling the intermediate files', required = False, default = -1, type = int)
	
	# Change performance of PheWAS
	parser.add_argument('--memory', help='Memory to use', required = True, type=int)
	parser.add_argument('--threads', help='Number of threads to use', required = False, default = 1, type=int)
	parser.add_argument('--regularization', help='Regularization parameter', required = False, default = False, action='store_true')
	parser.add_argument('--sample_threshold', help='Threshold for the minimum number of samples', default = 1000, required = False, type=int)

	# Reusing old data
	parser.add_argument('--reuse_phenotypes', help ='If you already have the phenotypes, pass the directory name. The phenotype files should start with the passed output prefixes.', required = False, default = '', type = str)
	parser.add_argument('--reuse_fixed_phenotypes', help ='If you already have the fixed phenotypes, pass the directory name. The phenotype files should start with the passed output prefixes.', required = False, default = '', type = str)
	parser.add_argument('--reuse_genos', help ='If you already have the extracted bfiles, pass the directory name', required = False, default = '', type = str)
	parser.add_argument('--reuse_raw_bfiles', help ='If you already have the raw bfiles, pass the directory name', required = False, default = '', type = str)
	parser.add_argument('--old_data_fields_dir', help = 'If you want to reuse phenotypes, you need to specify the directory where the old data fields are stored', required = False, default = '', type = str)
	
	# Miscellaneous arguments
	parser.add_argument('--delete_non_compressed_files', help = 'If you want to delete the non-compressed bed files', required = False, default = False, action = 'store_true')
	parser.add_argument('--covariate_name_map', help = 'File containing map from ukbb field name to custom covariate name [should be pandas dataframe]', required = False, default = None, type = str)
	parser.add_argument('--data_showcase_search_url', help = 'Url to search for data showcase', required = False, default = None, type = str)
	parser.add_argument('--gene_file', help = 'File containing list of genes and their positions. Used to annotate variants with nearby genes', required = False, default = None, type = str)
	parser.add_argument('--downstream', help = 'Number in KB to check for closest genes', required = False, default = 40, type = int)
	parser.add_argument('--upstream', help = 'Number in KB to check for closest genes', required = False, default = 40, type = int)
	parser.add_argument('--annotate', help = 'If you want to annotate the variants / genes with their function / consequences', required = False, default = False, action = 'store_true')

	args = parser.parse_args()
	
	# Arguments for PheWAS (either genotype based or phenotype)
	bfiles_directory = args.bfiles_directory
	variant_files = args.variant_file
	phenotypes_files = args.phenotypes_file
	phenotypes_lists = args.phenotypes_list
	ukbiobank_phenotype_file = args.ukbiobank_phenotype_file
	ukbiobank_url = args.ukbiobank_url
	dependent_phenotypes_file = args.dependent_phenotypes_file
	dependent_phenotypes_list = args.dependent_phenotypes_list
	sample_files = args.sample_file
	keep = args.keep
	alpha = args.alpha
	correction = args.correction
	covariates_file = args.covariates_file
	ethnicity_groups = args.ethnicity_groups
	output_prefixes = args.output_prefix

	# Data logging / location arguments
	directory_name = args.directory_name
	save_phenotypes = args.save_phenotypes
	pickle_intermediates = args.pickle_intermediates
	pickle_protocol = args.pickle_protocol

	# Change performance of PheWAS
	memory = args.memory
	threads = args.threads
	regularization = args.regularization
	sample_threshold = args.sample_threshold

	# Reusing old data
	reuse_phenotypes = args.reuse_phenotypes
	reuse_fixed_phenotypes = args.reuse_fixed_phenotypes
	reuse_genos = args.reuse_genos
	reuse_raw_bfiles = args.reuse_raw_bfiles
	old_data_fields_dir = args.old_data_fields_dir
	
	# Miscellaneous arguments
	delete_non_compressed_files = args.delete_non_compressed_files
	covariate_name_map = args.covariate_name_map
	data_showcase_search_url = args.data_showcase_search_url
	gene_file = args.gene_file
	downstream = args.downstream
	upstream = args.upstream
	annotate = args.annotate

	# ---------------------------------------VERIFY ARGS--------------------------------------- #

	genotype_phewas = False
	phenotype_phewas = False

	if variant_files and not phenotypes_files and not phenotypes_lists:
		genotype_phewas = True

		if len(variant_files) != len(sample_files):
			print('Error: You must specify the same number of phenotype files and sample files.')
			exit()
	elif phenotypes_files and not variant_files and not phenotypes_lists:
		phenotype_phewas = True

		if len(phenotypes_files) != len(sample_files):
			print('Error: You must specify the same number of phenotype files and sample files.')
			exit()
	elif phenotypes_lists and not variant_files and not phenotypes_files:
		phenotype_phewas = True

		if len(phenotypes_lists) != len(sample_files):
			print('Error: You must specify the same number of phenotype lists and sample files.')
			exit()
	else:
		print(("Error: You must specify either:\n"
				"\t1) A set of files containing phenotypes.\n" 
				"\t2) A set of files containing genetic variants.\n"
				"\t3) A set of lists containing phenotypes.\n"
				"Choose one and try again."))
		exit()

	if ((ukbiobank_url is not None and dependent_phenotypes_file is not None and not dependent_phenotypes_list) or 
		(ukbiobank_url is None and dependent_phenotypes_file is None and dependent_phenotypes_list)):
		print(("Error: You must specify either:\n"
				"\t1) A url from UKBB to grab phenotypes from (data showcase).\n" 
				"\t2) A set of files containing phenotypes we will regress against.\n"
				"\t3) A set of lists containing phenotypes we will regress against.\n"
				"Choose one and try again."))
		exit()

	if dependent_phenotypes_file is not None and dependent_phenotypes_list:
		print('Warning: You have specified both a file containing dependent phenotypes and a list of dependent phenotypes. The list will be used.')
		dependent_phenotypes_file = None

	if reuse_phenotypes != '' and old_data_fields_dir == '':
		print('Error: You must specify the directory where the old data fields are stored if you want to reuse phenotypes')
		exit()
	
	if reuse_phenotypes != '' and reuse_fixed_phenotypes != '':
		print('Warning: You provided both a directory to reuse phenotypes and a directory to reuse covariate fixed phenotypes. Only the fixed phenotypes will be reused.')
		reuse_fixed_phenotypes = ''

	if len(output_prefixes) != len(sample_files):
		print('Error: You must specify the same number of output prefixes as sample files.')
		exit()

	if phenotype_phewas and gene_file is not None:
		print('Warning: You have specified a gene file but are running a phenotype phewas. The gene file will be ignored.')
		gene_file = None

	if phenotype_phewas and annotate:
		print('Warning: You have specified to annotate variants / genes but are running a phenotype phewas. This argument will be ignored.')
		annotate = False

	if not phenotype_phewas and (downstream < 0 or upstream < 0):
		print('Error: You must specify a positive number for the number of KB downstream and upstream.')
		exit()

	if data_showcase_search_url is not None:
		global DEFAULT_SHOWCASE_URL 
		DEFAULT_SHOWCASE_URL = data_showcase_search_url

	# ----------------------------------------------------------------------------------------- #

	# ---------------------------------------CREATE DIRS--------------------------------------- #
	# make new directory for current run
	if directory_name == '':
		directory_name = '.'

	dir_path = directory_name 

	p = Path(dir_path)

	# if the directory already exists, don't delete it (we might be using files from inside it)
	if not p.exists():
		print('Creating directory: ' + dir_path)
		os.makedirs(dir_path)

	# make directories for intermediate files
	geno_dir = dir_path + '/geno'
	raw_dir = dir_path + '/raw'
	pheno_dir = dir_path + '/pheno'
	data_fields_dir = dir_path + '/data_fields'
	pheWAS_dir = dir_path + '/pheWAS'

	# compare the directories we want to create to the passed directory, if they are the same, don't delete the directory
	if genotype_phewas and geno_dir != reuse_genos:
		
		if reuse_genos == '': # only remake geno_dir if we aren't using another directory for the geno files
			print('deleting geno_dir since we aren\'t reusing it')
			shutil.rmtree(geno_dir, ignore_errors=True)
			os.makedirs(geno_dir)

	if genotype_phewas and raw_dir != reuse_raw_bfiles:
		
		if reuse_raw_bfiles == '': # only remake raw_bfile dir if we aren't using another directory for the raw files
			print('deleting raw_dir since we aren\'t reusing it')
			shutil.rmtree(raw_dir, ignore_errors=True)
			os.makedirs(raw_dir)

	if pheno_dir != reuse_phenotypes and pheno_dir != reuse_fixed_phenotypes:
		print('deleting pheno_dir since we aren\'t reusing it')
		shutil.rmtree(pheno_dir, ignore_errors=True)
		os.makedirs(pheno_dir)
	
	if data_fields_dir != old_data_fields_dir:
		print('deleting data_fields_dir since we aren\'t reusing it')
		shutil.rmtree(data_fields_dir, ignore_errors=True)
		os.makedirs(data_fields_dir)
	
	print('deleting and recreating old pheWAS_dir (if it exists)')
	shutil.rmtree(pheWAS_dir, ignore_errors=True)
	os.makedirs(pheWAS_dir)
	# ----------------------------------------------------------------------------------------- #

	# ---------------------------GET FIELD DESCRIPTIONS / CATEGORIES--------------------------- #
	field_dict = {}
	if ukbiobank_url is not None:
		field_dict = grab_data_fields(url = ukbiobank_url, out_dir = data_fields_dir)
	elif dependent_phenotypes_file is not None:
		# load the file containing the dependent_phenotypes into a list
		
		with open(dependent_phenotypes_file, 'r') as f:
			dependent_phenotypes = f.read().splitlines()

		field_dict = get_descriptions_and_categories(fields = dependent_phenotypes)
	else: # dependent_phenotypes_list is not None
		field_dict = get_descriptions_and_categories(fields = dependent_phenotypes_list)
	# ----------------------------------------------------------------------------------------- #

	covariates = None
	if covariates_file is not None:
		# copy covariates file to the new directory to keep all data field files together
		shutil.copy(covariates_file, data_fields_dir + '/cov_fields.txt')
		
		# read covariates into a list
		covariates = pd.read_csv(covariates_file, sep = '\t', usecols=[0], dtype = str).iloc[:,0].values.tolist()

	# -----------------------------------MAP SAMPLES 2 INDEP----------------------------------- #

	samples_2_ind = {}
	field_dict_for_independent_phenos = {}

	if variant_files:

		# sort the two lists (requires that the prefix of the files/strings match between the lists)
		sample_files, variant_files, output_prefixes =  (sorted(x) for x in [sample_files, variant_files, output_prefixes])

		# map them to each other
		samples_2_ind = {sfile:(indfile, oprefix) for sfile, indfile, oprefix in zip(sample_files, variant_files, output_prefixes)}

	elif phenotypes_files:

		# sort the two lists (requires that the prefix of the files/strings match between the lists)
		sample_files, phenotypes_files, output_prefixes =  (sorted(x) for x in [sample_files, phenotypes_files, output_prefixes])
		
		# then read in the phenotypes and then map 
		for pfile, sfile, oprefix in zip(phenotypes_files, sample_files, output_prefixes):
			# read the phenotype file into a list
			with open(pfile, 'r') as f:
				phenotypes = f.read().splitlines()

			for p in phenotypes:
				if covariates is not None and p in covariates:
					print('Error: The phenotype ' + p + ' is in the list of covariates. Please remove it from the covariates file or from the phenotypes and re-run.')
					exit()
				if p in field_dict.keys():
					print('Warning: The phenotype ' + p + ' is in the list of dependent variables. It will be removed from the list of dependent variables.')
					field_dict.pop(p)
			
			temp_dict = get_descriptions_and_categories(fields = phenotypes)

			# combine the two dictionaries
			field_dict_for_independent_phenos = {**field_dict_for_independent_phenos, **temp_dict}

			samples_2_ind[sfile] = (phenotypes, oprefix)

	else:
		samples_2_ind = {sfile:(indfile, oprefix) for sfile, indfile, oprefix in zip(sorted(sample_files), phenotypes_lists, sorted(output_prefixes))}

		for plist in phenotypes_lists:
		
			for p in plist:
				if covariates is not None and p in covariates:
					print('Error: The phenotype ' + p + ' is in the list of covariates. Please remove it from the covariates file or from the phenotypes and re-run.')
					exit()
				if p in field_dict.keys():
					print('Warning: The phenotype ' + p + ' is in the list of dependent variables. It will be removed from the list of dependent variables.')
					field_dict.pop(p)

			temp_dict = get_descriptions_and_categories(fields = plist)

			# combine the two dictionaries
			field_dict_for_independent_phenos = {**field_dict_for_independent_phenos, **temp_dict}
	
	# in case we have phenotypes as independent_variables as well, combine these dictionaries
	field_dict = {**field_dict, **field_dict_for_independent_phenos}
	# ----------------------------------------------------------------------------------------- #

	if reuse_genos == '' and genotype_phewas:
		# create bfiles with associated samples removed or kept and target variants extracted
		extract_variants_and_samples(sample_2_independent_dict = samples_2_ind, 
									bfiles_dir = bfiles_directory, 
									geno_dir = geno_dir, 
									threads = threads, 
									memory = memory, 
									keep = keep,
									delete_non_compressed = delete_non_compressed_files)
	elif genotype_phewas:
		print('Reusing genotype files')
		geno_dir = reuse_genos
	
	if reuse_raw_bfiles == '' and genotype_phewas:
		# create raw genotype files for the pheWAS
		create_raw_geno_files(geno_dir = geno_dir, 
								raw_dir = raw_dir)
	elif genotype_phewas:
		print('Reusing raw genotype files')
		raw_dir = reuse_raw_bfiles

	data_dict = {}
	phenof = ''
	if reuse_phenotypes == '' and reuse_fixed_phenotypes == '':

		data_fields = list(field_dict.keys())

		if phenotype_phewas:
			flatten = lambda l: [item for sublist in l for item in sublist]
			data_fields += flatten([pheno_list for pheno_list in list(samples_2_ind.values())[0]])

		# get subsetted phenotype data + covariate data field names
		data_dict, phenof = make_out_of_sample_pheno_file(sample_2_independent_dict = samples_2_ind, 
												pheno_file = ukbiobank_phenotype_file, 
												data_fields = data_fields, 
												cov = covariates, 
												keep = keep,
												out_dir = pheno_dir,
												data_fields_dir = data_fields_dir,
												save_phenos = save_phenotypes,
												pickle_intermediates = pickle_intermediates,
												pickle_protocol = pickle_protocol,
												phenotype_phewas = phenotype_phewas)	

	elif reuse_phenotypes != '' and reuse_fixed_phenotypes == '':
		# read in phenotype data
		# all the files in the reuse_phenotypes directory are assumed to be phenotype files
		for file in os.listdir(reuse_phenotypes):
			
			if '_covariates_fixed' in file:
				continue

			print('Reusing phenotype file: ' + file)

			if file.endswith('.pkl'):
				f = pd.read_pickle(reuse_phenotypes + '/' + file)
			else:
				f = pd.read_csv(reuse_phenotypes + '/' + file, sep = '\t', low_memory = False)

			for sfile, (_, output_prefix) in samples_2_ind.items():

				if file.startswith(output_prefix + '_'):
					data_dict[sfile] = f
					break
	
		phenof = pd.read_csv(old_data_fields_dir + '/' + 'ukb_data_fields.txt', sep = '\t', usecols = [0], dtype = str).iloc[:,0].values.tolist()

	pheWAS_ready_phenos_dict = {}
	if reuse_fixed_phenotypes == '':

		if covariates is not None:

			if covariate_name_map is not None:
				covariate_name_map = pd.read_csv(covariate_name_map, sep = '\t')
				covariate_name_map = covariate_name_map.set_index('Data_field')['Name'].to_dict()
			else:
				covariate_name_map = dict_UKB_fields_to_names
			# fix the covariate data (one hot encoding, calculating age, etc)
			pheWAS_ready_phenos_dict = fix_covariate_data(pheno_dict = data_dict, 
													phenotype_fields = phenof, 
													sample_2_independent_dict = samples_2_ind,
													out_dir = pheno_dir, 
													data_fields_dir = data_fields_dir,
													pickle_intermediates = pickle_intermediates,
													pickle_protocol = pickle_protocol,
													covariate_name_map = covariate_name_map)
		else:
			pheWAS_ready_phenos_dict = data_dict
	else:
		# read in the fixed phenotype data
		for file in os.listdir(reuse_fixed_phenotypes):

			if '_covariates_fixed' not in file:
				continue

			print('Reusing fixed phenotype file: ' + file)

			if file.endswith('.pkl'):
				f = pd.read_pickle(reuse_fixed_phenotypes + '/' + file)
			else:
				f = pd.read_csv(reuse_fixed_phenotypes + '/' + file, sep = '\t', low_memory = False)

			for sfile, (_, output_prefix) in samples_2_ind.items():

				if file.startswith(output_prefix + '_'):
					pheWAS_ready_phenos_dict[sfile] = f
					break

	if genotype_phewas:
		# run the pheWAS across all the variants 
		run_genotype_pheWAS(sample_2_independent_dict = samples_2_ind, 
							pheno_dict = pheWAS_ready_phenos_dict, 
							out_dir = pheWAS_dir, 
							raw_dir = raw_dir, 
							data_fields_dir = data_fields_dir,
							reg = regularization,
							thresh = sample_threshold,
							threads = threads,
							add_covariates = covariates,
							ethnicity_groups = ethnicity_groups)
	else:
		# run the pheWAS across the all the phenotypes
		run_phenotype_pheWAS(sample_2_independent_dict = samples_2_ind,
							pheno_dict = pheWAS_ready_phenos_dict,
							out_dir = pheWAS_dir,
							data_fields_dir = data_fields_dir,
							reg = regularization,
							thresh = sample_threshold,
							threads = threads,
							add_covariates = covariates, 
							ethnicity_groups = ethnicity_groups)


	# create the final output file
	top_results = save_top_variants(out_dir = pheWAS_dir, 
					field_dict = field_dict,
					correction = correction,
					alpha = alpha,
					gene_file = gene_file,
					sample_2_independent_dict = samples_2_ind,
					downstream=downstream,
					upstream=upstream)

	if annotate:
		functional_annotation_dir = dir_path + '/functional_annotation'
		
		print('deleting and recreating old functional_annotation_dir (if it exists)')
		shutil.rmtree(functional_annotation_dir, ignore_errors=True)
		os.makedirs(functional_annotation_dir)

		annotate_function(top_results = top_results, 
							out_dir = functional_annotation_dir)

	print('PheWAS finished, find all results and intermediate files in: ' + directory_name)


if __name__ == '__main__':
	main()