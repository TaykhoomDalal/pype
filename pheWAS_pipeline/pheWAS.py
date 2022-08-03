import statsmodels.formula.api as smf
import statsmodels.genmod.generalized_linear_model as sm_glm
from tqdm import tqdm
import numpy as np
import pandas as pd
import math 
import re
import traceback

# Code in this file was inspired from the pyPheWAS package, specifically from the pyPhewasCorev2.py file.

out_cols = ['Data_Field', 'rsID', 'Samples', '\"-log(p)\"','p-val', 'beta', 'std_error']

def run_regresssion(variant, phenotype, covariate_data = None, covariates = None, regularization = False) -> list:

	# rename variant column to genotype
	variant = variant.rename(columns={variant.columns[0]: 'genotype'})

	# add genotype and phenotype to dataframe
	data = variant.copy(deep = True)
	data['phenotype'] = phenotype

	# append '+' to covariates (if there are any) -> makes definition of 'f' more elegant
	if covariates != None:
		for cov in covariates:
			data[cov] = covariate_data[cov]

		# remove all non-alphanumeric characters from covariate names
		covariates = [re.sub(r'[^a-zA-Z0-9]+','', cov) for cov in covariates]
		covariates = " ".join(covariates)
	
	# if there are any missing values in the phenotype vector, drop the row from the data
	data.dropna(axis = 0, how = 'any', subset = ['phenotype'], inplace = True)

	# create formula for regression
	covariates = covariates.replace(" ", " + ")
	f = 'phenotype ~ genotype + ' + covariates

	# replace all non-alphanumerics with empty space to deal with statsmodels issues
	data.columns = [re.sub(r'[^a-zA-Z0-9]+','', x) for x in data.columns]

	try:
		# need to spend some time figuring out w
		results = smf.ols(formula = f, data = data).fit()

		# get results
		p = results.pvalues[1]
		beta = results.params[1]
		stderr = results.bse[1]
		reg_result = [-np.log10(p), p, beta, stderr]  # collect results

	except Exception as e:
		traceback.print_exc()
		print('Error when computing regression for phenotype %s' %(phenotype.name))
		exit()
		# reg_result = [np.nan, np.nan, np.nan, np.nan, np.nan]

	return reg_result

def run_phewas(phenotypes, genotypes, non_cov_pheno_list, reg, covariates = None, thresh=1000) -> pd.DataFrame:

	num_pheno = len(non_cov_pheno_list)
	num_variants = genotypes.shape[1]
	total_samples = '/' + str(genotypes.shape[0])

	if num_pheno == 0:
		print('No phenotypes to regress')
		return None

	regressions = pd.DataFrame(columns=out_cols)

	index = 0

	# grab the columns of the phenotype matrix corresponding to the covariate phenotypes
	covariate_data = phenotypes[covariates]

	# drop all the covariates from the phenotype data
	phenotypes.drop(covariates, axis=1, inplace=True)

	# keep all columns with more than 80% non-nan values
	covariate_data = covariate_data.dropna(thresh = len(covariate_data)*0.2, axis = 1)

	# update the covariates column names with the new names
	covariates = covariate_data.columns.tolist()

	# iterate over all phenotypes, dropping any that are of type object, that have fewer than thresh samples, or that have fewer than the # of covariates
	phenos_to_drop = []
	for i in range(num_pheno):
		phenotype_i = phenotypes[non_cov_pheno_list[i]]
		
		non_nan_samples_num = phenotype_i.shape[0] - phenotype_i.isnull().sum()
		
		# if we have a lot fewer observations (samples) than variables (covariates + our genotype), 
		# then there is no unique solution to OLS (unless we use regularization - user parameter)
		if non_nan_samples_num < len(covariates) + 1 and reg == False:
			print('Phenotype %s dropped - it has %d samples, less than the number of covariates (%d)' %(phenotype_i.name, non_nan_samples_num, len(covariates) + 1))
			phenos_to_drop.append(phenotype_i.name)
			continue
		
		# if our phenotype is an object and not a number, then we can't do regression
		if phenotype_i.dtype == 'object':
			print('Phenotype %s dropped - it is of type object and cannot be used for regression' %(phenotype_i.name))
			phenos_to_drop.append(phenotype_i.name)
			continue

		# to prevent false positives, only run regressions if more than thresh records have positive values
		# if our phenotype is a number, but has less than 1000 samples, then we can't do regression
		if non_nan_samples_num < thresh:
			print('Phenotype %s dropped - it has %d samples, less than the threshold (%d)' %(phenotype_i.name, non_nan_samples_num, thresh))
			phenos_to_drop.append(phenotype_i.name)
			continue
	
	# drop the phenotypes that we don't want to regress on and update the list of phenotypes
	phenotypes.drop(phenos_to_drop, axis=1, inplace=True)
	non_cov_pheno_list = [x for x in non_cov_pheno_list if x not in phenos_to_drop]
	
	# update the number of phenotypes
	num_pheno = len(non_cov_pheno_list)

	for v in range(num_variants):
		variant_i = genotypes.iloc[:, [v]]

		# its possible that there are NaN values in the variant data, so we need to drop them
		variant_non_nan_samples_num = variant_i.shape[0] - variant_i.isnull().sum()[0]

		# if we have a lot fewer observations (samples) than variables (covariates + genotype),
		# then there is no unique solution to OLS (unless we use regularization - user parameter)
		if variant_non_nan_samples_num < len(covariates) + 1 and reg == False:
			print('Variant %s dropped - it has %d samples, less than the number of covariates (%d)' %(variant_i.columns[0], variant_non_nan_samples_num, len(covariates) + 1))
			continue
		
		# to prevent false positives, only run regressions if more than thresh records have positive values
		# if our phenotype is a number, but has less than 1000 samples, then we can't do regression
		if variant_non_nan_samples_num < thresh:
			print('Variant %s dropped - it has %d samples, less than the threshold (%d)' %(variant_i.columns[0], variant_non_nan_samples_num, thresh))
			continue

		# drop the rows with NaN values
		variant_i = variant_i.dropna(axis = 0)

		for p in tqdm(range(num_pheno), desc="Running Regressions for Variant %s" %(variant_i.columns[0])):
			phenotype_i = phenotypes[non_cov_pheno_list[p]]

			pheno_non_nan_samples_num = phenotype_i.shape[0] - phenotype_i.isnull().sum()

			# run the regression
			stat_info = run_regresssion(variant_i, phenotype_i, covariate_data, covariates, reg)

			# whichever number is smaller, that is the number of samples we used in the regression
			used_samples = min(pheno_non_nan_samples_num, variant_non_nan_samples_num)

			# save regression data
			info = [non_cov_pheno_list[p], genotypes.columns[v], str(used_samples) + total_samples] + stat_info 

			regressions.loc[index] = info
			index +=1

	return regressions.dropna(subset=['p-val']).sort_values(by='p-val')  # sort by significance