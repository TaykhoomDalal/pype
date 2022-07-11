import statsmodels.api as smf
from tqdm import tqdm
import numpy as np
import pandas as pd
import math 
import re

# Code in this file was adapted from the pyPheWAS package, specifically from
# the pyPhewasCorev2.py file. I modified the code to run the pheWAS with non ICD-9/10
# codes as well as added functionality that was required for my project.
# The original code can be found at https://github.com/MASILab/pyPheWAS/blob/master/pyPheWAS/pyPhewasCorev2.py

out_cols = ['Data_Field', 'rsID', 'Samples', '\"-log(p)\"','p-val', 'beta', 'Conf-interval beta', 'std_error']

def fit_pheno_model(variant, phenotype, covariate_data = None, covariates = None, reg = False) -> list:
	"""
	Runs a linear regression for a specific phenotype vector

	The returned results list consists of (in order) the -log\ :sub:`10`\ (p-value), p-value, beta, beta's confidence interval,
	and beta's standard error, estimated from the lin model for the ``phenotype`` variable.

	:param variant: variant data for each individual
	:param phenotype: the aggregate phenotype vector
	:param covariate_data: covariate data for each individual
	:param covariates: *[optional]* covariates to include in the regressions separated by '+' (e.g. 'sex+ageAtDx')
	:param reg: *[optional]* regularized maximum likelihood optimization flag (False [default] = not regularized; True = regularized)

	:type variant: pandas DataFrame
	:type phenotype: pandas Series
	:type covariate_data: pandas DataFrame
	:type covariates: list of covariate names
	:type reg: bool

	:returns: regression results
	:rtype: list
	"""
	# rename variant column name by removing all non-alphanumeric characters from variant name
	variant.columns = [re.sub('[^0-9a-zA-Z]+', '', x) for x in variant.columns]

	data = variant.copy()
	data['y'] = phenotype

	# append '+' to covariates (if there are any) -> makes definition of 'f' more elegant
	if covariates != None:
		for cov in covariates:
			data[cov] = covariate_data[cov]

		# remove all non-alphanumeric characters from covariate names
		covariates = [re.sub(r'[^a-zA-Z0-9]+','', cov) for cov in covariates]
		covariates = " ".join(covariates)
	
	# if there are any missing values in the phenotype vector, drop the row from the data
	data.dropna(axis = 0, how = 'any', subset = ['y'], inplace = True)

	# create formula for regression
	predictors = covariates.replace(" ", " + ")
	f = variant.columns.tolist()[0] + ' ~ y + ' + predictors

	# replace all non-alphanumerics with empty space to deal with statsmodels issues
	data.columns = [re.sub(r'[^a-zA-Z0-9]+','', x) for x in data.columns]

	try:
		if reg == False:
			# fit logit without regulatization
			model = smf.OLS.from_formula(formula = f, data = data).fit()
		else:
			# fit with regularization
			variant_name = variant.columns[0]
			model = smf.OLS(data[variant_name], data[[phenotype.name] + covariates]).fit_regularized(method='l1', alpha=0.1)

		# get results
		p = model.pvalues.y
		beta = model.params.y
		conf = model.conf_int()
		conf_int = '[%s,%s]' % (conf[0]['y'], conf[1]['y'])
		stderr = model.bse.y
		reg_result = [-np.log10(p), p, beta, conf_int, stderr]  # collect results

	except Exception as e:
		print('ERROR (%s) computing regression for phenotype %s' %(e, phenotype.name))
		print(data)
		print(phenotype)
		print(covariate_data)
		reg_result = [np.nan, np.nan, np.nan, np.nan, np.nan]

	return reg_result

def run_phewas(phenotypes, genotypes, non_cov_pheno_list, reg, covariates = None, thresh=1000) -> pd.DataFrame:
	"""
	Iterate over all phenotypes, running a linear regression against each variant in the genotype
	
	The returned DataFrame includes the phenotype data field, variant in question, number of samples,
	-log\ :sub:`10`\ (p-value), p-value, beta, beta's confidence interval, standard error

	:param phenotypes: phenotype feature matrix
	:param genotypes: variant data
	:param non_cov_pheno_list: list of non-covariate phenotypes
	:param reg: whether to include regularization in the regression
	:param covariates: *[optional]* covariates to include in the regressions
	:param thresh: *[optional]* threshold for the minimum number of samples in a phenotype
	
	:type phenotypes: pandas DataFrame
	:type genotypes: pandas DataFrame
	:type non_cov_pheno_list: list
	:type reg: bool
	:type covariates: list
	:type thresh: int

	:returns: regression results for each phenotype
	:rtype: pandas DataFrame
	"""

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

	phenos_to_drop = []
	for i in range(num_pheno):
		phenotype_i = phenotypes[non_cov_pheno_list[i]]
		
		non_nan_samples_num = phenotype_i.shape[0] - phenotype_i.isnull().sum()
		
		# if we have a lot fewer observations (samples) than variables (covariates + our phenotype), 
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

		for p in tqdm(range(num_pheno), desc="Running Regressions for Variant %s" %(variant_i.columns[0])):
			phenotype_i = phenotypes[non_cov_pheno_list[p]]
			non_nan_samples_num = phenotype_i.shape[0] - phenotype_i.isnull().sum()
			# run the regression
			stat_info = fit_pheno_model(variant_i, phenotype_i, covariate_data, covariates, reg)

			# save regression data
			info = [non_cov_pheno_list[p], genotypes.columns[v], str(non_nan_samples_num) + total_samples] + stat_info 

			regressions.loc[index] = info
			index +=1

	return regressions.dropna(subset=['p-val']).sort_values(by='p-val')  # sort by significance