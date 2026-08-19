import numpy as np
import pandas as pd
import statsmodels.api as sm


RESULT_COLUMNS = [
	"Data_Field",
	"Independent_Var",
	"Samples",
	"-log(p)",
	"p-val",
	"beta",
	"std_error",
]


def run_regression(predictor, outcome, covariate_data=None, covariates=None):
	"""Run one complete-case linear regression."""
	covariates = list(covariates or [])
	data = pd.DataFrame({
		"predictor": predictor.iloc[:, 0],
		"outcome": outcome,
	})

	# Use safe internal names so user supplied column names cannot break the model matrix.
	if covariates:
		covariate_frame = covariate_data[covariates].copy()
		covariate_frame.columns = [f"covariate_{index}" for index in range(len(covariates))]
		data = pd.concat([data, covariate_frame], axis=1)

	data = data.apply(pd.to_numeric, errors="coerce").dropna()
	sample_count = len(data)
	parameter_count = 2 + len(covariates)

	if sample_count <= parameter_count:
		return [np.nan, np.nan, np.nan, np.nan, sample_count]
	if data["predictor"].nunique() < 2 or data["outcome"].nunique() < 2:
		return [np.nan, np.nan, np.nan, np.nan, sample_count]

	design = sm.add_constant(data.drop(columns="outcome"), has_constant="add")
	model = sm.OLS(data["outcome"], design).fit()
	pvalue = model.pvalues["predictor"]
	beta = model.params["predictor"]
	standard_error = model.bse["predictor"]

	return [-np.log10(pvalue), pvalue, beta, standard_error, int(model.nobs)]


def run_associations(phenotypes, predictors, outcomes, covariates=None, min_samples=1000):
	"""Run each predictor against each requested outcome."""
	if min_samples < 1:
		raise ValueError("min_samples must be at least 1.")
	if not phenotypes.index.equals(predictors.index):
		raise ValueError("Phenotype and predictor rows must use the same index.")

	covariates = list(covariates or [])
	missing_covariates = [column for column in covariates if column not in phenotypes]
	if missing_covariates:
		raise ValueError("Missing covariate columns: " + ", ".join(missing_covariates))

	missing_outcomes = [name for name in outcomes if name not in phenotypes]
	if missing_outcomes:
		raise ValueError("Missing outcome columns: " + ", ".join(missing_outcomes))
	outcome_names = [name for name in outcomes if name not in covariates]
	covariate_data = phenotypes[covariates] if covariates else None
	rows = []

	for predictor_name in predictors.columns:
		predictor = predictors[[predictor_name]]

		for outcome_name in outcome_names:
			log_pvalue, pvalue, beta, standard_error, sample_count = run_regression(
				predictor,
				phenotypes[outcome_name],
				covariate_data,
				covariates,
			)

			if sample_count < min_samples or not np.isfinite(pvalue):
				continue

			rows.append([
				outcome_name,
				predictor_name,
				f"{sample_count}/{len(predictors)}",
				log_pvalue,
				pvalue,
				beta,
				standard_error,
			])

	return pd.DataFrame(rows, columns=RESULT_COLUMNS).sort_values("p-val").reset_index(drop=True)


def phenome_wide_association(phenotypes, predictors, outcomes=None, covariates=None, min_samples=1000):
	"""Run a PheWAS from aligned phenotype and predictor dataframes."""
	covariates = list(covariates or [])
	if outcomes is None:
		outcomes = [
			column for column in phenotypes.columns if column not in covariates
		]
	else:
		outcomes = list(outcomes)
	return run_associations(
		phenotypes,
		predictors,
		outcomes,
		covariates=covariates,
		min_samples=min_samples,
	)
