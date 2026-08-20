import numpy as np
import pandas as pd
import statsmodels.api as sm


RESULT_COLUMNS = [
	"outcome",
	"predictor",
	"sample_count",
	"total_sample_count",
	"negative_log10_pvalue",
	"pvalue",
	"beta",
	"standard_error",
]


def _fit_linear_regression(predictor, outcome, covariate_data=None, covariates=None):
	covariates = list(covariates or [])
	data = pd.DataFrame({
		"predictor": predictor,
		"outcome": outcome,
	})

	# Safe internal names allow arbitrary user-facing covariate names.
	if covariates:
		covariate_frame = covariate_data[covariates].copy()
		covariate_frame.columns = [
			f"covariate_{index}" for index in range(len(covariates))
		]
		data = pd.concat([data, covariate_frame], axis=1)

	data = data.apply(pd.to_numeric, errors="coerce").dropna()
	sample_count = len(data)
	parameter_count = 2 + len(covariates)

	if sample_count <= parameter_count:
		return None
	if data["predictor"].nunique() < 2 or data["outcome"].nunique() < 2:
		return None

	design_matrix = sm.add_constant(
		data.drop(columns="outcome"),
		has_constant="add",
	)
	model = sm.OLS(data["outcome"], design_matrix).fit()
	pvalue = float(model.pvalues["predictor"])

	return {
		"sample_count": int(model.nobs),
		"negative_log10_pvalue": float(-np.log10(pvalue)),
		"pvalue": pvalue,
		"beta": float(model.params["predictor"]),
		"standard_error": float(model.bse["predictor"]),
	}


def phenome_wide_association(
	phenotypes,
	predictors,
	outcomes=None,
	covariates=None,
	min_sample_count=1000,
):
	"""Run each predictor against each selected outcome."""
	if min_sample_count < 1:
		raise ValueError("min_sample_count must be at least 1.")
	if not phenotypes.index.equals(predictors.index):
		raise ValueError("Phenotype and predictor rows must use the same index.")

	covariates = list(covariates or [])
	missing_covariates = [
		column for column in covariates if column not in phenotypes
	]
	if missing_covariates:
		raise ValueError(
			"Missing covariate columns: " + ", ".join(missing_covariates)
		)

	if outcomes is None:
		outcomes = [
			column for column in phenotypes.columns if column not in covariates
		]
	else:
		outcomes = list(outcomes)

	missing_outcomes = [column for column in outcomes if column not in phenotypes]
	if missing_outcomes:
		raise ValueError(
			"Missing outcome columns: " + ", ".join(missing_outcomes)
		)

	outcomes = [column for column in outcomes if column not in covariates]
	covariate_data = phenotypes[covariates] if covariates else None
	total_sample_count = len(phenotypes)
	rows = []

	for predictor_name in predictors.columns:
		for outcome_name in outcomes:
			result = _fit_linear_regression(
				predictors[predictor_name],
				phenotypes[outcome_name],
				covariate_data,
				covariates,
			)
			if (
				result is None
				or result["sample_count"] < min_sample_count
				or not np.isfinite(
					[
						result["pvalue"],
						result["beta"],
						result["standard_error"],
					]
				).all()
			):
				continue

			rows.append({
				"outcome": outcome_name,
				"predictor": predictor_name,
				"total_sample_count": total_sample_count,
				**result,
			})

	return (
		pd.DataFrame(rows, columns=RESULT_COLUMNS)
		.sort_values("pvalue")
		.reset_index(drop=True)
	)
