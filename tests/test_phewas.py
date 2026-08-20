import numpy as np
import pandas as pd
import pytest

from pype.phewas import _fit_linear_regression, phenome_wide_association


def regression_data():
	random_generator = np.random.default_rng(3)
	predictor = pd.Series(
		np.linspace(-2, 2, 12),
		name="variant",
	)
	covariate = pd.Series([0, 1] * 6, name="covariate")
	outcome = (
		1.75 * predictor
		+ 0.8 * covariate
		+ random_generator.normal(0, 0.03, 12)
	)

	predictor.loc[1] = np.nan
	covariate.loc[4] = np.nan
	outcome.loc[9] = np.nan
	return predictor, outcome, covariate


def test_regression_uses_complete_cases_and_returns_expected_effect():
	predictor, outcome, covariate = regression_data()
	result = _fit_linear_regression(
		predictor,
		outcome,
		pd.DataFrame({"covariate": covariate}),
		["covariate"],
	)

	assert result["sample_count"] == 9
	assert result["beta"] == pytest.approx(1.75, abs=0.05)
	assert result["pvalue"] < 1e-6


def test_phewas_returns_consistent_result_columns():
	predictor, outcome, covariate = regression_data()
	phenotypes = pd.DataFrame({
		"outcome": outcome,
		"covariate": covariate,
	})
	predictors = pd.DataFrame({"variant": predictor})

	results = phenome_wide_association(
		phenotypes,
		predictors,
		outcomes=["outcome"],
		covariates=["covariate"],
		min_sample_count=9,
	)

	assert results.loc[0, "sample_count"] == 9
	assert results.loc[0, "total_sample_count"] == 12
	assert results.loc[0, "predictor"] == "variant"
	assert results.loc[0, "outcome"] == "outcome"
	assert results.loc[0, "beta"] == pytest.approx(1.75, abs=0.05)


def test_phewas_omits_degenerate_regressions():
	phenotypes = pd.DataFrame({"outcome": [1.0, 2.0, 3.0, 4.0]})
	predictors = pd.DataFrame({"constant": [1.0, 1.0, 1.0, 1.0]})

	results = phenome_wide_association(
		phenotypes,
		predictors,
		min_sample_count=4,
	)

	assert results.empty
