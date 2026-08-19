import numpy as np
import pandas as pd
import pytest

from pype.phewas import phenome_wide_association, run_associations, run_regression


def _regression_data():
	rng = np.random.default_rng(3)
	independent = pd.DataFrame({"variant": np.linspace(-2, 2, 12)})
	covariate = pd.Series([0, 1] * 6, name="covariate")
	phenotype = 1.75 * independent["variant"] + 0.8 * covariate + rng.normal(0, 0.03, 12)

	independent.loc[1, "variant"] = np.nan
	covariate.loc[4] = np.nan
	phenotype.loc[9] = np.nan
	return independent, phenotype, covariate


def test_regression_uses_complete_cases_and_returns_expected_effect():
	independent, phenotype, covariate = _regression_data()

	result = run_regression(
		independent,
		phenotype,
		pd.DataFrame({"covariate": covariate}),
		["covariate"],
	)

	assert result[4] == 9
	assert result[2] == pytest.approx(1.75, abs=0.05)
	assert result[1] < 1e-6


def test_phewas_reports_the_actual_complete_case_sample_count():
	independent, phenotype, covariate = _regression_data()
	phenotypes = pd.DataFrame({"outcome": phenotype, "covariate": covariate})

	results = run_associations(
		phenotypes,
		independent,
		["outcome"],
		covariates=["covariate"],
		min_samples=9,
	)

	assert results.loc[0, "Samples"] == "9/12"
	assert results.loc[0, "beta"] == pytest.approx(1.75, abs=0.05)


def test_public_phewas_api_uses_clear_argument_names():
	independent, phenotype, covariate = _regression_data()
	phenotypes = pd.DataFrame({"outcome": phenotype, "covariate": covariate})

	results = phenome_wide_association(
		phenotypes,
		independent,
		outcomes=["outcome"],
		covariates=["covariate"],
		min_samples=9,
	)

	assert results.loc[0, "Samples"] == "9/12"
