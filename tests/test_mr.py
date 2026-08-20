import numpy as np
import pandas as pd
import pytest

from pype.mr.harmonization import harmonize, standardize_variants
from pype.mr.methods import (
	_weighted_median,
	_weighted_median_bootstrap,
	inverse_variance_weighted,
	kernel_density,
	mr_egger,
)


def harmonized_data():
	return pd.DataFrame({
		"exposure_beta": [0.10, 0.20, 0.15, 0.25, 0.30, 0.18],
		"outcome_beta": [0.05, 0.11, 0.07, 0.13, 0.16, 0.08],
		"exposure_standard_error": [0.01] * 6,
		"outcome_standard_error": [0.04] * 6,
	})


def test_ivw_matches_twosamplemr_under_dispersion_standard_error():
	exposure_beta = np.array([0.10, 0.20, 0.15, 0.25, 0.30])
	outcome_standard_error = np.full(5, 0.04)
	data = pd.DataFrame({
		"exposure_beta": exposure_beta,
		"outcome_beta": [0.051, 0.099, 0.076, 0.124, 0.151],
		"exposure_standard_error": [0.01] * 5,
		"outcome_standard_error": outcome_standard_error,
	})

	result = inverse_variance_weighted(data)
	expected_standard_error = np.sqrt(
		1 / np.sum(exposure_beta**2 / outcome_standard_error**2)
	)

	assert result["standard_error"] == pytest.approx(
		expected_standard_error
	)
	assert result["pvalue"] > 0


def test_egger_intercept_matches_twosamplemr_reference_result():
	result = mr_egger(harmonized_data())

	assert result["intercept"] == pytest.approx(-0.0117894736842105)
	assert result["intercept_standard_error"] == pytest.approx(
		0.0520526,
		rel=1e-6,
	)
	assert result["intercept_pvalue"] == pytest.approx(
		0.8319227,
		rel=1e-6,
	)


def test_weighted_median_bootstrap_uses_sample_standard_deviation():
	exposure_beta = np.array([0.1, 0.2, 0.3])
	outcome_beta = np.array([0.05, 0.11, 0.16])
	exposure_standard_error = np.array([0.01, 0.01, 0.02])
	outcome_standard_error = np.array([0.02, 0.02, 0.03])
	weights = np.array([1.0, 2.0, 3.0])

	np.random.seed(9)
	bootstrap_estimates = []
	for _ in range(20):
		sampled_exposure = np.random.normal(
			exposure_beta,
			exposure_standard_error,
		)
		sampled_outcome = np.random.normal(
			outcome_beta,
			outcome_standard_error,
		)
		bootstrap_estimates.append(
			_weighted_median(
				sampled_outcome / sampled_exposure,
				weights,
			)
		)
	expected = np.std(bootstrap_estimates, ddof=1)

	np.random.seed(9)
	actual = _weighted_median_bootstrap(
		exposure_beta,
		outcome_beta,
		exposure_standard_error,
		outcome_standard_error,
		weights,
		20,
	)
	assert actual == pytest.approx(expected)


def test_density_keeps_weights_aligned_when_ratio_is_infinite():
	result = kernel_density(
		np.array([0.3, np.inf, 0.7, 0.9]),
		0.1,
		np.array([0.1, 0.2, 0.3, 0.4]),
	)

	assert len(result["x"]) == 512
	assert np.isfinite(result["y"]).all()


def test_harmonize_flips_reversed_alleles_and_drops_incompatible_rows():
	exposure = pd.DataFrame({
		"variant_id": ["rs1", "rs2", "rs3", "rs4"],
		"chromosome": [1, 1, 1, 1],
		"effect_allele": ["A", "C", "A", "A"],
		"other_allele": ["G", "T", "C", "T"],
		"beta": [0.1, 0.2, 0.3, 0.4],
		"pvalue": [0.01, 0.02, 0.03, 0.04],
		"standard_error": [0.01, 0.02, 0.03, 0.04],
	})
	outcome = pd.DataFrame({
		"variant_id": ["rs1", "rs2", "rs3", "rs4"],
		"chromosome": [1, 1, 1, 1],
		"effect_allele": ["G", "G", "A", "A"],
		"other_allele": ["A", "A", "G", "T"],
		"beta": [0.4, 0.5, 0.6, 0.7],
		"pvalue": [0.04, 0.05, 0.06, 0.07],
		"standard_error": [0.04, 0.05, 0.06, 0.07],
	})

	with pytest.warns(UserWarning):
		result = harmonize(exposure, outcome)

	assert result["variant_id"].tolist() == ["rs1", "rs2"]
	assert result.loc[
		result["variant_id"] == "rs1",
		"outcome_beta",
	].item() == pytest.approx(-0.4)
	assert result.loc[
		result["variant_id"] == "rs2",
		"outcome_beta",
	].item() == pytest.approx(0.5)


def test_variant_values_are_normalized():
	data = pd.DataFrame({
		"variant_id": ["rs1"],
		"chromosome": ["chr1"],
		"effect_allele": ["a"],
		"other_allele": ["g"],
		"beta": [0.1],
		"standard_error": [0.02],
		"pvalue": [0.01],
		"sample_size": [1000],
	})

	result = standardize_variants(data)

	assert {
		"variant_id",
		"chromosome",
		"effect_allele",
		"other_allele",
		"beta",
		"standard_error",
		"pvalue",
		"sample_size",
	} <= set(result.columns)
