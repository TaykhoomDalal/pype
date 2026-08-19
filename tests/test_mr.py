import numpy as np
import pandas as pd
import pytest

from pype.mr.harmonization import harmonize
from pype.mr.methods import (
	density,
	run_mr_egger,
	run_mr_ivw,
	weighted_median,
	weighted_median_bootstrap,
)


def test_ivw_matches_twosamplemr_under_dispersion_standard_error():
	exposure = np.array([0.10, 0.20, 0.15, 0.25, 0.30])
	outcome = np.array([0.051, 0.099, 0.076, 0.124, 0.151])
	outcome_se = np.full(5, 0.04)
	data = pd.DataFrame({
		"beta_exp": exposure,
		"beta_out": outcome,
		"se_exp": np.full(5, 0.01),
		"se_out": outcome_se,
	})

	result = run_mr_ivw(data, "beta_exp", "beta_out", "se_exp", "se_out")
	expected_standard_error = np.sqrt(1 / np.sum(exposure**2 / outcome_se**2))

	assert result[2] == pytest.approx(expected_standard_error)
	assert result[0] > 0


def test_egger_intercept_matches_twosamplemr_reference_result():
	data = pd.DataFrame({
		"beta_exp": [0.10, 0.20, 0.15, 0.25, 0.30, 0.18],
		"beta_out": [0.05, 0.11, 0.07, 0.13, 0.16, 0.08],
		"se_exp": [0.01] * 6,
		"se_out": [0.04] * 6,
	})

	result = run_mr_egger(data, "beta_exp", "beta_out", "se_exp", "se_out")

	assert result[4] == pytest.approx(-0.0117894736842105)
	assert result[5] == pytest.approx(0.0520526, rel=1e-6)
	assert result[3] == pytest.approx(0.8319227, rel=1e-6)


def test_weighted_median_bootstrap_uses_sample_standard_deviation():
	beta_exp = np.array([0.1, 0.2, 0.3])
	beta_out = np.array([0.05, 0.11, 0.16])
	se_exp = np.array([0.01, 0.01, 0.02])
	se_out = np.array([0.02, 0.02, 0.03])
	weights = np.array([1.0, 2.0, 3.0])

	np.random.seed(9)
	medians = []
	for _ in range(20):
		exp_boot = np.random.normal(beta_exp, se_exp)
		out_boot = np.random.normal(beta_out, se_out)
		medians.append(weighted_median(out_boot / exp_boot, weights))
	expected = np.std(medians, ddof=1)

	np.random.seed(9)
	actual = weighted_median_bootstrap(
		beta_exp,
		beta_out,
		se_exp,
		se_out,
		weights,
		20,
	)

	assert actual == pytest.approx(expected)


def test_density_keeps_weights_aligned_when_ratio_is_infinite():
	result = density(
		np.array([0.3, np.inf, 0.7, 0.9]),
		0.1,
		np.array([0.1, 0.2, 0.3, 0.4]),
	)

	assert len(result["x"]) == 512
	assert np.isfinite(result["y"]).all()


def test_harmonize_flips_reversed_alleles_and_drops_incompatible_rows():
	exposure = pd.DataFrame({
		"rsID": ["rs1", "rs2", "rs3", "rs4"],
		"CHR": [1, 1, 1, 1],
		"Effect_Allele": ["A", "C", "A", "A"],
		"Non_Effect_Allele": ["G", "T", "C", "T"],
		"BETA": [0.1, 0.2, 0.3, 0.4],
		"P": [0.01, 0.02, 0.03, 0.04],
		"SE": [0.01, 0.02, 0.03, 0.04],
	})
	outcome = pd.DataFrame({
		"rsID": ["rs1", "rs2", "rs3", "rs4"],
		"CHR": [1, 1, 1, 1],
		"Effect_Allele": ["G", "G", "A", "A"],
		"Non_Effect_Allele": ["A", "A", "G", "T"],
		"BETA": [0.4, 0.5, 0.6, 0.7],
		"P": [0.04, 0.05, 0.06, 0.07],
		"SE": [0.04, 0.05, 0.06, 0.07],
	})

	with pytest.warns(UserWarning):
		result = harmonize(exposure, outcome, "_exp", "_out")

	assert result["rsID"].tolist() == ["rs1", "rs2"]
	assert result.loc[result["rsID"] == "rs1", "BETA_out"].item() == pytest.approx(-0.4)
	assert result.loc[result["rsID"] == "rs2", "BETA_out"].item() == pytest.approx(0.5)
