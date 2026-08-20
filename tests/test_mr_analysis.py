import inspect

import pandas as pd
import pytest

from pype.mr import analyze, load_config
from pype.mr.methods import METHOD_REGISTRY


def variant_data(variant_count):
	return pd.DataFrame({
		"variant_id": [
			f"rs{index}" for index in range(1, variant_count + 1)
		],
		"chromosome": [1, 1, 2, 2, 3, 3][:variant_count],
		"effect_allele": ["A", "C", "G", "T", "A", "C"][:variant_count],
		"other_allele": ["G", "T", "A", "C", "C", "A"][:variant_count],
		"beta": [0.10, 0.20, 0.15, 0.25, 0.30, 0.18][:variant_count],
		"pvalue": [1e-8] * variant_count,
		"standard_error": [0.01] * variant_count,
	})


def test_dataframe_api_runs_multiple_mr_methods():
	exposure = variant_data(5)
	outcome = exposure.copy()
	outcome["beta"] = [0.05, 0.11, 0.07, 0.13, 0.16]
	outcome["standard_error"] = [0.04] * 5

	results, diagnostics = analyze(
		exposure,
		outcome,
		exposure_name="trait_a",
		outcome_name="trait_b",
		methods=("ivw", "egger"),
		seed=1,
	)

	assert results["method"].tolist() == [
		"Inverse variance weighted",
		"MR-Egger",
	]
	assert results["variant_count"].tolist() == [5, 5]
	assert results["exposure"].tolist() == ["trait_a", "trait_a"]
	assert results["outcome"].tolist() == ["trait_b", "trait_b"]
	assert "Inverse variance weighted" in diagnostics
	assert "MR-Egger" in diagnostics

	single_result, _ = analyze(exposure, outcome, methods="ivw")
	assert single_result["method"].tolist() == [
		"Inverse variance weighted"
	]

	with pytest.raises(ValueError, match="must differ"):
		analyze(
			exposure,
			outcome,
			exposure_name="same",
			outcome_name="same",
		)


def test_dataframe_api_runs_mr_presso_diagnostics():
	exposure = variant_data(6)
	outcome = exposure.copy()
	outcome["beta"] = [0.05, 0.11, 0.07, 0.13, 0.16, 0.08]
	outcome["standard_error"] = [0.04] * 6
	config = load_config()
	config["presso"]["simulation_count"] = 50

	with pytest.warns(UserWarning):
		results, diagnostics = analyze(
			exposure,
			outcome,
			exposure_name="Trait A",
			outcome_name="Body mass index",
			methods=("presso",),
			config=config,
			seed=0,
		)

	assert results["method"].str.startswith("MR-PRESSO").all()
	assert "global_test" in diagnostics["presso"]


def test_default_parameter_names_match_method_signatures():
	config = load_config()
	for method_key, parameters in config.items():
		method = METHOD_REGISTRY[method_key][1]
		signature = inspect.signature(method)
		assert set(parameters) <= set(signature.parameters)
