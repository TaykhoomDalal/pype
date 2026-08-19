import pandas as pd
import pytest

from pype.mr import analyze, load_config


def test_dataframe_api_runs_multiple_mr_methods():
	exposure = pd.DataFrame({
		"rsID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
		"CHR": [1, 1, 2, 2, 3],
		"Effect_Allele": ["A", "C", "G", "T", "A"],
		"Non_Effect_Allele": ["G", "T", "A", "C", "C"],
		"BETA": [0.1, 0.2, 0.15, 0.25, 0.3],
		"P": [1e-8] * 5,
		"SE": [0.01] * 5,
	})
	outcome = exposure.copy()
	outcome["BETA"] = [0.05, 0.11, 0.07, 0.13, 0.16]
	outcome["SE"] = [0.04] * 5

	results, diagnostics = analyze(
		exposure,
		outcome,
		exposure_name="trait_a",
		outcome_name="trait_b",
		methods=("ivw", "egger"),
		seed=1,
	)

	assert results["MR_Method"].tolist() == ["Inverse Variance Weighted", "Egger"]
	assert results["Number_SNPs"].tolist() == [5, 5]
	assert diagnostics == {}

	single_result, _ = analyze(exposure, outcome, methods="ivw")
	assert single_result["MR_Method"].tolist() == ["Inverse Variance Weighted"]

	with pytest.raises(ValueError, match="must differ"):
		analyze(exposure, outcome, exposure_name="same", outcome_name="same")


def test_dataframe_api_runs_mr_presso_diagnostics():
	exposure = pd.DataFrame({
		"rsID": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6"],
		"CHR": [1, 1, 2, 2, 3, 3],
		"Effect_Allele": ["A", "C", "G", "T", "A", "C"],
		"Non_Effect_Allele": ["G", "T", "A", "C", "C", "A"],
		"BETA": [0.10, 0.20, 0.15, 0.25, 0.30, 0.18],
		"P": [1e-8] * 6,
		"SE": [0.01] * 6,
	})
	outcome = exposure.copy()
	outcome["BETA"] = [0.05, 0.11, 0.07, 0.13, 0.16, 0.08]
	outcome["SE"] = [0.04] * 6
	config = load_config()
	config["PRESSO"]["nbDist"] = 50

	results, diagnostics = analyze(
		exposure,
		outcome,
		exposure_name="Trait A",
		outcome_name="Body mass index",
		methods=("presso",),
		config=config,
		seed=0,
	)

	assert results["MR_Method"].str.startswith("MR PRESSO").all()
	assert "Global Test" in diagnostics["PRESSO"]
