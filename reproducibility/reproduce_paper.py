import argparse
import copy
import csv
import gzip
import io
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.formula.api import wls

from pype.mr import load_config
from pype.mr.methods import inverse_variance_weighted, run_methods


LE_GOALLEC_URL = (
	"https://media.springernature.com/original/springer-static/esm/"
	"art%3A10.1038%2Fs41467-022-29525-9/MediaObjects/"
	"41467_2022_29525_MOESM4_ESM.xlsx"
)
NEALE_BASE_URL = (
	"https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/"
	"round2/additive-tsvs"
)
OUTCOME_FILES = {
	"HbA1C": "30750_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
	"BMI": "21001_raw.gwas.imputed_v3.both_sexes.tsv.bgz",
	"Glucose": "30740_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
	"Waist Circumference": "48_raw.gwas.imputed_v3.both_sexes.tsv.bgz",
}
METHOD_KEYS = [
	"ivw",
	"egger",
	"simple_median",
	"weighted_median",
	"penalized_weighted_median",
	"simple_mode",
	"weighted_mode",
	"penalized_mode",
	"simple_mode_nome",
	"weighted_mode_nome",
]
EXPOSURE_VARIANTS = {
	"Abdomen": {
		"phenotype": "AbdAge",
		"variants": ["rs201407787", "rs2216113", "rs932274"],
	},
	"Liver": {
		"phenotype": "Liver Age",
		"variants": [
			"rs552571374",
			"rs201407787",
			"rs3791675",
			"rs13107325",
			"rs12539772",
			"rs11111209",
			"rs77353655",
			"rs76652635",
			"rs45515493",
		],
	},
}
PUBLISHED_TABLE_1 = {
	("Abdomen", "HbA1C"): (0.02652062, 0.02043320, 0.19431561),
	("Abdomen", "BMI"): (0.00608696, 0.11579183, 0.95807602),
	("Abdomen", "Glucose"): (-0.01507537, 0.01321411, 0.25392083),
	("Abdomen", "Waist Circumference"): (0.35556079, 0.10857999, 0.00105795),
	("Liver", "HbA1C"): (0.01529253, 0.04408397, 0.72866982),
	("Liver", "BMI"): (-0.01221840, 0.12306700, 0.92090000),
	("Liver", "Glucose"): (0.00658953, 0.00717388, 0.35833315),
	("Liver", "Waist Circumference"): (0.25086723, 0.15703291, 0.11014415),
}


def _download(url, path):
	path.parent.mkdir(parents=True, exist_ok=True)
	if not path.exists():
		urllib.request.urlretrieve(url, path)
	return path


def _load_exposures(input_directory):
	workbook = _download(LE_GOALLEC_URL, input_directory / "le_goallec_gwas.xlsx")
	data = pd.read_excel(workbook)
	data = data[data["rsID"] == data["IndSigSNP"]]
	exposures = {}

	for name, specification in EXPOSURE_VARIANTS.items():
		rows = data[
			(data["phenotype"] == specification["phenotype"])
			& data["rsID"].isin(specification["variants"])
		].set_index("rsID").loc[specification["variants"]]
		exposures[name] = rows.rename(columns={
			"beta": "exposure_beta",
			"se": "exposure_standard_error",
		})[
			[
				"chr",
				"pos",
				"exposure_beta",
				"exposure_standard_error",
			]
		]

	return exposures


def _extract_outcome(
	outcome_name,
	summary_filename,
	variant_positions,
	input_directory,
):
	path = input_directory / (
		outcome_name.lower().replace(" ", "_") + ".tsv"
	)
	if path.exists():
		result = pd.read_csv(path, sep="\t")
		if "position" not in result:
			result = pd.DataFrame({
				"position": result["variant"].str.split(":").str[1].astype(int),
				"outcome_beta": result["beta"].astype(float),
				"outcome_standard_error": result["se"].astype(float),
			})
		required_columns = {
			"position",
			"outcome_beta",
			"outcome_standard_error",
		}
		missing_columns = required_columns - set(result.columns)
		if missing_columns:
			raise ValueError(
				"Cached outcome data are missing: "
				+ ", ".join(sorted(missing_columns))
			)
		return result

	url = NEALE_BASE_URL + "/" + summary_filename
	selected = []
	with urllib.request.urlopen(url) as response:
		with gzip.GzipFile(fileobj=response) as compressed:
			reader = csv.DictReader(io.TextIOWrapper(compressed), delimiter="\t")
			for row in reader:
				position = int(row["variant"].split(":")[1])
				if position in variant_positions:
					selected.append({
						"position": position,
						"outcome_beta": float(row["beta"]),
						"outcome_standard_error": float(row["se"]),
					})

	result = pd.DataFrame(selected)
	result.to_csv(path, sep="\t", index=False)
	return result


def _analysis_frame(exposure, outcome):
	return (
		exposure.reset_index()
		.merge(outcome, left_on="pos", right_on="position", validate="one_to_one")
		[
			[
				"exposure_beta",
				"outcome_beta",
				"exposure_standard_error",
				"outcome_standard_error",
			]
		]
	)


def _format_results(raw_results, predictor_name, outcome_name):
	rows = []
	for display_name, result in raw_results.items():
		if result is None:
			continue
		if display_name == "MR-PRESSO":
			for estimate in result["mr_results"]:
				rows.append({
					"predictor": predictor_name,
					"outcome": outcome_name,
					**estimate,
				})
		else:
			rows.append({
				"predictor": predictor_name,
				"outcome": outcome_name,
				"method": display_name,
				"pvalue": result["pvalue"],
				"beta": result["beta"],
				"standard_error": result["standard_error"],
				"variant_count": result["variant_count"],
			})
	return rows


def _legacy_ivw(data):
	model = wls(
		"outcome_beta ~ -1 + exposure_beta",
		data=data,
		weights=1 / data["outcome_standard_error"] ** 2,
	).fit()
	beta = float(model.params["exposure_beta"])
	standard_error = float(model.bse["exposure_beta"])
	pvalue = 2 * norm.sf(abs(beta / standard_error))
	return beta, standard_error, pvalue


def reproduce(
	output_directory,
	bootstrap_iterations,
	simulation_count,
	seed,
):
	output_directory.mkdir(parents=True, exist_ok=True)
	input_directory = output_directory / "inputs"
	exposures = _load_exposures(input_directory)
	variant_positions = {
		int(position)
		for exposure in exposures.values()
		for position in exposure["pos"]
	}
	outcomes = {
		outcome_name: _extract_outcome(
			outcome_name,
			summary_filename,
			variant_positions,
			input_directory,
		)
		for outcome_name, summary_filename in OUTCOME_FILES.items()
	}

	config = load_config()
	for settings in config.values():
		if "bootstrap_iterations" in settings:
			settings["bootstrap_iterations"] = bootstrap_iterations
	config["presso"]["simulation_count"] = simulation_count

	all_results = []
	comparisons = []
	for predictor, exposure in exposures.items():
		for outcome_name, outcome in outcomes.items():
			data = _analysis_frame(exposure, outcome)
			np.random.seed(seed)
			selected_methods = METHOD_KEYS + (
				["presso"] if len(data) >= 4 else []
			)
			raw_results = run_methods(
				data,
				selected_methods,
				copy.deepcopy(config),
			)
			all_results.extend(
				_format_results(raw_results, predictor, outcome_name)
			)

			current = inverse_variance_weighted(data)
			legacy = _legacy_ivw(data)
			published = PUBLISHED_TABLE_1[(predictor, outcome_name)]
			comparisons.append({
				"predictor": predictor,
				"outcome": outcome_name,
				"published_beta": published[0],
				"published_standard_error": published[1],
				"published_pvalue": published[2],
				"legacy_beta": legacy[0],
				"legacy_standard_error": legacy[1],
				"legacy_pvalue": legacy[2],
				"current_beta": current["beta"],
				"current_standard_error": current["standard_error"],
				"current_pvalue": current["pvalue"],
				"published_match": (
					abs(legacy[0] - published[0]) < 0.001
					and abs(legacy[1] - published[1]) < 0.001
					and abs(legacy[2] - published[2]) < 0.01
				),
			})

	results = pd.DataFrame(all_results)
	comparison = pd.DataFrame(comparisons)
	results.to_csv(output_directory / "mr_results.tsv", sep="\t", index=False)
	comparison.to_csv(output_directory / "ivw_comparison.tsv", sep="\t", index=False)
	return results, comparison


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--output", type=Path, default=Path("paper_reproduction"))
	parser.add_argument("--bootstrap-iterations", type=int, default=1000)
	parser.add_argument("--simulation-count", type=int, default=1000)
	parser.add_argument("--seed", type=int, default=0)
	args = parser.parse_args()

	results, comparison = reproduce(
		args.output,
		args.bootstrap_iterations,
		args.simulation_count,
		args.seed,
	)
	print(comparison.to_string(index=False))
	print(f"\nWrote {len(results)} MR results to {args.output}")


if __name__ == "__main__":
	main()
