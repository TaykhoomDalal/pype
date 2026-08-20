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

from .mr import load_config
from .mr import methods


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
METHODS = [
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


def download(url, path):
	path.parent.mkdir(parents=True, exist_ok=True)
	if not path.exists():
		urllib.request.urlretrieve(url, path)
	return path


def load_exposures(input_directory):
	workbook = download(LE_GOALLEC_URL, input_directory / "le_goallec_gwas.xlsx")
	data = pd.read_excel(workbook)
	data = data[data["rsID"] == data["IndSigSNP"]]
	exposures = {}

	for name, specification in EXPOSURE_VARIANTS.items():
		rows = data[
			(data["phenotype"] == specification["phenotype"])
			& data["rsID"].isin(specification["variants"])
		].set_index("rsID").loc[specification["variants"]]
		exposures[name] = rows.rename(columns={
			"beta": "BETA_EXP",
			"se": "SE_EXP",
		})[["chr", "pos", "BETA_EXP", "SE_EXP"]]

	return exposures


def extract_outcome(name, filename, positions, input_directory):
	path = input_directory / (name.lower().replace(" ", "_") + ".tsv")
	if path.exists():
		result = pd.read_csv(path, sep="\t")
		if "position" not in result:
			result = pd.DataFrame({
				"position": result["variant"].str.split(":").str[1].astype(int),
				"BETA_OUT": result["beta"].astype(float),
				"SE_OUT": result["se"].astype(float),
			})
		return result

	url = NEALE_BASE_URL + "/" + filename
	selected = []
	with urllib.request.urlopen(url) as response:
		with gzip.GzipFile(fileobj=response) as compressed:
			reader = csv.DictReader(io.TextIOWrapper(compressed), delimiter="\t")
			for row in reader:
				position = int(row["variant"].split(":")[1])
				if position in positions:
					selected.append({
						"position": position,
						"BETA_OUT": float(row["beta"]),
						"SE_OUT": float(row["se"]),
					})

	result = pd.DataFrame(selected)
	result.to_csv(path, sep="\t", index=False)
	return result


def analysis_frame(exposure, outcome):
	return (
		exposure.reset_index()
		.merge(outcome, left_on="pos", right_on="position", validate="one_to_one")
		[["BETA_EXP", "BETA_OUT", "SE_EXP", "SE_OUT"]]
	)


def scalar(value):
	if isinstance(value, pd.Series):
		return value.iloc[0]
	if isinstance(value, np.ndarray):
		return value.flat[0]
	return value


def format_results(raw_results, predictor, outcome, snps):
	rows = []
	for method, result in raw_results.items():
		if result is None:
			continue
		if method == "PRESSO":
			for presso_result in result["MR_RESULTS"]:
				rows.append({
					"Predictor": predictor,
					"Outcome": outcome,
					"Method": presso_result[0],
					"P_value": float(scalar(presso_result[1])),
					"Effect_Size": float(scalar(presso_result[2])),
					"Standard_Error": float(scalar(presso_result[3])),
					"Number_SNPs": snps,
				})
		else:
			pvalue, beta, standard_error, *_ = result
			rows.append({
				"Predictor": predictor,
				"Outcome": outcome,
				"Method": method,
				"P_value": float(pvalue),
				"Effect_Size": float(beta),
				"Standard_Error": float(standard_error),
				"Number_SNPs": snps,
			})
	return rows


def legacy_ivw(data):
	model = wls(
		"BETA_OUT ~ -1 + BETA_EXP",
		data=data,
		weights=1 / data["SE_OUT"] ** 2,
	).fit()
	beta = float(model.params["BETA_EXP"])
	standard_error = float(model.bse["BETA_EXP"])
	pvalue = 2 * norm.sf(abs(beta / standard_error))
	return beta, standard_error, pvalue


def reproduce(output_directory, bootstrap, presso_distributions, seed):
	output_directory.mkdir(parents=True, exist_ok=True)
	input_directory = output_directory / "inputs"
	exposures = load_exposures(input_directory)
	positions = {
		int(position)
		for exposure in exposures.values()
		for position in exposure["pos"]
	}
	outcomes = {
		name: extract_outcome(name, filename, positions, input_directory)
		for name, filename in OUTCOME_FILES.items()
	}

	config = load_config()
	for settings in config.values():
		if "nboot" in settings:
			settings["nboot"] = bootstrap
	config["PRESSO"]["nbDist"] = presso_distributions

	all_results = []
	comparisons = []
	for predictor, exposure in exposures.items():
		for outcome_name, outcome in outcomes.items():
			data = analysis_frame(exposure, outcome)
			np.random.seed(seed)
			selected_methods = METHODS + (["presso"] if len(data) >= 4 else [])
			raw_results = methods.run_mr(
				selected_methods,
				data,
				"BETA_EXP",
				"BETA_OUT",
				"SE_EXP",
				"SE_OUT",
				copy.deepcopy(config),
				run_all=False,
			)
			all_results.extend(
				format_results(raw_results, predictor, outcome_name, len(data))
			)

			current = methods.run_mr_ivw(
				data,
				"BETA_EXP",
				"BETA_OUT",
				"SE_EXP",
				"SE_OUT",
			)
			legacy = legacy_ivw(data)
			published = PUBLISHED_TABLE_1[(predictor, outcome_name)]
			comparisons.append({
				"Predictor": predictor,
				"Outcome": outcome_name,
				"Published_Beta": published[0],
				"Published_SE": published[1],
				"Published_P": published[2],
				"Legacy_Beta": legacy[0],
				"Legacy_SE": legacy[1],
				"Legacy_P": legacy[2],
				"Current_Beta": float(current[1]),
				"Current_SE": float(current[2]),
				"Current_P": float(current[0]),
				"Published_Match": (
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
	parser.add_argument("--bootstrap", type=int, default=1000)
	parser.add_argument("--presso-distributions", type=int, default=1000)
	parser.add_argument("--seed", type=int, default=0)
	args = parser.parse_args()

	results, comparison = reproduce(
		args.output,
		args.bootstrap,
		args.presso_distributions,
		args.seed,
	)
	print(comparison.to_string(index=False))
	print(f"\nWrote {len(results)} MR results to {args.output}")


if __name__ == "__main__":
	main()
