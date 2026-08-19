import json
from importlib.resources import files

import numpy as np
import pandas as pd

from . import methods as estimators
from .harmonization import harmonize


DEFAULT_METHODS = ("ivw", "egger", "weighted_median")
SUPPORTED_METHODS = {
	"ivw",
	"egger",
	"simple_median",
	"weighted_median",
	"penalized_weighted_median",
	"simple_mode",
	"simple_mode_nome",
	"weighted_mode",
	"penalized_mode",
	"weighted_mode_nome",
	"presso",
}


def load_config(path=None):
	"""Load the bootstrap and MR-PRESSO settings."""
	config_path = path or files("pype.mr").joinpath("defaults.json")
	with open(config_path, encoding="utf-8") as file:
		return json.load(file)


def _scalar(value):
	if isinstance(value, pd.Series):
		return value.iloc[0]
	if isinstance(value, np.ndarray):
		return value.flat[0]
	return value


def analyze(exposure, outcome, exposure_name="exposure", outcome_name="outcome", methods=DEFAULT_METHODS, config=None, seed=None):
	"""Run MR methods on exposure and outcome variant dataframes."""
	if str(exposure_name) == str(outcome_name):
		raise ValueError("Exposure and outcome names must differ.")

	run_all = methods == "all"
	if isinstance(methods, str) and not run_all:
		selected_methods = [methods]
	else:
		selected_methods = [] if run_all else list(methods)
	unknown = set(selected_methods) - SUPPORTED_METHODS
	if unknown:
		raise ValueError("Unsupported MR methods: " + ", ".join(sorted(unknown)))

	suffix_exp = "_" + str(exposure_name)
	suffix_out = "_" + str(outcome_name)
	harmonized = harmonize(exposure, outcome, suffix_exp, suffix_out)
	config = config or load_config()

	random_state = np.random.get_state()
	try:
		if seed is not None:
			np.random.seed(seed)
		raw_results = estimators.run_mr(
			selected_methods,
			harmonized,
			"BETA" + suffix_exp,
			"BETA" + suffix_out,
			"SE" + suffix_exp,
			"SE" + suffix_out,
			config,
			run_all,
		)
	finally:
		np.random.set_state(random_state)

	rows = []
	diagnostics = {}
	for method, result in raw_results.items():
		if result is None:
			continue

		if method == "PRESSO":
			diagnostics["PRESSO"] = {
				key: value for key, value in result.items() if key != "MR_RESULTS"
			}
			for presso_result in result["MR_RESULTS"]:
				rows.append([
					presso_result[0],
					outcome_name,
					float(_scalar(presso_result[1])),
					float(_scalar(presso_result[2])),
					float(_scalar(presso_result[3])),
					len(harmonized),
				])
			continue

		pvalue, effect_size, standard_error, *_ = result
		rows.append([
			method,
			outcome_name,
			float(pvalue),
			float(effect_size),
			float(standard_error),
			len(harmonized),
		])

	results = pd.DataFrame(
		rows,
		columns=[
			"MR_Method",
			"Outcome",
			"P_value",
			"Effect_Size",
			"Standard_Error",
			"Number_SNPs",
		],
	)
	return results, diagnostics
