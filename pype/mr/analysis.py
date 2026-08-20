import json
from importlib.resources import files

import numpy as np
import pandas as pd

from .harmonization import harmonize
from .methods import METHOD_REGISTRY, run_methods


DEFAULT_METHODS = ("ivw", "egger", "weighted_median")
SUPPORTED_METHODS = set(METHOD_REGISTRY)
RESULT_COLUMNS = [
	"method",
	"exposure",
	"outcome",
	"pvalue",
	"beta",
	"standard_error",
	"variant_count",
]


def load_config(path=None):
	"""Load MR bootstrap, bandwidth, penalty, and simulation settings."""
	config_path = path or files("pype.mr").joinpath("defaults.json")
	with open(config_path, encoding="utf-8") as file:
		return json.load(file)


def analyze(
	exposure,
	outcome,
	exposure_name="exposure",
	outcome_name="outcome",
	methods=DEFAULT_METHODS,
	config=None,
	seed=None,
):
	"""Harmonize summary statistics and run selected MR methods."""
	if str(exposure_name) == str(outcome_name):
		raise ValueError("Exposure and outcome names must differ.")

	if methods == "all":
		method_keys = list(METHOD_REGISTRY)
	elif isinstance(methods, str):
		method_keys = [methods]
	else:
		method_keys = list(methods)

	unsupported_methods = set(method_keys) - SUPPORTED_METHODS
	if unsupported_methods:
		raise ValueError(
			"Unsupported MR methods: "
			+ ", ".join(sorted(unsupported_methods))
		)

	harmonized_data = harmonize(exposure, outcome)
	parameters = config or load_config()
	random_state = np.random.get_state()
	try:
		if seed is not None:
			np.random.seed(seed)
		raw_results = run_methods(
			harmonized_data,
			method_keys,
			parameters,
		)
	finally:
		np.random.set_state(random_state)

	rows = []
	diagnostics = {}
	for display_name, result in raw_results.items():
		if result is None:
			continue

		if display_name == "MR-PRESSO":
			diagnostics["presso"] = {
				key: value
				for key, value in result.items()
				if key != "mr_results"
			}
			for estimate in result["mr_results"]:
				rows.append({
					"method": estimate["method"],
					"exposure": exposure_name,
					"outcome": outcome_name,
					"pvalue": estimate["pvalue"],
					"beta": estimate["beta"],
					"standard_error": estimate["standard_error"],
					"variant_count": estimate["variant_count"],
				})
			continue

		rows.append({
			"method": display_name,
			"exposure": exposure_name,
			"outcome": outcome_name,
			"pvalue": result["pvalue"],
			"beta": result["beta"],
			"standard_error": result["standard_error"],
			"variant_count": result["variant_count"],
		})
		diagnostic_values = {
			key: value
			for key, value in result.items()
			if key not in {
				"pvalue",
				"beta",
				"standard_error",
				"variant_count",
			}
		}
		if diagnostic_values:
			diagnostics[display_name] = diagnostic_values

	return pd.DataFrame(rows, columns=RESULT_COLUMNS), diagnostics
