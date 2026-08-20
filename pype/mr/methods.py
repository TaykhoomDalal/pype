import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.stats import chi2, median_abs_deviation, norm, t


MR_COLUMNS = [
	"exposure_beta",
	"outcome_beta",
	"exposure_standard_error",
	"outcome_standard_error",
]


def _analysis_data(harmonized_data):
	data = harmonized_data[MR_COLUMNS].apply(
		pd.to_numeric,
		errors="coerce",
	).dropna()
	return data[data["exposure_beta"] != 0].reset_index(drop=True)


def _result(beta, standard_error, pvalue, variant_count, **diagnostics):
	return {
		"beta": float(beta),
		"standard_error": float(standard_error),
		"pvalue": float(pvalue),
		"variant_count": int(variant_count),
		**diagnostics,
	}


def inverse_variance_weighted(harmonized_data):
	"""Calculate the multiplicative random-effects IVW estimate."""
	data = _analysis_data(harmonized_data)
	variant_count = len(data)
	if variant_count < 2:
		warnings.warn("IVW requires at least 2 variants.", stacklevel=2)
		return None

	model = sm.WLS(
		data["outcome_beta"],
		data[["exposure_beta"]],
		weights=1 / data["outcome_standard_error"] ** 2,
	).fit()
	beta = model.params["exposure_beta"]
	residual_standard_deviation = np.sqrt(model.scale)
	standard_error = (
		model.bse["exposure_beta"]
		/ min(1, residual_standard_deviation)
	)
	pvalue = 2 * norm.sf(abs(beta / standard_error))
	cochrans_q_degrees_freedom = variant_count - 1
	cochrans_q = model.scale * cochrans_q_degrees_freedom
	cochrans_q_pvalue = chi2.sf(
		cochrans_q,
		cochrans_q_degrees_freedom,
	)
	i_squared = (
		max(
			0,
			(cochrans_q - cochrans_q_degrees_freedom) / cochrans_q,
		)
		if cochrans_q > 0
		else 0
	)

	return _result(
		beta,
		standard_error,
		pvalue,
		variant_count,
		cochrans_q=float(cochrans_q),
		cochrans_q_degrees_freedom=cochrans_q_degrees_freedom,
		cochrans_q_pvalue=float(cochrans_q_pvalue),
		i_squared=float(i_squared),
	)


def mr_egger(harmonized_data):
	"""Calculate the MR-Egger slope and intercept."""
	data = _analysis_data(harmonized_data)
	variant_count = len(data)
	if variant_count < 3:
		warnings.warn("MR-Egger requires at least 3 variants.", stacklevel=2)
		return None

	exposure_sign = np.sign(data["exposure_beta"])
	oriented_exposure = data["exposure_beta"].abs()
	oriented_outcome = data["outcome_beta"] * exposure_sign
	design_matrix = sm.add_constant(
		oriented_exposure.rename("exposure_beta"),
		has_constant="add",
	)
	model = sm.WLS(
		oriented_outcome,
		design_matrix,
		weights=1 / data["outcome_standard_error"] ** 2,
	).fit()
	residual_standard_deviation = np.sqrt(model.scale)
	degrees_freedom = variant_count - 2

	beta = model.params["exposure_beta"]
	standard_error = (
		model.bse["exposure_beta"]
		/ min(1, residual_standard_deviation)
	)
	pvalue = 2 * t.sf(abs(beta / standard_error), degrees_freedom)
	intercept = model.params["const"]
	intercept_standard_error = (
		model.bse["const"]
		/ min(1, residual_standard_deviation)
	)
	intercept_pvalue = 2 * t.sf(
		abs(intercept / intercept_standard_error),
		degrees_freedom,
	)
	cochrans_q = model.scale * degrees_freedom
	cochrans_q_pvalue = chi2.sf(cochrans_q, degrees_freedom)
	i_squared = (
		max(0, (cochrans_q - degrees_freedom) / cochrans_q)
		if cochrans_q > 0
		else 0
	)

	return _result(
		beta,
		standard_error,
		pvalue,
		variant_count,
		intercept=float(intercept),
		intercept_standard_error=float(intercept_standard_error),
		intercept_pvalue=float(intercept_pvalue),
		cochrans_q=float(cochrans_q),
		cochrans_q_degrees_freedom=degrees_freedom,
		cochrans_q_pvalue=float(cochrans_q_pvalue),
		i_squared=float(i_squared),
	)


def _weighted_median(ratio_estimates, weights):
	sorted_indices = np.argsort(ratio_estimates)
	ordered_estimates = np.asarray(ratio_estimates)[sorted_indices]
	ordered_weights = np.asarray(weights)[sorted_indices]
	cumulative_weights = (
		np.cumsum(ordered_weights) - 0.5 * ordered_weights
	)
	cumulative_weights /= np.sum(ordered_weights)
	lower_index = np.max(np.where(cumulative_weights < 0.5))
	return (
		ordered_estimates[lower_index]
		+ (
			ordered_estimates[lower_index + 1]
			- ordered_estimates[lower_index]
		)
		* (0.5 - cumulative_weights[lower_index])
		/ (
			cumulative_weights[lower_index + 1]
			- cumulative_weights[lower_index]
		)
	)


def _weighted_median_bootstrap(
	exposure_beta,
	outcome_beta,
	exposure_standard_error,
	outcome_standard_error,
	weights,
	bootstrap_iterations,
):
	bootstrap_estimates = np.zeros(bootstrap_iterations)
	for iteration in range(bootstrap_iterations):
		sampled_exposure = np.random.normal(
			exposure_beta,
			exposure_standard_error,
		)
		sampled_outcome = np.random.normal(
			outcome_beta,
			outcome_standard_error,
		)
		bootstrap_estimates[iteration] = _weighted_median(
			sampled_outcome / sampled_exposure,
			weights,
		)
	return np.std(bootstrap_estimates, ddof=1)


def _ratio_variance(data):
	return (
		data["outcome_standard_error"] ** 2
		/ data["exposure_beta"] ** 2
		+ data["outcome_beta"] ** 2
		* data["exposure_standard_error"] ** 2
		/ data["exposure_beta"] ** 4
	).to_numpy()


def simple_median(harmonized_data, bootstrap_iterations=1000):
	data = _analysis_data(harmonized_data)
	variant_count = len(data)
	if variant_count < 3:
		warnings.warn(
			"Simple median requires at least 3 variants.",
			stacklevel=2,
		)
		return None

	ratio_estimates = (
		data["outcome_beta"] / data["exposure_beta"]
	).to_numpy()
	weights = np.repeat(1 / variant_count, variant_count)
	beta = _weighted_median(ratio_estimates, weights)
	standard_error = _weighted_median_bootstrap(
		data["exposure_beta"].to_numpy(),
		data["outcome_beta"].to_numpy(),
		data["exposure_standard_error"].to_numpy(),
		data["outcome_standard_error"].to_numpy(),
		weights,
		bootstrap_iterations,
	)
	pvalue = 2 * norm.sf(abs(beta / standard_error))
	return _result(beta, standard_error, pvalue, variant_count)


def weighted_median(harmonized_data, bootstrap_iterations=1000):
	data = _analysis_data(harmonized_data)
	variant_count = len(data)
	if variant_count < 3:
		warnings.warn(
			"Weighted median requires at least 3 variants.",
			stacklevel=2,
		)
		return None

	ratio_estimates = (
		data["outcome_beta"] / data["exposure_beta"]
	).to_numpy()
	weights = 1 / _ratio_variance(data)
	beta = _weighted_median(ratio_estimates, weights)
	standard_error = _weighted_median_bootstrap(
		data["exposure_beta"].to_numpy(),
		data["outcome_beta"].to_numpy(),
		data["exposure_standard_error"].to_numpy(),
		data["outcome_standard_error"].to_numpy(),
		weights,
		bootstrap_iterations,
	)
	pvalue = 2 * norm.sf(abs(beta / standard_error))
	return _result(beta, standard_error, pvalue, variant_count)


def penalized_weighted_median(
	harmonized_data,
	bootstrap_iterations=1000,
	penalty_constant=20,
):
	data = _analysis_data(harmonized_data)
	variant_count = len(data)
	if variant_count < 3:
		warnings.warn(
			"Penalized weighted median requires at least 3 variants.",
			stacklevel=2,
		)
		return None

	ratio_estimates = (
		data["outcome_beta"] / data["exposure_beta"]
	).to_numpy()
	weights = 1 / _ratio_variance(data)
	initial_beta = _weighted_median(ratio_estimates, weights)
	penalty_pvalues = chi2.sf(
		weights * (ratio_estimates - initial_beta) ** 2,
		df=1,
	)
	penalized_weights = weights * np.minimum(
		1,
		penalty_pvalues * penalty_constant,
	)
	beta = _weighted_median(ratio_estimates, penalized_weights)
	standard_error = _weighted_median_bootstrap(
		data["exposure_beta"].to_numpy(),
		data["outcome_beta"].to_numpy(),
		data["exposure_standard_error"].to_numpy(),
		data["outcome_standard_error"].to_numpy(),
		penalized_weights,
		bootstrap_iterations,
	)
	pvalue = 2 * norm.sf(abs(beta / standard_error))
	return _result(beta, standard_error, pvalue, variant_count)


def _distribute_bins(values, weights, lower_bound, upper_bound, bin_count):
	distributed = np.zeros(2 * bin_count)
	bin_width = (upper_bound - lower_bound) / (bin_count - 1)
	for value, weight in zip(values, weights):
		if not np.isfinite(value):
			continue
		position = (value - lower_bound) / bin_width
		lower_index = int(np.floor(position))
		fraction = position - lower_index
		if 0 <= lower_index < bin_count - 1:
			distributed[lower_index] += (1 - fraction) * weight
			distributed[lower_index + 1] += fraction * weight
		elif lower_index == -1:
			distributed[0] += fraction * weight
		elif lower_index == bin_count - 1:
			distributed[bin_count - 1] += (1 - fraction) * weight
	return distributed


def kernel_density(values, bandwidth, weights):
	"""Compute the weighted Gaussian density used by TwoSampleMR modes."""
	values = np.asarray(values)
	weights = np.asarray(weights)
	if not np.issubdtype(values.dtype, np.number):
		raise ValueError("values must be numeric")
	if len(weights) != len(values):
		raise ValueError("values and weights must have equal length")
	if np.any(np.isnan(values)):
		raise ValueError("values contain missing values")
	if not np.all(np.isfinite(weights)) or np.any(weights < 0):
		raise ValueError("weights must be finite and non-negative")
	if not np.isfinite(bandwidth) or bandwidth <= 0:
		raise ValueError("bandwidth must be positive and finite")

	weight_sum = np.sum(weights)
	total_mass = (
		1
		if np.allclose(weight_sum, 1)
		else np.sum(weights[np.isfinite(values)]) / weight_sum
	)
	finite_mask = np.isfinite(values)
	values = values[finite_mask]
	weights = weights[finite_mask]

	bin_count = 512
	cut = 3
	density_minimum = np.min(values) - cut * bandwidth
	density_maximum = np.max(values) + cut * bandwidth
	lower_bound = density_minimum - 4 * bandwidth
	upper_bound = density_maximum + 4 * bandwidth
	distributed = (
		_distribute_bins(
			values,
			weights,
			lower_bound,
			upper_bound,
			bin_count,
		)
		* total_mass
	)

	kernel_coordinates = np.linspace(
		0,
		2 * (upper_bound - lower_bound),
		num=2 * bin_count,
	)
	kernel_coordinates[bin_count:] = -np.flip(
		kernel_coordinates[1 : bin_count + 1]
	)
	kernel = norm.pdf(kernel_coordinates, scale=bandwidth)
	density_values = ifft(
		fft(distributed) * np.conj(fft(kernel))
	).real
	density_values = np.maximum(
		0,
		density_values[:bin_count] / len(distributed),
	)

	source_coordinates = np.linspace(
		density_minimum - 4 * bandwidth,
		density_maximum + 4 * bandwidth,
		bin_count,
	)
	output_coordinates = np.linspace(
		density_minimum,
		density_maximum,
		bin_count,
	)
	output_density = interp1d(
		source_coordinates,
		density_values,
		fill_value="extrapolate",
	)(output_coordinates)
	return {
		"x": output_coordinates,
		"y": output_density,
	}


def _mode_estimate(
	ratio_estimates,
	ratio_standard_errors,
	bandwidth_scale,
):
	variant_count = len(ratio_estimates)
	bandwidth_rule = (
		0.9
		* min(
			np.std(ratio_estimates, ddof=1),
			median_abs_deviation(ratio_estimates, scale="normal"),
		)
		/ variant_count ** (1 / 5)
	)
	weights = ratio_standard_errors**-2
	weights /= np.sum(weights)
	bandwidth = max(1e-8, bandwidth_rule * bandwidth_scale)
	density_result = kernel_density(
		ratio_estimates,
		bandwidth,
		weights,
	)
	return density_result["x"][np.argmax(density_result["y"])]


def _mode_method(
	harmonized_data,
	bootstrap_iterations,
	bandwidth_scale,
	weighted,
	nome,
	penalty_constant=None,
):
	data = _analysis_data(harmonized_data)
	variant_count = len(data)
	if variant_count < 3:
		return None

	ratio_estimates = (
		data["outcome_beta"] / data["exposure_beta"]
	).to_numpy()
	full_standard_errors = np.sqrt(_ratio_variance(data))
	nome_standard_errors = (
		data["outcome_standard_error"]
		/ data["exposure_beta"].abs()
	).to_numpy()
	sampling_standard_errors = (
		nome_standard_errors if nome else full_standard_errors
	)
	mode_standard_errors = (
		sampling_standard_errors
		if weighted
		else np.ones(variant_count)
	)

	if penalty_constant is None:
		beta = _mode_estimate(
			ratio_estimates,
			mode_standard_errors,
			bandwidth_scale,
		)
	else:
		initial_beta = _mode_estimate(
			ratio_estimates,
			full_standard_errors,
			bandwidth_scale,
		)
		initial_weights = 1 / full_standard_errors**2
		penalty_pvalues = chi2.sf(
			initial_weights * (ratio_estimates - initial_beta) ** 2,
			df=1,
		)
		penalized_weights = initial_weights * np.minimum(
			1,
			penalty_pvalues * penalty_constant,
		)
		beta = _mode_estimate(
			ratio_estimates,
			np.sqrt(1 / penalized_weights),
			bandwidth_scale,
		)

	bootstrap_estimates = np.zeros(bootstrap_iterations)
	for iteration in range(bootstrap_iterations):
		sampled_ratios = np.random.normal(
			ratio_estimates,
			sampling_standard_errors,
		)
		if penalty_constant is None:
			bootstrap_estimates[iteration] = _mode_estimate(
				sampled_ratios,
				mode_standard_errors,
				bandwidth_scale,
			)
			continue

		initial_beta = _mode_estimate(
			sampled_ratios,
			full_standard_errors,
			bandwidth_scale,
		)
		initial_weights = 1 / full_standard_errors**2
		penalty_pvalues = chi2.sf(
			initial_weights * (sampled_ratios - initial_beta) ** 2,
			df=1,
		)
		penalized_weights = initial_weights * np.minimum(
			1,
			penalty_pvalues * penalty_constant,
		)
		bootstrap_estimates[iteration] = _mode_estimate(
			sampled_ratios,
			np.sqrt(1 / penalized_weights),
			bandwidth_scale,
		)

	standard_error = median_abs_deviation(
		bootstrap_estimates,
		scale="normal",
	)
	pvalue = 2 * t.sf(
		abs(beta / standard_error),
		variant_count - 1,
	)
	return _result(beta, standard_error, pvalue, variant_count)


def simple_mode(
	harmonized_data,
	bootstrap_iterations=1000,
	bandwidth_scale=1,
):
	return _mode_method(
		harmonized_data,
		bootstrap_iterations,
		bandwidth_scale,
		weighted=False,
		nome=False,
	)


def weighted_mode(
	harmonized_data,
	bootstrap_iterations=1000,
	bandwidth_scale=1,
):
	return _mode_method(
		harmonized_data,
		bootstrap_iterations,
		bandwidth_scale,
		weighted=True,
		nome=False,
	)


def penalized_mode(
	harmonized_data,
	bootstrap_iterations=1000,
	penalty_constant=20,
	bandwidth_scale=1,
):
	return _mode_method(
		harmonized_data,
		bootstrap_iterations,
		bandwidth_scale,
		weighted=True,
		nome=False,
		penalty_constant=penalty_constant,
	)


def simple_mode_nome(
	harmonized_data,
	bootstrap_iterations=1000,
	bandwidth_scale=1,
):
	return _mode_method(
		harmonized_data,
		bootstrap_iterations,
		bandwidth_scale,
		weighted=False,
		nome=True,
	)


def weighted_mode_nome(
	harmonized_data,
	bootstrap_iterations=1000,
	bandwidth_scale=1,
):
	return _mode_method(
		harmonized_data,
		bootstrap_iterations,
		bandwidth_scale,
		weighted=True,
		nome=True,
	)


def _matrix_power(matrix, exponent):
	eigenvalues, eigenvectors = np.linalg.eigh(matrix)
	return (
		eigenvectors
		@ np.diag(eigenvalues**exponent)
		@ eigenvectors.T
	)


def _leave_one_out_rss(outcome, exposures, weights, return_estimates):
	weighted_exposures = exposures * np.sqrt(weights)[:, np.newaxis]
	weighted_outcome = outcome * np.sqrt(weights)
	leave_one_out_estimates = []

	for variant_index in range(len(outcome)):
		keep_mask = np.arange(len(outcome)) != variant_index
		exposure_subset = weighted_exposures[keep_mask]
		outcome_subset = weighted_outcome[keep_mask]
		inverse_cross_product = _matrix_power(
			exposure_subset.T @ exposure_subset,
			-1,
		)
		estimate = (
			inverse_cross_product
			@ exposure_subset.T
			@ outcome_subset
		)
		leave_one_out_estimates.append(estimate)

	leave_one_out_estimates = np.vstack(leave_one_out_estimates)
	residuals = weighted_outcome - np.sum(
		weighted_exposures * leave_one_out_estimates,
		axis=1,
	)
	rss = np.sum(residuals**2)
	if return_estimates:
		return rss, leave_one_out_estimates
	return rss


def _simulate_presso_data(
	outcome,
	exposures,
	outcome_standard_error,
	exposure_standard_error,
	weights,
):
	predictions = np.zeros(len(outcome))
	for variant_index in range(len(outcome)):
		keep_mask = np.arange(len(outcome)) != variant_index
		model = sm.WLS(
			outcome[keep_mask],
			exposures[keep_mask],
			weights=weights[keep_mask],
		).fit()
		predictions[variant_index] = model.predict(
			exposures[[variant_index]]
		)[0]

	simulated_exposures = np.random.normal(
		exposures,
		exposure_standard_error,
	)
	simulated_outcome = np.random.normal(
		predictions,
		outcome_standard_error,
	)
	return simulated_outcome, simulated_exposures


def _sample_distortion_beta(outcome, exposures, weights, outlier_indices):
	all_indices = set(range(len(outcome)))
	non_outlier_indices = list(all_indices - set(outlier_indices))
	sampled_indices = list(outlier_indices) + [
		np.random.choice(non_outlier_indices)
		for _ in range(len(outcome) - len(outlier_indices))
	]
	model_indices = sampled_indices[: len(non_outlier_indices)]
	model = sm.WLS(
		outcome[model_indices],
		exposures[model_indices],
		weights=weights[model_indices],
	).fit()
	return np.asarray(model.params)


def mr_presso(
	harmonized_data,
	run_outlier_test=False,
	run_distortion_test=False,
	significance_threshold=0.05,
	simulation_count=1000,
):
	"""Run the MR-PRESSO global, outlier, and distortion analyses."""
	data = _analysis_data(harmonized_data)
	exposure_columns = ["exposure_beta"]
	variant_count = len(data)
	exposure_count = len(exposure_columns)
	minimum_variant_count = exposure_count + 3

	if variant_count < minimum_variant_count:
		warnings.warn(
			f"MR-PRESSO requires at least {minimum_variant_count} variants.",
			stacklevel=2,
		)
		return None
	if simulation_count <= variant_count:
		raise ValueError(
			"simulation_count must exceed the number of variants."
		)

	exposures = data[exposure_columns].to_numpy(copy=True)
	outcome = data["outcome_beta"].to_numpy(copy=True)
	exposure_standard_error = data[
		["exposure_standard_error"]
	].to_numpy()
	outcome_standard_error = data[
		"outcome_standard_error"
	].to_numpy()

	exposure_sign = np.sign(exposures[:, 0])
	exposures *= exposure_sign[:, np.newaxis]
	outcome *= exposure_sign
	weights = 1 / outcome_standard_error**2

	observed_rss = _leave_one_out_rss(
		outcome,
		exposures,
		weights,
		run_outlier_test,
	)
	simulations = [
		_simulate_presso_data(
			outcome,
			exposures,
			outcome_standard_error,
			exposure_standard_error,
			weights,
		)
		for _ in range(simulation_count)
	]
	simulated_rss = [
		_leave_one_out_rss(
			simulated_outcome,
			simulated_exposures,
			weights,
			run_outlier_test,
		)
		for simulated_outcome, simulated_exposures in simulations
	]

	observed_rss_value = (
		observed_rss[0] if run_outlier_test else observed_rss
	)
	simulated_rss_values = [
		value[0] if run_outlier_test else value
		for value in simulated_rss
	]
	global_pvalue = (
		np.sum(
			np.asarray(simulated_rss_values) > observed_rss_value
		)
		/ simulation_count
	)
	diagnostics = {
		"global_test": {
			"observed_rss": float(observed_rss_value),
			"pvalue": float(global_pvalue),
			"pvalue_upper_bound": (
				1 / simulation_count if global_pvalue == 0 else None
			),
		}
	}

	outlier_results = None
	if global_pvalue < significance_threshold and run_outlier_test:
		leave_one_out_estimates = observed_rss[1]
		rows = []
		for variant_index in range(variant_count):
			observed_difference = (
				outcome[variant_index]
				- exposures[variant_index]
				@ leave_one_out_estimates[variant_index]
			)
			simulated_differences = np.asarray([
				simulated_outcome[variant_index]
				- simulated_exposures[variant_index]
				@ leave_one_out_estimates[variant_index]
				for simulated_outcome, simulated_exposures in simulations
			])
			raw_pvalue = np.mean(
				simulated_differences**2 > observed_difference**2
			)
			rows.append({
				"variant_index": variant_index,
				"observed_rss": float(observed_difference**2),
				"pvalue": float(
					min(raw_pvalue * variant_count, 1)
				),
				"pvalue_upper_bound": (
					variant_count / simulation_count
					if raw_pvalue == 0
					else None
				),
			})
		outlier_results = pd.DataFrame(rows).set_index("variant_index")
		diagnostics["outlier_test"] = outlier_results

	model = sm.WLS(outcome, exposures, weights=weights).fit()
	mr_results = [{
		"method": "MR-PRESSO raw",
		"beta": float(model.params[0]),
		"standard_error": float(model.bse[0]),
		"pvalue": float(model.pvalues[0]),
		"variant_count": variant_count,
	}]

	if (
		run_distortion_test
		and outlier_results is not None
	):
		outlier_indices = outlier_results.index[
			outlier_results["pvalue"] <= significance_threshold
		].tolist()
		if 0 < len(outlier_indices) < variant_count:
			keep_mask = np.ones(variant_count, dtype=bool)
			keep_mask[outlier_indices] = False
			corrected_model = sm.WLS(
				outcome[keep_mask],
				exposures[keep_mask],
				weights=weights[keep_mask],
			).fit()
			mr_results.append({
				"method": "MR-PRESSO outlier-corrected",
				"beta": float(corrected_model.params[0]),
				"standard_error": float(corrected_model.bse[0]),
				"pvalue": float(corrected_model.pvalues[0]),
				"variant_count": int(np.sum(keep_mask)),
			})

			observed_distortion = (
				(model.params - corrected_model.params)
				/ np.abs(corrected_model.params)
			)
			simulated_betas = np.vstack([
				_sample_distortion_beta(
					outcome,
					exposures,
					weights,
					outlier_indices,
				)
				for _ in range(simulation_count)
			])
			simulated_distortion = (
				(model.params - simulated_betas)
				/ np.abs(simulated_betas)
			)
			distortion_pvalue = np.mean(
				np.abs(simulated_distortion)
				> np.abs(observed_distortion),
				axis=0,
			)[0]
			diagnostics["distortion_test"] = {
				"outlier_indices": outlier_indices,
				"distortion_coefficient": float(
					100 * observed_distortion[0]
				),
				"pvalue": float(distortion_pvalue),
				"pvalue_upper_bound": (
					1 / simulation_count
					if distortion_pvalue == 0
					else None
				),
			}

	if variant_count / simulation_count > significance_threshold:
		warnings.warn(
			"MR-PRESSO outlier p-values are too coarse for the "
			"requested threshold; increase simulation_count.",
			stacklevel=2,
		)

	return {
		"mr_results": mr_results,
		**diagnostics,
	}


METHOD_REGISTRY = {
	"ivw": ("Inverse variance weighted", inverse_variance_weighted),
	"egger": ("MR-Egger", mr_egger),
	"simple_median": ("Simple median", simple_median),
	"weighted_median": ("Weighted median", weighted_median),
	"penalized_weighted_median": (
		"Penalized weighted median",
		penalized_weighted_median,
	),
	"simple_mode": ("Simple mode", simple_mode),
	"weighted_mode": ("Weighted mode", weighted_mode),
	"penalized_mode": ("Penalized mode", penalized_mode),
	"simple_mode_nome": ("Simple mode (NOME)", simple_mode_nome),
	"weighted_mode_nome": (
		"Weighted mode (NOME)",
		weighted_mode_nome,
	),
	"presso": ("MR-PRESSO", mr_presso),
}


def run_methods(harmonized_data, method_keys, parameters):
	"""Run selected MR methods using the canonical harmonized schema."""
	results = {}
	for method_key in method_keys:
		display_name, method = METHOD_REGISTRY[method_key]
		results[display_name] = method(
			harmonized_data,
			**parameters.get(method_key, {}),
		)
	return results
