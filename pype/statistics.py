import numpy as np


def multiple_testing_correction(pvalues, alpha, correction):
	pvalues = np.asarray(pvalues, dtype=float)
	pvalues = pvalues[np.isfinite(pvalues)]

	if pvalues.size == 0:
		raise ValueError("At least one finite p-value is required.")
	if not 0 < alpha <= 1:
		raise ValueError("Alpha must be in the interval (0, 1].")
	if np.any((pvalues < 0) | (pvalues > 1)):
		raise ValueError("P-values must be in the interval [0, 1].")

	if correction == "bonferroni":
		return alpha / pvalues.size
	if correction == "sidak":
		return 1 - np.power(1 - alpha, 1 / pvalues.size)
	if correction == "fdr_bh":
		sorted_pvalues = np.sort(pvalues)
		critical_values = alpha * np.arange(1, pvalues.size + 1) / pvalues.size
		passing = sorted_pvalues <= critical_values
		return sorted_pvalues[passing][-1] if passing.any() else 0.0
	if correction == "no_correction":
		return alpha

	raise ValueError("Invalid multiple testing correction method.")
