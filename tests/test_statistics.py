import pytest

from pype.statistics import multiple_testing_correction


def test_benjamini_hochberg_uses_largest_passing_value():
	pvalues = [0.001, 0.01, 0.03, 0.2]

	assert multiple_testing_correction(pvalues, 0.05, "fdr_bh") == pytest.approx(0.03)


def test_benjamini_hochberg_returns_zero_without_discoveries():
	assert multiple_testing_correction([0.2, 0.4, 0.8], 0.05, "fdr_bh") == 0


def test_multiple_testing_rejects_invalid_pvalues():
	with pytest.raises(ValueError, match="P-values"):
		multiple_testing_correction([0.1, 1.1], 0.05, "bonferroni")
