import pandas as pd
import pytest

from pype.reproduce import PUBLISHED_TABLE_1, legacy_ivw


def test_paper_reference_contains_all_exposure_outcome_pairs():
	assert len(PUBLISHED_TABLE_1) == 8


def test_legacy_ivw_returns_the_original_statsmodels_standard_error():
	data = pd.DataFrame({
		"BETA_EXP": [0.1, 0.2, 0.3],
		"BETA_OUT": [0.05, 0.11, 0.16],
		"SE_OUT": [0.04, 0.04, 0.04],
	})

	beta, standard_error, pvalue = legacy_ivw(data)

	assert beta == pytest.approx(0.535714285714286)
	assert standard_error == pytest.approx(0.008748177652797067)
	assert pvalue == 0
