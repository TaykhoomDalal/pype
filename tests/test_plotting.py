import numpy as np
import pandas as pd

from pype.plotting import category_enrichment, manhattan, volcano


def _results():
	pvalues = np.array([1e-6, 2e-4, 0.02, 0.3, 0.04, 0.8])
	return pd.DataFrame({
		"PheWAS_Category": ["group"] * 6,
		"Category": ["Blood", "Blood", "Body", "Body", "Heart", "Heart"],
		"Data_Field": ["field_" + str(index) for index in range(6)],
		"Description": ["Phenotype " + str(index) for index in range(6)],
		"Independent_Var": ["variant_a"] * 3 + ["variant_b"] * 3,
		"p-val": pvalues,
		"-log(p)": -np.log10(pvalues),
		"beta": [0.4, -0.3, 0.2, -0.1, 0.5, -0.2],
	})


def test_plotting_writes_all_plot_types_without_titles(tmp_path):
	results = _results()
	manhattan_path = tmp_path / "manhattan.png"
	enrichment_path = tmp_path / "enrichment.png"
	volcano_path = tmp_path / "volcano.png"

	manhattan(results, manhattan_path, seed=0)
	category_enrichment(results, enrichment_path)
	volcano(results, volcano_path)

	assert manhattan_path.is_file()
	assert (tmp_path / "manhattan_significant_results.tsv").is_file()
	assert enrichment_path.is_file()
	assert (tmp_path / "volcano_variant_a.png").is_file()
	assert (tmp_path / "volcano_variant_b.png").is_file()
