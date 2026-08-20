import numpy as np
import pandas as pd

from pype.plotting import _label, category_enrichment, manhattan, volcano


def results_data():
	pvalues = np.array([1e-6, 2e-4, 0.02, 0.3, 0.04, 0.8])
	return pd.DataFrame({
		"category": ["Blood", "Blood", "Body", "Body", "Heart", "Heart"],
		"outcome": [f"field_{index}" for index in range(6)],
		"description": [f"Phenotype {index}" for index in range(6)],
		"predictor": ["variant_a"] * 3 + ["variant_b"] * 3,
		"pvalue": pvalues,
		"negative_log10_pvalue": -np.log10(pvalues),
		"beta": [0.4, -0.3, 0.2, -0.1, 0.5, -0.2],
	})


def test_plotting_writes_all_plot_types_without_titles(tmp_path, capsys):
	results = results_data()
	manhattan_path = tmp_path / "manhattan.png"
	enrichment_path = tmp_path / "enrichment.png"
	volcano_path = tmp_path / "volcano.png"

	manhattan(results, manhattan_path, seed=0)
	category_enrichment(results, enrichment_path)
	volcano(results, volcano_path)
	assert capsys.readouterr().out == ""

	assert manhattan_path.is_file()
	assert (tmp_path / "manhattan_significant_results.tsv").is_file()
	assert enrichment_path.is_file()
	assert (tmp_path / "volcano_variant_a.png").is_file()
	assert (tmp_path / "volcano_variant_b.png").is_file()


def test_annotation_labels_wrap_without_breaking_words():
	label = _label(
		pd.Series({
			"description": "A long phenotype description for plotting",
			"predictor": "variant",
		}),
		annotation_width=18,
	)

	assert "\n" in label
	assert "phenotype" in label
