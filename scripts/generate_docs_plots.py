import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from pype.plotting import category_enrichment, manhattan, volcano


ROOT = Path(__file__).parents[1]
OUTPUT = ROOT / "docs" / "assets"


def synthetic_results():
	rng = np.random.default_rng(12)
	categories = ["Blood", "Body", "Heart", "Metabolism", "Neurology"]
	rows = []

	for predictor_index, predictor in enumerate(["variant_a", "variant_b"]):
		for index in range(50):
			category = categories[index % len(categories)]
			pvalue = 10 ** -rng.uniform(0.05, 4)
			if index in {2, 11, 23}:
				pvalue = 10 ** -(7 + predictor_index + index / 100)
			rows.append({
				"category": category,
				"outcome": f"trait_{index + 1}",
				"description": f"Synthetic phenotype {index + 1}",
				"predictor": predictor,
				"pvalue": pvalue,
				"negative_log10_pvalue": -np.log10(pvalue),
				"beta": rng.normal(0, 0.35),
			})

	return pd.DataFrame(rows)


def main():
	OUTPUT.mkdir(parents=True, exist_ok=True)
	results = synthetic_results()

	with tempfile.TemporaryDirectory() as temporary_directory:
		temporary = Path(temporary_directory)
		manhattan(results, temporary / "manhattan.png", seed=4, width=11, height=5, dpi=160)
		manhattan(
			results,
			temporary / "manhattan_annotated.png",
			annotate=True,
			annotation_count=1,
			color_map="viridis",
			seed=4,
			width=11,
			height=5,
			dpi=160,
		)
		category_enrichment(
			results,
			temporary / "category_enrichment.png",
			width=9,
			height=5,
			dpi=160,
		)
		volcano(
			results,
			temporary / "volcano.png",
			annotate=True,
			annotation_count=3,
			width=8,
			height=5,
			dpi=160,
		)

		for filename in [
			"manhattan.png",
			"manhattan_annotated.png",
			"category_enrichment.png",
		]:
			shutil.copy2(temporary / filename, OUTPUT / filename)
		shutil.copy2(
			temporary / "volcano_variant_a.png",
			OUTPUT / "volcano.png",
		)


if __name__ == "__main__":
	main()
