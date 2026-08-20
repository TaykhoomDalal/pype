import pandas as pd

from pype.metadata import add_phenotype_metadata


def test_metadata_join_adds_descriptions_and_categories():
	results = pd.DataFrame({
		"outcome": ["body_mass_index", "sex"],
		"pvalue": [0.01, 0.02],
	})
	metadata = pd.DataFrame({
		"field_id": ["body_mass_index", "sex"],
		"label": ["Body mass index", "Sex"],
		"group": ["Body measurements", "Demographics"],
	})

	annotated = add_phenotype_metadata(
		results,
		metadata,
		metadata_field="field_id",
		description_column="label",
		category_column="group",
	)

	assert annotated["description"].tolist() == [
		"Body mass index",
		"Sex",
	]
	assert annotated["category"].tolist() == [
		"Body measurements",
		"Demographics",
	]
