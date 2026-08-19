import pandas as pd

from pype.metadata import add_phenotype_metadata


def test_metadata_join_adds_descriptions_and_categories():
	results = pd.DataFrame({
		"Data_Field": ["body_mass_index", "sex"],
		"p-val": [0.01, 0.02],
	})
	metadata = pd.DataFrame({
		"FieldID": ["body_mass_index", "sex"],
		"Field": ["Body mass index", "Sex"],
		"Category": ["Body measurements", "Demographics"],
	})

	annotated = add_phenotype_metadata(
		results,
		metadata,
		metadata_field="FieldID",
		description_column="Field",
	)

	assert annotated["Description"].tolist() == ["Body mass index", "Sex"]
	assert annotated["Category"].tolist() == ["Body measurements", "Demographics"]
