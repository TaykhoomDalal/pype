from pathlib import Path

import pandas as pd


def add_phenotype_metadata(
	results,
	metadata,
	result_field="outcome",
	metadata_field="outcome",
	description_column="description",
	category_column="category",
):
	"""Add phenotype descriptions and categories to PheWAS results."""
	if isinstance(metadata, (str, Path)):
		metadata = pd.read_csv(metadata, sep=None, engine="python")
	else:
		metadata = metadata.copy()

	required_columns = {metadata_field, category_column}
	missing_columns = required_columns - set(metadata.columns)
	if missing_columns:
		raise ValueError(
			"Missing metadata columns: " + ", ".join(sorted(missing_columns))
		)

	metadata_columns = [metadata_field, category_column]
	if description_column in metadata:
		metadata_columns.append(description_column)
	metadata = metadata[metadata_columns].copy()
	metadata["_metadata_key"] = metadata[metadata_field].astype(str)
	if metadata["_metadata_key"].duplicated().any():
		raise ValueError("Metadata contains duplicate phenotype fields.")

	data = results.copy()
	if result_field not in data:
		raise ValueError("Results are missing the field column: " + result_field)
	data["_metadata_key"] = data[result_field].astype(str)
	data = data.drop(
		columns=["category", "description"],
		errors="ignore",
	)

	renamed_columns = {category_column: "category"}
	if description_column in metadata:
		renamed_columns[description_column] = "description"
	metadata = metadata.rename(columns=renamed_columns)

	output_columns = ["_metadata_key", "category"]
	if "description" in metadata:
		output_columns.append("description")

	return (
		data.merge(
			metadata[output_columns],
			on="_metadata_key",
			how="left",
			validate="many_to_one",
		)
		.drop(columns="_metadata_key")
	)
