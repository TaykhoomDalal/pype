from pathlib import Path

import pandas as pd


def add_phenotype_metadata(results, metadata, result_field="Data_Field", metadata_field="Data_Field", description_column="Description", category_column="Category"):
	"""Add phenotype descriptions and categories from a metadata table."""
	if isinstance(metadata, (str, Path)):
		metadata = pd.read_csv(metadata, sep=None, engine="python")
	else:
		metadata = metadata.copy()

	required = {metadata_field, category_column}
	missing = required - set(metadata.columns)
	if missing:
		raise ValueError("Missing metadata columns: " + ", ".join(sorted(missing)))

	columns = [metadata_field, category_column]
	if description_column in metadata:
		columns.append(description_column)
	metadata = metadata[columns].copy()
	metadata["_metadata_key"] = metadata[metadata_field].astype(str)
	if metadata["_metadata_key"].duplicated().any():
		raise ValueError("Metadata contains duplicate phenotype fields.")

	data = results.copy()
	if result_field not in data:
		raise ValueError("Results are missing the field column: " + result_field)
	data["_metadata_key"] = data[result_field].astype(str)
	data = data.drop(columns=[category_column, description_column], errors="ignore")

	rename = {category_column: "Category"}
	if description_column in metadata:
		rename[description_column] = "Description"
	metadata = metadata.rename(columns=rename)
	keep = ["_metadata_key", "Category"]
	if "Description" in metadata:
		keep.append("Description")

	return data.merge(metadata[keep], on="_metadata_key", how="left", validate="many_to_one").drop(columns="_metadata_key")
