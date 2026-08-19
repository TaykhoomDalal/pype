import warnings

import pandas as pd


ALLELE_COMPLEMENTS = str.maketrans("ACGT", "TGCA")


def standardize_variants(data):
	"""Validate variant columns and normalize common allele column names."""
	data = data.copy()
	aliases = {
		"Effect_Allele": ["Effect_Allele", "Effect", "EA"],
		"Non_Effect_Allele": ["Non_Effect_Allele", "Non_Effect", "NEA"],
	}
	for canonical, choices in aliases.items():
		column = next((choice for choice in choices if choice in data), None)
		if column and column != canonical:
			data = data.rename(columns={column: canonical})

	required = {"rsID", "CHR", "Effect_Allele", "BETA", "P", "SE"}
	missing = required - set(data.columns)
	if missing:
		raise ValueError("Missing variant columns: " + ", ".join(sorted(missing)))

	data["rsID"] = data["rsID"].astype(str)
	data["CHR"] = (
		data["CHR"]
		.astype(str)
		.str.replace(r"^chr", "", regex=True)
		.replace({"X": "23", "Y": "24"})
	)
	data["CHR"] = pd.to_numeric(data["CHR"], errors="raise").astype("Int64")
	data["Effect_Allele"] = data["Effect_Allele"].astype(str).str.upper()
	if "Non_Effect_Allele" in data:
		data["Non_Effect_Allele"] = data["Non_Effect_Allele"].astype(str).str.upper()

	for column in ["BETA", "P", "SE", "N"]:
		if column in data:
			data[column] = pd.to_numeric(data[column], errors="coerce")

	data = data.dropna(subset=list(required)).drop_duplicates()
	keys = ["rsID", "CHR"]
	if data.duplicated(keys).any():
		if "N" not in data:
			duplicates = ", ".join(
				data.loc[data.duplicated(keys, keep=False), "rsID"].unique()[:10]
			)
			raise ValueError(
				"Duplicate variants require an N column to choose one row: " + duplicates
			)
		data = data.sort_values("N", ascending=False).drop_duplicates(keys)

	return data


def harmonize(exposure, outcome, suffix_exp, suffix_out):
	"""Align exposure and outcome alleles and return matched effect estimates."""
	exposure = standardize_variants(exposure)
	outcome = standardize_variants(outcome)
	keys = ["rsID", "CHR"]

	def prepare(data, suffix, side):
		columns = keys + ["Effect_Allele", "BETA", "P", "SE"]
		for optional in ["Non_Effect_Allele", "N"]:
			if optional in data:
				columns.append(optional)
		return data[columns].rename(columns={
			"Effect_Allele": "Effect_Allele_" + side,
			"Non_Effect_Allele": "Non_Effect_Allele_" + side,
			"BETA": "BETA" + suffix,
			"P": "P" + suffix,
			"SE": "SE" + suffix,
			"N": "N" + suffix,
		})

	merged = prepare(exposure, suffix_exp, "exp").merge(
		prepare(outcome, suffix_out, "out"),
		on=keys,
		how="inner",
		validate="one_to_one",
	)
	if merged.empty:
		raise ValueError("No variants overlap between the exposure and outcome data.")

	effect_exp = merged["Effect_Allele_exp"]
	effect_out = merged["Effect_Allele_out"]
	complement_out = effect_out.str.translate(ALLELE_COMPLEMENTS)

	if {"Non_Effect_Allele_exp", "Non_Effect_Allele_out"} <= set(merged.columns):
		other_exp = merged["Non_Effect_Allele_exp"]
		other_out = merged["Non_Effect_Allele_out"]
		complement_other_out = other_out.str.translate(ALLELE_COMPLEMENTS)
		same = ((effect_exp == effect_out) & (other_exp == other_out)) | (
			(effect_exp == complement_out) & (other_exp == complement_other_out)
		)
		reversed_alleles = ((effect_exp == other_out) & (other_exp == effect_out)) | (
			(effect_exp == complement_other_out) & (other_exp == complement_out)
		)
		ambiguous = same & reversed_alleles
	else:
		same = (effect_exp == effect_out) | (effect_exp == complement_out)
		reversed_alleles = pd.Series(False, index=merged.index)
		ambiguous = pd.Series(False, index=merged.index)

	if ambiguous.any():
		warnings.warn(
			"Dropping palindromic variants without allele frequencies: "
			+ ", ".join(merged.loc[ambiguous, "rsID"].astype(str).head(10)),
			stacklevel=2,
		)

	reversed_alleles &= ~same
	incompatible = ~(same | reversed_alleles)
	compatible = (same | reversed_alleles) & ~ambiguous
	if incompatible.any():
		warnings.warn(
			"Dropping variants with incompatible alleles: "
			+ ", ".join(merged.loc[incompatible, "rsID"].astype(str).head(10)),
			stacklevel=2,
		)

	merged = merged[compatible].copy()
	merged.loc[reversed_alleles[compatible], "BETA" + suffix_out] *= -1
	if merged.empty:
		raise ValueError("No variants remain after allele harmonization.")

	output_columns = keys + [
		"BETA" + suffix_exp,
		"BETA" + suffix_out,
		"P" + suffix_exp,
		"P" + suffix_out,
		"SE" + suffix_exp,
		"SE" + suffix_out,
	]
	for column in ["N" + suffix_exp, "N" + suffix_out]:
		if column in merged:
			output_columns.append(column)

	return merged[output_columns].reset_index(drop=True)
