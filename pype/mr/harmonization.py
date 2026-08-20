import warnings

import numpy as np
import pandas as pd


ALLELE_COMPLEMENTS = str.maketrans("ACGT", "TGCA")
REQUIRED_COLUMNS = {
	"variant_id",
	"chromosome",
	"effect_allele",
	"beta",
	"standard_error",
	"pvalue",
}


def standardize_variants(data):
	"""Validate and normalize summary statistics."""
	data = data.copy()
	missing_columns = REQUIRED_COLUMNS - set(data.columns)
	if missing_columns:
		raise ValueError(
			"Missing variant columns: " + ", ".join(sorted(missing_columns))
		)

	data = data.dropna(subset=list(REQUIRED_COLUMNS))
	data["variant_id"] = data["variant_id"].astype(str)
	data["chromosome"] = (
		data["chromosome"]
		.astype(str)
		.str.replace(r"^chr", "", regex=True)
		.replace({"X": "23", "Y": "24"})
	)
	data["chromosome"] = pd.to_numeric(
		data["chromosome"],
		errors="raise",
	).astype("Int64")
	data["effect_allele"] = data["effect_allele"].astype(str).str.upper()
	if "other_allele" in data:
		if data["other_allele"].isna().any():
			raise ValueError("other_allele contains missing values.")
		data["other_allele"] = data["other_allele"].astype(str).str.upper()

	for column_name in [
		"beta",
		"standard_error",
		"pvalue",
		"sample_size",
	]:
		if column_name in data:
			data[column_name] = pd.to_numeric(
				data[column_name],
				errors="coerce",
			)

	data = data.dropna(subset=list(REQUIRED_COLUMNS)).drop_duplicates()
	if not np.isfinite(
		data[["beta", "standard_error", "pvalue"]]
	).all().all():
		raise ValueError(
			"beta, standard_error, and pvalue must be finite."
		)
	if (data["standard_error"] <= 0).any():
		raise ValueError("standard_error must be greater than 0.")
	if ((data["pvalue"] < 0) | (data["pvalue"] > 1)).any():
		raise ValueError("pvalue must be in the interval [0, 1].")
	if not data["effect_allele"].str.fullmatch(r"[ACGT]+").all():
		raise ValueError("effect_allele must contain only A, C, G, and T.")
	if (
		"other_allele" in data
		and not data["other_allele"].str.fullmatch(r"[ACGT]+").all()
	):
		raise ValueError("other_allele must contain only A, C, G, and T.")
	if (
		"sample_size" in data
		and (data["sample_size"].dropna() <= 0).any()
	):
		raise ValueError("sample_size must be greater than 0.")
	key_columns = ["variant_id", "chromosome"]
	if data.duplicated(key_columns).any():
		if "sample_size" not in data:
			duplicate_ids = ", ".join(
				data.loc[
					data.duplicated(key_columns, keep=False),
					"variant_id",
				]
				.unique()[:10]
			)
			raise ValueError(
				"Duplicate variants require sample_size: " + duplicate_ids
			)
		data = (
			data.sort_values("sample_size", ascending=False)
			.drop_duplicates(key_columns)
		)

	return data


def _prepare_side(data, prefix):
	column_names = [
		"variant_id",
		"chromosome",
		"effect_allele",
		"beta",
		"standard_error",
		"pvalue",
	]
	for optional_column in ["other_allele", "sample_size"]:
		if optional_column in data:
			column_names.append(optional_column)

	renamed_columns = {
		"effect_allele": f"{prefix}_effect_allele",
		"other_allele": f"{prefix}_other_allele",
		"beta": f"{prefix}_beta",
		"standard_error": f"{prefix}_standard_error",
		"pvalue": f"{prefix}_pvalue",
		"sample_size": f"{prefix}_sample_size",
	}
	return data[column_names].rename(columns=renamed_columns)


def harmonize(exposure, outcome):
	"""Align exposure and outcome alleles and return matched statistics."""
	exposure = standardize_variants(exposure)
	outcome = standardize_variants(outcome)
	key_columns = ["variant_id", "chromosome"]

	harmonized = _prepare_side(exposure, "exposure").merge(
		_prepare_side(outcome, "outcome"),
		on=key_columns,
		how="inner",
		validate="one_to_one",
	)
	if harmonized.empty:
		raise ValueError(
			"No variants overlap between the exposure and outcome data."
		)

	exposure_effect = harmonized["exposure_effect_allele"]
	outcome_effect = harmonized["outcome_effect_allele"]
	outcome_effect_complement = outcome_effect.str.translate(
		ALLELE_COMPLEMENTS
	)

	if {
		"exposure_other_allele",
		"outcome_other_allele",
	} <= set(harmonized.columns):
		exposure_other = harmonized["exposure_other_allele"]
		outcome_other = harmonized["outcome_other_allele"]
		outcome_other_complement = outcome_other.str.translate(
			ALLELE_COMPLEMENTS
		)
		aligned = (
			(exposure_effect == outcome_effect)
			& (exposure_other == outcome_other)
		) | (
			(exposure_effect == outcome_effect_complement)
			& (exposure_other == outcome_other_complement)
		)
		reversed_alleles = (
			(exposure_effect == outcome_other)
			& (exposure_other == outcome_effect)
		) | (
			(exposure_effect == outcome_other_complement)
			& (exposure_other == outcome_effect_complement)
		)
		ambiguous_palindromes = aligned & reversed_alleles
	else:
		aligned = (
			(exposure_effect == outcome_effect)
			| (exposure_effect == outcome_effect_complement)
		)
		reversed_alleles = pd.Series(False, index=harmonized.index)
		ambiguous_palindromes = pd.Series(
			False,
			index=harmonized.index,
		)

	if ambiguous_palindromes.any():
		warnings.warn(
			"Dropping palindromic variants without allele frequencies: "
			+ ", ".join(
				harmonized.loc[
					ambiguous_palindromes,
					"variant_id",
				]
				.astype(str)
				.head(10)
			),
			stacklevel=2,
		)

	reversed_alleles &= ~aligned
	incompatible = ~(aligned | reversed_alleles)
	compatible = (aligned | reversed_alleles) & ~ambiguous_palindromes
	if incompatible.any():
		warnings.warn(
			"Dropping variants with incompatible alleles: "
			+ ", ".join(
				harmonized.loc[incompatible, "variant_id"]
				.astype(str)
				.head(10)
			),
			stacklevel=2,
		)

	harmonized = harmonized[compatible].copy()
	harmonized.loc[
		reversed_alleles[compatible],
		"outcome_beta",
	] *= -1
	if harmonized.empty:
		raise ValueError("No variants remain after allele harmonization.")

	output_columns = key_columns + [
		"exposure_beta",
		"outcome_beta",
		"exposure_standard_error",
		"outcome_standard_error",
		"exposure_pvalue",
		"outcome_pvalue",
	]
	for optional_column in [
		"exposure_sample_size",
		"outcome_sample_size",
	]:
		if optional_column in harmonized:
			output_columns.append(optional_column)

	return harmonized[output_columns].reset_index(drop=True)
