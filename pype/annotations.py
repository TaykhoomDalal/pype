import pprint
from collections import defaultdict
from pathlib import Path

import pandas as pd


def get_closest_genes(variants, genes, upstream_kb, downstream_kb):
	"""Return gene intervals near each variant."""
	required_variant_columns = {
		"variant_id",
		"chromosome",
		"position",
	}
	required_gene_columns = {
		"chromosome",
		"start",
		"end",
		"gene",
	}
	missing_variant_columns = required_variant_columns - set(variants.columns)
	missing_gene_columns = required_gene_columns - set(genes.columns)
	if missing_variant_columns:
		raise ValueError(
			"Missing variant columns: "
			+ ", ".join(sorted(missing_variant_columns))
		)
	if missing_gene_columns:
		raise ValueError(
			"Missing gene columns: "
			+ ", ".join(sorted(missing_gene_columns))
		)

	variants = variants.copy()
	genes = genes.copy()
	variants["variant_id"] = variants["variant_id"].astype(str)
	variants["chromosome"] = (
		variants["chromosome"].astype(str).str.removeprefix("chr")
	)
	variants["position"] = pd.to_numeric(
		variants["position"],
		errors="raise",
	)
	genes["chromosome"] = (
		genes["chromosome"].astype(str).str.removeprefix("chr")
	)
	genes[["start", "end"]] = genes[["start", "end"]].apply(
		pd.to_numeric,
		errors="raise",
	)

	# Combine transcript rows so each gene occupies one interval per chromosome.
	genes = (
		genes.groupby(["gene", "chromosome"], as_index=False)
		.agg(start=("start", "min"), end=("end", "max"))
	)
	merged = variants.merge(genes, on="chromosome", how="inner")
	return merged[
		(
			merged["position"]
			>= merged["start"] - upstream_kb * 1000
		)
		& (
			merged["position"]
			<= merged["end"] + downstream_kb * 1000
		)
	]


def annotate_genes(
	gene_intervals_file,
	variants,
	results,
	upstream_kb=10,
	downstream_kb=10,
	output_directory=None,
	output_prefix="pype",
):
	"""Add nearby gene symbols to variant association results."""
	gene_intervals = pd.read_csv(gene_intervals_file, sep="\t")
	nearby_genes = get_closest_genes(
		variants,
		gene_intervals,
		upstream_kb,
		downstream_kb,
	).sort_values("variant_id")

	genes_by_variant = defaultdict(list)
	for row in nearby_genes.itertuples():
		genes_by_variant[row.variant_id].append(row.gene)

	annotated_results = results.copy()
	result_variant_column = (
		"variant_id" if "variant_id" in annotated_results else "predictor"
	)
	if result_variant_column not in annotated_results:
		raise ValueError(
			"Results require a variant_id or predictor column."
		)
	annotated_results["gene"] = annotated_results[result_variant_column].map(
		lambda variant_id: ", ".join(genes_by_variant[str(variant_id)])
	)

	if output_directory is not None:
		output_directory = Path(output_directory)
		output_directory.mkdir(parents=True, exist_ok=True)
		annotated_results.to_csv(
			output_directory / f"{output_prefix}_results_with_genes.tsv",
			sep="\t",
			index=False,
		)
		nearby_genes.to_csv(
			output_directory / f"{output_prefix}_variant_gene_map.tsv",
			sep="\t",
			index=False,
		)

	return annotated_results, nearby_genes


def annotate_variants_and_genes(
	results,
	output_directory,
	variant_fields=("dbsnp",),
	gene_fields=("summary",),
):
	"""Write BioThings summaries for significant variants and nearby genes."""
	import mygene
	import myvariant

	output_directory = Path(output_directory)
	output_directory.mkdir(parents=True, exist_ok=True)
	results = results.copy()
	variant_column = "variant_id" if "variant_id" in results else "predictor"
	if variant_column not in results:
		raise ValueError("Results require a variant_id or predictor column.")
	results[variant_column] = (
		results[variant_column].astype(str).str.split("_").str[0]
	)

	variant_client = myvariant.MyVariantInfo()
	gene_client = mygene.MyGeneInfo()

	for variant_id in results[variant_column].dropna().unique():
		path = output_directory / f"{variant_id}.summary"
		with path.open("w", encoding="utf-8") as output:
			output.write("Variant: " + variant_id + "\n\n")
			variant = variant_client.getvariant(
				variant_id,
				fields=list(variant_fields),
			)
			if variant is not None:
				output.write("Variant information:\n\n")
				pprint.pprint(variant, stream=output)
			else:
				output.write("Variant information was unavailable.\n")

			if "gene" not in results:
				continue

			gene_values = results.loc[
				results[variant_column] == variant_id,
				"gene",
			].dropna()
			if gene_values.empty:
				continue

			genes = str(gene_values.iloc[0]).split(", ")
			output.write("\nNearby genes: " + ", ".join(genes) + "\n")
			for gene in genes:
				gene_result = gene_client.query(
					gene,
					fields=list(gene_fields),
				)
				hits = gene_result.get("hits", []) if gene_result else []
				if hits and "summary" in hits[0]:
					output.write(
						"\n" + gene + ":\n" + hits[0]["summary"] + "\n"
					)
