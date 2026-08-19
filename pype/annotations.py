import pprint
from collections import defaultdict
from pathlib import Path

import pandas as pd


def get_closest_genes(variants, genes, upstream, downstream):
	genes = genes.copy()
	genes.columns = ["CHR", "START", "END", "GENE"]
	genes["CHR"] = genes["CHR"].astype(str).str.removeprefix("chr")

	# Combine transcript rows so each gene occupies one interval per chromosome.
	genes = (
		genes.groupby(["GENE", "CHR"], as_index=False)
		.agg(START=("START", "min"), END=("END", "max"))
	)
	merged = variants.merge(genes, on="CHR", how="inner")
	return merged[
		(merged["POS"] >= merged["START"] - upstream * 1000)
		& (merged["POS"] <= merged["END"] + downstream * 1000)
	]


def annotate_genes(gene_file, rsid_df, down, up, regressions, output_dir, pheno_name, save=True):
	variants = rsid_df.copy()
	variants["CHR"] = variants["CHR"].astype(str)
	genes = pd.read_csv(gene_file, sep="\t")
	genes["#chrom"] = genes["#chrom"].astype(str).str.removeprefix("chr")
	valid_chromosomes = {str(value) for value in range(1, 23)} | {"X", "Y", "XY"}
	genes = genes[genes["#chrom"].isin(valid_chromosomes)]

	nearby = get_closest_genes(variants, genes, up, down).sort_values("rsID")
	gene_map = defaultdict(list)
	for row in nearby.itertuples():
		gene_map[row.rsID].append(row.GENE)

	results = regressions.copy()
	results["Gene"] = results["rsID"].map(
		lambda rsid: ", ".join(gene_map[rsid])
	)

	if save:
		output_directory = Path(output_dir)
		output_directory.mkdir(parents=True, exist_ok=True)
		results.to_csv(
			output_directory / f"{pheno_name}_pheWAS_results_with_nearby_genes.tab",
			sep="\t",
			index=False,
		)
		nearby.to_csv(
			output_directory / f"{pheno_name}_rsID_gene_map.tab",
			sep="\t",
			index=False,
		)

	return results, nearby


def annotate_variants_and_genes(top_results, variant_fields, gene_fields, out_dir):
	import mygene
	import myvariant

	output_directory = Path(out_dir)
	output_directory.mkdir(parents=True, exist_ok=True)
	results = top_results.copy()
	results["Independent_Var"] = results["Independent_Var"].str.split("_").str[0]

	variant_client = myvariant.MyVariantInfo()
	gene_client = mygene.MyGeneInfo()

	for rsid in results["Independent_Var"].dropna().unique():
		path = output_directory / f"{rsid}.summary"
		with path.open("w", encoding="utf-8") as output:
			output.write("The variant is: " + rsid + "\n\n")
			variant = variant_client.getvariant(rsid, fields=variant_fields)
			if variant is not None:
				output.write("Variant information:\n\n")
				pprint.pprint(variant, stream=output)
			else:
				output.write("No variant information was found.\n")

			if "Gene" not in results:
				continue

			gene_values = results.loc[
				results["Independent_Var"] == rsid,
				"Gene",
			].dropna()
			if gene_values.empty:
				continue

			genes = str(gene_values.iloc[0]).split(", ")
			output.write("\nNearby genes: " + ", ".join(genes) + "\n")
			for gene in genes:
				gene_result = gene_client.query(gene, fields=gene_fields)
				hits = gene_result.get("hits", []) if gene_result else []
				if hits and "summary" in hits[0]:
					output.write("\n" + gene + ":\n" + hits[0]["summary"] + "\n")
