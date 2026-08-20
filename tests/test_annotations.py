import pandas as pd

from pype.annotations import get_closest_genes


def test_gene_annotation_uses_canonical_columns():
	variants = pd.DataFrame({
		"variant_id": ["rs1"],
		"chromosome": [1],
		"position": [1050],
	})
	genes = pd.DataFrame({
		"chromosome": ["chr1"],
		"start": [1000],
		"end": [1100],
		"gene": ["GENE1"],
	})

	result = get_closest_genes(
		variants,
		genes,
		upstream_kb=1,
		downstream_kb=1,
	)

	assert result["variant_id"].tolist() == ["rs1"]
	assert result["gene"].tolist() == ["GENE1"]
