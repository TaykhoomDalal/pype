# PYPE: Python pipeline for PheWAS <img src="assets/pype_logo.png" width="15" alt="">

PYPE runs phenome-wide association studies (PheWAS), creates PheWAS plots, and performs Mendelian randomization (MR) in Python.

[Documentation](https://taykhoomdalal.github.io/pype/) | [Reference paper](https://doi.org/10.1016/j.patter.2024.100982)

## Features

- PheWAS with genotype or phenotype predictors
- Covariate-adjusted linear regression
- Bonferroni, Sidak, and Benjamini-Hochberg correction
- Manhattan, volcano, and category enrichment plots
- Variant-to-gene mapping and optional BioThings annotations
- Eleven MR methods, including MR-Egger, weighted median, weighted mode, and MR-PRESSO

## Installation

```bash
python -m pip install .
```

Optional annotation dependencies:

```bash
python -m pip install ".[annotations]"
```

## PheWAS

Pass one dataframe containing phenotypes and covariates, and another containing the predictors. Both dataframes must have the same sample index.

```python
import pandas as pd
import pype

phenotypes = pd.DataFrame({
    "trait_a": [1.2, 2.0, 1.5, 3.1],
    "trait_b": [0.0, 1.0, 0.0, 1.0],
    "age": [50, 61, 47, 70],
})

predictors = pd.DataFrame({
    "variant_1": [0, 1, 1, 2],
})

results = pype.phenome_wide_association(
    phenotypes,
    predictors,
    outcomes=["trait_a", "trait_b"],
    covariates=["age"],
    min_samples=3,
)
```

The result is a dataframe containing the phenotype, predictor, sample count, p-value, effect estimate, and standard error for each regression.

## Plotting

Plots use the PheWAS result dataframe plus a small metadata table containing phenotype descriptions and categories.

```python
import pandas as pd
import pype
from pype.plotting import category_enrichment, manhattan, volcano

metadata = pd.DataFrame({
    "Data_Field": ["trait_a", "trait_b"],
    "Description": ["Trait A", "Trait B"],
    "Category": ["Measurements", "Diagnoses"],
})
results = pype.add_phenotype_metadata(results, metadata)

manhattan(results, "manhattan.png")
category_enrichment(results, "category_enrichment.png")
volcano(results, "volcano.png")
```

## Mendelian randomization

Exposure and outcome dataframes must contain:

| Column | Description |
| --- | --- |
| `rsID` | Variant identifier |
| `CHR` | Chromosome |
| `Effect_Allele` | Effect allele |
| `Non_Effect_Allele` | Other allele |
| `BETA` | Effect estimate |
| `SE` | Standard error |
| `P` | P-value |
| `N` | Sample size, optional |

```python
import pandas as pd
import pype

exposure = pd.read_csv("exposure.tsv", sep="\t")
outcome = pd.read_csv("outcome.tsv", sep="\t")

results, diagnostics = pype.mendelian_randomization(
    exposure,
    outcome,
    exposure_name="exposure_trait",
    outcome_name="outcome_trait",
    methods=("ivw", "egger", "weighted_median"),
    seed=0,
)
```

Use `methods="all"` to run every available method. PYPE aligns effect alleles, flips reversed effects, and removes incompatible or ambiguous palindromic variants before analysis.

Available methods:

- Inverse variance weighted
- MR-Egger
- Simple, weighted, and penalized weighted median
- Simple, weighted, penalized, and NOME mode estimators
- MR-PRESSO

## Data

Tests use synthetic data. Study data are supplied by users under the access terms of their data source.

## Citation

Dalal T, Patel CJ. PYPE: A pipeline for phenome-wide association and Mendelian randomization in investigator-driven biobank scale analysis. Patterns. 2024. https://doi.org/10.1016/j.patter.2024.100982

## Development

```bash
python -m pytest
python -m build
```

## License

Apache License 2.0. See `LICENSE.md`.
