# Paper reproduction

`python reproducibility/reproduce_paper.py` reproduces the Mendelian randomization portion of the 2024 PYPE paper from public summary statistics.

It downloads:

- [Le Goallec et al. Supplementary Data 1](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29525-9/MediaObjects/41467_2022_29525_MOESM4_ESM.xlsx) for the abdomen and liver aging variants.
- Neale Lab raw GWAS summary statistics for [glycated haemoglobin](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30750_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz), [body mass index](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_raw.gwas.imputed_v3.both_sexes.tsv.bgz), [glucose](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30740_raw.gwas.imputed_v3.both_sexes.varorder.tsv.bgz), and [waist circumference](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/48_raw.gwas.imputed_v3.both_sexes.tsv.bgz).

The first run streams about 2 GB of compressed public GWAS files and saves only the selected variant rows. Later runs reuse those small extracted files.

## Run

From a GitHub source checkout:

```bash
python -m pip install -e . openpyxl
python reproducibility/reproduce_paper.py
```

Outputs are written to `paper_reproduction/`:

- `mr_results.tsv`: current PYPE MR results for every paper exposure/outcome pair.
- `ivw_comparison.tsv`: published Table 1 values, the legacy PYPE IVW calculation, and the corrected current IVW calculation.
- `inputs/`: the small public summary-statistic subsets used for the run.

## Interpretation

The effect estimates reproduce the published MR tables. The current IVW implementation follows TwoSampleMR's under-dispersion correction, so its standard errors and p-values match the R reference tables rather than the legacy PYPE values in cases where the original implementation underestimated the standard error.

Bootstrap standard errors for median and mode estimators vary with the random seed. Use `--seed` and `--bootstrap-iterations` for a repeatable run.

The participant-level PheWAS, Table S20, and the main PheWAS figures require authorized study data and are outside this public reproduction.
