# Data Generation

The BED/BIM/FAM files in this folder came from the [plink tutorial](https://zzz.bwh.harvard.edu/plink/binary.shtml)

To generate the phenotype data located in [input_ex](input_ex), I took the existing UKBB anonymized data from project 52887, choose the first 980 rows, and randomly shuffled the values in each column of the data to ensure that there was no way the original values could be recovered.

```python
for col in pheno.columns:
    pheno[col] = np.random.permutation(pheno[col].values)
```

From the BIM file, I picked a random variant and deleted the rest, and then converted the genotype encoding from 0 1 2 into A C T
