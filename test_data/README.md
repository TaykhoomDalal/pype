# Data Generation

The BED/BIM/FAM files in this folder came from the [plink tutorial](https://zzz.bwh.harvard.edu/plink/binary.shtml)

 

```python
for col in pheno.columns:
    pheno[col] = np.random.permutation(pheno[col].values)
```

From the BIM file, I picked a random variant and deleted the rest, and then converted the genotype encoding from 0 1 2 into A C T
