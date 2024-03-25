# Data Generation

The BED/BIM/FAM files located in [genotype_files](genotype_files) this folder came from the [plink tutorial](https://zzz.bwh.harvard.edu/plink/binary.shtml)

From the BIM file, I picked a random variant and deleted the rest, and then converted the genotype encoding from 0 1 2 into A C T

To generate the phenotype data located in [phenotype_files](phenotype_files), I took the UKBB data generated in [/example/UKBB](../../UKBB/) and converted the columns to a human readable format. I also reduced the number of columns by taking the first instance of each data field. The specifics are in [convert_UKBB_to_general.py](convert_UKBB_to_general.py)