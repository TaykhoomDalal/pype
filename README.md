# PYPE: Python pipeline for PheWAS <img src="assets/pype_logo.png" width="15" />
A Python Package for PheWAS Execution, Visualization, and Analysis

## Features
PYPE is a command line tool that can be used to run Phenome Wide Association Studies (PheWAS) on data produced by the UK Biobank (Version 1.0). 

Currently, PYPE takes as input BED/BIM/FAM files for the genotype data (only if running PheWAS with a set of variants), and TAB files for the phenotype data and outputs a set of TSV files for each chromosome (if using genotype data) or phenotype (if only using phenotype data) with the results of the PheWAS analysis (including p-values, beta values, and standard error). The user can then run the visualization scripts to generate manhattan plots across all categories with annotations:

<p align="center">
  <img src="assets/abd_manhattan_annotated.png" width="800" />
</p>

without annotations:

<p align="center">
  <img src="assets/abd_manhattan_not_annotated.png" width="800" />
</p>

individual categories of phenotypes with annotations:

<p align="center">
  <img src="assets/abd_liver_Circulating biochemistry_annotated.png" width="800" />
</p>

without annotations:

<p align="center">
  <img src="assets/abd_liver_Circulating biochemistry_not_annotated.png" width="800" />
</p>

volcano plots for each genetic variant tested with annotations:

<p align="center">
  <img src="assets/abd_liver_rs13107325_annotated.png" width="800" />
</p>

without annotations:

<p align="center">
  <img src="assets/abd_liver_rs13107325_not_annotated.png" width="800" />
</p>

and variant enrichment plots that show the number of significant variant-phenotype associations per phenotype category:  
<p align="center">
  <img src="assets/abd_cat_enrichment.png" width="800" />
</p>

PYPE also annotates each significant variant with the upstream and downstream genes close to the variant, providing output files that give descriptions of the variants' functional consequences, the genes' role in biological processes, as well as a variety of other annotations that the user can request. It does this using the lightweight REST API from the [BioThings API](https://biothings.io/). Lastly, using a port of the TwoSampleMR code from [MRC Integrative Epidemiology Unit at the University of Bristol](https://github.com/MRCIEU/TwoSampleMR), PYPE can be used to run Mendelian Randomization (MR) analysis on the variants from the PheWAS analysis to help researchers investigate the causal relationships between genetic variants and phenotypes.

## Installation

PYPE is currently a command line tool and can be installed by cloning the repository and installing the requirements.txt file. It is also extremely important to note that this package requires the installation of plink2, which can be found [here](https://www.cog-genomics.org/plink/2.0/), on either the user's local machine or on their High Performance Cluster (preferred), as a program on the PATH environment variable or accessible through some sort of environment [modules](https://modules.readthedocs.io/en/latest/).

## Usage

PYPE is a command line tool, so the best way of using it is to create a simple bash script where you can specify all the parameters that you want to use.

To see minimal reproducible examples of how to use PYPE, please see the [test_data](test_data/test_scripts/) folder.

### PheWAS

Here we will use the example data provided with PYPE to run a PheWAS analysis on fake UK Biobank phenotype and genotype data using [category 1006](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=1006), Physical measure summary, as the phenotypes of interest. 

```
./geno_phewas.sh 1006
```

There is also an example with running a PheWAS analysis with phenotypes as your independent variables, specifically using phenotypes [90012](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=90012) (Overall Acceleration Average) and [21001](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=21001) (Body Mass Index) against [category 1019](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=1019), Linked Health Outcomes, as the phenotypes of interest. 

```
./pheno_phewas.sh 1019
```

### Visualization

By running the script below, we can generate manhattan plots, volcano plots, and category enrichment plots for the PheWAS analysis that we ran above. 

```
./visualize.sh
```

### MR

By running the script below, we can run an MR analysis using a separate set of variants (not from the example data above as you need at least 2-3 variants to run MR):

```
./run_mr.sh
```

## Download HG19 Gene List

Go to [UCSD Genome Browser Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) and select the following parameters:

Dataset
* clade: Mammal
* genome: Human
* assembly: Feb. 2009 (GRCh37/hg19)
* group: Genes and Gene Predictions
* track: NCBI RefSeq Genes
* table: UCSC RefSeq (refGene)

Define region of interest
* region: genome

Retrieve and display data
* output format: selected fields from primary and related tables
* output file name: hg19_gene_list.txt
* output file format: tab-separated values

<p align="center">
  <img src="assets/UCSC_genome_browser_img1.png" width="800" />
</p>

Then click on get output and select the following fields from hg19.refGene:
* chrom
* cdsStart
* cdsEnd
* name2

<p align="center">
  <img src="assets/UCSC_genome_browser_img2.png" width="800" />
</p>

Finally, click on get output, downloading the file, and then save this file to the cached folder under the pype directory.

## Citation

TBD

## Future Plans

Support for PheWAS on data from a greater variety of studies, faster PheWAS execution, and more visualization options.

Move from annotating nearby genes using manual method using HG19 gene list to using the BioThings API.

## License
This project is licensed under the Apache 2.0 license.

## Author
Taykhoom Dalal (Undergraduate Student at the University of California, Los Angeles)

## Mentor
[Dr. Chirag Patel (Harvard Medical School)](https://dbmi.hms.harvard.edu/people/chirag-patel)
