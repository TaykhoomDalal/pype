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

Here we will use the example data provided with PYPE to run a PheWAS analysis on fake UK Biobank phenotype and genotype data using [category 1006](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=1006), Physical measure summary, as the phenotypes of interest. Both PheWAS scripts are run using the run_phewas.py file.

[geno_phewas.sh](test_data/test_scripts/geno_phewas.sh) 1006

There is also an example with running a PheWAS analysis with phenotypes as your independent variables, specifically using phenotypes [90012](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=90012) (Overall Acceleration Average) and [21001](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=21001) (Body Mass Index) against [category 1019](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=1019), Linked Health Outcomes, as the phenotypes of interest. 

[pheno_phewas.sh](test_data/test_scripts/pheno_phewas.sh) 1019

### Visualization

After running the PheWAS analysis, you need to use annotate_categories.py file to annotate the results with the category names. This file will create the file that can be used to create the various visualizations and is especially useful if you used the same independent variables across multiple analyses, as it will aggregate the results into one file that can be used for visualization. 

[annotate_categories.sh](test_data/test_scripts/annotate_categories.sh)

By running the script below, we can generate manhattan plots, volcano plots, and category enrichment plots for the PheWAS analysis that we ran above, as it calls the visualize.py file. 

[visualize.sh](test_data/test_scripts/visualize.sh)

### MR

By running the script below, we can run an MR analysis using a separate set of variants (not from the example data above as you need at least 2-3 variants to run MR). This script uses the run_mr.py file to do so.

[run_mr.sh](test_data/test_scripts/run_mr.sh)

## File Formats

#### PheWAS file formats

--ukbiobank_phenotype_file /path/to/ukbb/phenotype/file

This file should have individuals on the rows (listing individual IDs under the 'f.eid' column) and phenotypes on the columns. This data should come directly from the UK Biobank. Row one should contain the column names (as the header)

--bfiles_directory /path/to/bfiles/directory

This is the directory where the BIM/BED/FAM files are located. This data should come directly from the UK Biobank. If the data is not in this format (for example BGEN), it can be converted using the plink2 command line tool using the follow script for each bgen file for each chromosome:

```bash
plink2 --bgen /path/to/ukb_chr$1.bgen ref-last \
			--sample /path/to/ukb_chr$1.sample \
			--threads THREADS \
			--memory MEMORY \
			--make-bed \
			--out /path/to/output/with/prefix/chr$1
```

--covariates_file /path/to/covariates/file

This file should be a tab separated file with a header row with column one being "Field" and the second column being "Description". The Field column should contain fields from the [UKBB Data Showcase](https://biobank.ndph.ox.ac.uk/showcase/search.cgi). The Description column should contain the description of the field. This file is used as covariates for the PheWAS analysis. Ensure that your UKB phenotype file has these fields as columns in the data. Listed below are some common covariates that can be used:

```
Field	Description
31	Sex
34	Year of Birth
52	Month of Birth
53	Date of attending assessment center
54	Assessment Center
21000	Ethnic background 
22001	Genetic Sex 
22009	Genetic Principal Components
```

--sample_file /path/to/sample/file

This file should contain the sample IDs of individuals that you want to include or exclude from the PheWAS analysis (if you choose to include/exclude). This file should have no header row, simply should be a list of individual IDs that match the individual IDs in the phenotype file. For example:

```	
3816366
5752377
2272489
1983956
1646136
2343348
4934200
5778231
```

--variant_file /path/to/variant/file

This file should contain the genetic variants that you want to use as the independent variables in your analysis. This file should be tab separated with a header row with the following columns: rsID, CHR, POS, Non_Effect, Effect, BETA, SE, P. Here is an example of what the file should look like (Effect ==> Effect Allele and Non_Effect ==> Non Effect Allele):

```
rsID	CHR	POS	Non_Effect	Effect	BETA	SE	P
rs2042995	2	179558366	C	T	0.172193	0.0242183	1.2e-12
rs56078722	3	124420735	T	C	-0.161154	0.0279223	7.9e-09
rs572118015	3	196527497	C	T	2.21501	0.355689	4.7e-10
rs1892252	6	25772639	C	G	0.15306	0.0277552	3.5e-08
rs149897	6	28006650	G	A	-0.123681	0.022055	2e-08
rs111846618	7	73554836	C	T	-0.216041	0.0370781	5.7e-09
rs55963900	7	120977734	A	G	-0.13178	0.022519	4.9e-09
rs76871782	7	145200948	G	A	8.73297	1.47709	3.4e-09
rs1991860	8	75737818	T	C	-0.16939	0.0206412	2.3e-16
rs2202	10	30168699	T	C	-0.129675	0.0217878	2.7e-09
rs183837045	19	41024118	T	C	-0.554014	0.0942024	4.1e-09

```

--phenotypes_files /path/to/phenotypes/file

This file is used when running a PheWAS with phenotypes as both the independent variables and as the covariates. Currently these phenotypes must be from the UK Biobank. This file should have no header row, simply should be a list of phenotype data fields from the UK Biobank Data Showcase. For example:

```
90012
90013
21001
```

*Note that the argument --phenotypes_list can also be specified, where the same information is provided in a space separated list. For example:

```
--phenotypes_list 90012 90013 21001
```

### Visualization 

--phewas_results /path/to/phewas/results

This file should be the output of the annotate_categories.py file. This file should be tab separated with a header row with the following columns (there are some optional columns in here depending on the PheWAS analysis run): PheWAS_Category, Ethnicity (optional), Description, Data_Field, Independent_Var, Samples, -log(p), p-val, beta, std_err, Category, Gene (optional).

```
PheWAS_Category	Ethnicity	Description	Data_Field	Independent_Var	Samples	-log(p)	p-val	beta	std_error	Category	Gene
test	White	Sitting height	f.20015.0.0	rs6426192_2	910/969	114.17750887429588	6.644940946154206e-115	33.75599739085548	0.7741746311755939	Body size measures	DDX11L1, DDX11L17
```

--mapping /path/to/mapping/file

This is an optional tab separated file which contains a custom user mapping from the categories for each phenotype grabbed from the UK Biobank Data Showcase to a custom category. This file should have a header row with the following columns: Old_categories,	New_categories. An example of what the file should look like is shown below(note that the old categories in this example might have more than one word, but overall there are only two columns of values) and can also be seen at [category_mappings.txt](pype/cached/category_mappings.txt)

```
Old_categories		New_categories
Summary Maternity	Maternity
Stroke outcomes		Medical Conditions
Death register		Death
Infectious Disease Antigens	Disease
Blood biochemistry	Circulating biochemistry
Blood count	Blood Parameters
```

--variant_files /path/to/variant/file

These files should be the ones used in the PheWAS analysis and thus have the same format.

--gene_file /path/to/gene/file

This file should be a tab separated file (from the UCSC genome browser) with a header row with the following columns: #chrom, cdsStart, cdsEnd, name2. An example of what the file should look like is shown below, and can also be seen at [hg19_gene_list](pype/cached/hg19_gene_list):

```
#chrom	cdsStart	cdsEnd	name2
chr1	67000041	67208778	SGIP1
chr1	67000041	67208778	SGIP1
chr1	67000041	67208778	SGIP1
```

### Mendelian Randomization

--exposure_variants / --outcome_variants 

These files should be the ones used in the PheWAS and thus have the same format.

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
