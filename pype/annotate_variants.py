import pandas as pd
import argparse
import utility_funcs

def main():

	parser = argparse.ArgumentParser(description = 'Annotate variants and genes from a PheWAS significant SNPs files.')
	
	# Variant/Gene Annotation Arguments
	parser.add_argument('--significant_phewas_results', help = 'The significant PheWAS results file from run_phewas.py', required = True, type = str)
	parser.add_argument('--variant_fields', help = 'List of variant databases and fields to grab. The full list of choices can be found at: https://docs.myvariant.info/en/latest/doc/data.html#available-fields', required = False, default = ['dbsnp'], nargs = '+')
	parser.add_argument('--gene_fields', help = 'List of gene databases and fields to grab. The full list of choices can be found at: https://docs.mygene.info/en/latest/doc/data.html#available-fields', required = False, default = ['summary'], nargs = '+')
	parser.add_argument('--out_dir', help = 'The directory to write the output files to', required = True, type = str)
	
	# SNP-GENE Mapping Arguments
	parser.add_argument('--gene_file', help = 'File containing list of genes and their positions. Used to annotate variants with nearby genes', required = False, default = None, type = str)
	parser.add_argument('--downstream', help = 'Number in KB to check for closest genes', required = False, default = 10, type = int)
	parser.add_argument('--upstream', help = 'Number in KB to check for closest genes', required = False, default = 10, type = int)
	parser.add_argument('--variant_files', help='List of files with genetic variants to use as independent variables', required = False, action='append')
	parser.add_argument('--output_prefix', help = 'Prefix for SNP-Gene annotated output files', required = False, type = str, action = 'append')

	args = parser.parse_args()
	
	# Variant/Gene Annotation Arguments
	significant_phewas_results = args.significant_phewas_results
	variant_fields = args.variant_fields
	gene_fields = args.gene_fields
	out_dir = args.out_dir

	# SNP-GENE Mapping Arguments
	gene_file = args.gene_file
	downstream = args.downstream
	upstream = args.upstream
	variant_files = args.variant_files
	output_prefix = args.output_prefix


	# ---------------------------------------VERIFY ARGS--------------------------------------- #

	if upstream < 0 or downstream < 0:
		print('Error: The upstream and downstream values for annotating the genes must be non-negative.')
		exit()

	if gene_file is not None and len(variant_files) == 0:
		print('Error: You must specify the variants files when performing SNP-Gene mapping.')
		exit()

	if gene_file is not None and len(variant_files) != len(output_prefix):
		print('Error: You must specify a prefix for each variant file.')
		exit()

	# ----------------------------------------------------------------------------------------- #

	# ---------------------------------------RUN ANNOTATIONS--------------------------------------- #
		
	# read in the significant PheWAS results
	top_results = pd.read_csv(significant_phewas_results, sep = '\t')

	if gene_file is not None and 'Gene' not in top_results.columns:
		# create a aggregate dataframe of all the snps
		rsid_df = pd.DataFrame()

		# rename Independent_Var to rsID
		top_results = top_results.rename(columns = {'Independent_Var': 'rsID'})

		for rsid_f in variant_files:
			# read f and then concat to rsid_df
			rsid_df = pd.concat([rsid_df, pd.read_csv(rsid_f, sep = '\t')])
			rsid_df = rsid_df.drop_duplicates()
		
		rsid_df['rsID'] = rsid_df['rsID'].apply(lambda x: x.split('_')[0])
		old_results_rsIDs = top_results['rsID'].tolist()
		top_results['rsID'] = top_results['rsID'].apply(lambda x: x.split('_')[0])
		
		for prefix in output_prefix:
			top_results, _ = annotate_genes(gene_file = gene_file, 
												rsid_df = rsid_df, 
												down = downstream, 
												up = upstream, 
												regressions = top_results, 
												output_dir = out_dir, 
												pheno_name = output_prefix)
		
		top_results = top_results.rename(columns = {'rsID': 'Independent_Var'})
		top_results['Independent_Var'] = old_results_rsIDs

	# perform variant and gene annotation
	utility_funcs.annotateVariantsAndGenes(top_results = top_results, variant_fields = variant_fields, gene_fields = gene_fields, out_dir = out_dir)


if __name__ == '__main__':
	main()