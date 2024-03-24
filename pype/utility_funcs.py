import pickle
import mygene
import pprint
import myvariant
import numpy as np
import pandas as pd
from collections import defaultdict

def pickle_object(obj, filename, protocol=pickle.HIGHEST_PROTOCOL) -> None:
	"""
	Pickle an object to a file.
	
	Parameters
	----------
	obj : object
		Object to pickle.
	filename : str	
		Filename of the pickle file.
	protocol : int
		Protocol to use for pickling.

	Returns	
	-------
	None
	"""

	with open(filename, 'wb') as f:
		pickle.dump(obj, f, protocol = protocol)
	return

def unpickle_object(filename) -> object:
	"""
	Unpickle an object from a file.

	Parameters
	----------
	filename : str
		Filename of the pickle file.
	
	Returns
	-------
	obj : object
		Object that was pickled.
	"""
	with open(filename, 'rb') as f:
		obj = pickle.load(f)
	return obj

def multiple_testing_correction(pvalues, alpha, method) -> float:
	"""
	Perform multiple testing correction on pvalues, returning threshold to use.

	Parameters
	----------
	pvalues : array_like
		Array of pvalues.
	alpha : float
		Alpha value for the test.
	method : str
		Method for multiple testing correction.
		'bonferroni' for Bonferroni correction.
		'sidak' for Sidak correction.
		'fdr_bh' for Benjamini-Hochberg correction.
		'no_correction' for no correction.

	Returns
	-------
	threshold : float
		Threshold value for the test.
	"""

	pvalues = np.array(pvalues)
	n_pvalues = len(pvalues)

	if method == 'bonferroni':
		threshold = alpha / n_pvalues
	elif method == 'sidak':
		threshold = 1 - np.power(1. - alpha, 1. / n_pvalues)
	elif method == 'fdr_bh':
		# sort the pvalues in ascending order
		pvalues.sort()

		for rank, pval in enumerate(pvalues, 1):

			# benjamini-holchberg critical value
			critical_val = alpha * rank / n_pvalues

			# the largest pvaue that is less than the critical value
			if pval > critical_val:
				threshold = pval # the significant pvalues should be strictly less than this
				break
	elif method == 'no_correction':
		threshold = alpha
	else:
		raise ValueError('Invalid method for multiple testing correction.')
	
	print('Threshold for {} multiple testing correction: {}'.format(method, threshold))

	return threshold


def get_closest_genes(rsIDs, genes, upstream, downstream):

	# assume gene file always have 4 columns CHR, START, END, GENE
	genes.columns = ['CHR', 'START', 'END', 'GENE']
	
	# remove the 'chr' prefix if present
	genes['CHR'] = genes['CHR'].apply(lambda x: x.replace('chr', ''))
	
	# create dataframe to store results
	res = pd.DataFrame(columns = ['CHR', 'START', 'END', 'GENE'])
	
	# group by the gene and the associated CHR (seems to be some instances of gene isoforms on multiple chromosomes)
	gene_groups = genes.groupby(['GENE', 'CHR'])
	
	# for each grouping of gene, and then chromosome, get the tuple of the gene/chr, and save it in the results dataframe
	gene_chr_tuples = gene_groups.apply(lambda x: x.name).reset_index(drop = True)
	res['GENE'], res['CHR'] = gene_chr_tuples.apply(lambda x:x[0]).tolist(), gene_chr_tuples.apply(lambda x:x[1]).tolist()

	# we will combine all the isoforms of the genes on the same chromosome into a "single gene"
	# for the start, take the minimum of all the isoforms starting spot and for the end, take the maximum of all the isoforms ending spot
	res['START'] = gene_groups['START'].min().values.tolist()
	res['END'] = gene_groups['END'].max().values.tolist()
	
	# merge the rsIDs with the results dataframe based on the chromosome
	merge = pd.merge(rsIDs, res, how = 'inner', on = 'CHR')
	
	# get the rows that are within the upstream and downstream regions
	return merge.loc[(merge.POS >= merge.START - downstream*1000) & (merge.POS <= merge.END + upstream*1000)] 

def annotate_genes(gene_file, rsid_df, down, up, regressions, output_dir, pheno_name, save = True):

	rsid_df['CHR'] = rsid_df['CHR'].astype(str)
	
	#read in the gene file
	gene_df = pd.read_csv(gene_file, sep = '\t')
	
	# remove all the chr prefixes from the chromosome columns
	gene_df['#chrom'] = gene_df['#chrom'].apply(lambda x: x.replace('chr', ''))
	
	# only retain the the chromosomes with properly formatted chromosome names
	chroms = [str(x) for x in list(range(23)) + ['X', 'Y', 'XY']]
	gene_df = gene_df[gene_df['#chrom'].isin(chroms)]

	print('Annotating variants with nearby genes...')
	# get closest genes to each variant
	res = get_closest_genes(rsid_df, gene_df, down, up).sort_values(by = 'rsID')

	print('Annotating variants with nearby genes...done')

	gene_map = defaultdict(list)

	for _, row in res.iterrows():
		gene_map[row['rsID']].append(row['GENE'])

	regressions['Gene'] = regressions.apply(lambda x: ', '.join(gene_map[x['rsID']]), axis = 1)
	
	if save:
		regressions.to_csv(output_dir + '/' + pheno_name + '_pheWAS_results_with_nearby_genes.tab', sep = '\t', index = False)
		res.to_csv(output_dir + '/' + pheno_name + '_rsID_gene_map.tab', sep = '\t', index = False)

	return regressions, res

def annotateVariantsAndGenes(top_results, variant_fields, gene_fields, out_dir):
	
	mv = myvariant.MyVariantInfo()
	mg = mygene.MyGeneInfo()

	# for all variants in the Independent_Var column, strip everything after the first _
	top_results['Independent_Var'] = top_results['Independent_Var'].str.split('_').str[0]

	rsIDs = top_results['Independent_Var'].unique().tolist()
	for rsID in rsIDs:
		print('Finding information online about variant {} and any nearby genes'.format(rsID))
		# write the pretty print jsons to a file
		with open(out_dir + '/' + rsID + '.summary', 'w') as outfile:
			# write the rsID to the file
			outfile.write('The variant is: ' + rsID + '\n\n')

			try:
				# for this variant, get all the rows in the top_results dataframe, and then grab the genes
				genes = top_results[top_results['Independent_Var'] == rsID]['Gene'].unique()[0]
			except:
				continue # get here if gene_file was not provided and therefore there is no Gene column

			if not pd.isnull(genes):

				outfile.write('The genes found associated with this variant are: ' + genes + '\n\n')
				
				variant_json = mv.getvariant(rsID, fields = variant_fields)

				if variant_json is not None:
					outfile.write('Here is some information about the variant that you requested:\n\n')
					pprint.pprint(variant_json, stream = outfile)
				else:
					outfile.write('Unfortunately, no information could be found in the specified databases for this variant.\n')

				genes = genes.split(', ')

				for gene in genes:

					gene_json = mg.query(gene, field = gene_fields)

					if gene_json is not None:

						# Extract the top hit
						top_hit = gene_json['hits'][0]

						# Print the result
						if 'summary' in top_hit:
							outfile.write('\nHere is some information about gene {}, which was found to be close to this variant:\n\n'.format(gene))
							outfile.write(top_hit['summary'] + '\n')
						else:
							outfile.write('\nUnfortunately, no summary information could be found for gene {}.'.format(gene))
					else:
						outfile.write('\nUnfortunately, no summary information could be found for gene {}.'.format(gene))
			else:
				variant_json = mv.getvariant(rsID, fields = variant_fields)

				if variant_json is not None:
					outfile.write('Here is some information about the variant that you requested:\n\n')
					pprint.pprint(variant_json, stream = outfile)
				else:
					outfile.write('Unfortunately, no information could be found in the specified databases for this variant.\n')

				outfile.write('\nUnfortunately, no genes were found to be close to this variant.')