import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import textwrap 
from collections import defaultdict
import os
from adjustText import adjust_text
from matplotlib.patches import Patch
import matplotlib.ticker as mtick
from matplotlib.lines import Line2D
from utility_funcs import multiple_testing_correction
import shutil

def plot_significant_categories(idx_sig_cats, pheno, out_file, title, regressions, sig_thresh, annotated_regressions, N, cmp_orig_betas, transparency, color_map):
	# for each category in the significant categories
	for cat, _ in idx_sig_cats:
		
		if pheno == '': # if we are plotting aggregate results
			output_file_cat = out_file.split('.')[0] + '_' + cat + '.png'
		else: # if we are plotting results for a category within the results for a certain phenotype
			output_file_cat = out_file.split('.')[0] + '_' + cat + '_' + pheno + '.png'

		title_cat = title + ' - ' + cat
		data_cat = regressions.loc[regressions['Category'] == cat]

		manhattan_plot(data_cat, sig_thresh, 0, 1, title_cat, output_file_cat, cat, annotated_regressions, plt_top_cat = False, N = N, compare_orig_betas = cmp_orig_betas, transparency = transparency, pheno_color = color_map[cat])

def manhattan_plot(regressions, sig, low, high, title, output_file, pheno, annotated_regressions, plt_top_cat, N, compare_orig_betas, transparency, pheno_color = None):
	"""
	Function to plot manhattan plot for a given phenotype

	:param regressions: pandas dataframe containing the regression results
	:param thresh: p-value significance threshold
	:param title: the title of the file for the manhattan plot
	:param output_file: name of the outputfile along with extension
	
	:type regressions: pandas DataFrame
	:type significance level: float
	:type title: string
	:type output_file: string
	"""

	# save the top variants to a file (before filtering out the super low p-values)
	all_top = regressions.loc[regressions['"-log(p)"'] >= sig].groupby('Category').apply(lambda x : x.sort_values(by = '"-log(p)"', ascending = False).reset_index(drop = True))
	all_top.to_csv(output_file.split('.')[0] + '_top_variants.tsv', sep = '\t', index = False)

	# remove extreme outliers in order to be able to visualize the data more coherently
	regressions = regressions[regressions['"-log(p)"'].between(regressions['"-log(p)"'].quantile(low), regressions['"-log(p)"'].quantile(high))]
	
	# remove all variants with p-value above the significance level
	# groupby category and find median of all these significant variants
	# reset the index and sort by the log pvalues so largest will be at the right end of the plot
	medians = regressions.loc[regressions['"-log(p)"'] >= sig].groupby('Category')['"-log(p)"'].mean().reset_index().sort_values(by='"-log(p)"').reset_index(drop=True)

	# medians = significant_vars.groupby('Category')['"-log(p)"'].nlargest(1).reset_index(level=1, drop=True).sort_values()
	# print(medians.index)

	sig_cats = medians['Category'].values

	# don't forget to readd the categories that have no significant variants
	# calculate their median pvalues and add them to the dataframe
	# sorting by the log pvalues so largest will still be at the right end of the plot
	for cat in regressions['Category'].unique():
		if cat not in medians['Category'].values:
			cat_median = regressions.groupby('Category')['"-log(p)"'].mean()[cat]
			medians = medians.append(pd.DataFrame({'Category': cat, '"-log(p)"': cat_median}, index = [0]))
	
	medians = medians.sort_values(by='"-log(p)"').reset_index(drop=True)

	# based on the category ordering, sort the dataframe rows accordingly
	sorter = lambda x: x.map({cat:order for order, cat in enumerate(medians['Category'])})
	
	if pheno == '':
		print('All')
	else:
		print(pheno.capitalize())

	print(medians)
	print('\n')
	
	indices_of_significant_cats = [(cat, medians['Category'].tolist().index(cat)) for cat in sig_cats]

	# sort the categories by the med list
	regressions = regressions.sort_values(by='Category', key = sorter)

	# get the list of categories, to map the colors to the categories
	categories = regressions['Category'].unique()
	n_categories = len(categories)

	color_map = {}
	if n_categories > 1:
		# List of RGB triplets
		colors = sns.color_palette("rainbow", n_categories)

		# Map label to RGB
		color_map = dict(zip(categories, colors))
	else:
		color_map[pheno] = pheno_color


	if plt_top_cat:
		plot_significant_categories(indices_of_significant_cats, pheno, output_file, title, regressions, sig, annotated_regressions, N, compare_orig_betas, transparency, color_map)
	if regressions.empty:
		print('No significant variants found for this phenotype')
		return

	# rename the columns to be more readable
	regressions.rename(columns={'"-log(p)"': '-log(p)'}, inplace=True)
	
	# get ready to plot
	plt.figure(figsize=(15, 8), dpi=500)

	# data for plotting
	significant_vars = regressions.loc[regressions['-log(p)'] >= sig]
	non_significant_vars = regressions.loc[regressions['-log(p)'] < sig]
	significant_original_pos_dir = significant_vars.loc[significant_vars['beta'] >= 0].copy()
	significant_original_neg_dir = significant_vars.loc[significant_vars['beta'] < 0].copy()

	# the index for the collection for each category is simply the position of the category in the list of categories
	# when you use sns.stripplot with order = something, it will plot for all the categories, creating empty collections
	# for those that don't exist, thats why we add the initial index, then index + num categories
	# Note: if there are no significant positive, or negative beta variants, we need to check if our index is out of bounds
	cat_to_collection = defaultdict(list)
	for idx, cat in enumerate(categories):
		cat_to_collection[cat].append(idx) # accounts for the significant positive beta plots
		cat_to_collection[cat].append(idx + n_categories) # accounts for the significant negative beta plots

	# custom legend entries
	handles = []

	# if we are plotting 1 category, we have bigger points + different colors
	size = 7
	color = pheno_color

	# if we plotting more than 1 category, save space with points and the colors of the points will be different
	if n_categories > 1: 
		size = 5
		color = color

	if not significant_original_pos_dir.empty:
		if compare_orig_betas and n_categories == 1:
			pos_dir_palette ={"+": "red", "-": "orange"}
			significant_original_pos_dir.loc[:, 'original_beta_direction'] = significant_original_pos_dir.apply(lambda x: '+' if x['Original_beta'] >= 0 else '-', axis = 1)
			ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_original_pos_dir, hue="original_beta_direction", palette = pos_dir_palette, jitter=0.45, size = size, order = categories, linewidth=0.3, **{'marker': '^', 'alpha': transparency})
			handles.extend([Line2D([0], [0], color = 'red', marker = '^', linestyle='None',label = 'Pos/Pos', alpha = transparency),
							Line2D([0], [0], color = 'orange', marker = '^', linestyle='None',label = 'Pos/Neg', alpha = transparency)])
		else:
			# Plot the -log(p) values against the category values, with the colors mapped to the categories, with up arrow == the direction of the beta value is positive
			ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_original_pos_dir, palette = color_map, jitter=0.45, size = size, order = categories, linewidth=0.2, **{'marker': '^', 'alpha': transparency})
			handles.append(Line2D([0], [0], color = color, marker = '^', linestyle='None',label = 'Positive beta', alpha = transparency))
	
	if not significant_original_neg_dir.empty:
		if compare_orig_betas and n_categories == 1:
			neg_dir_palette ={"-": "blue", "+": "purple"}
			significant_original_neg_dir.loc[:, 'original_beta_direction'] = significant_original_neg_dir.apply(lambda x: '+' if x['Original_beta'] >= 0 else '-', axis = 1)
			ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_original_neg_dir, hue="original_beta_direction", palette = neg_dir_palette, jitter=0.45, size = size, order = categories, linewidth=0.3, **{'marker': 'v', 'alpha': transparency})
			handles.extend([Line2D([0], [0], color = 'blue', marker = 'v', linestyle='None',label = 'Neg/Neg', alpha = transparency),  
							Line2D([0], [0], color = 'purple', marker = 'v', linestyle='None', label = 'Neg/Pos', alpha = transparency)])
		else:
			# Plot the -log(p) values against the category values, with the colors mapped to the categories, with down arrow == the direction of the beta value is negative
			ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_original_neg_dir, palette = color_map, jitter=0.45, size = size, order = categories, linewidth=0.2, **{'marker': 'v', 'alpha': transparency})
			handles.append(Line2D([0], [0], color = color, marker = 'v', linestyle='None', label = 'Negative beta', alpha = transparency))
	
	# if we need to add annotations to the plot
	if annotated_regressions is not None:
		
		# if we are plotting more than 1 category, only annotate top 2 points
		annotate_top = 2

		# else, annotate the top N points
		if n_categories == 1:
			annotate_top = N

		# get top N variants for each category
		top_N = significant_vars.groupby('Category').apply(lambda x : x.sort_values(by = '-log(p)', ascending = False).head(annotate_top).reset_index(drop = True))
		
		collections = ax.collections

		# when we only have 1 category, there are some erroneous empty collections, so get rid of them
		# need to figure out why this is happening and make a better solution
		if n_categories == 1:
			collections = [ax.collections[i] for i in range(len(ax.collections)) if ax.collections[i].get_offsets().data.size != 0]

		texts = []
		for cat, _ in indices_of_significant_cats:

			# get the corresponding collections for the plot of this category
			collections_for_cat = []
			for idx in cat_to_collection[cat]:
				if idx < len(collections):

					collections_for_cat.append(collections[idx].get_offsets().data)

			collections_for_cat = np.concatenate(collections_for_cat)

			# only annotate top variants for each category (to avoid cluttering the plot, use only 2 for the big plots)
			top_xy_coords = sorted(collections_for_cat, key = lambda x: x[1], reverse = True)

			# # I DONT THINK WE NEED THIS BECAUSE OF HOW I CHANGED THE CODE -- CHECK IF WE CAN REMOVE
			# # choose all the points that are significant for this category
			# top_xy_coords = [(x,y) for x,y in top_xy_coords if y >= sig]

			# if we want to annotate less points than those that exist, only choose the top X points
			if annotate_top < len(top_xy_coords):
				top_xy_coords = top_xy_coords[:annotate_top]

			index = 0
			top = top_N.loc[top_N['Category'] == cat]
			for xy_coords in top_xy_coords:
				x,y = xy_coords
				
				# get rsID from the top variants dataframe
				rsID = top.iloc[index]['rsID'].split('_')[0]

				# only annotate first 3 genes, the rest can be found in the results file (for clarity)
				genes = ', '.join(annotated_regressions.loc[annotated_regressions['rsID'] == rsID, 'GENE'].unique().tolist()[:3])
				genes = textwrap.fill(genes, 25, break_long_words=False)

				# if we have more than 1 category, don't add the description to the plot
				if n_categories > 1:
					annotation = (rsID + ', ' + genes).strip()

					texts.append(ax.text(x, y, annotation, ha = 'center', va = 'center', fontsize = 3))

				else:
					desc = top.iloc[index]['Description']
					desc = textwrap.fill(desc, 35, break_long_words=False)
					annotation = (rsID + ',' + genes + '\n' + desc ).strip()
					
					texts.append(ax.text(x, y, annotation, fontsize = 6, ha = 'center', va = 'center'))

				index +=1
	
	if not non_significant_vars.empty:
		# plot for the values below the significance level (should be circles)
		ax = sns.stripplot(x = 'Category', y = '-log(p)', data = non_significant_vars, palette = color_map, jitter=0.45, size = size, order = categories, linewidth=0.2, **{'alpha': transparency})

	# add the legends (we do it here to avoid adding more collections to ax.collections - simplifies logic)
	if len(handles) > 0:
		ax.legend(handles = handles, loc = 'best')

	plt.title(title)
	plt.axhline(y=sig, color='red', ls='--', lw = 0.5)  # plot threshold

	# make it so that the labels wrap to the next line if they are too long
	labels = [ textwrap.fill(l, 12, break_long_words=False) for l in categories ]
	ax.set_xticks(ax.get_xticks().tolist())
	ax.set_xticklabels(labels)


	if n_categories > 5:
		# and rotate/shrink the labels
		plt.xticks(rotation = 90)
		ax.tick_params(axis='x', which='major', labelsize=5)

	ax.set_xmargin(0.1)
	# ax.set_ymargin(0.05)
	ax.autoscale_view()

	ax.set_ylim([None, ax.get_ylim()[1]*1.15])

	if annotated_regressions is not None:
		adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=0.2))

	else:
		output_file = output_file.split('.')[0] + '_no_annotations.png'

	# finally save the figure
	plt.savefig(output_file, bbox_inches ='tight', dpi = 300)
	plt.close()

def significant_bar_plot(regressions, sig,title, output_file, save = True):
	
	print('Plotting bar plot for %s' % output_file)
	regressions['Significant'] = regressions['"-log(p)"'] >= sig
	_, ax = plt.subplots(figsize = (12, 8))

	total_vars = {}
	significant_vars = {}

	for cat in regressions['Category'].unique():
		total_vars[cat] = regressions[regressions['Category'] == cat].shape[0]
		significant_vars[cat] = regressions.loc[regressions['Category'] == cat, 'Significant'].sum()

	categories = dict(sorted(significant_vars.items(), key=lambda item: item[1]/total_vars[item[0]])).keys()
	significant = [100*significant_vars[cat]/total_vars[cat] for cat in categories]

	p1 = ax.bar(categories, significant,  color='tab:blue', label = 'Significant associations')

	plt.xlabel('Categories', labelpad=12)
	ax.set_ylim([None, ax.get_ylim()[1]*1.05])
	plt.xticks(rotation=45, ha='right')
	plt.ylabel('Percentage of significant associations out of total', labelpad=12)

	ax.yaxis.set_major_formatter(mtick.PercentFormatter())

	index = 0
	categories = list(categories)
	sig_max = np.max(significant)
	for p in p1.patches:
		width, height = p.get_width(), p.get_height()
		x, y = p.get_xy() 
		width_addition = 0.2*width
		height_addition = sig_max*0.02

		if total_vars[categories[index]] < 100:
			width_addition = 0.3*width
		elif total_vars[categories[index]] > 1000:
			width_addition = 0.05*width

		ax.annotate(total_vars[categories[index]], (p.get_x() + width_addition, p.get_y()+height + height_addition), 
			fontsize=10, color='black', bbox=dict(facecolor='none', edgecolor='red', boxstyle='round'))
		index += 1
	
	legend_elements = [Patch(facecolor='none', edgecolor='red',label='Total associations')]
	
	ax.legend(handles=legend_elements, loc='upper left')
	# plt.yscale("log")
	plt.tight_layout()
	plt.title(title)
	
	if save:
		plt.savefig(output_file, bbox_inches ='tight', dpi = 300)
		plt.close()

# def volcano_plot(regressions, sig, title, output_file, save = True):

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

def annotate_genes(gene_file, rsid_df, down, up, regressions, output_dir, pheno_name):

	rsid_df['CHR'] = rsid_df['CHR'].astype(str)
	#read in the gene file
	gene_df = pd.read_csv(gene_file, sep = '\t')
	
	# remove all the chr prefixes from the chromosome columns
	gene_df['#chrom'] = gene_df['#chrom'].apply(lambda x: x.replace('chr', ''))
	
	# only retain the the chromosomes with properly formatted chromosome names
	chroms = [str(x) for x in list(range(23)) + ['X', 'Y', 'XY']]
	gene_df = gene_df[gene_df['#chrom'].isin(chroms)]

	# get closest genes to eaach variant
	res = get_closest_genes(rsid_df, gene_df, down, up).sort_values(by = 'rsID')

	gene_map = defaultdict(list)

	for _, row in res.iterrows():
		gene_map[row['rsID']].append(row['GENE'])

	regressions['Gene'] = regressions.apply(lambda x: ', '.join(gene_map[x['rsID']]), axis = 1)
	regressions.to_csv(output_dir + '/' + pheno_name + '_pheWAS_results_with_nearby_genes.tab', sep = '\t', index = False)

	return regressions, res

def main():
	parser = argparse.ArgumentParser(description='Plot significant phenotype data on a Manhattan Plot.')
	parser.add_argument('--input', help='Input file', required=True)
	parser.add_argument('--output', help='Output png', required=True)
	parser.add_argument('--phenotype_name', help='Phenotype name', required=True)
	parser.add_argument('--aggregate_title', help='Title of the plot of the aggregated data', required=False, default = None, nargs='+')
	parser.add_argument('--mapping', help='Mapping of old categories to new ones', required=False, default = None)
	parser.add_argument('--lower_outlier', help='Lower outlier threshold', required=False, default = 0, type = float)
	parser.add_argument('--upper_outlier', help='Upper outlier threshold', required=False, default = 0.995, type = float)
	parser.add_argument('--seed', help = 'Set seed to make the adjustText library deterministic', required = False, default = 0, type = int)
	parser.add_argument('--downstream', help = 'Number in KB to check for closest genes', required = False, default = 40, type = int)
	parser.add_argument('--upstream', help = 'Number in KB to check for closest genes', required = False, default = 40, type = int)
	parser.add_argument('--gene_file', help = 'File containing a list of genes use for finding closest genes', required = False, default = None)
	parser.add_argument('--rsid_files', help = 'Files containing a list of rsIDs to use for finding closest genes', required = True, action = 'append')
	parser.add_argument('--no_annotations', help = 'Don\'t add annotations', required = False, default = False, action = 'store_true')
	parser.add_argument('--plot_top_categories', help = 'Plot the top categories as well', required = False, default = False, action = 'store_true')
	parser.add_argument('--number_of_top_results', help = 'Number of top results to save', required = False, default = 10, type = int)
	parser.add_argument('--plot_manhattan', help = 'Plot the manhattan plot', required = False, default = False, action = 'store_true')
	parser.add_argument('--plot_bar', help = 'Plot the bar plot', required = False, default = False, action = 'store_true')
	parser.add_argument('--clear_old_files', help = 'Clear old files', required = False, default = False, action = 'store_true')
	parser.add_argument('--alpha', help = 'Significance threshold', required = False, default = 0.05, type = float)
	parser.add_argument('--correction', help = 'Correction method', required = False, default = 'bonferroni', choices = ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky', 'no_correction'])
	parser.add_argument('--compare_original_betas', help = 'Compare original betas with corrected betas', required = False, default = True, action = 'store_true')
	parser.add_argument('--transparency', help = 'Transparency of points', required = False, default = 0.75, type = float)

	args = parser.parse_args()
	input_file = args.input
	output_file = args.output
	phenotype_name = args.phenotype_name
	aggregate_title = args.aggregate_title
	mapping = args.mapping
	lower_outlier = args.lower_outlier
	upper_outlier = args.upper_outlier
	seed = args.seed
	downstream = args.downstream
	upstream = args.upstream
	gene_file = args.gene_file
	rsid_files = args.rsid_files
	no_annotations = args.no_annotations
	plt_top_cat = args.plot_top_categories
	N = args.number_of_top_results
	plt_manhattan = args.plot_manhattan
	plt_bar = args.plot_bar
	clear_old_files = args.clear_old_files
	alpha = args.alpha
	correction = args.correction
	compare_original_betas = args.compare_original_betas
	transparency = args.transparency

	# set a seed to make results reproducible (for the adjustText library)
	np.random.seed(seed)

	# Load data
	regressions = pd.read_csv(input_file, sep='\t')

	# get the directory of the output file
	output_dir = os.path.dirname(output_file)
	new_output_dir = output_dir + '/All' 
	output_file = new_output_dir + '/' + os.path.basename(output_file)

	# make directory for the aggregated results 
	if clear_old_files:
		shutil.rmtree(new_output_dir, ignore_errors=True)
	
	if not os.path.exists(new_output_dir):
		os.mkdir(new_output_dir)

	if aggregate_title is None:
		# if no title is provided, name the file based on the type of plot
		if plt_manhattan:
			aggregate_title = 'PheWAS Results for %s (ALL) Variants' % (phenotype_name.capitalize())
			output_file = output_file.split('.')[0] + '_manhattan.png'
		else:
			aggregate_title = '%s (ALL) Variant Enrichment per Category' % (phenotype_name.capitalize())
			output_file = output_file.split('.')[0] + '_bar.png'
	else:
		aggregate_title = ' '.join(aggregate_title)

	# make directory for each of the phenotypes
	phenos = regressions['Predictor'].unique()
	pheno_map = {}
	output_ext = '.' + output_file.rsplit('.', 1)[-1]

	if len(phenos) == 1:
		for p in phenos:
			new_output_dir_for_pheno = output_dir + '/' + p

			if clear_old_files:
				shutil.rmtree(new_output_dir_for_pheno, ignore_errors=True)

			if not os.path.exists(new_output_dir_for_pheno):
				os.mkdir(new_output_dir_for_pheno)

			title_and_output = ''
			if plt_manhattan:
				title_and_output = ('PheWAS Results for %s Variants' % (p.capitalize()), new_output_dir_for_pheno + '/'  + p.lower() + '_manhattan' + output_ext)
			else:
				title_and_output = ('%s Variant Enrichment per Category' % (p.capitalize()), new_output_dir_for_pheno + '/' + p.lower() + '_bar' + output_ext)

			pheno_map[p] = title_and_output

	if mapping is not None:
		# load mapping file - two column file, first with old categories, second with new categories
		mapping = dict(pd.read_csv(mapping, sep='\t').values)

		# map old categories to new ones
		regressions['Category'] = regressions['Category'].map(mapping)

	# create a aggregate dataframe of all the snps
	rsid_df = pd.DataFrame()

	for f in rsid_files:
		# read f and then concat to rsid_df
		rsid_df = pd.concat([rsid_df, pd.read_csv(f, sep = '\t')])
		rsid_df = rsid_df.drop_duplicates()
	
	# case where we have a variant like 3:100928901_CTT_C
	rsid_df['rsID'] = rsid_df['rsID'].apply(lambda x: x.split('_')[0])
	regressions['rsID'] = regressions['rsID'].apply(lambda x: x.split('_')[0])

	regressions['Original_beta'] = regressions.apply(lambda x: rsid_df.loc[rsid_df['rsID'] == x['rsID']]['BETA'].unique()[0], axis = 1)
	regressions['Original_pval'] = regressions.apply(lambda x: rsid_df.loc[rsid_df['rsID'] == x['rsID']]['P'].unique()[0], axis = 1)
	regressions['Direction'] = np.where(regressions['beta'] * regressions['Original_beta'] > 0, 'S', 'D')
	
	annotated_regressions = None
	if not no_annotations:
		if gene_file is not None:
			regressions, annotated_regressions = annotate_genes(gene_file = gene_file, 
											rsid_df = rsid_df, 
											down = downstream, 
											up = upstream, 
											regressions = regressions, 
											output_dir = output_dir, 
											pheno_name = phenotype_name)
		else:
			print('No gene file provided.')
			exit()

	if no_annotations:
		annotated_regressions = None	

	# correct the pvalues using the method of choice
	print('Correcting pvalues using %s correction' % (correction))
	sig = -np.log10(multiple_testing_correction(pvalues = regressions['p-val'].dropna(), alpha = alpha, method = correction))

	if plt_manhattan:
		manhattan_plot(regressions, sig, lower_outlier, upper_outlier, aggregate_title, output_file, '', annotated_regressions, plt_top_cat, N, compare_original_betas, transparency)

	if plt_bar:
		outlier_removed = regressions.copy(deep = True)
		significant_bar_plot(outlier_removed, sig, aggregate_title, output_file, True)

	if len(phenos) > 1:
		for pheno, (pheno_title, pheno_output_file) in pheno_map.items():
			# plot the data for each phenotype
			predictor_specific = regressions[regressions['Predictor'] == pheno].copy(deep = True)

			if plt_manhattan:
				manhattan_plot(predictor_specific, sig, lower_outlier, upper_outlier, pheno_title, pheno_output_file, pheno, annotated_regressions, plt_top_cat, N, compare_original_betas, transparency)
			
			if plt_bar:

				significant_bar_plot(predictor_specific, sig, pheno_title, pheno_output_file, True)


if __name__ == '__main__':
	main()