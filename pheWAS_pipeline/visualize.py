from pickle import NONE
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
	
def manhattan_plot(regressions, sig, lower_outlier, upper_outlier, title, output_file, pheno, res, plt_top_cat, N):
	"""
	Plots significant phenotype data on a Manhattan Plot.
	The significance of each phenotype (represented by :math:`-log_{10}(p)`\ ) is plotted along the
	y-axis, with phenotypes plotted along the x-axis.
	
	:param regressions: pandas dataframe containing the regression results
	:param thresh: p-value significance threshold
	:param title: the title of the file for the manhattan plot
	:param output_file: name of the outputfile along with extension
	
	:type regressions: pandas DataFrame
	:type significance level: float
	:type title: string
	:type output_file: string
	"""
	# print(np.isfinite(regressions['"-log(p)"']))

	# # convert the inf pvalues to 
	# regressions['"-log(p)"'] = np.where(np.isfinite(regressions['"-log(p)"']), regressions['"-log(p)"'],  np.ma.masked_invalid(regressions['"-log(p)"']).max() + 10)

	# print(np.isfinite(regressions['"-log(p)"']))
	# exit()
	# # save the top variants to a file (before filtering out the super low p-values)
	all_top = regressions.loc[regressions['"-log(p)"'] >= sig].groupby('Category').apply(lambda x : x.sort_values(by = '"-log(p)"', ascending = False).reset_index(drop = True))
	all_top.to_csv(output_file.split('.')[0] + '_top_variants.tsv', sep = '\t', index = False)

	# remove extreme outliers in order to be able to visualize the data more coherently
	regressions = regressions[regressions['"-log(p)"'].between(regressions['"-log(p)"'].quantile(lower_outlier), regressions['"-log(p)"'].quantile(upper_outlier))]
	
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

	if plt_top_cat:
		for cat, _ in indices_of_significant_cats:
			
			if pheno == '':
				output_file_cat = output_file.split('.')[0] + '_' + cat + '.png'
			else:
				output_file_cat = output_file.split('.')[0] + '_' + cat + '_' + pheno + '.png'

			title_cat = title + ' - ' + cat
			reg_cat = regressions.loc[regressions['Category'] == cat]

			manhattan_plot(reg_cat, sig, 0, 1, title_cat, output_file_cat, cat, res, plt_top_cat = False, N = N)
	
	if regressions.empty:
		return

	# rename the columns to be more readable
	regressions.rename(columns={'"-log(p)"': '-log(p)'}, inplace=True)
	
	# get the list of categories, to map the colors to the categories
	categories = regressions['Category'].unique()
	n_categories = len(categories)

	# List of RGB triplets
	colors = sns.color_palette("rainbow", n_categories)

	# Map label to RGB
	color_map = dict(zip(categories, colors))

	# get the axes and add a horizontal line at the threshold value
	# _, ax = plt.subplots(figsize=(12,8))
	plt.figure(figsize=(15, 8), dpi=500)

	significant_vars = regressions.loc[regressions['-log(p)'] >= sig]

	if n_categories > 1:
		# Plot the -log(p) values against the category values, with the colors mapped to the categories
		ax = sns.stripplot(x = 'Category', y = '-log(p)', data = regressions, palette = color_map, jitter=0.45, size = 3,  linewidth=0.2, **{'alpha': 0.5})
	else:
		non_significant_vars = regressions.loc[regressions['-log(p)'] < sig]
		# plot for the values below the significance level (should be circles)
		ax = sns.stripplot(x = 'Category', y = '-log(p)', data = non_significant_vars, color = 'green', jitter=0.45, size = 7,  linewidth=0.3, **{'alpha': 0.5})


		significant_original_pos_dir = significant_vars.loc[significant_vars['Original_beta'] >= 0].copy()
		significant_original_neg_dir = significant_vars.loc[significant_vars['Original_beta'] < 0].copy()
	

		if not significant_original_pos_dir.empty:
			pos_dir_palette ={"+": "red", "-": "orange"}
			significant_original_pos_dir.loc[:, 'new_beta_direction'] = significant_original_pos_dir.apply(lambda x: '+' if x['beta'] >= 0 else '-', axis = 1)
			ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_original_pos_dir, hue="new_beta_direction", palette = pos_dir_palette, jitter=0.45, size = 7,  linewidth=0.3, **{'marker': '^', 'alpha': 0.5})
		
		if not significant_original_neg_dir.empty:
			neg_dir_palette ={"-": "blue", "+": "purple"}
			significant_original_neg_dir.loc[:, 'new_beta_direction'] = significant_original_neg_dir.apply(lambda x: '+' if x['beta'] >= 0 else '-', axis = 1)
			ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_original_neg_dir, hue="new_beta_direction", palette = neg_dir_palette, jitter=0.45, size = 7,  linewidth=0.3, **{'marker': 'v', 'alpha': 0.5})
	

		# # get the rows of the significant values with the direction being S
		# significant_vars_same_dir = significant_vars.loc[significant_vars['Direction'] == 'S']

		# # get the rows of the significant values with the direction being D
		# significant_vars_diff_dir = significant_vars.loc[significant_vars['Direction'] == 'D']

		# # plot for the values above the significance level (should be ^ or v depending on the direction of the regression)
		# if not significant_vars_same_dir.empty:
		# 	ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_vars_same_dir, color = 'red', jitter=0.45, size = 7,  linewidth=0.3, **{'marker': '^', 'alpha': 0.5})
		
		# if not significant_vars_diff_dir.empty:
		# 	ax = sns.stripplot(x = 'Category', y = '-log(p)', data = significant_vars_diff_dir, color = 'blue', jitter=0.45, size = 7,  linewidth=0.3, **{'marker': 'v', 'alpha': 0.5})
	
		ax.legend(handles = [Line2D([0], [0], color = 'red', marker = '^', linestyle='None',label = 'Pos/Pos', alpha = 0.5),
							Line2D([0], [0], color = 'orange', marker = '^', linestyle='None',label = 'Pos/Neg', alpha = 0.5), 
							Line2D([0], [0], color = 'blue', marker = 'v', linestyle='None',label = 'Neg/Neg', alpha = 0.5),  
							Line2D([0], [0], color = 'purple', marker = 'v', linestyle='None', label = 'Neg/Pos', alpha = 0.5)], loc = 'best')
	
	# # save all top variants to a file
	# all_top = significant_vars.groupby('Category').apply(lambda x : x.sort_values(by = '-log(p)', ascending = False).reset_index(drop = True))
	# all_top.to_csv(output_file.split('.')[0] + '_top_variants.tsv', sep = '\t', index = False)
	
	if res is not None:

		xy = defaultdict(list)

		top = 2

		if n_categories == 1:
			top = N

		# get top N variants for each category
		top_N = significant_vars.groupby('Category').apply(lambda x : x.sort_values(by = '-log(p)', ascending = False).head(top).reset_index(drop = True))

		collections = [ax.collections[i].get_offsets().data for i in range(len(ax.collections)) if ax.collections[i].get_offsets().data.size != 0]

		if n_categories == 1:

			if len(collections) == 3:
				collections = [np.concatenate((collections[1], collections[2]))]
			elif len(collections) == 2:
				collections = [collections[1]]
			else:
				collections = [collections[0]]

		for cat, ind in indices_of_significant_cats:

			# only annotate top variants for each category (to avoid cluttering the plot, use only 2 for the big plots)
			top_xy_coords = sorted(collections[ind], key = lambda x: x[1], reverse = True)

			if top < len(top_xy_coords):
				top_xy_coords = top_xy_coords[:top]

			top_xy_coords = [(x,y) for x,y in top_xy_coords if y >= sig]
			xy[cat].extend(top_xy_coords)
			# print(cat, xy[cat])
	
		texts = []
		for cat, _ in indices_of_significant_cats:
			index = 0
			top = top_N.loc[top_N['Category'] == cat]
			for xy_coords in xy[cat]:
				x,y = xy_coords
				
				# get rsID from the top variants dataframe
				rsID = top.iloc[index]['rsID'].split('_')[0]

				# only annotate first 3 genes, the rest can be found in the results file (for clarity)
				genes = ', '.join(res.loc[res['rsID'] == rsID, 'GENE'].unique().tolist()[:3])
				genes = textwrap.fill(genes, 25, break_long_words=False)

				# # using the beta values, annotate the directionality of betas for each variant [expected, observed]
				# direction = 'Dir=%s,%s' % ('-' if top.iloc[index]['beta'] < 0 else '+', '-' if top.iloc[index]['Original_beta'] < 0 else '+')
				
				# if we have more than 1 category, don't add the description to the plot
				if n_categories > 1:
					annotation = (rsID + ', ' + genes).strip()

					# texts.append(ax.annotate(annotation,
					# 	xy=(x,y), xycoords='data',
					# 	xytext=(x, y+1), textcoords='data',
					# 	horizontalalignment="right",
					# 	arrowprops=dict(arrowstyle="->", connectionstyle="arc", alpha = 0.3),
					# 	fontsize = 3, ha = 'center', va = 'center',))
					texts.append(ax.text(x, y, annotation, ha = 'center', va = 'center', fontsize = 3))

				else:
					desc = top.iloc[index]['Description']
					desc = textwrap.fill(desc, 35, break_long_words=False)
					annotation = (rsID + ',' + genes + '\n' + desc ).strip()
					
					texts.append(ax.text(x, y, annotation, fontsize = 6, ha = 'center', va = 'center'))

				index +=1

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

	if res is not None:
		adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=0.2))

	else:
		output_file = output_file.split('.')[0] + '_no_annotations.png'

	# finally save the figure
	plt.savefig(output_file, bbox_inches ='tight', dpi = 300)
	plt.close()

def significant_bar_plot(regressions, sig,title, output_file, save = True):
	
	print('Plotting bar plot for %s' % output_file)
	regressions['Significant'] = regressions['"-log(p)"'] >= sig
	fig, ax = plt.subplots(figsize = (12, 8))

	total_vars = {}
	significant_vars = {}

	for cat in regressions['Category'].unique():
		total_vars[cat] = regressions[regressions['Category'] == cat].shape[0]
		significant_vars[cat] = regressions.loc[regressions['Category'] == cat, 'Significant'].sum()

	categories = dict(sorted(significant_vars.items(), key=lambda item: item[1]/total_vars[item[0]])).keys()
	not_significant = [(total_vars[cat] - significant_vars[cat]) for cat in categories]
	significant = [100*significant_vars[cat]/total_vars[cat] for cat in categories]

	p1 = ax.bar(categories, significant,  color='tab:blue', label = 'Significant associations')

	plt.xlabel('Categories', labelpad=12)
	ax.set_ylim([None, ax.get_ylim()[1]*1.05])
	plt.xticks(rotation=45, ha='right')
	plt.ylabel('Percentage of significant associations out of total', labelpad=12)

	ax.yaxis.set_major_formatter(mtick.PercentFormatter())

	index = 0
	categories = list(categories)
	for p in p1.patches:
		width, height = p.get_width(), p.get_height()
		x, y = p.get_xy() 
		# ax.annotate(str(round(height*100, 2)) + '%', (p.get_x() + 0.1*width, p.get_y()+height + 0.005), 
		#     fontsize=10, color='black')
		width_addition = 0.2*width

		if total_vars[categories[index]] < 100:
			width_addition = 0.3*width
		elif total_vars[categories[index]] > 1000:
			width_addition = 0.1*width

		ax.annotate(total_vars[categories[index]], (p.get_x() + width_addition, p.get_y()+height + 0.3), 
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


def main():
	parser = argparse.ArgumentParser(description='Plot significant phenotype data on a Manhattan Plot.')
	parser.add_argument('--input', help='Input file', required=True)
	parser.add_argument('--output', help='Output png', required=True)
	parser.add_argument('--aggregate_title', help='Title of the plot of the aggregated data', required=False, default = None, nargs='+')
	parser.add_argument('--mapping', help='Mapping of old categories to new ones', required=False, default = None)
	parser.add_argument('--lower_outlier', help='Lower outlier threshold', required=False, default = 0, type = float)
	parser.add_argument('--upper_outlier', help='Upper outlier threshold', required=False, default = 0.995, type = float)
	parser.add_argument('--seed', help = 'Set seed to make the adjustText library deterministic', required = False, default = 0, type = int)
	parser.add_argument('--downstream', help = 'Number in KB to check for closest genes', required = False, default = 40, type = int)
	parser.add_argument('--upstream', help = 'Number in KB to check for closest genes', required = False, default = 40, type = int)
	parser.add_argument('--gene_file', help = 'File containing a list of genes use for finding closest genes', required = False, default = None)
	parser.add_argument('--rsid_files', help = 'Files containing a list of rsIDs to use for finding closest genes', required = False, action = 'append')
	parser.add_argument('--no_annotations', help = 'Don\'t add annotations', required = False, default = False, action = 'store_true')
	parser.add_argument('--plot_top_categories', help = 'Plot the top categories as well', required = False, default = False, action = 'store_true')
	parser.add_argument('--number_of_top_results', help = 'Number of top results to save', required = False, default = 10, type = int)
	parser.add_argument('--plot_manhattan', help = 'Plot the manhattan plot', required = False, default = False, action = 'store_true')
	parser.add_argument('--plot_bar', help = 'Plot the bar plot', required = False, default = False, action = 'store_true')

	args = parser.parse_args()
	input_file = args.input
	output_file = args.output
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


	np.random.seed(seed)

	# Load data
	regressions = pd.read_csv(input_file, sep='\t')

	if aggregate_title is None:
		pheno_cat = input_file.split('/')[-1].split('_')[-1].split('.')[0]

		if plt_manhattan:
			aggregate_title = 'PheWAS Results for %s (ALL) Variants' % (pheno_cat.capitalize())
			output_file = output_file.split('.')[0] + '_manhattan.png'
		else:
			aggregate_title = '%s (All Predictors) Variant Enrichment per Category' % (pheno_cat.capitalize())
			output_file = output_file.split('.')[0] + '_bar.png'
	else:
		aggregate_title = ' '.join(aggregate_title)

	phenos = regressions['Predictor'].unique()

	pheno_map = {}

	output_ext = os.path.basename(output_file).split('.')[1]
	output_file_dir = os.path.dirname(output_file)
	
	if output_file_dir == '':
		output_file_dir = '.'

	for p in phenos:
		title_and_output = ''
		if plt_manhattan:
			title_and_output = ('PheWAS Results for %s Variants' % (p.capitalize()), output_file_dir + '/'  + p.lower() + '_manhattan.' + output_ext)
		else:
			title_and_output = ('%s Variant Enrichment per Category' % (p.capitalize()), output_file_dir + '/' + p.lower() + '_bar.' + output_ext)

		pheno_map[p] = title_and_output


	if mapping is not None:
		# load mapping file - two column file, first with old categories, second with new categories
		mapping = dict(pd.read_csv(mapping, sep='\t').values)

		# map old categories to new ones
		regressions['Category'] = regressions['Category'].map(mapping)

	gene_df = None
	rsid_df = None
	res = None
	
	if gene_file is not None:
		gene_df = pd.read_csv(gene_file, sep = '\t')
		
		gene_df['#chrom'] = gene_df['#chrom'].apply(lambda x: x.replace('chr', ''))
		chroms = [str(x) for x in list(range(23)) + ['X', 'Y', 'XY']]
		gene_df = gene_df[gene_df['#chrom'].isin(chroms)]

		rsid_df = pd.DataFrame()

		if rsid_files == []:
			print('Need to provide a list of rsID files if you want to use a gene file')
			exit()
		else:
			for f in rsid_files:
				# read f and then concat to rsid_df
				rsid_df = pd.concat([rsid_df, pd.read_csv(f, sep = '\t')])
				rsid_df = rsid_df.drop_duplicates()
			rsid_df['CHR'] = rsid_df['CHR'].astype(str)
		# get closest genes
		res = get_closest_genes(rsid_df, gene_df, downstream, upstream).sort_values(by = 'rsID')


		gene_map = defaultdict(list)

		for _, row in res.iterrows():
			gene_map[row['rsID']].append(row['GENE'])

		regressions['rsID'] = regressions['rsID'].apply(lambda x: x.split('_')[0])
		
		# case where we have a variant like 3:100928901_CTT_C
		rsid_df['rsID'] = rsid_df['rsID'].apply(lambda x: x.split('_')[0])

		regressions['Gene'] = regressions.apply(lambda x: ', '.join(gene_map[x['rsID']]), axis = 1)

		regressions['Original_beta'] = regressions.apply(lambda x: rsid_df.loc[rsid_df['rsID'] == x['rsID']]['BETA'].unique()[0], axis = 1)
		regressions['Original_pval'] = regressions.apply(lambda x: rsid_df.loc[rsid_df['rsID'] == x['rsID']]['P'].unique()[0], axis = 1)
		regressions['Direction'] = np.where(regressions['beta'] * regressions['Original_beta'] > 0, 'S', 'D')
		
	if no_annotations:
		res = None	

	# plot the aggregate data
	alpha = 0.05
	# Plot data
	sig = -np.log10(alpha/regressions.shape[0])

	if plt_manhattan:
		manhattan_plot(regressions, sig, lower_outlier, upper_outlier, aggregate_title, output_file, '', res, plt_top_cat, N)

	if plt_bar:
		# remove extreme outliers to be coherent with manhattan plot
		# outlier_removed = regressions[regressions['"-log(p)"'].between(regressions['"-log(p)"'].quantile(lower_outlier), regressions['"-log(p)"'].quantile(upper_outlier))]
		outlier_removed = regressions.copy(deep = True)
		significant_bar_plot(outlier_removed, sig, aggregate_title, output_file, True)

	for pheno, (pheno_title, pheno_output_file) in pheno_map.items():
		# plot the data for each phenotype
		predictor_specific = regressions[regressions['Predictor'] == pheno].copy(deep = True)

		if plt_manhattan:
			manhattan_plot(predictor_specific, sig,lower_outlier, upper_outlier, pheno_title, pheno_output_file, pheno, res, plt_top_cat, N)
		
		if plt_bar:

			# remove extreme outliers to be coherent with manhattan plot
			# predictor_specific = predictor_specific[predictor_specific['"-log(p)"'].between(predictor_specific['"-log(p)"'].quantile(lower_outlier), predictor_specific['"-log(p)"'].quantile(upper_outlier))]

			significant_bar_plot(predictor_specific, sig, pheno_title, pheno_output_file, True)


if __name__ == '__main__':
	main()