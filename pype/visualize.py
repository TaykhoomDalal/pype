import os
import shutil
import argparse
import textwrap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from adjustText import adjust_text
from collections import defaultdict
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from utility_funcs import multiple_testing_correction, annotate_genes

def plot_significant_categories(idx_sig_cats, out_file, title, phewas_results, sig_thresh, annotate, N, cmp_orig_betas, transparency, color_map):
	# for each category in the significant categories
	for cat, _ in idx_sig_cats:
		
		output_file_cat = out_file.rsplit('.', 1)[0] + '_' + cat + '.png'

		title_cat = title + ' - ' + cat
		data_cat = phewas_results.loc[phewas_results['Category'] == cat]

		manhattan_plot(data_cat, sig_thresh, 0, 1, title_cat, output_file_cat, cat, annotate, plt_top_cat = False, N = N, compare_orig_betas = cmp_orig_betas, transparency = transparency, pheno_color = color_map[cat])

def manhattan_plot(phewas_results, sig, low, high, title, output_file, pheno, annotate, plt_top_cat, N, compare_orig_betas, transparency, pheno_color):
	"""
	Function to plot manhattan plots from PheWAS results
	"""

	# save the significant results to a file (before filtering out the super low p-values)
	all_top = phewas_results.loc[phewas_results["-log(p)"] >= sig].groupby('Category').apply(lambda x : x.sort_values(by = "-log(p)", ascending = False).reset_index(drop = True))
	all_top.to_csv(output_file.rsplit('.', 1)[0] + '_significant_results.tsv', sep = '\t', index = False)

	# remove extreme outliers in order to be able to visualize the data more coherently
	outlier_removed_phewas_results = phewas_results[phewas_results["-log(p)"].between(phewas_results["-log(p)"].quantile(low), phewas_results["-log(p)"].quantile(high))]
	
	# remove all results with p-value above the significance level
	# groupby category and find median of all these significant results
	# reset the index and sort by the log pvalues so largest will be at the right end of the plot
	medians = outlier_removed_phewas_results.loc[outlier_removed_phewas_results["-log(p)"] >= sig].groupby('Category')["-log(p)"].mean().reset_index().sort_values(by="-log(p)").reset_index(drop=True)


	sig_cats = medians['Category'].values

	# don't forget to re-add the categories that have no significant results
	# calculate their median pvalues and add them to the dataframe
	# sorting by the log pvalues so largest will still be at the right end of the plot
	for cat in outlier_removed_phewas_results['Category'].unique():
		if cat not in medians['Category'].values:
			cat_median = outlier_removed_phewas_results.groupby('Category')["-log(p)"].mean()[cat]
			medians = pd.concat([medians, pd.DataFrame({'Category': cat, "-log(p)": cat_median}, index = [0])])
	
	medians = medians.sort_values(by="-log(p)").reset_index(drop=True)

	# based on the category ordering, sort the dataframe rows accordingly
	sorter = lambda x: x.map({cat:order for order, cat in enumerate(medians['Category'])})
	
	print(pheno.capitalize())
	print(medians)
	print('\n')
	
	# get the indices of the significant categories
	indices_of_significant_cats = [(cat, medians['Category'].tolist().index(cat)) for cat in sig_cats]

	# if there are no significant categories, then don't plot anything
	if len(indices_of_significant_cats) == 0:
		print('No significant categories found')
		return

	# sort the categories by the median list so we can plot in the right order
	outlier_removed_phewas_results = outlier_removed_phewas_results.sort_values(by='Category', key = sorter)

	# get the list of categories, to map the colors to the categories
	categories = outlier_removed_phewas_results['Category'].unique()
	n_categories = len(categories)

	color_map = {}
	if n_categories > 1:
		# List of RGB triplets
		colors = sns.color_palette(pheno_color, n_categories)

		# Map label to RGB
		color_map = dict(zip(categories, colors))
	else:
		# otherwise there is only one category, so just use the color passed through the pheno_color argument
		color_map[pheno] = pheno_color

	if plt_top_cat:
		# sort the 

		plot_significant_categories(indices_of_significant_cats, output_file, title, phewas_results.sort_values(by='Category', key = sorter), sig, annotate, N, compare_orig_betas, transparency, color_map)
	
	if outlier_removed_phewas_results.empty:
		print('No significant results found for this phenotype')
		return
	
	# get ready to plot
	plt.figure(figsize=(15, 8), dpi=500)

	# data for plotting
	significant_results = outlier_removed_phewas_results.loc[outlier_removed_phewas_results["-log(p)"] >= sig]
	non_significant_results = outlier_removed_phewas_results.loc[outlier_removed_phewas_results["-log(p)"] < sig]

	# get the significant results which have beta values that are positive or negative (to plot them differently)
	significant_phewas_pos_beta = significant_results.loc[significant_results['beta'] >= 0].copy()
	significant_phewas_neg_beta = significant_results.loc[significant_results['beta'] < 0].copy()

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
	if n_categories == 1:
		size = 7
		color = pheno_color
	# if we plotting more than 1 category, save space with points and the colors of the points will be different
	elif n_categories > 1: 
		size = 5
		color = 'black' # the color of the points will be determined by the category, but the legend will be black (can't have multiple colors in the legend for the same entry)
	else:
		raise ValueError('Number of categories must be greater than 0')
	if not significant_phewas_pos_beta.empty:
		if compare_orig_betas and n_categories == 1:
			pos_dir_palette ={"+": "red", "-": "orange"}
			significant_phewas_pos_beta.loc[:, 'original_beta_direction'] = significant_phewas_pos_beta.apply(lambda x: '+' if x['Original_beta'] >= 0 else '-', axis = 1)
			ax = sns.stripplot(x = 'Category', y = "-log(p)", data = significant_phewas_pos_beta, hue="original_beta_direction", palette = pos_dir_palette, jitter=0.45, size = size, order = categories, linewidth=0.3, **{'marker': '^', 'alpha': transparency})
			handles.extend([Line2D([0], [0], color = 'red', marker = '^', linestyle='None',label = 'Pos/Pos', alpha = transparency),
							Line2D([0], [0], color = 'orange', marker = '^', linestyle='None',label = 'Pos/Neg', alpha = transparency)])
		else:
			# Plot the -log(p) values against the category values, with the colors mapped to the categories, with up arrow == the direction of the beta value is positive
			ax = sns.stripplot(x = 'Category', y = "-log(p)", data = significant_phewas_pos_beta, palette = color_map, jitter=0.45, size = size, order = categories, linewidth=0.2, **{'marker': '^', 'alpha': transparency})
			handles.append(Line2D([0], [0], color = color, marker = '^', linestyle='None',label = 'Positive beta', alpha = transparency))
	
	if not significant_phewas_neg_beta.empty:
		if compare_orig_betas and n_categories == 1:
			neg_dir_palette ={"-": "blue", "+": "purple"}
			significant_phewas_neg_beta.loc[:, 'original_beta_direction'] = significant_phewas_neg_beta.apply(lambda x: '+' if x['Original_beta'] >= 0 else '-', axis = 1)
			ax = sns.stripplot(x = 'Category', y = "-log(p)", data = significant_phewas_neg_beta, hue="original_beta_direction", palette = neg_dir_palette, jitter=0.45, size = size, order = categories, linewidth=0.3, **{'marker': 'v', 'alpha': transparency})
			handles.extend([Line2D([0], [0], color = 'blue', marker = 'v', linestyle='None',label = 'Neg/Neg', alpha = transparency),  
							Line2D([0], [0], color = 'purple', marker = 'v', linestyle='None', label = 'Neg/Pos', alpha = transparency)])
		else:
			# Plot the -log(p) values against the category values, with the colors mapped to the categories, with down arrow == the direction of the beta value is negative
			ax = sns.stripplot(x = 'Category', y = "-log(p)", data = significant_phewas_neg_beta, palette = color_map, jitter=0.45, size = size, order = categories, linewidth=0.2, **{'marker': 'v', 'alpha': transparency})
			handles.append(Line2D([0], [0], color = color, marker = 'v', linestyle='None', label = 'Negative beta', alpha = transparency))

	# if we need to add annotations to the plot
	if annotate is not None:
		
		
		if n_categories > 1: # if we are plotting more than 1 category, only annotate top 2 points
			annotate_top = 2
		# else, annotate the top N points
		elif n_categories == 1:
			annotate_top = N
		else:
			raise ValueError("n_categories must be greater than 0")

		# get top N variants for each category
		top_N = significant_results.groupby('Category').apply(lambda x : x.sort_values(by = "-log(p)", ascending = False).head(annotate_top).reset_index(drop = True))
		
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

			# only annotate top results for each category (to avoid cluttering the plot, use only 2 for the big plots)
			top_xy_coords = sorted(collections_for_cat, key = lambda x: x[1], reverse = True)

			# if we want to annotate less points than those that exist, only choose the top X points
			if annotate_top < len(top_xy_coords):
				top_xy_coords = top_xy_coords[:annotate_top]

			index = 0

			# for each of the top points for this category, annotate the plot
			top = top_N.loc[top_N['Category'] == cat]
			for xy_coords in top_xy_coords:
				x,y = xy_coords
				
				# get the independent variable from the top dataframe (either variants or phenotypes as independent variables)
				independent_variable_name = top.iloc[index]['Independent_Var'].split('_')[0]

				annotation = ''

				# if we are annotating phenotypes or if we don't have SNP-Gene mapping
				if annotate == 'phenotype' or (annotate == 'genotype' and 'Gene' not in top.columns):
					# we simply annotate the independent variable name
					annotation = independent_variable_name.strip()

				elif annotate == 'genotype' and 'Gene' in top.columns:
					# only annotate first 3 genes, the rest can be found in the results file (for clarity)
					genes = ', '.join(top.loc[top['Independent_Var'] == independent_variable_name, 'Gene'].unique().tolist()[:3])
					genes = textwrap.fill(genes, 25, break_long_words=False)

					annotation = (independent_variable_name + ', ' + genes).strip()

				# if we have more than 1 category, don't add the description to the plot and make the font smaller
				if n_categories > 1:
					texts.append(ax.text(x, y, annotation, ha = 'center', va = 'center', fontsize = 3))
				else:
					desc = top.iloc[index]['Description']
					desc = textwrap.fill(desc, 35, break_long_words=False)
					annotation = (annotation + '\n' + desc ).strip()
					
					texts.append(ax.text(x, y, annotation, fontsize = 6, ha = 'center', va = 'center'))

				index +=1
	
	if not non_significant_results.empty:
		# plot for the values below the significance level (should be circles)
		ax = sns.stripplot(x = 'Category', y = "-log(p)", data = non_significant_results, palette = color_map, jitter=0.45, size = size, order = categories, linewidth=0.2, **{'alpha': transparency})

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

	if annotate is not None:
		adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=0.2))

	# finally save the figure
	plt.savefig(output_file, bbox_inches ='tight', dpi = 300)
	plt.close()

def category_enrichment_plot(phewas_results, sig,title, output_file):
	
	print('Plotting bar plot for %s' % output_file)

	# annotate each row with whether it is significant or not
	phewas_results['Significant'] = np.where(phewas_results["-log(p)"] >= sig, True, False)
	_, ax = plt.subplots(figsize = (12, 8))

	total_vars = {}
	significant_results = {}

	# for each category get the total number of variants and the number of significant variants
	for cat in phewas_results['Category'].unique():
		total_vars[cat] = phewas_results[phewas_results['Category'] == cat].shape[0]
		significant_results[cat] = phewas_results.loc[phewas_results['Category'] == cat, 'Significant'].sum()

	# sort the categories by the ratio of significant variants to total variants
	categories = dict(sorted(significant_results.items(), key=lambda item: item[1]/total_vars[item[0]])).keys()
	
	# get the percentages for each category (significant variants / total variants)
	significant = [100*significant_results[cat]/total_vars[cat] for cat in categories]

	# plot the results as a bar plot
	p1 = ax.bar(categories, significant,  color='tab:blue', label = 'Significant associations')

	plt.xlabel('Categories', labelpad=12)
	ax.set_ylim([None, ax.get_ylim()[1]*1.05])
	plt.xticks(rotation=45, ha='right')
	plt.ylabel('Percentage of significant associations out of total', labelpad=12)

	ax.yaxis.set_major_formatter(mtick.PercentFormatter())

	# annotate each bar with the total number of variants
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
	
	# create a manual legend for these annotations
	legend_elements = [Patch(facecolor='none', edgecolor='red',label='Total associations')]
	
	ax.legend(handles=legend_elements, loc='upper left')
	plt.tight_layout()
	plt.title(title)
	
	plt.savefig(output_file, bbox_inches ='tight', dpi = 300)
	plt.close()

def volcano_plot(phewas_results, sig, title, output_file, phewas_type, N, max_genes, annotate, transparency, compare_orig_betas, logp_low, logp_high, beta_low, beta_high):

	# make a volcano plot for each independent variable that was tested
	ind_vars = phewas_results['Independent_Var'].unique()

	output_file_split = output_file.rsplit('.', 1)
	ext = '.' + output_file_split[-1]
	output_base_name = output_file_split[0]

	for ind_var in ind_vars:
	
		ind_var_i = phewas_results.loc[phewas_results['Independent_Var'] == ind_var].copy(deep = True)

		# save the significant results to a file (before filtering out the super low p-values)
		all_top = ind_var_i.loc[ind_var_i["-log(p)"] >= sig].sort_values(by = "-log(p)", ascending = False).reset_index(drop = True)
		all_top.to_csv(output_base_name + '_' + ind_var +'.tsv', sep = '\t', index = False)

		# remove extreme outliers in order to be able to visualize the data more coherently
		outlier_removed_ind_var_i = ind_var_i[ind_var_i["-log(p)"].between(ind_var_i["-log(p)"].quantile(logp_low), ind_var_i["-log(p)"].quantile(logp_high))]
	
		# for visualization purposes, in my examples, found to be the best, but can be changed based on the data
		outlier_removed_ind_var_i = outlier_removed_ind_var_i.loc[(outlier_removed_ind_var_i['beta'] < beta_high) & (outlier_removed_ind_var_i['beta'] > beta_low)].reset_index()

		# data for plotting
		significant_results = outlier_removed_ind_var_i.loc[outlier_removed_ind_var_i["-log(p)"] >= sig]
		non_significant_results = outlier_removed_ind_var_i.loc[outlier_removed_ind_var_i["-log(p)"] < sig]

		# get the significant results which have beta values that are positive or negative (to plot them differently)
		significant_phewas_pos_beta = significant_results.loc[significant_results['beta'] >= 0].copy()
		significant_phewas_neg_beta = significant_results.loc[significant_results['beta'] < 0].copy()
		
		# custom legend entries
		handles = []

		# get ready to plot
		plt.figure(figsize=(15, 8), dpi=500)

		# plot the results
		if not significant_phewas_pos_beta.empty:
			if compare_orig_betas: # if the arrow is pointing up, the PheWAS beta is positive, and if the color is red, the original beta was positive, else the original beta was negative
				pos_dir_palette ={"+": "red", "-": "orange"}
				significant_phewas_pos_beta.loc[:, 'original_beta_direction'] = significant_phewas_pos_beta.apply(lambda x: '+' if x['Original_beta'] >= 0 else '-', axis = 1)
				ax = sns.scatterplot(x = 'beta', y = "-log(p)", data = significant_phewas_pos_beta, hue="original_beta_direction", palette = pos_dir_palette, linewidth=0.3, **{'marker': '^', 'alpha': transparency})
				handles.extend([Line2D([0], [0], color = 'red', marker = '^', linestyle='None',label = 'Pos/Pos', alpha = transparency),
								Line2D([0], [0], color = 'orange', marker = '^', linestyle='None',label = 'Pos/Neg', alpha = transparency)])
			else:
				# up arrow == the direction of the beta value is positive
				ax = sns.scatterplot(x = 'beta', y = "-log(p)", data = significant_phewas_pos_beta, color = 'red', linewidth=0.2, **{'marker': '^', 'alpha': transparency})
				handles.append(Line2D([0], [0], color = 'red', marker = '^', linestyle='None', label = 'Positive PheWAS beta', alpha = transparency))
		
		if not significant_phewas_neg_beta.empty:
			if compare_orig_betas: # if the arrow is pointing down, the PheWAS beta is negative, and if the color is blue, the original beta was negative, else the original beta was positive
				neg_dir_palette ={"-": "blue", "+": "purple"}
				significant_phewas_neg_beta.loc[:, 'original_beta_direction'] = significant_phewas_neg_beta.apply(lambda x: '+' if x['Original_beta'] >= 0 else '-', axis = 1)
				ax = sns.scatterplot(x = 'beta', y = "-log(p)", data = significant_phewas_neg_beta, hue="original_beta_direction", palette = neg_dir_palette, linewidth=0.3, **{'marker': 'v', 'alpha': transparency})
				handles.extend([Line2D([0], [0], color = 'blue', marker = 'v', linestyle='None',label = 'Neg/Neg', alpha = transparency),  
								Line2D([0], [0], color = 'purple', marker = 'v', linestyle='None', label = 'Neg/Pos', alpha = transparency)])
			else:
				# down arrow == the direction of the beta value is negative
				ax = sns.scatterplot(x = 'beta', y = "-log(p)", data = significant_phewas_neg_beta, color = 'blue', linewidth=0.2, **{'marker': 'v', 'alpha': transparency})
				handles.append(Line2D([0], [0], color = 'blue', marker = 'v', linestyle='None', label = 'Negative PheWAS beta', alpha = transparency))
		if not non_significant_results.empty:
			# plot for the values below the significance level (should be circles)
			ax = sns.scatterplot(x = 'beta', y = "-log(p)", data = non_significant_results, color = 'black', linewidth=0.2, **{'alpha': transparency})
		
		# add the legends
		if len(handles) > 0:
			ax.legend(handles = handles, loc = 'best')

		# draw a line at the threshold
		plt.axhline(sig, color = 'red')
		
		ind_title = title + ' : ' + ind_var

		if annotate:
			texts = []
			
			# get top N most significant results
			outlier_removed_ind_var_i = outlier_removed_ind_var_i.sort_values(by = '-log(p)', ascending = False)

			descriptions = outlier_removed_ind_var_i['Description']

			results_annotated = N
			# annotate the top N most significant results with their description
			for i, desc in enumerate(descriptions):
				if outlier_removed_ind_var_i.iloc[i]["-log(p)"] > sig:
					
					# check if either the x or y value is not a finite number
					if outlier_removed_ind_var_i.iloc[i]["-log(p)"] == np.inf or outlier_removed_ind_var_i.iloc[i]["beta"] == np.inf:
						continue
					else:
						texts.append(plt.text(outlier_removed_ind_var_i.iloc[i]['beta'], outlier_removed_ind_var_i.iloc[i]["-log(p)"], desc, fontsize = 9, ha = 'center', va = 'center'))
						results_annotated -= 1

						if results_annotated == 0:
							break
			
			if phewas_type == 'genotype' and 'Gene' in outlier_removed_ind_var_i.columns:
				genes = outlier_removed_ind_var_i['Gene'].unique()
			
				if genes.any():
					genes = ''
				else:
					if len(genes) > max_genes:
						genes = ", ".join(genes[:max_genes])
					genes = '(' + genes + ')'

					ind_title = ind_title + ' ' + ", ".join(genes)
			
			adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=0.2))

		plt.title(ind_title)
		plt.savefig(output_base_name + '_' + ind_var + ext, bbox_inches ='tight', dpi = 300)
		plt.close()

def createDirectory(directory, clear):
	# if clear_old_files selected, delete the directory
	if clear:
		print("Clearing old files from the directory located at: " + directory)
		shutil.rmtree(directory, ignore_errors=True)

	# check to see if directory exists, creating it if not
	if not os.path.exists(directory):
		print("Creating directory at: " + directory)
		os.mkdir(directory)

def createPlotSpecificDirectories(main_directory, output_prefix, clear, plot_type, separate_plots):
	
	# remove backslash from end of directory if present
	if main_directory[-1] == '/':
		main_directory = main_directory[:-1]

	# overall directory
	plot_type_directory = main_directory + '/' + output_prefix + '_' + plot_type
	createDirectory(plot_type_directory, clear)

	# directory for each plot at the top level
	plot_type_agg_directory = plot_type_directory + '/' + output_prefix + '_aggregated'
	createDirectory(plot_type_agg_directory, clear)

	# directory for each plot by category
	if separate_plots:
		plot_type_separate_directory = plot_type_directory + '/' + output_prefix + '_separate'
		createDirectory(plot_type_separate_directory, clear)

		return plot_type_agg_directory, plot_type_separate_directory
	
	return plot_type_agg_directory, None

def main():
	parser = argparse.ArgumentParser(description='Visualize PheWAS results on a variety of plot types.')

	# Required Arguments
	parser.add_argument('--phewas_results', help = "Input PheWAS results file", required=True, type = str)
	parser.add_argument('--directory_name', help = "Name of the directory to save the plots in", required=True, type = str)
	parser.add_argument('--output_prefix', help = "Prefix for output files", required = True, type = str)
	parser.add_argument('--output_extension', help = "Extension for output files", required = True, choices = [".png", ".jpg", ".jpeg"], type = str)
	parser.add_argument('--phewas_type', help = "Whether the PheWAS was run with genotypes or phenotypes as the independent variables", required = True, choices = ["genotype", "phenotype"], type = str)

	# Optional Miscellaneous Arguments
	parser.add_argument('--mapping', help = "File containing mapping of old categories to new ones (should be pandas DataFrame)", required=False, default = None, type = str)
	parser.add_argument('--clear_old_files', help = "Clear old files", required = False, default = False, action = 'store_true')
	parser.add_argument('--plot_phewas_categories_separately', help = "Plot PheWAS results for each category separately as well as in aggregate", required = False, default = False, action = 'store_true')
	parser.add_argument('--compare_original_betas', help = "Whether to annotate significant results with their original beta values and PheWAS beta values", required = False, default = False, action = 'store_true')
	parser.add_argument('--single_category', help = "Whether the PheWAS results file contains results across just one independent variable category", required = False, default = False, action = 'store_true')

	# Gene Annotation
	parser.add_argument('--variant_files', help = "Files containing the list of variants (from the PheWAS) to use for SNP-to-Gene mapping", required = False, action = 'append')
	parser.add_argument('--gene_file', help = "File containing the list of genes to use for SNP-to-Gene mapping", required = False, default = None)
	parser.add_argument('--downstream', help = "Number in KB to check for closest genes", required = False, default = 10, type = int)
	parser.add_argument('--upstream', help = "Number in KB to check for closest genes", required = False, default = 10, type = int)
	parser.add_argument('--save_gene_annotation', help = "Save the gene annotation file", required = False, default = False, action = 'store_true')

	# Controlling Plot Aesthetics
	parser.add_argument('--seed', help = "Seed value to make the adjustText library (annotations) deterministic", required = False, default = None, type = int)
	parser.add_argument('--annotate', help = "Add annotations to the manhattan or volcano plots", required = False, default = False, action = 'store_true')
	parser.add_argument('--transparency', help = "Transparency of points", required = False, default = 0.75, type = float)
	parser.add_argument('--color_map', help = "Seaborn color map to use", required = False, default = 'rainbow', type = str)

	# Manhattan Plot Specific
	parser.add_argument('--plot_manhattan', help = "Plot the manhattan plot", required = False, default = False, action = 'store_true')
	parser.add_argument('--plot_top_data_field_categories', help = "Plot each data field category in the manhattan plot with significant results", required = False, default = False, action = 'store_true')
	parser.add_argument('--annotate_top_N_manhattan', help = "Number of how many (max) significant results to annotate in each manhattan plot", required = False, default = 10, type = int)
	parser.add_argument('--manhattan_title', help = "Title of the top level manhattan plot of the data", required=False, default = None, nargs='+')
	parser.add_argument('--logp_lower_outlier_thresh_manhattan', help = "Lower outlier threshold (ratio) of points to remove for the manhattan plot (for improving how the plot looks)", required=False, default = 0, type = float)
	parser.add_argument('--logp_upper_outlier_thresh_manhattan', help = "Upper outlier threshold (ratio) of points to remove for the manhattan plot (for improving how the plot looks)", required=False, default = 0.995, type = float)

	# Category Enrichment Plot Specific
	parser.add_argument('--plot_category_enrichment', help = "Plot category enrichment bar plot", required = False, default = False, action = 'store_true')
	parser.add_argument('--category_enrichment_title', help = "Title of the category enrichment plot", required=False, default = None, nargs='+')

	# Volcano Plot Specific
	parser.add_argument('--plot_volcano', help = "Plot the volcano plot", required = False, default = False, action = 'store_true')
	parser.add_argument('--volcano_title', help = "Title of the volcano plot (will have the SNP-Gene mapping appended if annotations and gene file specified)", required=False, default = None, nargs='+')
	parser.add_argument('--max_genes', help = "Maximum number of genes to include in the title for the volcano plot", required=False, default = 3, type = int)
	parser.add_argument('--annotate_top_N_volcano', help = "Number of how many (max) significant results to annotate in each volcano plot", required = False, default = 5, type = int)
	parser.add_argument('--beta_lower_outlier_val', help = "Lower outlier values of betas to remove for the volcano plot (for improving how the plot looks)", required=False, default = -10, type = float)
	parser.add_argument('--beta_upper_outlier_val', help = "Upper outlier values of betas to remove for the volcano plot (for improving how the plot looks)", required=False, default = 10, type = float)
	parser.add_argument('--logp_lower_outlier_thresh_volcano', help = "Lower outlier threshold (ratio) of points to remove for the volcano plot (for improving how the plot looks)", required=False, default = 0, type = float)
	parser.add_argument('--logp_upper_outlier_thresh_volcano', help = "Upper outlier threshold (ratio) of points to remove for the volcano plot (for improving how the plot looks)", required=False, default = 1, type = float)


	# Significance Level Information
	parser.add_argument('--alpha', help = "Significance threshold", required = False, default = 0.05, type = float)
	parser.add_argument('--correction', help = "Correction method", required = False, default = 'bonferroni', choices = ['bonferroni', 'sidak','fdr_bh', 'no_correction'], type = str)
	
	args = parser.parse_args()

	# Required Arguments
	phewas_results = args.phewas_results
	directory_name = args.directory_name
	output_prefix = args.output_prefix
	output_extension = args.output_extension
	phewas_type = args.phewas_type

	# Optional Miscellaneous Arguments
	mapping = args.mapping
	clear_old_files = args.clear_old_files
	plot_phewas_categories_separately = args.plot_phewas_categories_separately
	compare_original_betas = args.compare_original_betas
	single_category = args.single_category

	# Gene Annotation
	variant_files = args.variant_files
	gene_file = args.gene_file
	downstream = args.downstream
	upstream = args.upstream
	save_gene_annotation = args.save_gene_annotation

	# Controlling Plot Aesthetics
	logp_lower_outlier_thresh_manhattan = args.logp_lower_outlier_thresh_manhattan
	logp_upper_outlier_thresh_manhattan = args.logp_upper_outlier_thresh_manhattan
	seed = args.seed
	annotate = args.annotate
	transparency = args.transparency
	color_map = args.color_map

	# Manhattan Plot Specific
	plot_manhattan = args.plot_manhattan
	plot_top_data_field_categories = args.plot_top_data_field_categories
	annotate_top_N_manhattan = args.annotate_top_N_manhattan
	manhattan_title = args.manhattan_title

	# Category Enrichment Plot Specific
	plot_category_enrichment = args.plot_category_enrichment
	category_enrichment_title = args.category_enrichment_title

	# Volcano Plot Specific
	plot_volcano = args.plot_volcano
	volcano_title = args.volcano_title
	max_genes = args.max_genes
	annotate_top_N_volcano = args.annotate_top_N_volcano
	beta_lower_outlier_val = args.beta_lower_outlier_val
	beta_upper_outlier_val = args.beta_upper_outlier_val
	logp_lower_outlier_thresh_volcano = args.logp_lower_outlier_thresh_volcano
	logp_upper_outlier_thresh_volcano = args.logp_upper_outlier_thresh_volcano
	
	# Significance Level Information
	alpha = args.alpha
	correction = args.correction

	# ---------------------------------------VERIFY ARGS--------------------------------------- #
	
	annotation_type = None
	
	if compare_original_betas and len(variant_files) == 0:
		print("Error: You must specify the variant files to plot the original betas versus the PheWAS betas")
		exit()
	
	if phewas_type == 'phenotype':
		if len(variant_files) != 0:
			print("Warning: You have specified variant files for annotation but have indicated that this is a Phenotype PheWAS. The variant files will be ignored.")
			variant_files = []

		if gene_file is not None:
			print("Warning: You have specified a gene file for annotation but have indicated that this is a Phenotype PheWAS. The gene file will be ignored.")
			gene_file = None

		if compare_original_betas is not False:
			print("Warning: You have specified to compare the original betas to the PheWAS betas but have indicated that this is a Phenotype PheWAS. No comparison will occur.")
			compare_original_betas = False
		
		if annotate:
			annotation_type = 'phenotype'
	
	elif phewas_type == 'genotype':
		if annotate:
			annotation_type = 'genotype'
	
	if plot_phewas_categories_separately and single_category:
		print("Error: You have specified to plot the PheWAS categories separately but also specified that there is only one category to plot. Please only specify one of these options.")
		exit()
	
	if plot_manhattan or plot_volcano:
		if logp_lower_outlier_thresh_manhattan < 0 or logp_lower_outlier_thresh_manhattan > 1:
			print("Error: The log(p) lower outlier must be between 0 and 1.")
			exit()
		
		if logp_upper_outlier_thresh_manhattan < 0 or logp_upper_outlier_thresh_manhattan > 1:
			print("Error: The log(p) upper outlier must be between 0 and 1.")
			exit()
		
		if logp_lower_outlier_thresh_manhattan >= logp_upper_outlier_thresh_manhattan:
			print("Error: The log(p) lower outlier cannot be greater than or equal to the log(p) upper outlier.")
			exit()
	
	if plot_volcano:
		
		if logp_lower_outlier_thresh_volcano < 0 or logp_lower_outlier_thresh_volcano > 1:
			print("Error: The log(p) lower outlier must be between 0 and 1.")
			exit()
		
		if logp_upper_outlier_thresh_volcano < 0 or logp_upper_outlier_thresh_volcano > 1:
			print("Error: The log(p) upper outlier must be between 0 and 1.")
			exit()
		
		if logp_lower_outlier_thresh_volcano >= logp_upper_outlier_thresh_volcano:
			print("Error: The log(p) lower outlier cannot be greater than or equal to the log(p) upper outlier.")
			exit()

		if beta_lower_outlier_val >= beta_upper_outlier_val:
			print("Error: The beta lower outlier cannot be greater than or equal to the beta upper outlier.")
			exit()
	
	if transparency < 0 or transparency > 1:
		print("Error: The transparency must be between 0 and 1.")
		exit()

	# ----------------------------------------------------------------------------------------- #

	# -------------------------------------SETUP DATA CODE------------------------------------- #
	if seed is not None:
		# set a seed to make results reproducible (for the adjustText library)
		np.random.seed(seed)

	# load the data to be plotted
	phewas_results = pd.read_csv(phewas_results, sep = '\t')

	# any rows with missing values, positive or negative infinity, set them to the max value for that column
	phewas_results = phewas_results.replace([np.inf, -np.inf], np.nan)
	phewas_results = phewas_results.fillna(phewas_results.max(numeric_only=True))
	

	# if there is a user specified mapping, then use it
	if mapping is not None:
		# load mapping file - two column file, first with old categories, second with new categories
		mapping = dict(pd.read_csv(mapping, sep='\t').values)

		# map old categories to new ones
		phewas_results['Category'] = phewas_results['Category'].map(mapping)

	# correct the pvalues using the method of choice
	print('Correcting pvalues using %s correction' % (correction))
	try:
		significance_level = -np.log10(multiple_testing_correction(pvalues = phewas_results['p-val'].dropna(), alpha = alpha, method = correction))
	except:
		print("Error: %s correction failed due to taking the log of a negative number or zero. Try using a different correction method." % (correction))
		exit()

	# ----------------------------------------------------------------------------------------- #

	# ----------------------------------GENOTYPE BASED PHEWAS---------------------------------- #

	# if the user wants to annotate and has specified variant files, then load them
	if len(variant_files) != 0:

		phewas_results = phewas_results.rename(columns = {'Independent_Var': 'rsID'})

		variant_file_list = []

		# load each variant file and append to a list
		for variant_file in variant_files:
			variant_file_list.append(pd.read_csv(variant_file, sep = '\t'))
		
		# concatenate all the variant files into one dataframe, drop duplicates, and reset the index
		variant_data = pd.concat(variant_file_list, axis = 0).drop_duplicates().reset_index(drop = True)

		# before moving on, fix the rsID column by removing the allele information if it exists
		variant_data['rsID'] = variant_data['rsID'].apply(lambda x: x.split('_')[0])
		phewas_results['rsID'] = phewas_results['rsID'].apply(lambda x: x.split('_')[0])

		# for each row in our phewas results, add a column for the variant specified in the row with the value of the original beta & original p-value from the variant files
		phewas_results['Original_beta'] = phewas_results.apply(lambda x: variant_data.loc[variant_data['rsID'] == x['rsID']]['BETA'].unique()[0], axis = 1)
		phewas_results['Original_pval'] = phewas_results.apply(lambda x: variant_data.loc[variant_data['rsID'] == x['rsID']]['P'].unique()[0], axis = 1)

		# for each row in our phewas results, add a column specifying whether the direction of the original beta is the same as the direction of the phewas beta
		phewas_results['Direction'] = np.where(phewas_results['beta'] * phewas_results['Original_beta'] > 0, 'S', 'D')
	
		if gene_file is not None:
			phewas_results, _ = annotate_genes(gene_file = gene_file, 
												rsid_df = variant_data, 
												down = downstream, 
												up = upstream, 
												regressions = phewas_results, 
												output_dir = directory_name, 
												pheno_name = output_prefix,
												save = save_gene_annotation)
		phewas_results = phewas_results.rename(columns = {'rsID': 'Independent_Var'})

	# ----------------------------------------------------------------------------------------- #

	# ---------------------------------------CREATE DIRS--------------------------------------- #

	# create original directory if it doesn't exist
	if not os.path.exists(directory_name):
		print("Creating directory at: " + directory_name)
		os.mkdir(directory_name)

	# for each plot type, create a directory for it, the aggregated plots, and the separate plots if present, deleting each if they already exist
	if plot_manhattan:
		if not annotate:
			name = 'manhattan_not_annotated'
		else:
			name = 'manhattan_annotated'

		manhattan_agg_dir, manhattan_sep_dir = createPlotSpecificDirectories(main_directory = directory_name, 
																			output_prefix = output_prefix,
																			clear = clear_old_files, 
																			plot_type = name, 
																			separate_plots = plot_phewas_categories_separately)
	if plot_category_enrichment:
		category_enrichment_agg_dir, category_enrichment_sep_dir = createPlotSpecificDirectories(main_directory = directory_name,
																								output_prefix = output_prefix,
																								clear = clear_old_files,
																								plot_type = 'category_enrichment',
																								separate_plots = plot_phewas_categories_separately)
	if plot_volcano:
		if not annotate:
			name = 'volcano_not_annotated'
		else:
			name = 'volcano_annotated'

		volcano_agg_dir, volcano_sep_dir = createPlotSpecificDirectories(main_directory = directory_name,
																		output_prefix = output_prefix,
																		clear = clear_old_files,
																		plot_type = name,
																		separate_plots = plot_phewas_categories_separately)
	
	MANHATTAN = 0
	CATEGORY_ENRICHMENT = 1
	VOLCANO = 2

	FILE = 0
	TITLE = 1

	if plot_phewas_categories_separately:
		phewas_categories = phewas_results['PheWAS_Category'].unique()

		phewas_categories_map = {}

		for category in phewas_categories:

			phewas_categories_map[category] = [(None, None), (None, None), (None, None)]

			if plot_manhattan:
				category_directory = manhattan_sep_dir + '/' + category
				createDirectory(category_directory, clear_old_files)
				
				out_file = category_directory + '/' + output_prefix + '_' + category.lower() + output_extension
				out_title = 'PheWAS Separated Manhattan Plot (%s)' % (category.capitalize())
				phewas_categories_map[category][MANHATTAN] = (out_file, out_title)
			
			if plot_category_enrichment:
				category_directory = category_enrichment_sep_dir + '/' + category
				createDirectory(category_directory, clear_old_files)
				
				out_file = category_directory + '/' + output_prefix + '_' + category.lower() + output_extension
				out_title = 'PheWAS Separated Category Enrichment Plot (%s)' % (category.capitalize())
				phewas_categories_map[category][CATEGORY_ENRICHMENT] = (out_file, out_title)
			
			if plot_volcano:
				category_directory = volcano_sep_dir + '/' + category
				createDirectory(category_directory, clear_old_files)

				out_file = category_directory + '/' + output_prefix + '_' + category.lower() + output_extension
				out_title = 'PheWAS Separated Volcano Plot (%s)' % (category.capitalize())
				phewas_categories_map[category][VOLCANO] = (out_file, out_title)


	# ----------------------------------------------------------------------------------------- #
	
	# ---------------------------------------PLOT RESULTS-------------------------------------- #


	if plot_manhattan:

		manhattan_title = ' '.join(manhattan_title) if manhattan_title is not None else 'PheWAS Aggregate Manhattan Plot (Aggregate)'

		manhattan_plot(phewas_results  = phewas_results, 
						sig = significance_level, 
						low = logp_lower_outlier_thresh_manhattan,
						high = logp_upper_outlier_thresh_manhattan,
						title = manhattan_title,
						output_file = manhattan_agg_dir + '/' + output_prefix + output_extension,
						pheno = 'Aggregate',
						annotate = annotation_type,
						plt_top_cat = plot_top_data_field_categories, 
						N = annotate_top_N_manhattan, 
						compare_orig_betas = compare_original_betas, 
						transparency = transparency,
						pheno_color = color_map)

	if plot_category_enrichment:

		category_enrichment_title = ' '.join(category_enrichment_title) if category_enrichment_title is not None else 'PheWAS Aggregate Category Enrichment Plot (Aggregate)'

		outlier_removed = phewas_results.copy(deep = True)
		category_enrichment_plot(phewas_results = outlier_removed, 
								sig = significance_level, 
								title = category_enrichment_title, 
								output_file = category_enrichment_agg_dir + '/' + output_prefix + output_extension)


	if plot_volcano:

		volcano_title = ' '.join(volcano_title) if volcano_title is not None else 'PheWAS Volcano Plot (Aggregate)'

		volcano_plot(phewas_results = phewas_results,
						sig = significance_level,
						title = volcano_title,
						output_file = volcano_agg_dir + '/' + output_prefix + output_extension,
						phewas_type = phewas_type,
						N = annotate_top_N_volcano,
						max_genes = max_genes,
						annotate = annotate,
						transparency = transparency,
						compare_orig_betas = compare_original_betas,
						logp_low = logp_lower_outlier_thresh_volcano,
						logp_high = logp_upper_outlier_thresh_volcano,
						beta_low = beta_lower_outlier_val,
						beta_high = beta_upper_outlier_val)
						
	if plot_phewas_categories_separately:
		for category, plot_info in phewas_categories_map.items():

			# plot the data for each phenotype
			predictor_specific = phewas_results[phewas_results['PheWAS_Category'] == category].copy(deep = True)

			if plot_manhattan:
				print('Plotting %s Manhattan plot' % (category))
				manhattan_plot(phewas_results  = predictor_specific, 
								sig = significance_level, 
								low = logp_lower_outlier_thresh_manhattan,
								high = logp_upper_outlier_thresh_manhattan,
								title = plot_info[MANHATTAN][TITLE],
								output_file = plot_info[MANHATTAN][FILE],
								pheno = category,
								annotate = annotation_type,
								plt_top_cat = plot_top_data_field_categories, 
								N = annotate_top_N_manhattan, 
								compare_orig_betas = compare_original_betas, 
								transparency = transparency,
								pheno_color = color_map)

			if plot_category_enrichment:
				print('Plotting %s bar plot' % (category))
				category_enrichment_plot(phewas_results = predictor_specific, 
								sig = significance_level, 
								title = plot_info[CATEGORY_ENRICHMENT][TITLE], 
								output_file = plot_info[CATEGORY_ENRICHMENT][FILE])

			if plot_volcano:
				print('Plotting %s volcano plot' % (category))
				volcano_plot(phewas_results = predictor_specific, 
							sig = significance_level,
							title = plot_info[VOLCANO][TITLE], 
							output_file = plot_info[VOLCANO][FILE],
							phewas_type = phewas_type,
							N = annotate_top_N_volcano,
							max_genes = max_genes,
							annotate = annotate,
							transparency = transparency,
							compare_orig_betas = compare_original_betas,
							logp_low = logp_lower_outlier_thresh_volcano,
							logp_high = logp_upper_outlier_thresh_volcano,
							beta_low = beta_lower_outlier_val,
							beta_high = beta_upper_outlier_val)


if __name__ == '__main__':
	main()