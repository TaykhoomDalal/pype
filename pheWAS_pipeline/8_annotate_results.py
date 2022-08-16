import argparse
import pandas as pd
import numpy as np
import os
import requests
from bs4 import BeautifulSoup
import re
from html import unescape

def main():
	parser = argparse.ArgumentParser(description='Annotate results file with category labels')
	parser.add_argument('-i', '--input_files', help='Input results file', required=True, action = 'append')
	parser.add_argument('-c', '--categories', help='Category to parse', required=True, action = 'append')
	parser.add_argument('-s', '--singulars', help = 'Singular category name', required = True, action = 'append')
	parser.add_argument('-o', '--output_file', help='Output file', required=True)

	args = parser.parse_args()
	input_file = args.input_files
	categories_list = args.categories
	singular_list = args.singulars
	aggregate_file = args.output_file

	aggregate_results = pd.DataFrame()

	if len(input_file) != len(categories_list) or len(input_file) != len(singular_list):
		print('Number of input files and categories must match')
		exit()

	for input_file, category, singular in zip(input_file, categories_list, singular_list):

		# make the output file name
		output_file = input_file.rsplit('.')[0] + '_with_cat.tsv'

		print('Parsing url for category ' + category + ' information')

		# read the input file
		df = pd.read_csv(input_file, sep='\t', index_col = None)

		# Get the html file
		html = requests.get('https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=' + str(category)).text

		# Create beautiful soup object which can be used to parse the html file
		soup = BeautifulSoup(html, 'html.parser')
		
		if singular == "True":
			ARROW = unescape("&#9204;")

			# this is bad but it works - basically find the section of the html document where the category exists
			# based on observation, it always comes right after a \n, so split on that, and take the second element,
			# since the first element contains the Catgory X words, and the second element contains the category
			# then split on the arrow (if it exists) and take the first element, removing any leading or trailing whitespace
			category_name = soup.find("div", class_="main").text.split('\n')[1].split(ARROW)[0].strip()

			df['Category'] = category_name
		else:
			# Get all tags which correspond to the data field descriptions
			data_html_tags = soup.find_all("a", class_="subtle", href=re.compile("field.cgi"))

			# parse and grab the actual data-fields descriptions
			data_desc = [tag.string for tag in data_html_tags]

			# Get all tags which correspond to the data field categories
			category_html_tags = soup.find_all("a", class_="subtle", href=re.compile("label.cgi"))

			# parse and grab the actual data-fields categories
			data_cat = [tag.string for tag in category_html_tags]

			# Create a dictionary with the data field descriptions as keys and the data field categories as values
			data_dict = dict(zip(data_desc, data_cat))

			df['Category'] = df['Description'].replace(data_dict)

		df.to_csv(output_file, sep='\t', index=False)

		print('Finished parsing url for ' + category + ' information')
		print('Output file: ' + output_file + '\n')

		aggregate_results = aggregate_results.append(df)

	aggregate_results = aggregate_results.drop(columns = ['Bonferroni_correction'], axis = 1, errors='ignore')
	aggregate_results.to_csv(aggregate_file, sep='\t', index=False)

	print('Finished writing aggregate results file')
	print('Output file: ' + aggregate_file + '\n')

if __name__ == '__main__':
	main()