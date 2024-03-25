import os
import argparse
import pandas as pd

def main():
	parser = argparse.ArgumentParser(description='Annotate results file with custom category labels')
	parser.add_argument('-i', '--input_file', help='Aggregate pheWAS results file (output of run_pheWAS for a specific category of phenotypes).', required=True, action = 'append')
	parser.add_argument('-m', '--mapping_file', help='File containing mapping of phenotypes to categories. This must be a tsv with phenotypes in the first column and categories in the second.', required=True)
	parser.add_argument('-o', '--output_file', help='Output aggregate file containg pheWAS results from every specified input_file, annotated with category information.', required=True)

	args = parser.parse_args()
	input_files = args.input_file
	mapping_file = args.mapping_file
	aggregate_file = args.output_file

	mapping_df = pd.read_csv(mapping_file, sep='\t', index_col = None)

	if mapping_df.shape[1] != 2:
		print('Mapping file must have two columns. The first column should contain the phenotypes and the second column should contain the categories.')
		exit()

	mapping_df.columns = ['Phenotype', 'Category']

	aggregate_results = pd.DataFrame()
	temp = []

	for input_file in input_files:

		output_file = input_file.rsplit('.', 1)[0] + '_with_cat.tsv'

		# read the input file
		df = pd.read_csv(input_file, sep='\t', index_col = None)

		# for each phenotype in the mapping file, get the category and add it to the df
		merged_df = pd.merge(df, mapping_df, how = 'left', left_on = 'Data_Field', right_on = 'Phenotype').drop('Phenotype', axis = 1)

		# write the output file
		merged_df.to_csv(output_file, sep = '\t', index = False)

		print('Wrote output file: ' + output_file + '\n')

		temp.append(merged_df)
	
	aggregate_results = pd.concat(temp)

	# check to see if the directory for the aggregate file exists
	if not os.path.exists(os.path.dirname(aggregate_file)):

		print('Creating directory for aggregate results file since it does not exist')

		os.makedirs(os.path.dirname(aggregate_file))

		print('Finished creating directory for aggregate results file\n')
	
	aggregate_results.to_csv(aggregate_file, sep='\t', index=False)

	print('Finished writing aggregate results file')
	print('Wrote output file: ' + aggregate_file + '\n')
	

if __name__ == '__main__':
	main()