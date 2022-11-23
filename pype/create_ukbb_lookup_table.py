import argparse
import pandas as pd
import utility_funcs as uf
from bs4 import BeautifulSoup

dict_type_to_pandas_type = {'Integer' : 'Int64', 
							'Sequence' : 'Int64', 
							'Categorical (single)' : 'Int64', 
							'Boolean' : 'boolean',
							'Categorical (multiple)' : 'category',
							'Text' : 'string',
							'Date': 'string',
							'Time': 'string',
							'Curve': 'string',
							'Continuous': 'Float64',
							'Curve' : 'str'}

# https://stackoverflow.com/questions/2935658/beautifulsoup-get-the-contents-of-a-specific-table
def main():
	parser = argparse.ArgumentParser(description='Create lookup table for pheWAS')
	parser.add_argument('--input', help='Input UKBB html file associated with phenotype.tab file', required=True)
	parser.add_argument('--output', help='Output lookup table', required=True)
	parser.add_argument('--pickle', help='Pickle lookup table', required = False, action='store_true')
	parser.add_argument('--pickle_protocol', help='Pickle protocol', required = False, default=-1)
	parser.add_argument('--make_dictionary', help='Make dictionary rather than dataframe', required = False, action='store_true')

	args = parser.parse_args()
	input_file = args.input
	output_file = args.output
	pickle_flag = args.pickle
	pickle_protocol = args.pickle_protocol
	make_dictionary = args.make_dictionary

	# Create beautiful soup object which can be used to parse the html file
	with open(input_file) as f:
		soup = BeautifulSoup(f, 'html.parser')

	# the second table is the one we want
	htmltable = soup.findAll('table')[1]

	rows_in_table = htmltable.find_all('tr')

	# there has to be a header row in the table
	columns = ['UDI', 'Count', 'Type', 'Description']

	# never iteratively append to a pandas dataframe, use 
	data_rows = []
	old_type_and_description = []
	for row in rows_in_table[1:]: # for every table row
		fixed_row = [list(col.stripped_strings)[0] for col in row.find_all('td')]
		
		try:
			
			# assuming that there are 5 columns in the row (Column, UDI (Unique Data Identifier), Count, Type, Description)
			# we don't need the Column column, so drop it

			fixed_row = fixed_row[1:]

			# parse the UDI column --> should be a string, format for standard data fields is field_id-instance_index.array_index
			# change it to how the column is named in the .tab file
			fixed_row[0] = 'f.' + fixed_row[0].replace('-', '.')

			# parse the count column --> should be an integer
			fixed_row[1] = int(fixed_row[1])

			if len(fixed_row) < len(columns):
				fixed_row += old_type_and_description
			else:
				
				# parse the type column --> should be a string, then map to pandas type using dict_type_to_pandas_type
				fixed_row[2] = dict_type_to_pandas_type[fixed_row[2]]

				# Note: the description column --> should be a string (and is already parsed properly)

				old_type_and_description = fixed_row[2:4]
		except:
			print('Error: row has less than 4 columns, and no old type and description')
			exit(1)
		
		data_rows.append(fixed_row)
	

	# create a pandas dataframe
	ukbb_df = pd.DataFrame(data_rows, columns=columns)

	if make_dictionary:
		ukbb_df.set_index('UDI',inplace=True)
		ukbb_df = ukbb_df[['Type', 'Description']].to_dict('index')

		# if this option is selected, then write the dictionary to a pickle file
		uf.pickle_object(ukbb_df, output_file.split('.')[0] + '.pkl', protocol=pickle_protocol)
		
		return

	if pickle_flag:
		# pickle the dataframe --> -1 is the HIGHEST_PROTOCOL --> least backwards compatiblility
		ukbb_df.to_pickle(output_file, protocol = pickle_protocol)
	else:
		# save the dataframe to a tsv file
		ukbb_df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
	main()