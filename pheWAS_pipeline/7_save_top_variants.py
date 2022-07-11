import argparse
import os
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup
import re
import requests
from os import listdir
from os.path import isfile, join

def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Pass the input html file')
    parser.add_argument('-u', '--url', help='Url to grab description of fields', required=True)
    parser.add_argument('-o', '--out_dir', help='Directory of the PheWAS Results (and where to store the results)', required=True)

    args = parser.parse_args()
    url = args.url	
    out_dir = args.out_dir

    output = out_dir + '/aggregate_pheWAS_results.tsv'

    # if the output file already exists, delete it
    if os.path.exists(output):
        os.remove(output)

    out_dir = os.path.dirname(output)
    out_file = os.path.basename(output)
    top_SNPS_file = out_dir + '/TOP_SNPS_' + out_file

    # if the top SNPS output file already exists, delete it
    if os.path.exists(top_SNPS_file):
        os.remove(top_SNPS_file)

    # Get the html file
    html = requests.get(url).text

    # Create beautiful soup object which can be used to parse the html file
    soup = BeautifulSoup(html, 'html.parser')
    
    # Get all tags which correspond to the Data_Fields for ukbb
    desc_and_field_tags = soup.find_all("a", class_="subtle", href=re.compile("field.cgi"))
    category_tags = soup.find_all("a", class_="subtle", href=re.compile("label.cgi"))

    # parse and grab the actual data-fields
    field_descriptions = [tag.string for tag in desc_and_field_tags]
    fields = [tag.get('href').split("=")[1] for tag in desc_and_field_tags]
    
    field_dict = {field: description for field, description in zip(fields, field_descriptions)}

     # create a list of the files in the out_dir (these are our pheWAS results)
    results_files = [f for f in listdir(out_dir) if isfile(join(out_dir, f))]

    # Create a dataframe to store the results
    # get the headers from the first file in the list 
    cols = ['Predictor','Ethnicity', 'Description'] + pd.read_csv(out_dir + '/' + results_files[0], sep='\t', nrows = 0).columns.tolist()
    results = pd.DataFrame(columns = cols)
    top_SNPS = pd.DataFrame(columns = cols)

    # Loop through all the results files and add them to the dataframe
    for file in results_files:
        # Read the file
        df = pd.read_csv(out_dir + '/' + file, sep = '\t', header = 0)
        
        data_fields = df['Data_Field'].tolist()

        # get the ethnic group that was included in the covariates 
        # (naming is like: Organ_chromosome_version_bfile_Ethnicity_Group_results.tsv)
        # first split gets us: 'Organ_chromosome_version_bfile' 'Group_results.tsv']
        # second split gets us: 'Group', ''
        group = file.split('_Ethnicity_')[1].split('_results.tsv')[0]

        # get age predictor prefix from file name # I.E. Pancreas_chr10_v3_bfile_Ethnicity_Asian_results.tsv --> Pancreas
        age_predictor_prefix = file.split('_')[0].split('Age')[0]

        for field in data_fields:
            # split the field at the period twice, once to get ['f', 'X.X.X'], and 
            # the second time to get ['f', 'X', 'X.X'], the first X represents the field
            field_without_periods = field.split(".", 2)[1]
            
            # Get the description, category, and category description for the field
            description = field_dict[field_without_periods]

            # Add the description, predictor, ethnicity to the dataframe
            df.loc[df['Data_Field'] == field, 'Description'] = description
            df.loc[df['Data_Field'] == field, 'Predictor'] = age_predictor_prefix
            df.loc[df['Data_Field'] == field, 'Ethnicity'] = group

        results = results.append(df)

        # Get the top SNPS based on bonferroni corrected significance level 
        # (only regressing on one SNP at a time so the number of tests is the same as the number of p vals)
        try:
            alpha = 0.05
            bonferroni_correction = alpha / sum(np.isfinite(df['p-val']))
            df = df[df['p-val'] <= bonferroni_correction].sort_values(by='p-val')

            top_SNPS = top_SNPS.append(df)

            # add the significance threshold used
            results['Bonferroni_correction'] = bonferroni_correction
            top_SNPS['Bonferroni_correction'] = bonferroni_correction

        except:
            print('Error: No p-vals found for ' + file)
            print(df['p-val'])
            continue

    # Save the results to a file
    results.to_csv(output, sep = '\t', index=False)

    # Save the top SNPs (based on pval) to a file
    top_SNPS.to_csv(top_SNPS_file, sep = '\t', index=False)


if __name__ == "__main__":
    main()