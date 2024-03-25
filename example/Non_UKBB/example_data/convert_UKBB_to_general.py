import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
import re

def query_ukbb(url):
    # Get the html file
    html = requests.get(url).text

    # Create beautiful soup object which can be used to parse the html file
    soup = BeautifulSoup(html, 'html.parser')

    # Find the <td> tag that contains the text "Category:"
    category_td = soup.find("td", text="Category:")
    
    # Initialize variables for description and last category
    description_text = None
    last_category_text = None

    # Assuming the structure is consistent and the category value is always in the next <td>
    if category_td:
        category_value_td = category_td.find_next_sibling("td")
        if category_value_td:
            # Find all <a> tags within the category_value_td
            a_tags = category_value_td.find_all("a")
            if a_tags:
                # Extract text of the last <a> tag
                last_category_text = a_tags[-1].text
                
    # Find the <td> tag that contains the text "Description:"
    description_td = soup.find("td", text="Description:")
    if description_td:
        value_td = description_td.find_next_sibling("td")
        if value_td:
            description_text = value_td.text

    # Return both the description and last category text
    return description_text, last_category_text

def main():
    UKBB_DATA = pd.read_csv("../../UKBB/example_data/phenotype_files/example.tab", sep="\t", low_memory=False)
    
    url = "https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=X"
    # Get the data fields from the UKBB website
    new_columns = {}
    seen_fields = set()
    pheno_to_category = {}

    for full_field in UKBB_DATA.columns:
        # split the field at the period twice, once to get ['f', 'X.X.X'], and 
        # the second time to get ['f', 'X', 'X.X'], the first X represents the field, the 
		# second X represents the instance index, and the third X represents the array index
        field = full_field.split(".", 2)[1]
	
        if field in seen_fields:
            continue
        else:
            seen_fields.add(field)

        if field == 'eid':
            new_columns[full_field] = 'ID'
            continue
		
        field_description, category = query_ukbb(url.replace('id=X', 'id=' + field))

        if field_description is None:
            print("Field description not found for field: ", field)
        else:
            print(field_description, "\t", category)
            # get the field description from the UKBB website
            new_columns[full_field]  = field_description
            pheno_to_category[field_description] = category
    
    # Rename the columns
    UKBB_DATA.rename(columns=new_columns, inplace=True)

    # remove all columns that start with 'f.' since they are not useful
    UKBB_DATA = UKBB_DATA[[col for col in UKBB_DATA.columns if not col.startswith('f.')]]

    # Save the data to a new file
    UKBB_DATA.to_csv("phenotype_files/example.tsv", index=False, sep = '\t')

    # Save the phenotype to category mapping to a new file as a dataframe with two columns 'Phenotype' and 'Category'
    pheno_to_category_df = pd.DataFrame(list(pheno_to_category.items()), columns=['Phenotype', 'Category'])

    pheno_to_category_df.to_csv("phenotype_files/example_mapper.tsv", index=False, sep = '\t')
        
if __name__ == '__main__':
	main()