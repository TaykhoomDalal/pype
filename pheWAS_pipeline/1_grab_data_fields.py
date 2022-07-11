import argparse
import os
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup
import re
import requests

def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Pass the input html file')
    parser.add_argument('-u', '--url', help='Url of file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)

    args = parser.parse_args()
    url = args.url
    output = args.output

    # Get the html file
    html = requests.get(url).text

    # Create beautiful soup object which can be used to parse the html file
    soup = BeautifulSoup(html, 'html.parser')
    
    # Get all tags which correspond to the data fields for ukbb
    html_tags = soup.find_all("a", class_="basic", href=re.compile("field.cgi"))

    # parse and grab the actual data-fields
    fields = [tag.get('href').split("=")[1] for tag in html_tags]

    # write the data-fields to the output file
    with open(output, 'w') as f:
        f.write("\n".join(fields))

if __name__ == '__main__':
    main()