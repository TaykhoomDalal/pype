import pandas as pd
import os
import argparse

def extract_ukb_data_fields(col_names, fields_to_keep, file, write_file = False):
    """
    This function takes a list of headers from a phenotype file and a list of fields to keep and returns a list with the UKBB data fields.
    It also writes the list to a file for future reference
    """
    # get the column names that are in the list of columns to keep
    
    exact_col_names = []
    for field in fields_to_keep:

        # data fields are in the form of "f.field.[etc]"
        field = "f." + str(field) + "." 
        for col_name in col_names:
            if field in col_name:
                exact_col_names.append(col_name)
    
    # ensure columns are unique (use dict.fromkeys to remove duplicates - keeps order as compared to set())
    exact_col_names = list(dict.fromkeys(exact_col_names))
    
    if write_file:
        # write the list of columns to a file
        with open(file, "w") as data_file:
            data_file.write("\n".join(exact_col_names))

    # return the list with the fields to keep
    return exact_col_names
    
def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Extract the individuals who had an MRI')
    parser.add_argument('-s', '--sample_files', help='List of samples', action='append', required=True)
    parser.add_argument('-p', '--pheno_file', help='PheWAS phenotype file', required=True)
    parser.add_argument('-f', '--fields', help='List of data fields to keep', required=True)
    parser.add_argument('-c', '--covariates_file', help='Covariate file', required=True)
    parser.add_argument('-o', '--out_dir', help='Directory of where to store the phenotype data', required=True)
    parser.add_argument('-k', '--keep', help ='Whether the samples are in / out of sample', action='store_true')
    parser.add_argument('-d', '--data_fields_dir', help ='Directory of where to store the data fields', required=True)
    args = parser.parse_args()

    # Parse the command line arguments
    sample_files = args.sample_files
    pheno_file = args.pheno_file
    fields = args.fields
    covariates_file = args.covariates_file    
    out_dir = args.out_dir
    keep = args.keep
    data_fields_dir = args.data_fields_dir
    
    # read covariates into a list
    with open(covariates_file, "r") as f:
        cov = f.read().splitlines()

    # read the data fields into a list
    with open(fields, "r") as f:
        data_fields = f.read().splitlines()

    data_fields = sorted(data_fields)
    data_fields_with_cov = sorted(data_fields + cov)
    exact_col_names = []
    
    # get the column names
    col_names = pd.read_csv(pheno_file, sep ='\t', nrows = 0).columns.tolist()

    # we only need to keep the phenotypes we want to run the PheWAS on
    exact_col_names = ['f.eid'] + extract_ukb_data_fields(col_names, data_fields_with_cov, "", False)

    # create file with covariate data field names
    extract_ukb_data_fields(col_names, cov, data_fields_dir + "/ukb_cov_fields.txt", True)

    # create file with data field names
    extract_ukb_data_fields(col_names, data_fields, data_fields_dir + "/ukb_data_fields.txt", True)
    
    # full_data is the dataframe with all the data
    index = 0
    predictor_dict = {}
    sample_file_data_list = {}

    for file in sample_files:
        predictor_dict[file] = pd.DataFrame(columns = exact_col_names)
        sample_file_data_list[file] = pd.read_csv(file, sep='\t', header=None).iloc[:,0].tolist()

    # can't read entire .tab file at once, very slow, so read it in chunk by chunk
    chunksize = 1000 # number or rows we see at a time
    with pd.read_csv(pheno_file, sep = '\t', chunksize=chunksize, low_memory = False, usecols = exact_col_names) as reader:
        
        # for each chunk of data
        for chunk in reader:
            
            for file in sample_files:
                
                if keep:
                    data = chunk[chunk['f.eid'].isin(sample_file_data_list[file])] # get the data for the chunk that is in the sample_IDs
                else:
                    data = chunk[~chunk['f.eid'].isin(sample_file_data_list[file])] # get the data for the chunk that is not in sample_IDs  

                prev_data = predictor_dict[file]

                prev_data = pd.concat([prev_data, data]) # append the data to the full data frame

                predictor_dict[file] = prev_data

                print('chunk ' + str(index) +' for file '+ os.path.basename(file) +' done')
            index += 1

    for file in sample_files:
        # write the subsetted phenotype files to a file
        predictor_dict[file].to_csv(out_dir + '/' + os.path.basename(file).split('_')[0] + '_subsetted_pheno.tab', sep='\t', header=True, index=False)



if __name__ == '__main__':
    main()