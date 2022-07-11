import argparse
import pandas as pd
import numpy as np
import math
from datetime import datetime
import os

# adapted code from Alan's MI_Classes.py file
dict_UKB_fields_to_names =  {'f.31.0.0': 'Sex', 
                            'f.34.0.0': 'Year_of_birth', 
                            'f.52.0.0': 'Month_of_birth',
                            'f.53.0.0': 'Date_attended_center_0', 
                            'f.53.1.0': 'Date_attended_center_1',
                            'f.53.2.0': 'Date_attended_center_2', 
                            'f.53.3.0': 'Date_attended_center_3',
                            'f.21000.0.0': 'Ethnicity', 
                            'f.21000.1.0': 'Ethnicity_1', 
                            'f.21000.2.0': 'Ethnicity_2',
                            'f.22001.0.0': 'Sex_genetic'}

def compute_sex(pheno):
    # Use genetic sex when available
    pheno['Sex_genetic'].fillna(pheno['Sex'], inplace = True)

    pheno.drop(['Sex'], axis=1, inplace=True)
    pheno.rename(columns={'Sex_genetic': 'Sex'}, inplace=True)
    pheno.dropna(subset=['Sex'], inplace=True)

    return pheno

def compute_age(pheno):
    # Recompute age with greater precision by leveraging the month of birth
    pheno['Year_of_birth'] = pheno['Year_of_birth'].astype(int)
    pheno['Month_of_birth'] = pheno['Month_of_birth'].astype(int)
    pheno['Date_of_birth'] = pheno.apply(
        lambda row: datetime(row.Year_of_birth, row.Month_of_birth, 15), axis=1)
    for i in range(4):
        i = str(i)
        pheno['Date_attended_center_' + i] = \
            pheno['Date_attended_center_' + i].apply(
                lambda x: pd.NaT if pd.isna(x) else datetime.strptime(x, '%Y-%m-%d'))
        pheno['Age_' + i] = pheno['Date_attended_center_' + i] - pheno['Date_of_birth']
        pheno['Age_' + i] = pheno['Age_' + i].dt.days / 365.25
        pheno.drop(['Date_attended_center_' + i], axis=1, inplace=True)
        
        # Save age as a float32 instead of float64
        pheno['Age_' + i] = np.float32(pheno['Age_' + i])

    pheno.drop(['Year_of_birth', 'Month_of_birth', 'Date_of_birth'], axis=1, inplace=True)
    pheno.dropna(how='all', subset=['Age_0', 'Age_1','Age_2', 'Age_3'], inplace=True)

    return pheno

def encode_ethnicity(pheno):
    # Fill NAs for ethnicity on instance 0 if available in other instances
    eids_missing_ethnicity = pheno['f.eid'][pheno['Ethnicity'].isna()]
    for eid in eids_missing_ethnicity:
        sample = pheno.loc[eid, :]
        if not math.isnan(sample['Ethnicity_1']):
            pheno.loc[eid, 'Ethnicity'] = pheno.loc[eid, 'Ethnicity_1']
        elif not math.isnan(sample['Ethnicity_2']):
            pheno.loc[eid, 'Ethnicity'] = pheno.loc[eid, 'Ethnicity_2']
    pheno.drop(['Ethnicity_1', 'Ethnicity_2'], axis=1, inplace=True)
    
    # One hot encode ethnicity
    dict_ethnicity_codes = {'1': 'Ethnicity.White', '1001': 'Ethnicity.British', '1002': 'Ethnicity.Irish',
                            '1003': 'Ethnicity.White_Other',
                            '2': 'Ethnicity.Mixed', '2001': 'Ethnicity.White_and_Black_Caribbean',
                            '2002': 'Ethnicity.White_and_Black_African',
                            '2003': 'Ethnicity.White_and_Asian', '2004': 'Ethnicity.Mixed_Other',
                            '3': 'Ethnicity.Asian', '3001': 'Ethnicity.Indian', '3002': 'Ethnicity.Pakistani',
                            '3003': 'Ethnicity.Bangladeshi', '3004': 'Ethnicity.Asian_Other',
                            '4': 'Ethnicity.Black', '4001': 'Ethnicity.Caribbean', '4002': 'Ethnicity.African',
                            '4003': 'Ethnicity.Black_Other',
                            '5': 'Ethnicity.Chinese',
                            '6': 'Ethnicity.Other_ethnicity',
                            '-1': 'Ethnicity.Do_not_know',
                            '-3': 'Ethnicity.Prefer_not_to_answer',
                            '-5': 'Ethnicity.NA'}
    pheno['Ethnicity'] = pheno['Ethnicity'].fillna(-5).astype(int).astype(str)
    ethnicities = pd.get_dummies(pheno['Ethnicity'])
    pheno.drop(['Ethnicity'], axis=1, inplace=True)
    ethnicities.rename(columns=dict_ethnicity_codes, inplace=True)
    ethnicities['Ethnicity.White'] = ethnicities['Ethnicity.White'] + ethnicities['Ethnicity.British'] + \
                                    ethnicities['Ethnicity.Irish'] + ethnicities['Ethnicity.White_Other']
    ethnicities['Ethnicity.Mixed'] = ethnicities['Ethnicity.Mixed'] + \
                                    ethnicities['Ethnicity.White_and_Black_Caribbean'] + \
                                    ethnicities['Ethnicity.White_and_Black_African'] + \
                                    ethnicities['Ethnicity.White_and_Asian'] + \
                                    ethnicities['Ethnicity.Mixed_Other']
    ethnicities['Ethnicity.Asian'] = ethnicities['Ethnicity.Asian'] + ethnicities['Ethnicity.Indian'] + \
                                    ethnicities['Ethnicity.Pakistani'] + ethnicities['Ethnicity.Bangladeshi'] + \
                                    ethnicities['Ethnicity.Asian_Other']
    ethnicities['Ethnicity.Black'] = ethnicities['Ethnicity.Black'] + ethnicities['Ethnicity.Caribbean'] + \
                                    ethnicities['Ethnicity.African'] + ethnicities['Ethnicity.Black_Other']
    ethnicities['Ethnicity.Other'] = ethnicities['Ethnicity.Other_ethnicity'] + \
                                    ethnicities['Ethnicity.Do_not_know'] + \
                                    ethnicities['Ethnicity.Prefer_not_to_answer'] + \
                                    ethnicities['Ethnicity.NA']
    pheno = pheno.join(ethnicities)
    return pheno

def main():
    parser = argparse.ArgumentParser(description='Fix covariates')
    parser.add_argument('-p', '--pheno_files', required=True, help='Phenotype files', action = 'append')
    parser.add_argument('-f', '--fields_file', required=True, help = 'Phenotype fields')
    parser.add_argument('-o', '--out_dir', required=True, help = 'Phenotype directory')
    parser.add_argument('-d', '--data_fields_dir', required=True, help = 'Data fields directory')
    args = parser.parse_args()

    # parse arguments
    pheno_files = args.pheno_files
    phenotype_fields_file = args.fields_file
    out_dir = args.out_dir
    data_fields_dir = args.data_fields_dir

    with open(phenotype_fields_file) as f:
        phenotype_fields = f.read().splitlines()

    for file in pheno_files:
        file_prefix = os.path.basename(file).split('_')[0]

        pheno = pd.read_csv(file, sep='\t', low_memory = False, header = 0)

        # Formatting
        pheno.rename(columns=dict_UKB_fields_to_names, inplace=True)
        pheno.set_index('f.eid', drop=False, inplace=True)
        pheno.index.name = 'column_names'

        pheno = compute_sex(pheno)
        pheno = compute_age(pheno)
        pheno = encode_ethnicity(pheno)

        # set of all covariates, including the new covariate related columns added
        all_covariates = [col for col in pheno.columns if col not in phenotype_fields + ['f.eid']]

        # keep all columns with more than 80% non-nan values
        remove_columns = pheno.loc[:, all_covariates].dropna(thresh = len(pheno[all_covariates])*0.2, axis = 1).columns
        dropped_cols = [col for col in all_covariates if col not in remove_columns]
        cols_to_keep = [col for col in pheno if col not in dropped_cols]

        print(file_prefix, cols_to_keep)

        remaining_covariates = [col for col in all_covariates if col not in dropped_cols] 
        # fixed covariate file path
        covariate_file = os.path.join(data_fields_dir, file_prefix + '_fixed_cov_fields.txt')
        with open(covariate_file, 'w') as f:
            f.write('\n'.join(remaining_covariates))
        
        pheno = pheno[cols_to_keep].reset_index().drop('column_names', axis = 1)

        pheno.to_csv(out_dir + '/' + file_prefix + '_subsetted_pheno_covariates_fixed.tab', index = False, sep = '\t')

if __name__ == "__main__":
    main()