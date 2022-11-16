import mr
import pandas as pd
import mr_utils
import utility_funcs
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--variants_file', help = 'File containing variants from PheWAS', required = False, type = str)
    parser.add_argument('--traits', help = 'List of traits to search for in OpenGWAS to use for MR', required = False, type = str, action ='append', nargs = '+')
    parser.add_argument('--output', help ='Output file', required = False, type = str)
    parser.add_argument('--mr_type', help = 'MR test to run', choices = ['ivw', 'egger', 'simple-median'], required = True, type = str)
    parser.add_argument('--similarity_func', help = 'Similarity function to use for matching trait strings in OpenGWAS', choices = ['jaccard', 'levenshtein'], required = False, type = str)
    parser.add_argument('--cache_all_studies', help = 'Cache all studies in OpenGWAS', required = False, action = 'store_true')
    parser.add_argument('--batches', help = 'GWAS summary dataset batches from Open GWAS', choices=mr_utils.open_gwas_batches.keys(), required = False, type = str, nargs = '+')
    parser.add_argument('--strip_names', help = 'If trying to match traits is not working well try standardizing to strip all non alpha characters from the traits', required = False, action = 'store_true')

    args = parser.parse_args()
    variants_file = args.variants_file
    traits = args.traits
    output = args.output
    mr_type = args.mr_type
    similarity_function = args.similarity_func
    cache_all_studies = args.cache_all_studies
    batches = args.batches
    strip_names = args.strip_names

    if cache_all_studies:
        mr_utils.cache_all_studies()
    
    if similarity_function != None:
        similarity_func = mr_utils.get_similarity_func(similarity_function)
    else:
        similarity_func = mr_utils.get_similarity_func('levenshtein')

    variants = pd.read_csv(variants_file, sep = '\t')

    variants['CHR'] = variants['CHR'].replace('X', 23).replace('Y', 24)

    if len(traits) >= 1:
        for trait_family in traits:
            for trait in trait_family:

                traits_to_study_mapping = mr_utils.find_studies_based_on_traits([[trait]], similarity = similarity_func, batch_list = batches, strip = strip_names)

                external_gwas_data = mr_utils.extract_snps_from_outcomes(variants['rsID'].to_list(), list(traits_to_study_mapping.values())[0], proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
                
                data = {}

                for phenotype in external_gwas_data['phenotype'].unique():
                    external_gwas_data_i = external_gwas_data.loc[external_gwas_data['phenotype'] == phenotype]
                    data[phenotype] = mr_utils.harmonize(variants, external_gwas_data_i, '_heart', '_' + phenotype)
                    print('Running MR for the {} phenotype'.format(phenotype))
                    pval, beta, std_err, *additional = mr.run_mr(mr_type, data[phenotype], 'BETA_heart', 'BETA_' + phenotype, 'SE_heart', 'SE_' + phenotype)

                    print('MR {} p-value: {}, beta: {}, std_err: {}\n'.format(mr_type.capitalize(), pval, beta, std_err))


    

    


if __name__ == '__main__':
    main()