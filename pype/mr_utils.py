import re
import numpy as np
import pandas as pd
import utility_funcs
from collections import defaultdict
from difflib import SequenceMatcher
from ieugwaspy.backwards import legacy_ids
from ieugwaspy.query import gwasinfo, tophits, associations

open_gwas_batches = {"bbj-a" : "Biobank Japan release of disease traits",
                    "ebi-a" : "Datasets that satisfy minimum requirements imported from the EBI database of complete GWAS summary data",
                    "eqtl-a" : "eQTLGen 2019 results (comprising all cis and some trans regions of gene expression in whole blood)",
                    "finn-b" : "FinnGen biobank analysis round 5",
                    "ieu-a" : "GWAS summary datasets generated by many different consortia that have been manually collected and curated (initially developed for MR-Base)",
                    "ieu-b" : "GWAS summary datasets generated by many different consortia that have been manually collected and curated (initially developed for MR-Base (round 2))",
                    "met-a" : "Human blood metabolites analysed by Shin et al 2014",
                    "met-b" : "Human immune system traits analysed by Roederer et al 2015",
                    "met-c" : "Circulating metabolites analysed by Kettunen et al 2016",
                    "met-d" : "Metabolic biomarkers in the UK Biobank measured by Nightingale Health 2020",
                    "prot-a" : "Complete GWAS summary data on protein levels as described by Sun et al 2018",
                    "prot-b" : "Complete GWAS summary data on protein levels as described by Folkersen et al 2017",
                    "prot-c" : "Complete GWAS summary data on protein levels as described by Suhre et al 2017",
                    "ubm-a" : "Complete GWAS summary data on brain region volumes as described by Elliott et al 2018",
                    "ukb-a" : "Neale lab analysis of UK Biobank phenotypes (round 1)",
                    "ukb-b" : "IEU analysis of UK Biobank phenotypes",
                    "ukb-d" : "Neale lab analysis of UK Biobank phenotypes (round 2)",
                    "ukb-e" : "Pan-ancestry genetic analysis of the UK Biobank performed at the Broad Institute"}

def findWord(word, phrase):
    return re.compile(r'\b({0})\b'.format(word), flags=re.IGNORECASE).search(phrase)

def setAllStudiesPath(path):
    global ALL_STUDIES_PATH
    ALL_STUDIES_PATH = path

def cache_all_studies():
    all_studies = gwasinfo()
    utility_funcs.pickle_object(all_studies, ALL_STUDIES_PATH)

def get_all_studies():
    try:
        return utility_funcs.unpickle_object(ALL_STUDIES_PATH)
    except:
        raise ValueError('No cached studies found. Please run_mr with --cache_all_studies flag')

def get_similarity_func(func):
    if func == 'jaccard':
        return lambda x, y: len(set.intersection(*[set(x), set(y)])) / len(set.union(*[set(x), set(y)]))
    elif func == 'levenshtein':
        return lambda x, y: SequenceMatcher(None, x, y).ratio()
    else:
        raise ValueError('Invalid similarity function')

def extract_instruments(outcome_dict, p1 = 5e-8, clump = 1, r2 = 0.001, kb = 10, access_token = "NULL", force_server = False):
    
    data_list = []
    
    for phenotype, studies in outcome_dict.items():
    
        outcomes = legacy_ids(list(set(studies)))
    
        snp_top_hits = tophits(outcomes, pval = p1, clump = clump, r2=r2, kb = kb*1000, force_server = force_server, access_token = access_token)

        if len(snp_top_hits) == 0:
            continue

        for snp in snp_top_hits:
            data_list.append([phenotype,
                            snp['rsid'], 
                            snp['chr'],
                            snp['position'],
                            snp['nea'],
                            snp['ea'],
                            snp['beta'],
                            snp['p'],
                            snp['se'],
                            snp['n'],
                            snp['eaf'],
                            snp['id']])
            
    data = pd.DataFrame(data_list, columns = ['phenotype', 'rsID', 'CHR', 'POS', 'Non_Effect', 'Effect', 'BETA', 'P', 'SE', 'N', 'Effect_Allele_Frequency', 'Study'])
    
    data = data.replace('', np.nan)
    
    data['CHR'] = data['CHR'].replace('X', 23).replace('Y', 24)
    
    data[['CHR', 'POS', 'N']] = data[['CHR', 'POS', 'N']].astype('float').astype("Int64")
    data[['BETA', 'P', 'SE', 'Effect_Allele_Frequency']] = data[['BETA', 'P', 'SE', 'Effect_Allele_Frequency']].astype("float64")
    return data

def find_studies_based_on_traits(traits, similarity, batch_list = None, strip = False):
    
    trait_to_study = defaultdict(list)

    all_studies = get_all_studies()

    regex = re.compile('[^a-zA-Z]')
    strip_str = lambda string: regex.sub('', string).lower()
    
    study_found_in_batch = lambda study: len([1 for batch in batch_list if batch in study]) != 0
    
    for study_name, study_data in all_studies.items():
        if strip:
            study_trait = strip_str(study_data['trait'])
        else:
            study_trait = study_data['trait']
        if batch_list != None:
            if not study_found_in_batch(study_name):
                continue
        
        for trait_list in traits:
            
            trait = trait_list[0]
            
            for alternate_trait_name in trait_list: 
            
                if strip:
                    fixed_trait = strip_str(alternate_trait_name)
                else:
                    fixed_trait = alternate_trait_name
                    
                if similarity(fixed_trait, study_trait) > 0.8 or findWord(fixed_trait, study_trait):
                    trait_to_study[trait].append(study_name)
                
    return trait_to_study

def harmonize(sample1, sample2, suffix1, suffix2):
    
    # TODO: add more in depth harmonization capabilties

    sample1['CHR'] = sample1['CHR'].replace('X', 23).replace('Y', 24)
    sample2['CHR'] = sample2['CHR'].replace('X', 23).replace('Y', 24)

    # drop uneeded chr and pos cols
    s1 = sample1.drop(['POS'], axis = 1, errors = 'ignore')
    s2 = sample2.drop(['POS'], axis = 1, errors = 'ignore')
    
    # merge the two dataframes, ensuring the chromosomes/position as well as alleles are the same
    merged = pd.merge(s1, s2, on=['rsID', 'CHR'], suffixes = [suffix1, suffix2]) # on=['rsID', 'CHR', 'Non_Effect', 'Effect']
    
    # rename the columns (for use with ldsc)
    # merged = merged.rename(columns = {'Effect': 'A1',
                                     # 'Non_Effect': 'A2'})
    
    # retain only the columns we care about in this order
    if 'N' in sample1.columns or 'N' in sample2.columns:
        merged = merged[['rsID', 'CHR', 'BETA' + suffix1, 'BETA' + suffix2, 'P' + suffix1, 'P' + suffix2, 'SE' + suffix1, 'SE' + suffix2, 'N']]
        
        merged['N'] = merged['N'].replace(np.nan, -1)

        # drop duplicates, keeping the values with greatest number of samples (if data is available)
        merged = merged.loc[merged.groupby('rsID')['N'].idxmax()].reset_index(drop = True)
    else:
        merged = merged[['rsID', 'CHR', 'BETA' + suffix1, 'BETA' + suffix2, 'P' + suffix1, 'P' + suffix2, 'SE' + suffix1, 'SE' + suffix2]]
   
    # if the dataframe is empty, we print an error message
    if merged.empty:
        print("The dataframe is empty")

    return merged


def extract_snps_from_outcomes(snps, outcome_studies, proxies = True, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = "NULL", splitsize=10000, proxy_splitsize = 500, verbose = False):    
    '''
    Extracts the target SNPs from the outcome studies and returns a dataframe with the SNPs and their proxies

    Parameters
    ----------
    snps : list
        List of SNPs to extract
    outcome_studies : list
        List of outcome studies to extract the SNPs from
    proxies : bool, optional
        Whether to extract proxies for the SNPs. The default is True.
    rsq : float, optional
        The minimum r2 to use for proxy extraction. The default is 0.8.
    align_alleles : bool, optional
        Whether to align the alleles of the proxies to the target SNPs. The default is 1.
    palindromes : bool, optional
        Whether to include palindromic SNPs in the proxy extraction. The default is 1.
    maf_threshold : float, optional
        The minimum MAF to use for proxy extraction. The default is 0.3.
    access_token : str, optional
        The access token to use for the API. The default is "NULL".
    splitsize : int, optional
        The number of SNPs to extract at a time. The default is 10000.
    proxy_splitsize : int, optional
        The number of SNPs to extract proxies for at a time. The default is 500.
    verbose : bool, optional
        Whether to print the progress of the extraction. The default is False.

    Returns
    -------
    snp_df : pd.DataFrame
        A dataframe containing the target SNPs from the outcome studies and their proxies
    '''

    snps = list(set(snps))
    
    firstpass = extract_snps_from_outcomes_internal(snps, outcome_studies, proxies = False, access_token=access_token, splitsize = splitsize, verbose = verbose)
    
    if proxies:
        for outcome in outcome_studies:
            if firstpass is not None:
                firstpass_outcomei_rows = firstpass.loc[firstpass['id'] == outcome]
                missed_snps = [snp for snp in snps if snp not in firstpass_outcomei_rows['rsid'].unique()]
            else:
                missed_snps = snps
            if len(missed_snps) > 0:
                if verbose:
                    print("Finding proxies for " + str(len(missed_snps)) + " SNPs in outcome " + outcome)
                temp = extract_snps_from_outcomes_internal(missed_snps, [outcome], True, rsq, align_alleles, palindromes, maf_threshold, access_token, proxy_splitsize, verbose)
                
                if temp is not None:
                    firstpass = pd.concat([firstpass, temp])
    
    # change the column names to match the format of our expected columns
    firstpass = firstpass.rename(columns = {'chr':'CHR', 'position': 'POS', 'beta': 'BETA', 'se': 'SE', 'p':'P', 'n': 'N', 'id': 'Study', 'rsid':'rsID', 'ea': 'Effect', 'nea':'Non_Effect', 'eaf':'Effect_Allele_Frequency', 'trait': 'phenotype'})

    firstpass = firstpass.replace('', np.nan)
    firstpass['CHR'] = firstpass['CHR'].replace('X', 23).replace('Y', 24)
    firstpass[['CHR', 'POS', 'N']] = firstpass[['CHR', 'POS', 'N']].astype('float').astype("Int64")
    firstpass[['BETA', 'P', 'SE', 'Effect_Allele_Frequency']] = firstpass[['BETA', 'P', 'SE', 'Effect_Allele_Frequency']].astype("float64")
    
    return firstpass

def extract_snps_from_outcomes_internal(snps, outcome_studies, proxies = True, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, access_token = "NULL", splitsize=10000, verbose = False): 
    
    snps = list(set(snps))
    outcomes = list(set(outcome_studies))
    
    if verbose:
        print("Extracting data for " + str(len(snps)) + " SNP(s) from " + str(len(outcome_studies)) + " GWAS(s)")
    
    if proxies:
        proxies = 1
    else:
        proxies = 0
    
    assoc = pd.DataFrame()
    if len(snps) < splitsize and len(outcomes) < splitsize:
        assoc = pd.DataFrame(associations(variantlist = snps,
                             id = outcomes,
                             proxies = proxies,
                             r2 = rsq,
                             align_alleles = align_alleles,
                             palindromes = palindromes,
                             maf_threshold = maf_threshold,
                             access_token = access_token))
    elif len(snps) > len(outcomes):
        snps_chunked = [snps[i:i + splitsize] for i in range(0, len(snps), splitsize)] 
        
        assoc = []
        
        for outcome in outcomes:
            if verbose:
                print("Outcome: " + outcome)
            for snp_chunk in snps_chunked:
                assoc.append(pd.DataFrame(associations(variantlist = snp_chunk,
                                                       id = outcome,
                                                       proxies = proxies,
                                                       r2 = rsq,
                                                       align_alleles = align_alleles,
                                                       palindromes = palindromes,
                                                       maf_threshold = maf_threshold,
                                                       access_token = access_token)))
        assoc = pd.concat(assoc)
    else: # len(outcomes) > len(snps)
        outcomes_chunked = [outcomes[i:i + splitsize] for i in range(0, len(outcomes), splitsize)] 
        
        assoc = []
        
        for snp in snps:
            if verbose:
                print("Outcome: " + outcome)
            for outcome_chunk in outcomes_chunked:
                assoc.append(pd.DataFrame(associations(variantlist = snp,
                                                       id = outcome_chunk,
                                                       proxies = proxies,
                                                       r2 = rsq,
                                                       align_alleles = align_alleles,
                                                       palindromes = palindromes,
                                                       maf_threshold = maf_threshold,
                                                       access_token = access_token)))
        assoc = pd.concat(assoc)
    if assoc.empty:
        return None
    else:
        return assoc

