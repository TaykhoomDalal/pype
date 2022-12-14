import numpy as np
from scipy.stats import norm, chi2, t
from statsmodels.formula.api import wls


def run_mr(mr_type, harmonized_data, beta_exp, beta_out, se_exp, se_out):
    '''
    Run MR using the specified methods. 

    Parameters
    ----------
    mr_type : list
        List of MR methods to run.
    harmonized_data : pandas.DataFrame
        Dataframe containing the harmonized data.
    beta_exp : str
        Name of the column containing the exposure effect sizes.
    beta_out : str
        Name of the column containing the outcome effect sizes.
    se_exp : str
        Name of the column containing the exposure standard errors.
    se_out : str    
        Name of the column containing the outcome standard errors.
    
    Returns
    -------
    mr_results : dict
        Dictionary containing the MR results.
    
    '''

    mr_results = {}
    
    for method in mr_type:
        if method == 'ivw':
            mr_results[method] = run_mr_ivw(harmonized_data, beta_exp, beta_out, se_exp, se_out)
        elif method == 'simple_median':
            mr_results[method] = run_mr_simple_median(harmonized_data, beta_exp, beta_out, se_exp, se_out)
        elif method == 'egger':
            mr_results[method] = run_mr_egger(harmonized_data, beta_exp, beta_out, se_exp, se_out)
        
    return mr_results

def run_mr_ivw(harmonized_data, beta_exp, beta_out, se_exp, se_out):
    '''
    Calculate the inverse variance weighted MR estimate.

    Parameters
    ----------
    harmonized_data : pandas.DataFrame
        Dataframe containing the harmonized data.
    beta_exp : str
        Name of the column containing the exposure effect sizes.
    beta_out : str
        Name of the column containing the outcome effect sizes.
    se_exp : str    
        Name of the column containing the exposure standard errors.
    se_out : str    
        Name of the column containing the outcome standard errors.
    
    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''



    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 2:
        print('Found {} variants. Need at least 2. Exiting MR.'.format(data[['BETA_EXP']].shape[0]))
        return None
    
    f = 'BETA_OUT ~ -1 + BETA_EXP'
    
    fitted_model = wls(formula = f, data = data, weights = 1/data['SE_OUT']**2).fit()
    
    beta = fitted_model.params.BETA_EXP
    standard_error = fitted_model.bse.BETA_EXP
    model_pvalue = 2*(1 - norm.cdf(abs(beta/standard_error))) #fitted_model.pvalues.BETA_EXP
    confidence_interval = fitted_model.conf_int()
    Cochrans_q_degrees_freedom = exp_beta_len - 1
    residual_standard_deviation = np.sqrt(fitted_model.scale)
    Cochrans_q = residual_standard_deviation**2 * Cochrans_q_degrees_freedom
    Cochrans_q_pval = 1 - chi2.cdf(Cochrans_q, Cochrans_q_degrees_freedom)
    i_squared = (Cochrans_q - Cochrans_q_degrees_freedom) / Cochrans_q
    
    return [model_pvalue, beta, standard_error, Cochrans_q, Cochrans_q_degrees_freedom, Cochrans_q_pval, i_squared]

def run_mr_egger(harmonized_data, beta_exp, beta_out, se_exp, se_out):
    '''
    Calculate the Egger MR estimate.

    Parameters
    ----------
    harmonized_data : pandas.DataFrame
        Dataframe containing the harmonized data.
    beta_exp : str
        Name of the column containing the exposure effect sizes.
    beta_out : str
        Name of the column containing the outcome effect sizes.
    se_exp : str    
        Name of the column containing the exposure standard errors.
    se_out : str    
        Name of the column containing the outcome standard errors.

    Returns
    ------- 
    list
        List containing the model-pvalue, beta, and standard error.
    '''


    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR.'.format(data[['BETA_EXP']].shape[0]))
        return None
    
    change_0_to_1 = lambda x: np.sign(x.replace(0,1))
    
    data['BETA_OUT'] = data['BETA_OUT'] * change_0_to_1(data['BETA_EXP'])
    data['BETA_EXP'] = abs(data['BETA_EXP'])
    
    f = 'BETA_OUT ~ BETA_EXP'
    
    fitted_model = wls(formula = f, data = data, weights = 1/data['SE_OUT']**2).fit()
    
    beta = fitted_model.params.BETA_EXP
    standard_error = fitted_model.bse.BETA_EXP
    model_pvalue = 2*(1 - t.cdf(abs(beta/standard_error), exp_beta_len - 2)) #fitted_model.pvalues.BETA_EXP
    confidence_interval = fitted_model.conf_int().T['BETA_EXP']
    
    beta_intercept = fitted_model.params.Intercept
    standard_error_intercept = fitted_model.bse.Intercept
    pvalue_intercept = 2*(1 - t.cdf(abs(beta_intercept/standard_error_intercept), exp_beta_len - 2)) #fitted_model.pvalues.Intercept
    confidence_interval_intercept = fitted_model.conf_int().T['Intercept']
    
    Cochrans_q_degrees_freedom = exp_beta_len - 2
    residual_standard_deviation = np.sqrt(fitted_model.scale)
    
    Cochrans_q = residual_standard_deviation**2 * Cochrans_q_degrees_freedom
    Cochrans_q_pval = 1 - chi2.cdf(Cochrans_q, Cochrans_q_degrees_freedom)
    i_squared = (Cochrans_q - Cochrans_q_degrees_freedom) / Cochrans_q
    
    return [model_pvalue, 
            beta, standard_error, 
            pvalue_intercept, beta_intercept,
            standard_error_intercept,
            Cochrans_q, Cochrans_q_degrees_freedom, 
            Cochrans_q_pval, i_squared]

def run_mr_simple_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000):
    '''
    Calculate the simple median MR estimate.

    Parameters
    ----------
    harmonized_data : pandas.DataFrame
        Dataframe containing the harmonized data.
    beta_exp : str
        Name of the column containing the exposure effect sizes.
    beta_out : str
        Name of the column containing the outcome effect sizes.
    se_exp : str    
        Name of the column containing the exposure standard errors.
    se_out : str    
        Name of the column containing the outcome standard errors.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.

    '''
    
    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR.'.format(data[['BETA_EXP']].shape[0]))
        return None
    
    beta_iv = np.array(data['BETA_OUT'] / data['BETA_EXP'])
    
    inv_rep_len = np.repeat(1/exp_beta_len, exp_beta_len)
    
    beta = weighted_median(beta_iv, inv_rep_len)
    standard_error = weighted_median_bootstrap(data['BETA_EXP'], data['BETA_OUT'], data['SE_EXP'], data['SE_OUT'], inv_rep_len, nboot)
    model_pvalue = 2*(1 - t.cdf(abs(beta/standard_error), exp_beta_len - 2))
    
    return [model_pvalue, beta, standard_error]

def weighted_median(beta_iv, weights):
    '''
    Calculate the weighted median of a list of numbers.
    credit:
    '''
    i = np.argsort(beta_iv)
    c = np.cumsum(weights[i])
    return beta_iv[i[np.searchsorted(c, 0.5 * c[-1])]]
    
def weighted_median_bootstrap(beta_exp, beta_out, se_exp, se_out, weights, nboot):
    '''
    Calculate standard error of the weighted median using bootstrap.

    Parameters
    ----------
    beta_exp : pandas.Series
        Series containing the exposure effect sizes.
    beta_out : pandas.Series
        Series containing the outcome effect sizes.
    se_exp : pandas.Series
        Series containing the exposure standard errors.
    se_out : pandas.Series
        Series containing the outcome standard errors.
    weights : numpy.array
        Array containing the weights for each variant.
    nboot : int
        Number of bootstrap iterations.

    Returns
    -------
    float
        Standard error of the weighted median.

    '''
    # make list of size nboot filled with 0
    medians = [0] * nboot
    
    for i in range(nboot):
        b_exp_bstrapped = np.random.normal(loc = beta_exp, scale = se_exp, size = beta_exp.shape[0])
        b_out_bstrapped = np.random.normal(loc = beta_out, scale = se_out, size = beta_out.shape[0])
        beta_iv_bstrapped = np.array(b_out_bstrapped) / np.array(b_exp_bstrapped)

        # set output of weighted_median to medians[i]
        medians[i] = weighted_median(np.array(beta_iv_bstrapped), np.array(weights))

    
    return np.std(medians)
    
def run_mr_median(harmonized_data, beta_exp, beta_out, se_exp, se_out):
    '''
    This function runs all the different MR-Median variations including the simple median, the weighted median, and the weighted median with penalization.

    Parameters
    ----------
    harmonized_data : pandas dataframe
        Dataframe containing the harmonized data.
    beta_exp : string
        Name of the column containing the effect sizes of the exposure.
    beta_out : string
        Name of the column containing the effect sizes of the outcome.
    se_exp : string
        Name of the column containing the standard errors of the exposure.
    se_out : string
        Name of the column containing the standard errors of the outcome.

    Returns
    -------
    '''
    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR.'.format(data[['BETA_EXP']].shape[0]))
        return None