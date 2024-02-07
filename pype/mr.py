import numpy as np
import pandas as pd
from statsmodels.formula.api import wls
from scipy.fft import fft, ifft
from scipy.interpolate import interp1d
from scipy.stats import norm, chi2, t, gaussian_kde, median_abs_deviation

def run_mr(mr_type, harmonized_data, beta_exp, beta_out, se_exp, se_out, mr_args, all = False):
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
    mr_args : dict
        Dictionary containing the extra arguments for the various MR functions.
    all : boolean
        Boolean indicating whether or not to run all the MR functions on the data.
        
    Returns
    -------
    mr_results : dict
        Dictionary containing the MR results.
    
    '''

    mr_results = {}
    
    if all == True:
        mr_results['Egger'] = run_mr_egger(harmonized_data, beta_exp, beta_out, se_exp, se_out)
        mr_results['Inverse Variance Weighted'] = run_mr_ivw(harmonized_data, beta_exp, beta_out, se_exp, se_out)
        mr_results['Simple Median'] = run_mr_simple_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Simple Median"])
        mr_results['Weighted Median'] = run_mr_weighted_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Weighted Median"])
        mr_results['Penalized Weighted Median'] = run_mr_penalized_weighted_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Penalized Weighted Median"])
        mr_results['Simple Mode'] = run_mr_simple_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Simple Mode"])
        mr_results['Weighted Mode'] = run_mr_weighted_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Weighted Mode"])
        mr_results['Penalized Mode'] = run_mr_penalized_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Penalized Mode"])
        mr_results['Simple Mode (NOME)'] = run_mr_simple_mode_nome(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Simple Mode (NOME)"])
        mr_results['Weighted Mode (NOME)'] = run_mr_weighted_mode_nome(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Weighted Mode (NOME)"])
        mr_results['PRESSO'] = run_mr_presso(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["PRESSO"])
        return mr_results

    for method in mr_type:
        if method == 'ivw':
            mr_results['Inverse Variance Weighted'] = run_mr_ivw(harmonized_data, beta_exp, beta_out, se_exp, se_out)
        elif method == 'egger':
            mr_results['Egger'] = run_mr_egger(harmonized_data, beta_exp, beta_out, se_exp, se_out)
        elif method == 'simple_median':
            mr_results['Simple Median'] = run_mr_simple_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Simple Median"])
        elif method == 'weighted_median':
            mr_results['Weighted Median'] = run_mr_weighted_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Weighted Median"])
        elif method == 'penalized_weighted_median':
            mr_results['Penalized Weighted Median'] = run_mr_penalized_weighted_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Penalized Weighted Median"])
        elif method == 'simple_mode':
            mr_results['Simple Mode'] = run_mr_simple_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Simple Mode"])
        elif method == 'simple_mode_nome':
            mr_results['Simple Mode (NOME)'] = run_mr_simple_mode_nome(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Simple Mode (NOME)"])
        elif method == 'weighted_mode':
            mr_results['Weighted Mode'] = run_mr_weighted_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Weighted Mode"])
        elif method == 'penalized_mode':
            mr_results['Penalized Mode'] = run_mr_penalized_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Penalized Mode"])
        elif method == 'weighted_mode_nome':
            mr_results['Weighted Mode (NOME)'] = run_mr_weighted_mode_nome(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["Weighted Mode (NOME)"])
        elif method == 'presso':
            mr_results['PRESSO'] = run_mr_presso(harmonized_data, beta_exp, beta_out, se_exp, se_out, **mr_args["PRESSO"])
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
    
    print("Running Inverse Variance Weighted (IVW) Regression")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 2:
        print('Found {} variants. Need at least 2. Exiting MR IVW.'.format(data[['BETA_EXP']].shape[0]))
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

    print("Running Egger's Regression")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Egger.'.format(data[['BETA_EXP']].shape[0]))
        return None
    
    to_flip = np.sign(data['BETA_EXP']) == -1
    data['BETA_OUT'] = data['BETA_OUT'] * np.sign(data['BETA_EXP'])
    data['BETA_EXP'] = abs(data['BETA_EXP'])
    
    f = 'BETA_OUT ~ BETA_EXP'
    
    fitted_model = wls(formula = f, data = data, weights = 1/data['SE_OUT']**2).fit()
    
    beta = fitted_model.params.BETA_EXP
    standard_error = fitted_model.bse.BETA_EXP / min(1, fitted_model.scale ** 0.5)
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
    nboot : int
        Number of bootstrap replications to calculate the standard error.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''
    
    print("Running Simple Median")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Simple Median.'.format(data[['BETA_EXP']].shape[0]))
        return None
    
    beta_iv = np.array(data['BETA_OUT'] / data['BETA_EXP'])
    
    inv_rep_len = np.repeat(1/exp_beta_len, exp_beta_len)
    
    beta = weighted_median(beta_iv, inv_rep_len)
    standard_error = weighted_median_bootstrap(data['BETA_EXP'].values, data['BETA_OUT'].values, data['SE_EXP'].values, data['SE_OUT'].values, inv_rep_len, nboot)
    model_pvalue = 2 * (1 - norm.cdf(abs(beta/standard_error)))
    return [model_pvalue, beta, standard_error]

def run_mr_weighted_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000):
    '''
    Calculate the weighted median MR estimate.

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
    nboot : int
        Number of bootstrap replications to calculate the standard error.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''

    print("Running Weighted Median")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Weighted Median.'.format(data[['BETA_EXP']].shape[0]))
        return None

    beta_iv = np.array(data['BETA_OUT'] / data['BETA_EXP'])
    VBj = np.array(((data['SE_OUT'])**2) / (data['BETA_EXP'])**2 + (data['BETA_OUT']**2) * ((data['SE_EXP']**2)) / (data['BETA_EXP'])**4)

    beta = weighted_median(beta_iv, 1/VBj)
    standard_error = weighted_median_bootstrap(data['BETA_EXP'].values, data['BETA_OUT'].values, data['SE_EXP'].values, data['SE_OUT'].values, 1/VBj, nboot)
    model_pvalue = 2 * (1 - norm.cdf(abs(beta/standard_error)))
    return [model_pvalue, beta, standard_error]

def run_mr_penalized_weighted_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000, penalty_con = 20):
    '''
    Calculate the penalized weighted median MR estimate.

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
    nboot : int
        Number of bootstrap replications to calculate the standard error.
    penalty_con : int
        Constant term in penalisation.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''

    print("Running Penalized Weighted Median")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Weighted Median.'.format(data[['BETA_EXP']].shape[0]))
        return None

    beta_iv = np.array(data['BETA_OUT'] / data['BETA_EXP'])
    beta_ivw = np.sum(data['BETA_OUT'] * data['BETA_EXP'] / data['SE_OUT']**2) / np.sum(data['BETA_EXP']**2 / data['SE_OUT']**2)
    VBj = np.array(((data['SE_OUT'])**2) / (data['BETA_EXP'])**2 + (data['BETA_OUT']**2) * ((data['SE_EXP']**2)) / (data['BETA_EXP'])**4)

    beta_weighted_median = run_mr_weighted_median(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot)
    penalty = chi2.sf((1/VBj) * (beta_iv - beta_weighted_median[1])**2, df=1)
    penalized_weights = (1/VBj) * np.minimum(1, penalty * penalty_con)

    beta = weighted_median(beta_iv, penalized_weights)
    standard_error = weighted_median_bootstrap(data['BETA_EXP'].values, data['BETA_OUT'].values, data['SE_EXP'].values, data['SE_OUT'].values, penalized_weights, nboot)
    model_pvalue = 2 * (1 - norm.cdf(abs(beta/standard_error)))
    return [model_pvalue, beta, standard_error]

def weighted_median(b_iv, weights):
    '''
    Calculate the weighted median of a list of numbers.
    credit: https://github.com/MRCIEU/TwoSampleMR/blob/master/R/mr.R
    '''
    # Sort b_iv and weights based on b_iv
    sorted_indices = np.argsort(b_iv)
    betaIV_order = np.array(b_iv)[sorted_indices]
    weights_order = np.array(weights)[sorted_indices]

    # Calculate cumulative sum of weights
    weights_sum = np.cumsum(weights_order) - 0.5 * weights_order
    weights_sum /= np.sum(weights_order)

    # Find index where cumulative sum is below 0.5
    below = np.max(np.where(weights_sum < 0.5))

    # Calculate weighted median
    b = betaIV_order[below] + (betaIV_order[below + 1] - betaIV_order[below]) * \
        (0.5 - weights_sum[below]) / (weights_sum[below + 1] - weights_sum[below])

    return b
    
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

def binDist(x, w, xlo, xhi, n):
    '''
    Distributes weights across bins based on the positions of the x values within the 
    specified range (from xlo to xhi) and the number of bins (n).

    credit: massdist.c from R source code
    '''
    # Initialize the output array
    ans = np.zeros(2 * n)
    
    # Calculate the width of each bin
    xdelta = (xhi - xlo) / (n - 1)
    
    # Iterate through each x value
    for i in range(len(x)):
        if np.isfinite(x[i]):
            # Calculate the position of x within the bins
            xpos = (x[i] - xlo) / xdelta
            
            # Determine the closest bin index and the fractional part
            ix = int(np.floor(xpos))
            fx = xpos - ix
            
            # Distribute the weight across the two closest bins
            if ix >= 0 and ix < (n - 1):
                ans[ix] += (1 - fx) * w[i]
                ans[ix + 1] += fx * w[i]
            elif ix == -1:
                ans[0] += fx * w[i]
            elif ix == (n - 1):
                ans[n - 1] += (1 - fx) * w[i]
    
    return ans

def density(x, bw, weights):
    '''
    Computes kernel density estimates with gaussian kernel. Trimmed version of R stats density (S3)
    credit: density.R from R source code
    '''
    if not np.issubdtype(np.array(x).dtype, np.number):
        raise ValueError("argument 'x' must be numeric")
    if len(weights) != len(x):
        raise ValueError("'x' and 'weights' have unequal length")
    if np.any(np.isnan(x)):
        raise ValueError("'x' contains missing values")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("'weights' must all be finite and not negative")
    
    wsum = np.sum(weights)
    totMass = 1 if np.allclose(wsum, 1) else np.sum(weights[~np.isnan(x) & np.isfinite(x)]) / wsum
    
    x = x[np.isfinite(x)]
    weights = weights[np.isfinite(x)]
    
    if not np.isfinite(bw) or bw <= 0:
        raise ValueError("'bw' is not positive or not finite")
    
    n = 512
    cut = 3
    from_ = np.min(x) - cut * bw
    to = np.max(x) + cut * bw
    lo = from_ - 4 * bw
    up = to + 4*bw
    
    y = binDist(x, weights, lo, up, n) * totMass

    kords = np.linspace(0, 2*(up -lo), num = 2 * n)
    kords[n:] = -np.flip(kords[1:n+1])
    kords = norm.pdf(kords, scale=bw)
    kords = ifft(fft(y) * np.conj(fft(kords))).real
    kords = np.maximum(0, kords[:n] / len(y))
    
    xords = np.linspace(from_ - 4 * bw, to + 4 * bw, n)
    x_interp = np.linspace(from_, to, n)
    y_interp = interp1d(xords, kords, fill_value="extrapolate")(x_interp)
    
    return {"x": x_interp, "y": y_interp}

def beta_MODE(BetaIV_in, seBetaIV_in, phi):
    
    # Bandwidth rule - modified Silverman's rule proposed by Bickel (2002)
    n = len(BetaIV_in)
    s = 0.9 * min(np.std(BetaIV_in, ddof=1), median_abs_deviation(BetaIV_in, scale = 'normal')) / (n**(1/5))
    
    # Standardized weights
    weights = seBetaIV_in**-2 / sum(seBetaIV_in**-2)
    
    # Define the actual bandwidth
    h = max(0.00000001, s * phi)

    density_result = density(BetaIV_in, h, weights)

    # Find the x-value where the y-value (density) is at its maximum
    max_density_x_value = density_result['x'][density_result['y'] == np.max(density_result['y'])][0]

    return max_density_x_value

def run_mr_simple_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000, penalty_con = 20, phi = 1, alpha = 0.05):
    '''
    Calculate the Simple Mode MR estimate.

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
    nboot : int
        Number of bootstrap replications to calculate the standard error.
    penalty_con : int
        Constant term in penalisation.
    phi : int
        Bandwidth parameter.
    alpha : float
        Width of confidence interval.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''   

    print("Running Simple Mode")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Simple Mode.'.format(data[['BETA_EXP']].shape[0]))
        return None

    BetaIV = data['BETA_OUT'] / data['BETA_EXP']
    seBetaIV = np.vstack([np.sqrt((data['SE_OUT']**2) / (data['BETA_EXP']**2) + ((data['BETA_OUT']**2) * (data['SE_EXP']**2)) / (data['BETA_EXP']**4)),
                            data['SE_OUT'] / np.abs(data['BETA_EXP'])]).T
    beta = beta_MODE(BetaIV, np.ones_like(BetaIV), phi)

    # calculate boostrapped betas 
    beta_boot = np.zeros(nboot)

    for i in range(nboot):
        # Resample each ratio estimate using SEs derived not assuming NOME
        BetaIV_boot = norm.rvs(size=len(BetaIV), loc=BetaIV, scale=seBetaIV[:, 0])
        
        # Simple mode, not assuming NOME
        beta_boot[i] = beta_MODE(BetaIV_in=BetaIV_boot, seBetaIV_in=np.ones(len(BetaIV)), phi=phi)
    
    standard_error = median_abs_deviation(beta_boot, scale='normal')

    CIlow_Mode = beta - norm.ppf(1 - alpha/2) * standard_error
    CIupp_Mode = beta + norm.ppf(1 - alpha/2) * standard_error

    model_pvalue = 2 * t.sf(abs(beta / standard_error), len(data['BETA_EXP']) - 1)

    return [model_pvalue, beta, standard_error]

def run_mr_simple_mode_nome(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000, penalty_con = 20, phi = 1, alpha = 0.05):
    '''
    Calculate the Simple Mode (NOME) MR estimate.

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
    nboot : int
        Number of bootstrap replications to calculate the standard error.
    penalty_con : int
        Constant term in penalisation.
    phi : int
        Bandwidth parameter.
    alpha : float
        Width of confidence interval.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''

    print("Running Simple Mode (NOME)")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Simple Mode (NOME).'.format(data[['BETA_EXP']].shape[0]))
        return None

    BetaIV = data['BETA_OUT'] / data['BETA_EXP']
    seBetaIV = np.vstack([np.sqrt((data['SE_OUT']**2) / (data['BETA_EXP']**2) + ((data['BETA_OUT']**2) * (data['SE_EXP']**2)) / (data['BETA_EXP']**4)),
                        data['SE_OUT'] / np.abs(data['BETA_EXP'])]).T
    beta = beta_MODE(BetaIV, np.ones_like(BetaIV), phi)

    ########################

    # calculate boostrapped betas 
    beta_boot = np.zeros(nboot)

    for i in range(nboot):
        # Resample each ratio estimate using SEs derived under NOME
        BetaIV_boot_NOME = norm.rvs(size=len(BetaIV), loc=BetaIV, scale=seBetaIV[:, 1])
        
        # Simple mode, assuming NOME
        beta_boot[i] = beta_MODE(BetaIV_in=BetaIV_boot_NOME, seBetaIV_in=np.ones(len(BetaIV)), phi=phi)
    
    standard_error = median_abs_deviation(beta_boot, scale='normal')

    CIlow_Mode = beta - norm.ppf(1 - alpha/2) * standard_error
    CIupp_Mode = beta + norm.ppf(1 - alpha/2) * standard_error

    model_pvalue = 2 * t.sf(abs(beta / standard_error), len(data['BETA_EXP']) - 1)

    return [model_pvalue, beta, standard_error]

def run_mr_weighted_mode_nome(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000, penalty_con = 20, phi = 1, alpha = 0.05):
    '''
    Calculate the Weighted Mode (NOME) MR estimate.

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
    nboot : int
        Number of bootstrap replications to calculate the standard error.
    penalty_con : int
        Constant term in penalisation.
    phi : int
        Bandwidth parameter.
    alpha : float
        Width of confidence interval.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''

    print("Running Weighted Mode (NOME)")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Weighted Mode (NOME).'.format(data[['BETA_EXP']].shape[0]))
        return None

    BetaIV = data['BETA_OUT'] / data['BETA_EXP']
    seBetaIV = np.vstack([np.sqrt((data['SE_OUT']**2) / (data['BETA_EXP']**2) + ((data['BETA_OUT']**2) * (data['SE_EXP']**2)) / (data['BETA_EXP']**4)),
                            data['SE_OUT'] / np.abs(data['BETA_EXP'])]).T

    beta = beta_MODE(BetaIV, seBetaIV[:, 1], phi)

    ########################

    # calculate boostrapped betas 
    beta_boot = np.zeros(nboot)

    for i in range(nboot):
        # Resample each ratio estimate using SEs derived under NOME
        BetaIV_boot_NOME = norm.rvs(size=len(BetaIV), loc=BetaIV, scale=seBetaIV[:, 1])

        # Weighted mode, assuming NOME
        beta_boot[i] = beta_MODE(BetaIV_in=BetaIV_boot_NOME, seBetaIV_in=seBetaIV[:, 1], phi=phi)
    
    standard_error = median_abs_deviation(beta_boot, scale='normal')

    CIlow_Mode = beta - norm.ppf(1 - alpha/2) * standard_error
    CIupp_Mode = beta + norm.ppf(1 - alpha/2) * standard_error

    model_pvalue = 2 * t.sf(abs(beta / standard_error), len(data['BETA_EXP']) - 1)

    return [model_pvalue, beta, standard_error]

def run_mr_penalized_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000, penalty_con = 20, phi = 1, alpha = 0.05):
    '''
    Calculate the Penalized Mode MR estimate.

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
    nboot : int
        Number of bootstrap replications to calculate the standard error.
    penalty_con : int
        Constant term in penalisation.
    phi : int
        Bandwidth parameter.
    alpha : float
        Width of confidence interval.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''

    print("Running Penalized Mode")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Penalized Weighted Mode.'.format(data[['BETA_EXP']].shape[0]))
        return None

    BetaIV = data['BETA_OUT'] / data['BETA_EXP']
    seBetaIV = np.vstack([np.sqrt((data['SE_OUT']**2) / (data['BETA_EXP']**2) + ((data['BETA_OUT']**2) * (data['SE_EXP']**2)) / (data['BETA_EXP']**4)),
                            data['SE_OUT'] / np.abs(data['BETA_EXP'])]).T

    beta_WeightedMode = beta_MODE(BetaIV, seBetaIV[:, 0], phi)

    weights = 1 / seBetaIV[:, 0]**2
    penalty = chi2.sf(weights * (BetaIV - beta_WeightedMode)**2, df=1)
    pen_weights = weights * np.minimum(1, penalty * penalty_con)
    beta = beta_MODE(BetaIV, np.sqrt(1/pen_weights), phi)

    # calculate boostrapped betas 
    beta_boot = np.zeros(nboot)

    for i in range(nboot):
        # Resample each ratio estimate using SEs derived not assuming NOME
        BetaIV_boot = norm.rvs(size=len(BetaIV), loc=BetaIV, scale=seBetaIV[:, 0])
        
        # Weighted mode, not assuming NOME
        beta_boot[i] = beta_MODE(BetaIV_in=BetaIV_boot, seBetaIV_in=seBetaIV[:, 0], phi=phi)

        # Penalized mode, not assuming NOME
        weights = 1 / seBetaIV[:, 0]**2
        penalty = chi2.sf(weights * (BetaIV_boot - beta_boot[i])**2, df=1)
        pen_weights = weights * np.minimum(1, penalty * penalty_con)
        beta_boot[i] = beta_MODE(BetaIV_in=BetaIV_boot, seBetaIV_in=np.sqrt(1/pen_weights), phi=phi)
    
    standard_error = median_abs_deviation(beta_boot, scale='normal')

    CIlow_Mode = beta - norm.ppf(1 - alpha/2) * standard_error
    CIupp_Mode = beta + norm.ppf(1 - alpha/2) * standard_error

    model_pvalue = 2 * t.sf(abs(beta / standard_error), len(data['BETA_EXP']) - 1)

    return [model_pvalue, beta, standard_error]

def run_mr_weighted_mode(harmonized_data, beta_exp, beta_out, se_exp, se_out, nboot = 1000, penalty_con = 20, phi = 1, alpha = 0.05):
    '''
    Calculate the Weighted Mode MR estimate.

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
    nboot : int
        Number of bootstrap replications to calculate the standard error.
    penalty_con : int
        Constant term in penalisation.
    phi : int
        Bandwidth parameter.
    alpha : float
        Width of confidence interval.

    Returns
    -------
    list
        List containing the model-pvalue, beta, and standard error.
    '''

    print("Running Weighted Mode")

    data = harmonized_data[[beta_exp, beta_out, se_exp, se_out]]
    
    data = data.rename(columns={beta_exp : 'BETA_EXP', 
                                beta_out : 'BETA_OUT', 
                                se_exp : 'SE_EXP', 
                                se_out : 'SE_OUT'})
    
    exp_beta_len = data['BETA_EXP'].shape[0]
    
    if exp_beta_len < 3:
        print('Found {} variants. Need at least 3. Exiting MR Weighted Mode.'.format(data[['BETA_EXP']].shape[0]))
        return None

    BetaIV = data['BETA_OUT'] / data['BETA_EXP']
    seBetaIV = np.vstack([np.sqrt((data['SE_OUT']**2) / (data['BETA_EXP']**2) + ((data['BETA_OUT']**2) * (data['SE_EXP']**2)) / (data['BETA_EXP']**4)),
                            data['SE_OUT'] / np.abs(data['BETA_EXP'])]).T

    beta = beta_MODE(BetaIV, seBetaIV[:, 0], phi)

    ############

    # calculate boostrapped betas 
    beta_boot = np.zeros(nboot)

    for i in range(nboot):
        # Resample each ratio estimate using SEs derived not assuming NOME
        BetaIV_boot = norm.rvs(size=len(BetaIV), loc=BetaIV, scale=seBetaIV[:, 0])
        
        # Weighted mode, not assuming NOME
        beta_boot[i] = beta_MODE(BetaIV_in=BetaIV_boot, seBetaIV_in=seBetaIV[:, 0], phi=phi)
    
    standard_error = median_abs_deviation(beta_boot, scale='normal')

    CIlow_Mode = beta - norm.ppf(1 - alpha/2) * standard_error
    CIupp_Mode = beta + norm.ppf(1 - alpha/2) * standard_error

    model_pvalue = 2 * t.sf(abs(beta / standard_error), len(data['BETA_EXP']) - 1)

    return [model_pvalue, beta, standard_error]

def matrix_power_eigen(x, n):
    values, vectors = np.linalg.eig(x)
    M_power_n = vectors @ np.diag(values**n) @ np.linalg.inv(vectors)
    return M_power_n

def getRSS_LOO(BetaOutcome, BetaExposure, data, returnIV):
    # Check if BetaExposure is a list (multiple exposures) or a single string (single exposure)
    
    is_multiple_exposures = len(BetaExposure) > 1
    
    # Adjust dataW construction based on the type of BetaExposure
    if is_multiple_exposures:
        dataW = data[[BetaOutcome] + BetaExposure].multiply(np.sqrt(data["Weights"].values[:, np.newaxis]), axis=0)
        X = dataW[BetaExposure].values
    else:
        BetaExposure = BetaExposure[0]
        dataW = data[[BetaOutcome, BetaExposure]].multiply(np.sqrt(data["Weights"].values[:, np.newaxis]), axis=0)
        X = dataW[[BetaExposure]].values
    
    Y = dataW[[BetaOutcome]].values

    CausalEstimate_LOO = []
    for i in range(X.shape[0]):
        X_loo = np.delete(X, i, axis=0)
        Y_loo = np.delete(Y, i, axis=0)
        # Ensure X_loo is 2D for matrix operations
        if X_loo.ndim == 1:
            X_loo = X_loo.reshape(-1, 1)
        if Y_loo.ndim == 1:
            Y_loo = Y_loo.reshape(-1, 1)
        inv_XTX = matrix_power_eigen(X_loo.T @ X_loo, -1)
        estimate = (inv_XTX @ X_loo.T) @ Y_loo

        CausalEstimate_LOO.append(estimate)

    CausalEstimate_LOO = np.squeeze(np.array(CausalEstimate_LOO))

    # Adjust the RSS calculation based on the dimensions of CausalEstimate_LOO
    if is_multiple_exposures:
        # Proper handling for multiple exposures
        RSS = np.sum((Y - np.dot(X, CausalEstimate_LOO.T)) ** 2)
    else:
        # For single exposure
        RSS = np.sum((Y.flatten() - CausalEstimate_LOO * X.flatten()) ** 2)

    if returnIV:
        return RSS, CausalEstimate_LOO
    else:
        return RSS

def getRandomData(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data):
    # Prepare the formula for linear regression
    formula = f"{BetaOutcome} ~ -1 + {' + '.join(BetaExposure)}"
    
    # Run leave-one-out regression and predict the outcome for the left-out observation
    predictions = []
    for i in range(len(data)):
        # Exclude the i-th row
        train_data = data.drop(index=i)
        test_data = data.iloc[[i]]
        
        # Fit the linear model
        model = wls(formula=formula, data=train_data, weights=train_data['Weights']).fit()
        
        # Predict the outcome for the left-out observation
        pred = model.predict(test_data)
        predictions.append(pred.iloc[0])
    
    # Generate random data based on the original data and predictions
    dataRandom = pd.DataFrame({
        exp: np.random.normal(loc=data[exp], scale=data[SdExp]) for exp, SdExp in zip(BetaExposure, SdExposure)
    })
    dataRandom[BetaOutcome] = [np.random.normal(loc=pred, scale=data.iloc[i][SdOutcome]) for i, pred in enumerate(predictions)]
    dataRandom['Weights'] = data['Weights']
    
    return dataRandom

def getRandomBias(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, refOutlier):
    indices = list(refOutlier) + [np.random.choice(list(set(range(len(data))) - set(refOutlier))) for _ in range(len(data) - len(refOutlier))]
    formula = f"{BetaOutcome} ~ -1 + {' + '.join(BetaExposure)}"
    data_subset = data.iloc[indices[:len(indices) - len(refOutlier)]]
    mod_random = wls(formula, data=data_subset, weights=data_subset['Weights']).fit()
    return mod_random.params[BetaExposure]

def run_mr_presso(harmonized_data, beta_exp, beta_out, se_exp, se_out, outlier_test = False, distortion_test = False, sigthresh = 0.05, nbDist = 1000):
    '''
    Calculate the MR PRESSO test.

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
    dict
        Dictionary containing the MR results, as well as the results of the Global, Outlier, and Distortion Tests.
        For more information check out https://github.com/rondolab/MR-PRESSO/tree/master
    '''

    print("Running Pleiotropy RESidual Sum and Outlier (PRESSO)")

    if isinstance(beta_exp, str):
        beta_exp = [beta_exp]
    if isinstance(se_exp, str):
        se_exp = [se_exp]

    selected_columns = [beta_out] + beta_exp + [se_out] + se_exp

    data = harmonized_data[selected_columns]

    # drop rows with any NA values
    data = data.dropna()

    # Adjust specified columns by the sign of the first exposure variable
    sign_of_first_exposure = np.sign(data[beta_exp[0]])

    # Loop through each column to be adjusted and multiply by the sign
    for column in [beta_out] + beta_exp:
        data[column] *= sign_of_first_exposure

    # Adjust weights based on the square of the inverse of se_out
    data['Weights'] = 1 / data[se_out] ** 2

    if data.shape[0] <= len(beta_exp) + 2:
        print('Found {} variants. Need at least {}. Exiting MR PRESSO.'.format(data.shape[0], len(beta_exp) + 3))
        return None

    if len(data) >= nbDist:
        print("Not enough elements to compute empirical P-values for MR PRESSO, increase nbDist from {} to above {}".format(nbDist, len(data)))
        return None

    RSSobs = getRSS_LOO(beta_out, beta_exp, data, outlier_test)

    random_data = [getRandomData(BetaOutcome=beta_out, BetaExposure=beta_exp,
                                SdOutcome=se_out, SdExposure=se_exp, data=data)
                                for _ in range(nbDist)]

    # Computing RSS for each random dataset
    RSSexp = [getRSS_LOO(beta_out, beta_exp, random_dataset, outlier_test)
                        for random_dataset in random_data]

    # Calculate p-value based on outlier_test
    if outlier_test:
        # Assuming RSSobs is a list with the observed RSS as its first element
        GlobalTest = {'RSSobs': RSSobs[0], 'P_value': sum(rss[0] > RSSobs[0] for rss in RSSexp) / nbDist}
    else:
        GlobalTest = {'RSSobs': RSSobs, 'P_value': sum(rss > RSSobs for rss in RSSexp) / nbDist}

    beta_outcome_index = random_data[0].columns.get_loc(beta_out)
    beta_exposure_indices = [random_data[0].columns.get_loc(col) for col in beta_exp]

    if GlobalTest['P_value'] < sigthresh and outlier_test:
        OutlierTest = []

        for SNV in range(len(data)):
            # Extracting the row across all random_data matrices for the current SNV
            randomSNP = np.vstack([mat.iloc[SNV] for mat in random_data])

            # Calculating the difference and expected values based on BetaExposure length
            if len(beta_exp) == 1:
                Dif = (data.iloc[SNV][beta_out] - data.iloc[SNV][beta_exp] * RSSobs[1][SNV])[0]
                Exp = randomSNP[:, beta_outcome_index] - randomSNP[:, beta_exposure_indices[0]] * RSSobs[1][SNV]
            else:
                Dif = (data.iloc[SNV][beta_out] - np.sum(data.iloc[SNV][beta_exp] * RSSobs[1][:, SNV])).values
                Exp = randomSNP[:, beta_outcome_index] - np.sum(randomSNP[:, beta_exposure_indices] * RSSobs[1][:, SNV], axis=1)

            # Calculating p-value for the current SNV
            pval = np.sum(Exp**2 > Dif**2) / len(random_data)
            OutlierTest.append([Dif**2, pval])

        # Converting the results to a DataFrame
        OutlierTest = pd.DataFrame(OutlierTest, columns=['RSSobs', 'P_value'], index=data.index)

        # Applying Bonferroni correction
        OutlierTest['P_value'] = np.minimum(OutlierTest['P_value'] * len(data), 1.0)
    else:
        outlier_test = False     

    formula_all = f"{beta_out} ~ -1 + {' + '.join(beta_exp)}"
    mod_all = wls(formula_all, data=data, weights=data['Weights']).fit()
    
    if distortion_test and outlier_test:
        refOutlier = OutlierTest[OutlierTest['P_value'] <= sigthresh].index.tolist()
    
        if len(refOutlier) > 0:
            if len(refOutlier) < len(data):
                BiasExp = [getRandomBias(beta_out, beta_exp, se_out, se_exp, data, refOutlier) for _ in range(nbDist)]
                BiasExp = pd.DataFrame(BiasExp)
                
                data_noOutliers = data.drop(refOutlier)
                mod_noOutliers = wls(formula_all, data=data_noOutliers, weights=data_noOutliers['Weights']).fit()
                
                BiasObs = (mod_all.params[beta_exp] - mod_noOutliers.params[beta_exp]) / abs(mod_noOutliers.params[beta_exp])
                BiasExp = (mod_all.params[beta_exp].values - BiasExp) / abs(BiasExp)
                
                BiasTest = {'Outliers Indices': refOutlier, 
                            'Distortion Coefficient': (100 * BiasObs)[0], 
                            'P_value': np.sum(np.abs(BiasExp) > np.abs(BiasObs)) / nbDist}
            else:
                BiasTest = {'Outliers Indices': "All SNPs considered as outliers", 
                            'Distortion Coefficient': np.nan, 
                            'P_value': np.nan}
        else:
            BiasTest = {'Outliers Indices': "No significant outliers", 
                        'Distortion Coefficient': np.nan, 
                        'P_value': np.nan}
    else:
        outlier_test = False
    
    if GlobalTest['P_value'] == 0:
        GlobalTest['P_value_significance'] = f"<{1/nbDist}"
    else:
        GlobalTest['P_value_significance'] = f"not <{1/nbDist}"
    
    res = {'Global Test': GlobalTest}

    if outlier_test:
        # Formatting the p-value in OutlierTest
        OutlierTest['P_value'] = OutlierTest['P_value'].apply(lambda p: f"<{len(data)/nbDist}" if p == 0 else p)
        res['Outlier Test'] = OutlierTest
        
        if distortion_test:
            # Formatting the p-value in BiasTest
            BiasTest['P_value'] = BiasTest['P_value'].apply(lambda p: f"<{1/nbDist}" if p == 0 else p)
            res['Distortion Test'] = BiasTest
        
        # Check for instability in outlier test
        if len(data) / nbDist > sigthresh:
            print("Outlier test unstable for MR PRESSO. The significance threshold of {} for the outlier test is not achievable with only {} distributions to compute the null distribution. The current precision is <{}. Increase nbDist.".format(sigthresh, nbDist, len(data)/nbDist))

    coefficients_all = mod_all.summary2().tables[1]
    OriginalMR = ["MR PRESSO Raw", coefficients_all['P>|t|'], coefficients_all['Coef.'], coefficients_all['Std.Err.'], coefficients_all['t']]

    # Checking if mod_noOutliers exists and constructing the OutlierCorrectedMR DataFrame
    try:
        coefficients_no_outliers = mod_noOutliers.summary2().tables[1]
        OutlierCorrectedMR = ["MR PRESSO Outlier-corrected", coefficients_no_outliers['P>|t|'], coefficients_no_outliers['Coef.'], coefficients_no_outliers['Std.Err.'],coefficients_no_outliers['t']]

        res['MR_RESULTS'] = [OriginalMR, OutlierCorrectedMR]

    except NameError:
        if outlier_test == True and distortion_test == True:
            print("No outliers were identified for MR PRESSO, therefore no results will be returned for the outlier-corrected MR PRESSO.")
        res['MR_RESULTS'] = [OriginalMR]

    # Return the results dictionary
    return res