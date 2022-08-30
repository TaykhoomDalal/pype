import pickle
import numpy as np

def pickle_object(obj, filename, protocol=pickle.HIGHEST_PROTOCOL) -> None:
	"""
	Pickle an object to a file.
	
	Parameters
	----------
	obj : object
		Object to pickle.
	filename : str	
		Filename of the pickle file.
	protocol : int
		Protocol to use for pickling.

	Returns	
	-------
	None
	"""

	with open(filename, 'wb') as f:
		pickle.dump(obj, f, protocol = protocol)
	return

def unpickle_object(filename) -> object:
	"""
	Unpickle an object from a file.

	Parameters
	----------
	filename : str
		Filename of the pickle file.
	
	Returns
	-------
	obj : object
		Object that was pickled.
	"""
	with open(filename, 'rb') as f:
		obj = pickle.load(f)
	return obj

def multiple_testing_correction(pvalues, alpha, method) -> float:
	"""
	Perform multiple testing correction on pvalues, returning threshold to use.

	Parameters
	----------
	pvalues : array_like
		Array of pvalues.
	alpha : float
		Alpha value for the test.
	method : str
		Method for multiple testing correction.
		'bonferroni' for Bonferroni correction.
		'sidak' for Sidak correction.
		'fdr_bh' for Benjamini-Hochberg correction.
		'no_correction' for no correction.

	Returns
	-------
	threshold : float
		Threshold value for the test.
	"""

	pvalues = np.array(pvalues)
	n_pvalues = len(pvalues)

	if method == 'bonferroni':
		threshold = alpha / n_pvalues
	elif method == 'sidak':
		threshold = 1 - np.power(1. - alpha, 1. / n_pvalues)
	elif method == 'fdr_bh':
		# sort the pvalues in ascending order
		pvalues.sort()

		for rank, pval in enumerate(pvalues, 1):

			# benjamini-holchberg critical value
			critical_val = alpha * rank / n_pvalues

			# the largest pvaue that is less than the critical value
			if pval > critical_val:
				threshold = pval # the significant pvalues should be strictly less than this
				break
	elif method == 'no_correction':
		threshold = alpha
	else:
		raise ValueError('Invalid method for multiple testing correction.')
	
	print('Threshold for {} multiple testing correction: {}'.format(method.capitalize(), threshold))

	return threshold
