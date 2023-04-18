import numpy as np
from scipy import stats
from scipy.spatial import distance


def __fix_zero_bins(data):
    """add a small value to bins of count 0 for divergense calculations
    
    Parameters
    ----------
    data : np.ndarray
        the counts or density of the bins
    
    Returns
    -------
    data : np.ndarray
        the counts or density of the bins, with very small values added to the ones that had value of 0
    """
    zero_bins = []
    non_zero_values = []
    for i, value in enumerate(data):
        if value==0:
            zero_bins.append(i)
        else:
            non_zero_values.append(value)
    zero_value = (min(non_zero_values)/1000)/len(zero_bins)
    for idx in zero_bins:
        data[idx] = zero_value
    
    return data


def calculate_divergence(data, reference_data, bins, div_type='JS', data_range=None):
    """Calculate divergence between two distributions. When calculating JS divergence, it does not matter which
    distribution is the reference one, as it is symmetrical. However, this becomes important when calculating
    KL divergence, so if 'data' and 'reference_data' is swapped around, the values will be different.
    
    Parameters
    ----------
    data : np.ndarray
        data for the distribution
    
    reference_data : np.ndarray
        data for the reference distribution
        
    bins : int
        number of bins
        
    div_type : ['JS', 'KL']
        divergence type
    
    data_range : (float, float) or (int, int)
        range of binned data. If none, the min/max of all data will be used
    
    Returns
    -------
    divergence : float
        divergence between two distributions
    """
    #calculate range if needed
    if data_range is None:
        combined_data = np.append(data, reference_data)
        data_range = (combined_data.min(), combined_data.max())
    
    #bin data
    data_bins = np.histogram(data, bins=bins, range=data_range, density=True)[0]
    data_bins = __fix_zero_bins(data_bins)
    reference_bins = np.histogram(reference_data, bins=bins, range=data_range, density=True)[0]
    reference_bins = __fix_zero_bins(reference_bins)

    #calculate divergence
    if div_type == 'KL':
        divergence = stats.entropy(pk=data_bins, qk=reference_bins)
    elif div_type == 'JS':
        divergence = distance.jensenshannon(data_bins, reference_bins)
    
    return divergence