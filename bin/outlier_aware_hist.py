import numpy as np
from matplotlib import pyplot as plt

def mad(data):
    median = np.median(data)
    diff = np.abs(data - median)
    mad = np.median(diff)
    return mad

def calculate_bounds(data, z_thresh=3.5):
    MAD = mad(data)
    median = np.median(data)
    const = z_thresh * MAD / 0.6745
    return (median - const, median + const)

def outlier_aware_hist(data, lower=None, upper=None, cdf=False):
    if not lower or lower < data.min():
        lower = data.min()
        lower_outliers = False
    else:
        lower_outliers = True
        
    if not upper or upper > data.max():
        upper = data.max()
        upper_outliers = False
    else:
        upper_outliers = True
    
    if cdf:
        n, bins, patches = plt.hist(np.clip(data, lower, upper), bins='auto', cumulative=True, density=True, histtype='step')
    else:
        n, bins, patches = plt.hist(data, range=(lower, upper), bins='auto')
        
        if lower_outliers:
            n_lower_outliers = (data < lower).sum()
            patches[0].set_height(patches[0].get_height() + n_lower_outliers)
            patches[0].set_facecolor('c')
            patches[0].set_label('Lower outliers: ({:.2f}, {:.2f})'.format(data.min(), lower))

        if upper_outliers:
            n_upper_outliers = (data > upper).sum()
            patches[-1].set_height(patches[-1].get_height() + n_upper_outliers)
            patches[-1].set_facecolor('m')
            patches[-1].set_label('Upper outliers: ({:.2f}, {:.2f})'.format(upper, data.max()))

        if lower_outliers or upper_outliers:
            plt.legend()
