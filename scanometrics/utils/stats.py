"""
stats.py : stats module for ScanOMetrics software
"""

import numpy as np


def fdr(p, q=0.05):
    """Compute False Discovery Rate (FDR)'s parametric and non-parametric thresholds. Returns 0.0 if no sample is found
     below the FDR.

    :param p: vector of input p-values
    :type p: numpy (N,) array
    :param q: False Discovery Rate level, defaults to 0.05
    :type q: float, optional
    :return: FDR parametric (based on independence or positive dependence) and non-parametric thresholds
    :rtype: tuple
    """

    p = p[np.isfinite(p)]  # ignore NaN's
    p = np.sort(p.reshape(-1))
    V = len(p)
    I = np.arange(1, V+1)  # rank of ith p-value

    cVID = 1
    cVN = np.sum(1. / np.arange(1, V+1))

    cond_idx = np.argwhere(p <= I / V * q / cVID)[:, 0]
    if len(cond_idx) >= 1:
        pID = p[cond_idx.max()]
    else:
        pID = 0.0
    cond_idx = np.argwhere(p <= I / V * q / cVN)[:, 0]
    if len(cond_idx) >= 1:
        pN = p[cond_idx.max()]
    else:
        pN = 0.0
    return pID, pN
