"""
Tests for python_module
"""


from scanometrics.utils.stats import fdr
import numpy as np

def test_fdr():
    """Tests scanometrics.utils.stats.fdr() function, compared to octave output"""
    p = np.array([0.01, np.nan, 0.003, 0.2, 0.3, 0.004, 0.2, 0.1, 0.02, np.nan, 0.04, 0.1, 0.003, 0.0004, 0.07])
    pID, pN = fdr(p, 0.05)
    assert pID == 0.02 and pN == 4e-3  # Check that output matches the octave implementation
    assert np.array_equal(p,np.array([0.01, np.nan, 0.003, 0.2, 0.3, 0.004, 0.2, 0.1, 0.02, np.nan, 0.04, 0.1, 0.003, 0.0004, 0.07]),equal_nan=True)  # Check that p was not modified
