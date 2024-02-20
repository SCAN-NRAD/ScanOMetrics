"""
Test possibility of formulating polynomial fit as a least-square fit, accounting for NaNs
"""

import numpy as np

# Mimics 500 subjects aged between 1 and 100 years old, for which 4 metrics are taken
np.random.seed(238237)
age = np.random.uniform(1, 100, (500, 1))
A = np.hstack((np.ones(age.shape), age, age**2, age**3))
b = np.hstack((2 -3*age + 2.5*age**3,
               age**2-0.5*age**3,
               1.2*age,
               -3.5*age**2+2*age**3))
np.random.seed(2323419)
b += np.random.normal(0, 0.05, b.shape)

x = np.linalg.pinv(A).dot(b).T  # gets close GT estimate

# Now sets 10 subjects to nan for 1st metric
outliers = np.zeros(len(age), dtype='bool')
outliers[np.array([5,10,53,86,176,246,321,378,432,478])] = 1
b[outliers, 0] = np.nan
np.linalg.pinv(A).dot(b).T

# The estimate we should get without these values
np.linalg.pinv(A[~outliers, :]).dot(b[~outliers, 0]).T

# TODO: test if fit changes with higher order A, or find out how to estimate An+1 from A
