from scanometrics.core import ScanOMetrics_project
from os import chdir

chdir('/home/david/Documents/DATA')
SOM = ScanOMetrics_project("test_precomputed_dataset", cov2float={'sex': {'M': 0, 'm':0, 'F':1, 'f':1}, 'scanner':{'Skyra_fit':0 , 'MAGNETOM_Vida':1}, 'sequence': {'t1_mpr_adni_sag_iso_1mm_ns':0, '1_mp2rage_sag_p3_iso':1, 'MDEFT':2, 't1_mpr_sag_iso_1mm_ns':3}})
SOM.load_subjects()
SOM.load_proc_metrics()
SOM.flag_outliers(1.5, 'measured_metrics')
SOM.set_normative_model('Polynomial')
# Sets uncertainty to one as computation of uncertainty is done afterwards.
import numpy as np

dummy_uncertainty = {}
for k in SOM.measured_metrics.keys():
    dummy_uncertainty[k] = np.ones(SOM.measured_metrics[k].shape[1])
SOM.normativeModel.fit(SOM.measured_metrics, SOM.metric_names, dummy_uncertainty, SOM.outliers, SOM.covariate_values,
                       SOM.covariate_names, flag_opt=1, N_cycl=1)
# Re-esimate outliers based on residuals
SOM.flag_outliers(1.5, 'residuals')
# Estimate uncertainties
SOM.compute_uncertainty()
# Fit again with clean data
SOM.normativeModel.fit(SOM.measured_metrics, SOM.metric_names, SOM.uncertainty,
                       SOM.outliers, SOM.covariate_values, SOM.covariate_names,
                       flag_opt=1,
                       deg_max=int(np.floor(len(SOM.covariate_values[:, SOM.covariate_names.index('age')]) / (2 * 20))),
                       N_cycl=100)  # reduce deg_max by a factor two, N_cycl=100
# Save fitted normative model with default name in default folder
SOM.save_normative_model()
