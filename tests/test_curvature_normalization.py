"""
Tests if sum of gaussian curvatures approximates 2*pi as the gaussian curvature of a closed surface should be.
"""

from scanometrics.core import ScanOMetrics_project
import numpy as np
import matplotlib.pyplot as plt

bids_database = '/home/david/Documents/DATA/test_precomputed_dataset'
proc_pipeline = 'freesurfer'
normative_model = 'Polynomial'


cov2float = cov2float={'sex': {'M': 0, 'm':0, 'F':1, 'f':1}, 'scanner':{'Skyra_fit':0 , 'MAGNETOM_Vida':1}, 'sequence': {'t1_mpr_adni_sag_iso_1mm_ns':0, '1_mp2rage_sag_p3_iso':1, 'MDEFT':2, 't1_mpr_sag_iso_1mm_ns':3}}

SOM = ScanOMetrics_project(bids_database, proc_pipeline, cov2float)

# Load all subjects using load_subject_options dictionary
load_subject_options = {'covariates_include':['age', 'sex', 'scanner', 'sequence'],
                        'sub_ses_exclude':['CH_BE_Insel_Healthy_HC000008_ses-2']}
for key in ['subjects_include', 'subjects_exclude', 'sub_ses_include', 'sub_ses_exclude',
            'covariates_include', 'covariates_exclude']:
    if key not in load_subject_options.keys():
        load_subject_options[key] = None
matching_covariates = ['sex', 'sequence', 'scanner']
SOM.load_subjects(load_subject_options["subjects_include"], load_subject_options["subjects_exclude"],
                  load_subject_options["sub_ses_include"], load_subject_options["sub_ses_exclude"],
                  load_subject_options["covariates_include"], load_subject_options["covariates_exclude"])
subject_IDmergedSes = []
for sub in SOM.subject:
    for ses in SOM.subject[sub]:
        subject_IDmergedSes.append('%s_%s' % (sub, ses))
SOM.load_proc_metrics(subject_IDmergedSes)
ROIs = ['caudalmiddlefrontal', 'lateralorbitofrontal', 'medialorbitofrontal',
        'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis',
        'precentral', 'rostralmiddlefrontal', 'superiorfrontal', 'frontalpole',
        'inferiorparietal', 'postcentral', 'precuneus', 'superiorparietal', 'supramarginal',
        'cuneus', 'lateraloccipital', 'lingual', 'pericalcarine',
        'bankssts', 'entorhinal', 'fusiform', 'inferiortemporal', 'middletemporal',
        'parahippocampal', 'superiortemporal', 'temporalpole', 'transversetemporal',
        'caudalanteriorcingulate', 'isthmuscingulate', 'posteriorcingulate',
        'rostralanteriorcingulate',
        'insula']

for hemi in ['lh', 'rh']:
    col_idx = np.array([SOM.metric_names.index('aparc_%s_%s_area' % (hemi, ctx_label)) for ctx_label in ROIs])
    norm_area = SOM.measured_metrics['orig'][:, col_idx].copy()
    norm_area /= norm_area.sum(axis=1, keepdims=True)
    col_idx = np.array([SOM.metric_names.index('aparc_%s_%s_gauscurv' % (hemi, ctx_label)) for ctx_label in ROIs])
    gauscurvs = SOM.measured_metrics['orig'][:, col_idx]
    norm_gauscurvs = (gauscurvs*norm_area).sum(axis=1)
    fig1 = plt.figure()
    ax1 = fig1.gca()
    ax1.plot(SOM.covariate_values[:, SOM.covariate_names.index('age')], norm_gauscurvs, '+')
    ax1.grid()
    fig1.show()
