from scanometrics.cross_validation.loocv import LOOCV
import os

bids_database = 'CH-First-Insel'
proc_pipeline = 'freesurfer'
normative_model = 'Polynomial'
load_subject_options = {'cov2float': {'sex': {'M': 0, 'F': 1},
                                      'scanner': {'Skyra_fit': 0, 'MAGNETOM_Vida': 1},
                                      'sequence': {'t1_mpr_adni_sag_iso_1mm_ns': 0, 't1_mp2rage_sag_p3_iso': 1, 'MDEFT': 2, 't1_mpr_sag_iso_1mm_ns': 3}},
                        'covariates_include': ['age', 'sex', 'scanner', 'sequence'],
                        'sub_ses_exclude': ['CH_BE_Insel_Healthy_HC000008_ses-2']}
matching_covariates = ['sex', 'sequence', 'scanner']
outputs, pIDs = LOOCV(bids_database, proc_pipeline, normative_model, load_subject_options, 1.5, matching_covariates)

