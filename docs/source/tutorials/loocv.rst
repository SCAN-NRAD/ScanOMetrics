Leave One Out Cross-Validation (tutorial)
=========================================

This tutorial shows how to use ScanOMetrics to cross-validate a normative model and dataset, using Leave-One-Out Cross Validation (LOOCV).

As an example, we use the OASIS3 dataset, organized following BIDS guidelines.

Filtering healthy subjects
--------------------------

The following code assumes you are in the 'code' folder from the OASIS bids directory. It loads the OASIS3 subjects
information, including their diagnosis-subtype, which we use to filter out healthy subjects for which all scans have a
CDR = 0:

::


    from scanometrics.utils.logging import set_verbose
    from scanometrics.core import ScanOMetrics_project
    import numpy as np
    from multiprocessing import Process, Pool
    import os
    import pickle
    from itertools import repeat
    from time import sleep

    # Define mapping of categorical covariates to numerical values
    cov2float = {'sex': {'F': 1, 'M': 0},
                 'diagnosis_subtype': {'CDR-0': 0, 'CDR-0.5': 1, 'CDR-1': 2, 'CDR-2': 3, 'CDR-3': 4},
                 'scanner_type': {'TrioTim': 0, 'Avanto': 1, 'Sonata': 2, 'Biograph_mMR': 3},
                 'scanner': {'OASIS-21926': 0, 'OASIS-26321': 1, 'OASIS-35177': 2, 'OASIS-35248': 3, 'OASIS-51010': 4, 'OASIS-NA': 5},
                 'sequence': {'': 0, 'MPRAGE_GRAPPA2': 1, 'Sagittal_3D_Accelerated_T1': 2, 't1_mpr_1mm_p2': 3, 't1_mpr_1mm_p2_pos50': 4, 't1_mpr_ns_sag': 5, 't1_mprage_sag_isoWU': 6}
                 }

    # Load subjects and extract list of healthy subject IDs
    bids_directory = '..'
    SOM = ScanOMetrics_project(bids_directory, proc_pipeline='dldirect', cov2float=cov2float.copy(), n_threads=1)
    SOM.load_subjects()
    HC_subj_ids = []
    for subj_id in SOM.subject.keys():
        is_HC = True
        for ses_id in SOM.subject[subj_id].keys():
            for acq_id in SOM.subject[subj_id][ses_id].keys():
                if SOM.subject[subj_id][ses_id][acq_id]['diagnosis_subtype'] != 0:
                    is_HC = False
        if is_HC:
            HC_subj_ids.append(subj_id)

Then the following code defines the LOOCV function:

::


    # Define LOOCV analysis for a single subject
    def run_loocv(subject_to_exclude, HC_list, bids_folder):
        SOM_in = ScanOMetrics_project(shared_bids_folder, 'dldirect', acq_pattern='*T1w', dataset_id='OASIS3_dldirect_LOOCV-exclude-%s' % subject_to_exclude,
                                      cov2float=cov2float.copy(),
                                      n_threads=1)
        LOO_HC_list = HC_list.copy()
        LOO_HC_list.remove(subject_to_exclude)
        SOM_in.load_subjects(subjects_include=LOO_HC_list)
        SOM_in.load_proc_metrics()
        SOM_in.set_normative_model('Polynomial')
        SOM_in.normativeModel.flag_outliers(1.5, 'measured_metrics')
        SOM_in.normativeModel.uncertainty = {}
        for k in SOM_in.measured_metrics.keys():
            SOM_in.normativeModel.uncertainty[k] = np.ones(SOM_in.measured_metrics[k].shape[1])
        SOM_in.normativeModel.fit(flag_opt=1, N_cycl=100, global_deg_max=22)  # Leave default deg_max, frac=20, alpha=0.01
        SOM_in.normativeModel.flag_outliers(1.5, 'residuals')
        SOM_in.normativeModel.compute_uncertainty()
        SOM_in.normativeModel.fit(flag_opt=1, N_cycl=100)  # deg_max is computed depending on number of samples for the fit
        output_folder = os.path.join(shared_bids_folder, 'derivatives', 'scanometrics', SOM_in.normativeModel.model_dataset_id)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        SOM_in.save_normative_model(output_filename=os.path.join(output_folder, SOM_in.normativeModel.model_dataset_id+'.pkl'))
        SOM_out = ScanOMetrics_project(shared_bids_folder, 'dldirect', dataset_id='OASIS3_dldirect_LOOCV-exclude-%s' % subject_to_exclude,
                                       acq_pattern='*T1w', cov2float=cov2float.copy(),
                                       n_threads=1)
        SOM_out.load_subjects(subjects_include=[subject_to_exclude])
        SOM_out.normativeModel = SOM_in.normativeModel
        SOM_out.load_proc_metrics(ref_metric_values=SOM_out.normativeModel.measured_metrics['orig'].copy(),
                                  ref_metric_names=SOM_out.normativeModel.metric_names.copy(),
                                  ref_covariate_values=SOM_out.normativeModel.covariate_values.copy(),
                                  ref_covariate_names=SOM_out.normativeModel.covariate_names.copy())
        evaluation_results = SOM_out.evaluate_singleSubject_allSes(subject_to_exclude, ['sex', 'sequence', 'scanner'], create_html_report=False)
        output_folder = os.path.join(shared_bids_folder, 'derivatives', 'scanometrics', SOM_out.normativeModel.model_dataset_id,
                                     subject_to_exclude)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        with open(os.path.join(output_folder, 'loocv_evaluation_results.pkl'), 'wb') as fout:
            pickle.dump(evaluation_results, fout, pickle.HIGHEST_PROTOCOL)

    # Run LOOCV over all healthy subjects, in parallel using n_threads
    n_threads = 20
    with Pool(n_threads) as pool:
        pool.starmap(run_loocv, zip(HC_subj_ids, repeat(HC_subj_ids, len(HC_subj_ids))), bids_directory)



Full script
***********

::

    from scanometrics.utils.logging import set_verbose
    from scanometrics.core import ScanOMetrics_project
    import numpy as np
    from multiprocessing import Process, Pool
    import os
    import pickle
    from itertools import repeat
    from time import sleep

    # Define LOOCV analysis for a single subject
    def run_loocv(subject_to_exclude, HC_list, bids_folder):
        SOM_in = ScanOMetrics_project(shared_bids_folder, 'dldirect', acq_pattern='*T1w', dataset_id='OASIS3_dldirect_LOOCV-exclude-%s' % subject_to_exclude,
                                      cov2float=cov2float.copy(),
                                      n_threads=1)
        LOO_HC_list = HC_list.copy()
        LOO_HC_list.remove(subject_to_exclude)
        SOM_in.load_subjects(subjects_include=LOO_HC_list)
        SOM_in.load_proc_metrics()
        SOM_in.set_normative_model('Polynomial')
        SOM_in.normativeModel.flag_outliers(1.5, 'measured_metrics')
        SOM_in.normativeModel.uncertainty = {}
        for k in SOM_in.measured_metrics.keys():
            SOM_in.normativeModel.uncertainty[k] = np.ones(SOM_in.measured_metrics[k].shape[1])
        SOM_in.normativeModel.fit(flag_opt=1, N_cycl=100, global_deg_max=22)  # Leave default deg_max, frac=20, alpha=0.01
        SOM_in.normativeModel.flag_outliers(1.5, 'residuals')
        SOM_in.normativeModel.compute_uncertainty()
        SOM_in.normativeModel.fit(flag_opt=1, N_cycl=100)  # deg_max is computed depending on number of samples for the fit
        output_folder = os.path.join(shared_bids_folder, 'derivatives', 'scanometrics', SOM_in.normativeModel.model_dataset_id)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        SOM_in.save_normative_model(output_filename=os.path.join(output_folder, SOM_in.normativeModel.model_dataset_id+'.pkl'))
        SOM_out = ScanOMetrics_project(shared_bids_folder, 'dldirect', dataset_id='OASIS3_dldirect_LOOCV-exclude-%s' % subject_to_exclude,
                                       acq_pattern='*T1w', cov2float=cov2float.copy(),
                                       n_threads=1)
        SOM_out.load_subjects(subjects_include=[subject_to_exclude])
        SOM_out.normativeModel = SOM_in.normativeModel
        SOM_out.load_proc_metrics(ref_metric_values=SOM_out.normativeModel.measured_metrics['orig'].copy(),
                                  ref_metric_names=SOM_out.normativeModel.metric_names.copy(),
                                  ref_covariate_values=SOM_out.normativeModel.covariate_values.copy(),
                                  ref_covariate_names=SOM_out.normativeModel.covariate_names.copy())
        evaluation_results = SOM_out.evaluate_singleSubject_allSes(subject_to_exclude, ['sex', 'sequence', 'scanner'], create_html_report=False)
        output_folder = os.path.join(shared_bids_folder, 'derivatives', 'scanometrics', SOM_out.normativeModel.model_dataset_id,
                                     subject_to_exclude)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        with open(os.path.join(output_folder, 'loocv_evaluation_results.pkl'), 'wb') as fout:
            pickle.dump(evaluation_results, fout, pickle.HIGHEST_PROTOCOL)

    # Define mapping of categorical covariates to numerical values
    cov2float = {'sex': {'F': 1, 'M': 0},
                 'diagnosis_subtype': {'CDR-0': 0, 'CDR-0.5': 1, 'CDR-1': 2, 'CDR-2': 3, 'CDR-3': 4},
                 'scanner_type': {'TrioTim': 0, 'Avanto': 1, 'Sonata': 2, 'Biograph_mMR': 3},
                 'scanner': {'OASIS-21926': 0, 'OASIS-26321': 1, 'OASIS-35177': 2, 'OASIS-35248': 3, 'OASIS-51010': 4, 'OASIS-NA': 5},
                 'sequence': {'': 0, 'MPRAGE_GRAPPA2': 1, 'Sagittal_3D_Accelerated_T1': 2, 't1_mpr_1mm_p2': 3, 't1_mpr_1mm_p2_pos50': 4, 't1_mpr_ns_sag': 5, 't1_mprage_sag_isoWU': 6}
                 }

    # Load subjects and extract list of healthy subject IDs
    bids_directory = '..'
    SOM = ScanOMetrics_project(bids_directory, proc_pipeline='dldirect', cov2float=cov2float.copy(), n_threads=1)
    SOM.load_subjects()
    HC_subj_ids = []
    for subj_id in SOM.subject.keys():
        is_HC = True
        for ses_id in SOM.subject[subj_id].keys():
            for acq_id in SOM.subject[subj_id][ses_id].keys():
                if SOM.subject[subj_id][ses_id][acq_id]['diagnosis_subtype'] != 0:
                    is_HC = False
        if is_HC:
            HC_subj_ids.append(subj_id)

    # Run LOOCV over all healthy subjects, in parallel using n_threads
    n_threads = 20
    with Pool(n_threads) as pool:
        pool.starmap(run_loocv, zip(HC_subj_ids, repeat(HC_subj_ids, len(HC_subj_ids))), bids_directory)
