"""
Leave-one out cross-validation script. Allows user to set a bids directory as input to train a 'normative_model' on data
processed with a 'proc_pipeline', and test the training on a leaved-one-out subject.
"""


from scanometrics.core import ScanOMetrics_project
import numpy as np


def LOOCV(bids_database, proc_pipeline='freesurfer', normative_model='Polynomial', cov2float=None, load_subject_options={},
          stats2table_folder=None, sigma_outliers=1.5, matching_covariates=['age', 'sex', 'sequence', 'scanner'], n_cycl=100,
          save_loocv_reports=False):
    """Cross-validate dataset and model of interest using Leave-One-Out Cross-Validation.

    :param bids_database: path to bids directory
    :type bids_database: string
    :param proc_pipeline: name of pipeline used to process bids database, to be used as input for the normative model.
                          e.g.: 'freesurfer', 'dldirect'
    :type proc_pipeline: string
    :param normative_model: name of normative model (e.g.: 'Polynomial')
    :type normative_model: string
    :param cov2float: dictionary mapping categorical covariates to float.
    :type cov2float: dictionary
    :param load_subject_options: dictionary to gather loading options passed as argument to the SOM initializer (i.e.
                                 'cov2float') and SOM.load_subjects() method (i.e. 'cov2float', 'subjects_include',
                                 'subjects_exclude', 'sub_ses_include', 'sub_ses_exclude', 'covariates_include', and
                                 'covariates_exclude')
    :type load_subject_options: dictionary
    :param stats2table_folder: path to folder containing stats2table files
    :type stats2table_folder: string
    :param sigma_outliers: standard deviation threshould above which to consider outliers
    :type sigma_outliers: float
    :param matching_covariates: list of covariate names to be used for matching selection.
    :type matching_covariates: list
    :param n_cycl: number of samplings to perform model estimation.
    :type n_cycl: int.
    :param save_loocv_reports: bool option specifying whether to save loocv html reports (time and memory consuming).
    :type save_loocv_reports: bool
    """


    # Create SOM instance with BIDS database folder, ran proc_pipeline and cov2float conversion dictionary
    # ====================================================================================================
    SOM_all = ScanOMetrics_project(bids_database, proc_pipeline, cov2float=cov2float)

    # Load all subjects using load_subject_options dictionary
    # =======================================================
    SOM_all.load_subjects(**load_subject_options)

    # Load metrics for all subjects and sessions available
    if stats2table_folder is None:
        SOM_all.metric_proc_pipeline.load_from_ID = True
        SOM_all.load_proc_metrics()
    else:
        SOM_all.metric_proc_pipeline.load_from_ID = False
        SOM_all.load_proc_metrics(stats2table_folder=stats2table_folder)

    # Get list of all subjects for the leave-one-out loop
    all_subjects = list(SOM_all.subject.keys())

    # Statistics to be gathered in outputs and pIDs
    outputs = []
    pIDs = np.zeros(len(all_subjects))
    for subj_i, subject_to_exclude in enumerate(all_subjects):
        # Create SOM instance excluding current subject
        SOM_in = ScanOMetrics_project(bids_database, proc_pipeline, dataset_id='%s_LOOCV-%03d' % (SOM_all.dataset_id, subj_i+1), cov2float=cov2float)
        # Change load_subject_options considering subject to exclude
        options_with_exclude = load_subject_options.copy()
        # Remove subject to exclude from subjects_include in case its defined there
        if options_with_exclude['subjects_include'] is not None and subject_to_exclude in options_with_exclude['subjects_include']:
            options_with_exclude['subjects_include'].remove(subject_to_exclude)
        # Same if subject to exclude is part of sub_ses_include
        if options_with_exclude['sub_ses_include'] is not None:
            for sub_ses in options_with_exclude['sub_ses_include']:
                if subject_to_exclude in sub_ses:
                    options_with_exclude['sub_ses_include'].remove(sub_ses)
        # Add subject to exclude to subject_exclude
        if options_with_exclude['subjects_exclude'] is not None:
            options_with_exclude['subjects_exclude'].append(subject_to_exclude)
        else:
            options_with_exclude['subjects_exclude'] = [subject_to_exclude]
        # Load included subjects
        SOM_in.load_subjects(**options_with_exclude)

        if stats2table_folder is None:
            SOM_in.metric_proc_pipeline.load_from_ID = True
            SOM_in.load_proc_metrics()
        else:
            SOM_in.metric_proc_pipeline.load_from_ID = False
            SOM_in.load_proc_metrics(stats2table_folder=stats2table_folder)

        # Define normative model to train on the include dataset.
        SOM_in.set_normative_model(normative_model)
        # Flag outliers based on sigma_outliers deviation from distribution of measured_metrics
        SOM_in.normativeModel.flag_outliers(sigma_outliers, 'measured_metrics')
        # Set uncertainty to 1.0 as initial value
        SOM_in.normativeModel.uncertainty = {}
        for k in SOM_in.measured_metrics.keys():
            SOM_in.normativeModel.uncertainty[k] = np.ones(SOM_in.measured_metrics[k].shape[1])
        # Fit normative model
        SOM_in.normativeModel.fit(flag_opt=1, N_cycl=n_cycl)  # Leave default deg_max, frac=20, alpha=0.01
        # Reflag outliers based on residuals
        SOM_in.normativeModel.flag_outliers(sigma_outliers, 'residuals')

        # Compute uncertainty
        SOM_in.normativeModel.compute_uncertainty()

        # Fit again with clean data
        SOM_in.normativeModel.fit(flag_opt=1,
                                  deg_max=int(np.floor(len(SOM_in.covariate_values[:, SOM_in.covariate_names.index('age')])/(2*20))),
                                  N_cycl=n_cycl)  # reduce deg_max by a factor two, N_cycl=100

        # Evaluate excluded subject against model trained on others
        SOM_all.normativeModel = SOM_in.normativeModel

        tmp_output, pIDs[subj_i] = SOM_all.evaluate_singleSubject_allSes(subject_to_exclude,
                                                                         ['sex', 'sequence', 'scanner'],
                                                                         create_html_report=save_loocv_reports)
        outputs.append(tmp_output)
    return outputs, pIDs
