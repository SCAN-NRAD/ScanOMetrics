Fit normative data (tutorial)
=============================

The following tutorial illustrates how to use the Scan-O-Metrics package to fit the polynomial normative model to
the OASIS3 dataset.

Make sure the Scan-O-Metrics package has been properly :ref:`installed <installation>`, and that you have the OASIS3 dataset
formated according to the BIDS recommendations.

Import the ScanOMetrics_project class, as well as `numpy` and `os` for further steps::

    from scanometrics.core import ScanOMetrics_project
    import numpy as np
    import os

Set the path to your BIDS dataset directory::

    bids_directory = "~/Downloads/OASIS3"  # <- example to be adapted

The initialisation of a Scan-O-Metrics project has an optional `cov2float` argument which can be used to set the mapping
from categorical covariates to numerical values. By default, it contains the mapping for a 'sex' covariate, using a
dictionary structure to convert males 'M' to 0 and females 'F' to 1. The OASIS3 dataset has other co-variates of interest,
which are mapped to numerical values using the following dictionary::

    cov2float = {'sex': {'F': 1, 'M': 0},
                 'diagnosis_subtype': {'CDR-0': 0, 'CDR-0.5': 1, 'CDR-1': 2, 'CDR-2': 3, 'CDR-3': 4},
                 'scanner_type': {'TrioTim': 0, 'Avanto': 1, 'Sonata': 2, 'Biograph_mMR': 3},
                 'scanner': {'OASIS-21926': 0, 'OASIS-26321': 1, 'OASIS-35177': 2, 'OASIS-35248': 3, 'OASIS-51010': 4, 'OASIS-NA': 5},
                 'sequence': {'': 0, 'MPRAGE_GRAPPA2': 1, 'Sagittal_3D_Accelerated_T1': 2, 't1_mpr_1mm_p2': 3, 't1_mpr_1mm_p2_pos50': 4, 't1_mpr_ns_sag': 5, 't1_mprage_sag_isoWU': 6}
                 }

By default, the proccesing pipeline is set to 'freesurfer'. Alternative pipelines include 'dldirect'. Both pipelines can
use either the `DesikanKilliany` or `Destrieux` parcellation atlases. To instanciate a Scan-O-Metrics project using the
`dldirect` processed data and the `DesikanKilliany` atlas, run the following::

    SOM = ScanOMetrics_project(bids_directory, proc_pipeline='dldirect', cov2float=cov2float, atlas='DesikanKilliany')

New subjects can then be loaded using the `load_subjects()` method. By default, `load_subjects()` loads all available
subjects and covariates it finds in `participants.tsv`, as well as all scans and covariates available in each subjects'
`<participant_id>_sessions.tsv` file. Here, we use the `diagnosis_subtype` covariate to identify healthy subjects that do
not have any scan with a CDR > 0::

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

We then reload the HC subjects only::

    SOM.load_subjects(subjects_include=HC_subj_ids)

The following code loads the subjects' measured metrics::

    SOM.load_proc_metrics()

The normative model (in this case a polynomial model) can be set with `set_normative_model()` method::

    SOM.set_normative_model('Polynomial')

Outliers based on standard deviation from samples with matching age can be flagged using `flag_outliers(sigma, reference)`
where `sigma` sets how many standard deviations defines an outlier, and `reference` allows switching between deviations in
'measured_metrics' data, or the 'residuals' after model fitting, which is used in a 2nd step (see below)::

    SOM.normativeModel.flag_outliers(1.5, 'measured_metrics')

Finally, the normative model can be fitted with `normativeModel.fit()`. As fitting relies on uncertainty estimates, the
uncertainty is first set to 1 for all variables::

    # Sets uncertainty to one as computation of uncertainty is done afterwards.
    SOM.normativeModel.uncertainty = {}
    for k in SOM.measured_metrics.keys():
        SOM.normativeModel.uncertainty[k] = np.ones(SOM.measured_metrics[k].shape[1])
    SOM.normativeModel.fit(flag_opt=1, N_cycl=100, global_deg_max=22)

Once there's a first estimate of the model, outliers can be flagged based on the deviation of residuals::

    SOM.normativeModel.flag_outliers(1.5, 'residuals')

This time the uncertainty can be computed::

    SOM.normativeModel.compute_uncertainty()

And the model recomputed with updated outliers and uncertainty::

    SOM.normativeModel.fit(flag_opt=1, N_cycl=100)  # deg_max is computed depending on number of samples for the fit

Finally, the estimated model can be saved using `save_normative_model()`. The function will
use the dataset and model names to create an ID name to save the model. Models available
in the default folder can be listed with `SOM.list_normative_models()`, and the model of
interest can be loaded with `SOM.load_normative_model(model_id)`.::

    output_folder = os.path.join(bids_directory, 'derivatives', 'scanometrics', SOM.normativeModel.model_dataset_id)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    SOM.save_normative_model(output_filename=os.path.join(output_folder, SOM.normativeModel.model_dataset_id + '.pkl'))


Full script
***********

::

    from scanometrics.core import ScanOMetrics_project
    import numpy as np
    import os

    # Set path to bids_directory
    bids_directory = "~/Downloads/OASIS3"  # <- example to be adapted
    # Map categorical covariates to numerical values
    cov2float = {'sex': {'F': 1, 'M': 0},
                 'diagnosis_subtype': {'CDR-0': 0, 'CDR-0.5': 1, 'CDR-1': 2, 'CDR-2': 3, 'CDR-3': 4},
                 'scanner_type': {'TrioTim': 0, 'Avanto': 1, 'Sonata': 2, 'Biograph_mMR': 3},
                 'scanner': {'OASIS-21926': 0, 'OASIS-26321': 1, 'OASIS-35177': 2, 'OASIS-35248': 3, 'OASIS-51010': 4, 'OASIS-NA': 5},
                 'sequence': {'': 0, 'MPRAGE_GRAPPA2': 1, 'Sagittal_3D_Accelerated_T1': 2, 't1_mpr_1mm_p2': 3, 't1_mpr_1mm_p2_pos50': 4, 't1_mpr_ns_sag': 5, 't1_mprage_sag_isoWU': 6}
                 }
    # Create Scan-O-metrics instance
    SOM = ScanOMetrics_project(bids_directory, proc_pipeline='dldirect', cov2float=cov2float)
    # Load all subjects in the bids 'participants.tsv' file
    SOM.load_subjects()
    # Filter out HC subjects (no scan with CDR > 0)
    HC_subj_ids = []
    for subj_id in SOM.subject.keys():
        is_HC = True
        for ses_id in SOM.subject[subj_id].keys():
            for acq_id in SOM.subject[subj_id][ses_id].keys():
                if SOM.subject[subj_id][ses_id][acq_id]['diagnosis_subtype'] != 0:
                    is_HC = False
        if is_HC:
            HC_subj_ids.append(subj_id)
    # Reload HC subjects only
    SOM.load_subjects(subjects_include=HC_subj_ids)
    # Load HC metrics
    SOM.load_proc_metrics()
    # Set normative model
    SOM.set_normative_model('Polynomial')
    # Flag outliers based on deviations from measurements
    SOM.normativeModel.flag_outliers(1.5, 'measured_metrics')
    # Sets uncertainty to one for an initial model fit
    SOM.normativeModel.uncertainty = {}
    for k in SOM.measured_metrics.keys():
        SOM.normativeModel.uncertainty[k] = np.ones(SOM.measured_metrics[k].shape[1])
    SOM.normativeModel.fit(flag_opt=1, N_cycl=100, global_deg_max=22)
    # Flat outliers based on deviations from residuals
    SOM.normativeModel.flag_outliers(1.5, 'residuals')
    # Compute uncertainty based on the 1st model fit and updated outliers
    SOM.normativeModel.compute_uncertainty()
    # Fit the model with updated outliers and uncertainties
    SOM.normativeModel.fit(flag_opt=1, N_cycl=100)  # deg_max is computed depending on number of samples for the fit
    # Save fitted model to output folder
    output_folder = os.path.join(bids_directory, 'derivatives', 'scanometrics', SOM.normativeModel.model_dataset_id)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    SOM.save_normative_model(output_filename=os.path.join(output_folder, SOM.normativeModel.model_dataset_id + '.pkl'))
