Evaluate single subject
=======================

This tutorial evaluates a single subject against a normative model.

Loading and setting up subject and normative model
**************************************************

First create a SOM instance with the bids directory containing your subject of interest::

    from scanometrics.core import ScanOMetrics_project
    bids_directory = "~/Downloads/OASIS3"  # <-- example to be edited
    SOM = ScanOMetrics_project(bids_directory)

To list the available models in the scanometrics default folder, you can run
`SOM.list_normative_models()`. To load a given model you can run::

    model_name = "Polynomial_dldirect_v1-0-3_DesikanKilliany_OASIS3_som-v0-1-0.pkl"  # <-- example to be edited
    SOM.load_normative_model(model_name)

where `<model_name>` is the identifier for the model you want to load. Here, `OASIS3` refers to
the OASIS3 healthy control dataset, after removing scans that had more than 7.6% of abnormal
regions in a Leave-one-out cross validation.

Normative models usually rely on age and other covariate values, some of which are categorical
and need conversion to numerical values for their effect to be accounted for. The mapping is loaded
along with the normative model, in a dictionary called `cov2float`. This makes sure that the
mapping done when loading new subjects is the same that was used for training the
normative model (e.g.: sex mapping to 0 for males and 1 for females). If your subject has new
covariates compared to the training data, and which should be mapped to a numerical value,
you should add the mapping function to the `SOM.cov2float` dictionary. For example, the sex
covariate is mapped to a numerical value using the following key/value pair::

    print(SOM.cov2float['sex'])
    # out[]: {'M':0, 'F':1}

You are now ready to load a new subject::

    subject_id = "sub-OAS0001"  # <-- example to be edited
    SOM.load_subjects(subjects_include=[subject_id])

To load the new subject data, we use the `load_proc_metrics()` method, which requires setting the normative
model as reference for metric normalization. We use the 'orig' data in the normative model to do so::

    SOM.load_proc_metrics(ref_metric_values=SOM.normativeModel.measured_metrics['orig'],
                          ref_metric_names=SOM.normativeModel.metric_names,
                          ref_covariate_values=SOM.normativeModel.covariate_values,
                          ref_covariate_names=SOM.normativeModel.covariate_names)

Comparison against normative curves
***********************************

To compare the new subject against normative curves, you can use the èvaluate_singleSubject_allSes()` method::

    output, pID = SOM.evaluate_singleSubject_allSes(subject_id, matching_covariates= ['sex', 'sequence', 'scanner'], create_html_report=True)

Full script
***********
::

    from scanometrics.core import ScanOMetrics_project
    # Editable variables: bids directory, model_name, and subject_id
    bids_directory = "~/Downloads/OASIS3"  # <-- to be edited
    model_name = "Polynomial_dldirect_v1-0-3_DesikanKilliany_OASIS3_som-v0-1-0.pkl"  # <-- to be edited. Use SOM.list_normative_models() to list available models
    subject_id = "sub-OAS0001"  # <-- to be edited
    # Evaluation of a subject against a normative model
    SOM = ScanOMetrics_project(bids_directory)
    SOM.load_normative_model(model_name)
    SOM.load_subjects(subjects_include=[subject_id])
    SOM.load_proc_metrics(ref_metric_values=SOM.normativeModel.measured_metrics['orig'],
                          ref_metric_names=SOM.normativeModel.metric_names,
                          ref_covariate_values=SOM.normativeModel.covariate_values,
                          ref_covariate_names=SOM.normativeModel.covariate_names)
    output, pID = SOM.evaluate_singleSubject_allSes(subject_id, matching_covariates= ['sex', 'sequence', 'scanner'], create_html_report=True)
