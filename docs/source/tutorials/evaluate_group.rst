Group comparison
================

This tutorial illustrates how to perform a statistical comparison between a group of scans and
matching subjects within the normative model.

Loading and setting up subjects and normative model
***************************************************

First create a SOM instance with the bids directory containing your subjects of interest::

    from scanometrics.core import ScanOMetrics_project
    bids_directory = "~/Downloads/OASIS3"  # <-- example to be edited
    SOM = ScanOMetrics_project(bids_directory)

To list the available models in the scanometrics default folder, you can run
`SOM.list_normative_models()`. To load a given model you can run::

    model_name = "Polynomial_dldirect_cleanOASIS3"  # <-- example to be edited
    SOM.load_normative_model(model_name)

where <model_name> is the identifier for the model you want to load. Here, `cleanOASIS3` refers to
the OASIS3 healthy control dataset, after removing scans that had more than 5% of abnormal
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

You are now ready to load your subjects of interest. The default is to load all subjects available
in the BIDS directory (filtering of subjects to be tested can be done afterwards, see below)::

    SOM.load_subjects()

To load subjects data, we use the `load_proc_metrics()` method, which requires setting the normative
model as reference for metric normalization. We use the 'orig' data in the normative model to do so::

    SOM.load_proc_metrics(ref_metric_values=SOM.normativeModel.measured_metrics['orig'].copy(),
                          ref_metric_names=SOM.normativeModel.metric_names.copy(),
                          ref_covariate_values=SOM.normativeModel.covariate_values.copy(),
                          ref_covariate_names=SOM.normativeModel.covariate_names.copy())

Statistical comparison
**********************

To compare the loaded subjects against the normative model, you can use the `test_group_differences()` method::

    group_analysis = SOM.test_group_differences(matching_covariates=['sex', 'scanner', 'sequence'])

As group comparisons is usually used for a certain condition, you can select the set of subjects
using the `group_covariate` and `group_label` optional arguments. The `group_covariate` variable
should be an array of integers, usually a column of SOM.covariate_values, with the same length
as the number of loaded subjects. The `group_label` should be the value for the subjects you are
interested into. For example, if you would like to test scans with severe AD in the OASIS3 dataset,
you could run the following::

    # Perform analysis on CDR=3 subjects (maps to 4 in cov2float)
    group_analysis = SOM.test_group_differences(group_covariate=SOM.covariate_values[:, SOM.covariate_names.index('diagnosis_subtype'),
                                                group_label=4,
                                                matching_covariates=['sex', 'scanner', 'sequence'])

The `group_analysis` variable is a dictionary with several keys:
    - 'gp_residuals': the residuals of each metric available in both the loaded data and the normative model, for each loaded subject.
    - 'gp-vs-gp_dfs': the degrees of freedom for each tested metric.
    - 'gp-vs-gp_ts': the value of the statistical variable for each metric.
    - 'gp-vs-gp_pvals': the p-value of each test.
    - 'gp-vs-gp_logps': the signed log10 of the p-value.
    - 'gp-vs-gp_cohen-d': the effect size of each tested metric.

Full script
***********
::

    from scanometrics.core import ScanOMetrics_project
    # Editable variables: bids directory, model_name, and subject_id
    bids_directory = "~/Downloads/OASIS3"  # <-- to be edited
    model_name = "Polynomial_dldirect_cleanOASIS3"  # <-- to be edited. Use SOM.list_normative_models() to list available models
    group_variable = 'diagnosis_subtype'  # <-- optional, to be edited. Name of the column in SOM.covariate_values to be used for selecting subjects to analyse
    group_label = 4  # <-- optional, to be edited. Selects subjects with CDR-3 as an example (corresponds to 4 in the covariate_values variable)
    # Evaluation of a subject against a normative model
    SOM = ScanOMetrics_project(bids_directory)
    SOM.load_normative_model(model_name)
    SOM.load_subjects()
    SOM.load_proc_metrics(ref_metric_values=SOM.normativeModel.measured_metrics['orig'].copy(),
                          ref_metric_names=SOM.normativeModel.metric_names.copy(),
                          ref_covariate_values=SOM.normativeModel.covariate_values.copy(),
                          ref_covariate_names=SOM.normativeModel.covariate_names.copy())
    group_analysis = SOM.test_group_differences(group_covariate=group_variable,
                                                group_label=group_label,
                                                matching_covariates=['sex', 'scanner', 'sequence'])
