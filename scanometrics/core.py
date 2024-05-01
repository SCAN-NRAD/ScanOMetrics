"""
Core ScanOMetrics classes and methods
"""


import pickle
import numpy as np
import os
from scanometrics import normative_models, processing, __version__
from scanometrics.utils import logging
from scanometrics.utils.stats import fdr
from csv import DictReader
from scipy.stats import norm, ttest_ind, chi2
from numpy.polynomial import Polynomial as poly_model
import marshal
import types
from shutil import copytree
import matplotlib.pyplot as plt
from glob import glob
import re

class ScanOMetrics_project:
    """
    Class defining a ScanOMetrics project.
    Subjects should be stored according to BIDS data structure for multiple sessions: https://bids.neuroimaging.io/
    Participants IDs and fixed co-variates should be saved in <bids_database>/participants.tsv
    Session specific co-variates (e.g. age, sequence, scanner) should be saved in the respective session.tsv files
    """

    ###############
    # Constructor #
    ###############

    def __init__(self, bids_database, proc_pipeline='dldirect', dataset_id=None, cov2float={'sex': {'M': 0, 'm': 0,
        'F': 1, 'f': 1}}, acq_pattern='*T1w', ses_delimiter="_", acq_delimiter="_", n_threads=-1):
        """
        ScanOMetrics_project constructor from bids_database path. BIDS database should at least contain a participants.tsv
        file to load subject names and fixed covariate_values. Pipeline for data processing can be set through the
        `proc_pipeline` argument (defaults to 'dldirect'). Dataset ID can be set through the dataset_id argument (defaults to
        <proc_pipeline>_<bids_database>). The cov2float argument is used to convert categorical covariate values to
        numerical values. It is defined here to be saved with the project and avoid changes during data loading, project
        saving, etc... The acq_pattern argument sets the sequence basename to use as main image to process (defaults to
        '*T1w', useful to process subjects with multiple acquisition sequences).

        :param bids_database: path to BIDS database, can be relative to pwd or absolute
        :type bids_database: string
        :param proc_pipeline: name of pipeline for processing MRI volumes (scanometrics/processing/<proc_pipeline>.py
                              should be an existing file)
        :type proc_pipeline: string
        :param dataset_id: unique string identifying the processed bids dataset. Defaults: <proc_pipeline>_<bids_database>
                           (e.g. dldirect_OASIS3)
        :type dataset_id: string
        :param cov2float: dictionary mapping categorical covariate_values to float (defaults to {'sex': {'M':0, 'm':0,
                          'F': 1, 'f': 1}})
        :type cov2float: dictionary
        :param acq_pattern: pattern to look for with glob() when looking for anatomical nifti files. Influences the acq
                            label assigned to the corresponding scans.
        :type acq_pattern: string
        :param n_threads: sets number of threads to run methods with multithreading options
        :type n_threads: int
        """
        # Set directory of bids database with processed subjects
        self.bids_database = bids_database
        # Reset list of subjects (loaded by user calling self.load_subjects())
        self.subject = {}
        self.ses_delimiter = ses_delimiter
        self.acq_delimiter = acq_delimiter
        # Reset metrics array and metric_names, to be filled when calling self.load_subjects()
        self.measured_metrics = {'orig': np.array([], dtype='float'), 'norm': np.array([], dtype='float')}
        self.covariate_values = np.array([], dtype='float')
        self.cov2float = cov2float.copy()
        self.metric_names = []
        self.metric_plotting_info = {}
        self.covariate_names = []
        # Define and set processing pipeline used to generate the measured metrics
        self.metric_proc_pipeline = None
        self.set_proc_pipeline(proc_pipeline)
        # Reset normative model
        self.normativeModel = None
        # Set dataset_id
        if dataset_id is None:
            dataset_id = '_'.join([proc_pipeline, os.path.basename(os.path.normpath(self.bids_database))])
        self.dataset_id = dataset_id
        self.acq_pattern = acq_pattern
        # Set number of threads to maximum if set to -1
        self.n_threads = n_threads

    ###################
    # SUBJECT LOADING #
    ###################

    def add_ses_row(self, ID, row, ses_id, ses_row, acqs=['T1w']):
        """
        Add a subject's session to self.covariate_values array and self.subject dictionary (appends to existing
        sessions or creates new subject). Missing covariates are set to np.nan. Sessions are combined with acq labels to
        allows subjects to have multiple scan inputs.

        :param ID: subject_id, usually taken from <bids_database>/participants.tsv table
        :type ID: string
        :param row: dictionary with subject's covariates taken from <bids_database>/participants.tsv table
        :type row: dict
        :param ses_id: current session_id as taken from <bids_database>/<subject_id>/sessions.tsv table
        :type ses_id: string
        :param ses_row: dictionary with subject's session specific covariate values, taken from sessions.tsv table
        :type ses_row: dict
        :param acqs: acquisition labels to use as input. A session is created for each session/acq label combination.
        :type acqs: string
        """
        # Transfer all session covariate values to row (handy as overwrites potential values there)
        for k in ses_row.keys():
            row[k] = ses_row[k]
        # Try converting everything to float, ValueErrors might rise from categorical variables, try cov2float and error exit if doesn't work
        for k in row.keys():
            try:
                row[k] = np.float64(row[k])
            except ValueError as ve:
                if k in self.cov2float.keys():
                    row[k] = self.cov2float[k][row[k]]
                else:
                    logging.ERROR('ValueError for %s key (' % (k) + str(ve) + ') when running load_subjects(). Try adding categorical to'
                                  'numerical mapping in the cov2float dictionary when creating the ScanOMetrics_project'
                                  'and run load_subjects() again. Check that input covariate files do not have "NA", '
                                  '"na", "N/A" etc... entries instead of "NaN", "nan", or "NAN" and replace them '
                                  'accordingly')
        # Fill missing values with nan
        for k in [k for k in self.covariate_names if k not in row.keys()]:
            row[k] = np.nan
        for acq in acqs:
            # Fill in self.covariate_values according to self.covariate_names order
            self.covariate_values.append(np.array([row[k] for k in self.covariate_names], dtype='float')[None, :])
            # Keep track of row subject ID, ses_id and acq in self.subject
            if ID in self.subject.keys():
                # Subject already exists, add new session to session dictionary. Copy `row` as values can get overwritten
                # in the next session, and would overwrite values in self.subject
                if ses_id in self.subject[ID]:
                    self.subject[ID][ses_id][acq] = row.copy()
                else:
                    self.subject[ID][ses_id] = {acq: row.copy()}
            else:
                # Subject is new, add a session dictionary
                self.subject[ID] = {ses_id: {acq: row.copy()}}

    def load_subjects(self, subjects_include=None, subjects_exclude=None, sub_ses_include=None, sub_ses_exclude=None,
                      sub_ses_acq_include=None, sub_ses_acq_exclude=None, covariates_include=None, covariates_exclude=None,
                      subjects_table=None, subj_ses_acq_pattern=r"^(?P<subjID>sub-.+?)(?=(_ses-|$))_(?P<sesID>ses-.+?)(?:(?=(_acq-|$))_(?P<acqID>acq-.+))?"):
        """Load subject data from <self.bids_database>/participants.tsv, with inclusion/exclusion lists based on subject,
        session and acq labels, and inclusion/exclusion of covariate_values based on covariate name. Also keeps track of
        session/acquisition labels, and adds relevant information to self.subject and self.covariate_values. Repeated
        measures can be ignored by including/excluding specific sessions at runtime. Assumes repeated measures are a
        separate scan in <self.bids_database>/sub-<ID>/ses-<label>/anat/sub-<ID>_ses-<label>_acq-<label>.nii.gz from
        which all metrics are extracted. Session info should be in <self.bids_database>/sub-<ID>/sub-<ID>_sessions.tsv
        (one row per session, with a mandatory column named 'session_id' and parameter(s) that change(s) between
        sessions). The function first loops through participants.tsv and session.tsv files to gather subject IDs and
        session IDs, then filters out elements based on exclusion/inclusion criteria. Current implementation assumes a
        single session if there are no session.tsv files, and a session_id is automatically added, following bids
        recommendation that only a single row per subject must appear in the participants.tsv file.
        NB: everything relies on participants.tsv and sessions.tsv files. participants.tsv MUST have only one line per
        subject, and repeated scans MUST be encoded through a <participant_id>_sessions.tsv file in its bids folder.

        :param subjects_include: list of subjects to include. Defaults to None to load all subjects in participants.tsv
        :param subjects_exclude: list of subjects to exclude. Defaults to None to load all subjects in subjects_include
        :param sub_ses_include: list of sub-<subject_id>_ses-<session_id> to include. Defaults to None to load all ses.
        :param sub_ses_exclude: list of sub-<subject_id>_ses-<session_id> to exclude. Defaults to None to laod all ses.
        :param sub_ses_acq_include: list of sub-<subject_id>_ses-<session_id>_<acq_label> to include. Defaults to None to load all acqs.
        :param sub_ses_acq_exclude: list of sub-<subject_id>_ses-<session_id>_<acq_label> to exclude. Defaults to None to laod all acqs.
        :param covariates_include: list of covariate_values to include. Defaults to None to load all covariate_values in participants.tsv
        :param covariates_exclude: list of covariate_values to exclude. Defaults to None to load all covariate_values in covariates_include
        :param subjects_table: path to single tsv table with all subject names and covariates, overwrites BIDS scrapping
        :param subj_ses_pattern: binary string used as regex expression to split subject and session IDs from scan_id
                                 when scrapping participants from a single table instead of BIDS structure.
        """
        # TODO: test what happens if the session.tsv exists because of repeats regarding modalities unrelated to
        #  anatomical scans ? In that case the sessions.tsv will contain covariates that change with the unrelated scans,
        #  but the subject might have a single anatomical scan... Might not be an issue as the acq list should be empty
        #  when adding such sessions, but it should be double checked.

        # Initialize subject and covariate_values to empty dict and array
        self.subject = {}
        self.covariate_values = []  # Start with a list to be converted to numpy array after looping through tsv files
        self.covariate_names = []
        # Load from a single table if subjects_table is not None
        if subjects_table is not None:
            with open(subjects_table, 'r') as f_in:
                dict_reader = DictReader(f_in, delimiter='\t')
                if covariates_include is not None:
                    self.covariate_names = covariates_include
                else:
                    self.covariate_names = dict_reader.fieldnames.copy()
                if 'participant_id' in self.covariate_names:
                    self.covariate_names.remove('participant_id')
                if 'session_id' in self.covariate_names:
                    self.covariate_names.remove('session_id')
                # Exclude covariate_names in covariates_exclude
                if covariates_exclude is not None:
                    for k in [t for t in covariates_exclude if t in self.covariate_names]:
                        print("Removing %s from 'covariate_names'" % k)
                        self.covariate_names.remove(k)
                if 'scan_id' in self.covariate_names:
                    self.covariate_names.remove('scan_id')
                if covariates_exclude is not None:
                    for covariate_name in [covariate_name for covariate_name in covariates_exclude if covariate_name in self.covariate_names]:
                        self.covariate_names.remove(covariate_name)
                for row in dict_reader:
                    scan_id = row['scan_id']
                    matches = re.match(subj_ses_acq_pattern, scan_id)
                    if matches:
                        subj_id = matches.group("subjID") or ""
                        ses_id = matches.group("sesID") or ""
                        acq_label = matches.group("acqID") or ""
                    else:
                        subj_id = scan_id
                        ses_id = ''
                        acq_label = ''
                    if subjects_include is not None and subj_id not in subjects_include:
                        continue
                    if subjects_exclude is not None and subj_id in subjects_exclude:
                        continue
                    if sub_ses_acq_include is not None and scan_id not in sub_ses_acq_include:
                        continue
                    if sub_ses_acq_exclude is not None and scan_id in sub_ses_acq_exclude:
                        continue
                    # Try converting everything to float, ValueErrors might rise from categorical variables, try cov2float and error exit if doesn't work
                    for k in self.covariate_names:
                        try:
                            row[k] = np.float64(row[k])
                        except ValueError as ve:
                            if k in self.cov2float.keys():
                                row[k] = self.cov2float[k][row[k]]
                            else:
                                logging.ERROR('ValueError for key %s, value "%s" and scan %s_%s_%s (error message was "' % (k, row[k], subj_id, ses_id, acq_label) + str(
                                    ve) + '") when running load_subjects(). Try adding categorical to'
                                          'numerical mapping in the cov2float dictionary when creating the ScanOMetrics_project'
                                          'and run load_subjects() again. Check that input covariate files do not have "NA", '
                                          '"na", "N/A" etc... entries instead of "NaN", "nan", or "NAN" and replace them '
                                          'accordingly')
                    self.covariate_values.append(np.array([row[k] for k in self.covariate_names], dtype='float')[None, :])
                    if subj_id not in self.subject.keys():
                        self.subject[subj_id] = {ses_id: {acq_label: {}}}
                        for covariate_name in self.covariate_names:
                            self.subject[subj_id][ses_id][acq_label][covariate_name] = row[covariate_name]
                    else:
                        if ses_id not in self.subject[subj_id].keys():
                            self.subject[subj_id][ses_id] = {acq_label: {}}
                            for covariate_name in self.covariate_names:
                                self.subject[subj_id][ses_id][acq_label][covariate_name] = row[covariate_name]
                        else:
                            if acq_label not in self.subject[subj_id][ses_id].keys():
                                self.subject[subj_id][ses_id][acq_label] = {}
                                for covariate_name in self.covariate_names:
                                    self.subject[subj_id][ses_id][acq_label][covariate_name] = row[covariate_name]
        else:
            # Build covariate_names
            if covariates_include is not None:
                # Set covariate_names to the list of covariates specified by the user
                self.covariate_names = covariates_include
            else:
                # If list is None, loop through participants and sub-<label>_sessions tsvs to collect unique covariate_names
                with open(os.path.join(self.bids_database, 'participants.tsv'), 'r') as f:
                    reader = DictReader(f, delimiter='\t')
                    self.covariate_names = reader.fieldnames.copy()
                    for row in reader:
                        ID = row.pop('participant_id')
                        ses_file = os.path.join(self.bids_database, ID, ID + '_sessions.tsv')
                        if os.path.exists(ses_file):
                            with open(ses_file, 'r') as f2:
                                reader2 = DictReader(f2, delimiter='\t')
                                self.covariate_names += [k for k in reader2.fieldnames if k not in self.covariate_names]
            if 'participant_id' in self.covariate_names:
                self.covariate_names.remove('participant_id')
            if 'session_id' in self.covariate_names:
                self.covariate_names.remove('session_id')
            # Exclude covariate_names in covariates_exclude
            if covariates_exclude is not None:
                for k in [t for t in covariates_exclude if t in self.covariate_names]:
                    logging.PRINT("Removing %s from 'covariate_names'" % k)
                    self.covariate_names.remove(k)
            logging.PRINT('Selected covariate_names are %s' % self.covariate_names)
            # Loop again, this time filling values for each subject
            with open(os.path.join(self.bids_database, 'participants.tsv'), 'r') as f:
                reader = DictReader(f, delimiter='\t')
                for row in reader:
                    logging.PRINT(row)
                    ID = row.pop('participant_id')
                    # Filter out subject not in subjects_include or in subjects_exclude lists
                    if ((subjects_include is not None) and (ID not in subjects_include)) or\
                       ((subjects_exclude is not None) and (ID in subjects_exclude)):
                        continue
                    # Keep row keys that are in self.covariate_names
                    row = {k: row[k] for k in self.covariate_names if k in row.keys()}
                    # Check if sessions.tsv file exists and add values to row
                    ses_file = os.path.join(self.bids_database, ID, ID + '_sessions.tsv')
                    logging.PRINT('ses_file=%s' % ses_file)
                    if os.path.exists(ses_file):
                        logging.PRINT('ses_file exists, opening...')
                        with open(ses_file, 'r') as f_ses:
                            ses_reader = DictReader(f_ses, delimiter='\t')
                            for ses_row in ses_reader:  # Loop through repeats, update row dictionary each time
                                ses_id = ses_row.pop('session_id')
                                logging.PRINT('ses_id=%s, sub_ses_exclude=%s, sub_ses_include=%s' % (ses_id, sub_ses_exclude, sub_ses_include))
                                if ((sub_ses_exclude is not None) and ((ID+'_'+ses_id) in sub_ses_exclude))\
                                or ((sub_ses_include is not None) and ((ID+'_'+ses_id) not in sub_ses_include)):
                                    continue
                                logging.PRINT('Adding subject %s (session %s)' % (ID, ses_id))
                                if sub_ses_acq_include is not None:
                                    sub_acqs = []
                                    for sub_ses_acq in sub_ses_acq_include:
                                        if '%s_%s_' % (ID, ses_id) in sub_ses_acq:
                                            sub_acqs.append(sub_ses_acq.removeprefix('%s_%s_' % (ID, ses_id)))
                                else:
                                    sub_acqs = [os.path.basename(f).removeprefix('%s_%s_' % (ID, ses_id)).removesuffix('.gz').removesuffix('.nii') for f in glob(os.path.join(self.bids_database, ID, ses_id, 'anat', self.acq_pattern + '.nii*'))]
                                    if sub_ses_acq_exclude is not None:
                                        for sub_ses_acq in sub_ses_acq_exclude:
                                            if '%s_%s_' % (ID, ses_id) in sub_ses_acq and sub_ses_acq.removeprefix('%s_%s_' % (ID, ses_id)) in sub_acqs:
                                                sub_acqs.remove(sub_ses_acq.removeprefix('%s_%s_' % (ID, ses_id)))
                                self.add_ses_row(ID, row, ses_id, ses_row, sub_acqs)
                    else:
                        # No session files were found, directory should already contain anat, dwi, func, etc... folders
                        logging.WARNING('Subject %s has no session.tsv file. We recommend adopting a folder structure with'
                                        ' sessions folders to keep track of longitudinal data.' % ID)
                        logging.PRINT('Adding subject %s (no session)' % ID)
                        # list .nii and .nii.gz files matching acqs (can use wildcards for multiple matches)
                        sub_acqs = [os.path.basename(f).removeprefix('%s_%s_' % (ID, ses_id)).removesuffix('.gz').removesuffix('.nii') for f in glob(os.path.join(self.bids_database, ID, 'anat', self.acq_pattern + '.nii*'))]
                        self.add_ses_row(ID, row, '', {}, sub_acqs)
        # Update self.covariate_values with np matrix
        self.covariate_values = np.vstack(self.covariate_values)

    ############################
    #   DATA (PRE)PROCESSING   #
    ############################

    def set_proc_pipeline(self, metric_proc_pipeline):
        """
        Sets the pipeline used to process MRI scans and compute morphometric values for each participant. Should match
        `scanometrics/processing/<metric_proc_pipeline>.py`.
        """
        if hasattr(processing, metric_proc_pipeline):
            self.metric_proc_pipeline = getattr(processing, metric_proc_pipeline).proc_pipeline(self.bids_database, ses_delimiter=self.ses_delimiter, acq_delimiter=self.acq_delimiter)
        else:
            logging.ERROR("""Processing pipeline %s not found in scanometrics.processing module. Make sure the module file
exists, that it has been added to the __init__.py file, and that scanometrics is up-to-date"
with 'pip install -U .'""" % metric_proc_pipeline)

    def run_proc_pipeline(self, n_threads=-1):
        """Runs metric_proc_pipeline(), which generates measured metric values based on n_threads. Can be used to process
        normative data, or a set of subjects to evaluate against a trained dataset.
        Proc pipelines work in a dedicated 'derivatives' folder. In the case of Freesurfer, it expects a single
        directory with a single folder for each combination of [subject_id, session_id, acq_label]. Each scan in copied
        into a '<subject_id>_<session_id>_<acq_label>' folder. This implies that subj_id in FS has to be compared to
        <subject_id>_<session_id>_<acq_label> when checking if the correct subject is being loaded. Also means that BIDS
        naming should respect <bids_directory>/<subject_id>/<session_id>/anat/<subject_id>_<session_id>_<acq_label>T1w.nii.gz.
        acq_label can be used to give a subject specific T1 acq_label array. Regarding BIDS guidelines, the file
        participants.tsv MUST have a single line per subject (each subject has to appear only once). This implies that
        repeated scans must be tracked with a sessions.tsv file encoding the variables changing across sessions (eg age,
        scanner or sequence), overwritting values from participatns.tsv if needed. load_subjects() takes care of loading
        covariates from the participants.tsv, and loop through the sessions.tsv to add session specific covariates. This
        is currently achieved by having self.subject as a dictionary with subject_id as keys, and each value is another
        dict with session_id as keys, with a last dictionary with covariate names and values as key/value pairs."""
        # Call pipeline specific run_pipeline()
        self.metric_proc_pipeline.run_pipeline(self.subject, self.n_threads)

    def proc2table(self, n_threads=-1):
        """Method to convert pipeline specific outputs to a common table format. Has to be added to be able to import
        data processed externally but following the folder naming and organization implemented in the method
        self.metric_proc_pipeline.run_pipeline(). Externally processed data with different folder organization should be
        gathered into tables saved in <bids_database>/derivatives/<proc_pipeline>/<subject_id>_<session_id> folders
        by the user. Tables should match naming expected by self.metric_proc_pipeline.proc2metric(). In the case of
        freesurfer, this assumes that data to import is in a set of tables for each subject, located in
        <bids_directory>/derivatives/freesurfer/<subject_id>_<session_id> and that subject IDs in the generated tables
        will be <subject_id>_<session_id>. If user processed freesurfer externally without using this format, they
        should change participants.tsv and folder structure to match it.
        TODO: implement a warning instead of error when stumbling uppon an ID mismatch when loading tables.

        :param n_threads: number of threads to use
        :type n_threads: int
        """
        # Split subjSesAcq_list and run parallel processing
        self.metric_proc_pipeline.proc2table(self.subject, self.n_threads)

    def load_proc_metrics(self, subjects=None, stats2table_folder=None, ref_rows=None, ref_metric_values=None,
                          ref_metric_names=None, ref_covariate_values=None, ref_covariate_names=None, metric_include=None,
                          metric_exclude=None):
        """Load metrics computed by processing pipeline. We let the processing class to implement the specific loading
        of metrics according to it. Fills with Nans the values that don't exist for a given subject.
        Can be saved with save\_proc\_metrics() (eg a 75x1258 matrix requires 760 kB). Covariates is a list of variable
        names in participants.tsv to keep and save in a covariate_values numpy array.

        :param subjects: subject dictionary following the SOM.subject structure, which takes participants IDs from the
                         participants.tsv table and session IDs from sessions.tsv tables, according to the following
                         structure: {'<participant_id>': {'<session_id>': {covariate_name[0]: covariate_value[0], ...}}}
        :type subjects: dict
        :param stats2table_folder: path to folder containing stats2table files, to be loaded instead of subject specific
                                   files.
        :type stats2table_folder: string
        :param ref_rows: numpy array of indexes of rows in covariate_values to be passed to
                         metric_prac_pipeline.prac2metric(). Indexes are used to specify which subject and session
                         combinations should be used for reference during the normalization step.
        :type ref_rows: numpy array
        :param ref_metric_values: numpy array with values to use as reference for normalization of metrics.
        :type ref_metric_values: numpy array.
        :param ref_metric_names: list of metric names to use as reference.
        :type ref_metric_names: list
        :param ref_covariate_values: list of covariate values to find matching scans to normalize with.
        :type ref_covariate_values: list
        :param ref_covariate_names: list of covariate names to base search for matching scans on
        :type ref_covariate_names: list
        :param metric_include: list of metrics to keep after loading
        :type metric_include: list
        :param metric_exclude: list of metrics to exclude after loading
        :type metric_exclude: list
        """
        if subjects is None:
            subjects = self.subject.copy()
        self.metric_names, self.measured_metrics, self.metric_plotting_info = \
                self.metric_proc_pipeline.proc2metric(subjects, self.covariate_values, self.covariate_names,
                                                      stats2table_folder=stats2table_folder, ref_rows=ref_rows,
                                                      ref_metric_values=ref_metric_values,
                                                      ref_metric_names=ref_metric_names,
                                                      ref_covariate_values=ref_covariate_values,
                                                      ref_covariate_names=ref_covariate_names,
                                                      metric_include=metric_include, metric_exclude=metric_exclude)

    ###################
    # NORMATIVE MODEL #
    ###################
    def list_normative_models(self):
        """
        Lists available models in scanometrics/resources/normative_models. This is the default location for normative
        models distributed with ScanOMetrics.
        """
        models_path = os.path.join(os.path.dirname(__file__), 'resources',  'normative_models')
        normative_model_files = glob(os.path.join(models_path, '*.pkl'))
        logging.PRINT('Available models in %s:' % models_path)
        for f in normative_model_files:
            logging.PRINT('\t- %s' % (os.path.basename(f)))
        logging.PRINT("Eg: to load the 1st model, run `SOM.load_normative('%s')`" % (os.path.splitext(os.path.basename(normative_model_files[0]))[0]))


    def set_normative_model(self, model_name='Polynomial'):
        """
        Sets normative model based on model name (defaults to 'Polynomial') and a training set ID (defaults to dldirect_OASIS3).
        Intended to initialize normModel dictionary for further fitting or loading already trained model.

        :param model_name: string identifying the model name. Should correspond to a <model_name>.py file in the folder
                           scanometrics/processing
        :type model_name: string
        """
        if hasattr(normative_models, model_name):
            self.normativeModel = getattr(normative_models, model_name)(self.measured_metrics, self.metric_names,
                                  self.covariate_values, self.covariate_names, self.dataset_id,
                                  self.metric_proc_pipeline.proc_pipeline_name, self.cov2float, self.subject)
        else:
            logging.ERROR('Model %s not available in scanometrics.normative_models' % model_name)

    def save_normative_model(self, output_filename=None):
        """
        Saves normative model to pkl file. By default, will save a <model_name>_<proc_pipeline>_<bids_directory>.pkl file
        in the `scanometrics/resources/normative_models` folder (e.g. Polynomial_dldirect_OASIS3). Provide output_filename
        to save to another file.

        :param output_filename: path to pkl file (including filename and pkl extension)
        :type output_filename: string
        """
        if output_filename is None:
            output_filename = os.path.join(os.path.dirname(__file__), 'resources',  'normative_models', self.normativeModel.model_dataset_id+'.pkl')
        # Save pickle file
        with open(output_filename, 'wb') as fout:
            pickle.dump(self.normativeModel, fout, pickle.HIGHEST_PROTOCOL)

    def load_normative_model(self, model_dataset_id):
        """
        Loads normative model from pkl file into SOM.normativeModel structure, and overwrites SOM processing pipeline
        and cov2float to match those in the normative model.

        :param model_dataset_id: name of model to load. Can be one of the outputs of `list_normative_models()` without
                                 pkl extension (the file should be in `scanometrics/resources/normative_models` folder).
                                 Can also be the path to a pkl file (should contain the .pkl extension).
        :type model_dataset_id: string
        """
        if os.path.exists(model_dataset_id):
            input_file = model_dataset_id
        else:
            input_file = os.path.join(os.path.dirname(__file__), 'resources', 'normative_models', model_dataset_id+'.pkl')
            logging.PRINT("%s was not found to be an existing file, setting input_file to %s" % (model_dataset_id, input_file))
        if os.path.exists(input_file):
            with open(input_file, 'rb') as fin:
                self.normativeModel = pickle.load(fin)
            # Set current processing pipeline to the same used to train the normative dataset
            self.set_proc_pipeline(self.normativeModel.proc_pipeline_id)
            self.cov2float = self.normativeModel.cov2float.copy()
        else:
            logging.ERROR("Normative model file not found (%s). Please specify model_dataset_id without the '.pkl' at "
                          "the end of the file. For a list of available normative models, you can execute the "
                          "scanometrics.core.list_normative_models() function. You can also specify the whole path to a"
                          " pkl file, including the extension." % input_file)
        if self.normativeModel.som_version != __version__:
            logging.WARNING("SOM version used for training of %s model was %s, but local version of ScanOMetrics is %s."
                            "You might want to review changes in ScanOMetrics to check for discrepancies in processing"
                            "pipelines." % (self.normativeModel.model_dataset_id, self.normativeModel.som_version,
                                            __version__))

    #######################
    # ACCESSORY FUNCTIONS #
    #######################

    def get_subjSesAcq_array(self):
        """
        Get the array concatenating all <subj_ID>_<ses_ID>_<acq_label> strings in the loaded dataset.
        """
        subject_IDmergedSes = []
        for sub_id in self.subject.keys():
            for ses_id in self.subject[sub_id].keys():
                for acq_label in self.subject[sub_id][ses_id].keys():
                    subjSesAcq_id = self.get_subjSesAcq_id(sub_id, ses_id, acq_label)
                    subject_IDmergedSes.append(subjSesAcq_id)
        return subject_IDmergedSes

    def get_subjSesAcq_T1s(self):
        """
        Get list of paths to all T1 files in the loaded dataset.
        """
        subjSesAcqs_T1s = []
        for subj_id in self.subject.keys():
            for ses_id in self.subject[subj_id].keys():
                for acq_label in self.subject[subj_id][ses_id].keys():
                    subjSesAcq_id = self.get_subjSesAcq_id(subj_id, ses_id, acq_label)
                    subjSesAcqs_T1s.append(os.path.join(self.bids_database, subj_id, ses_id, 'anat', '%s.nii.gz' % subjSesAcq_id))
        return subjSesAcqs_T1s


    def get_subjSesAcq_row(self, subject_id, session_id, acq_label):
        """
        Get row index for particular subject, session and acq_label (returns a single value).
        """
        subject_IDmergedSes = self.get_subjSesAcq_array()
        subjSesAcq_id = self.get_subjSesAcq_id(subject_id, session_id, acq_label)
        return subject_IDmergedSes.index(subjSesAcq_id)

    def get_subjSesAcq_id(self, subject_id, session_id, acq_label):
        """
        Get row index for particular subject, session and acq_label (returns a single value).
        """
        subjSesAcq_id = subject_id
        if session_id != "":
            subjSesAcq_id += self.ses_delimiter+session_id
        if acq_label != "":
            subjSesAcq_id += self.acq_delimiter+acq_label
        return subjSesAcq_id

    def get_subj_rows(self, subject_id):
        """
        Get list of row indexes for a particular subject (returns a list for all subject sessions and acq labels).
        """
        subSes_rows = []
        subject_IDmergedSes = self.get_subjSesAcq_array()
        for session_id in self.subject[subject_id].keys():
            for acq_label in self.subject[subject_id][session_id].keys():
                subjSesAcq_id = self.get_subjSesAcq_id(subject_id, session_id, acq_label)
                subSes_rows.append(subject_IDmergedSes.index(subjSesAcq_id))
        return np.array(subSes_rows)

    ######################
    # RUN CORE FUNCTIONS #
    ######################

    def evaluate_singleSubject_allSes(self, subject_id, matching_covariates, min_num_ctrl=5, alpha_uncorr=0.01, q_fdr=0.05,
                                      create_html_report=True):
        """
        Evaluation method to assess whether a new scan significantly differs from normative ranges. Loops through all
        available sessions and scans available in self.subject[subject_id].
        
        :param subject_id: ID of the participant to evaluate.
        :type subject_id: string
        :param matching_covariates: array of covariate names to filter set of normative controls to evaluate against.
        :type matching_covariates: list
        :param min_num_ctrl: minimum number of normative subject matches for the evaluation to be compted. If there are less matches, the evaluation is ignored.
        :type min_num_ctrl: int
        :param alpha_uncorr: statistical significance level for uncorrected p-values.
        :type alpha_uncorr: float
        :param q_fdr: statistical significance level for fdr-corrected p-values.
        :type q_fdr: float
        :param create_html_report: flag to activate/deactivate generation of html report with figures (time and memory consuming).
        :type create_html_report: bool
        """

        # NB: there can be SOM instances with dissociated subjects between normativeModel/self datasets (most usual case
        # (eg LOOCV were SOM_all gets its normativeModel from fit on dataset excluding subject of interest, or general
        # case were normativeModel is trained on a dataset and the SOM instance has only new subjects). However, SOM
        # should handle the case were new session for a subject in the training dataset comes in... Does this affect the
        # implementation of predict_singleSubject_residuals ? Probably yes, as residuals computed on normativeModel trained
        # with subject included in the training dataset. User might need to retrain the model excluding the whole subject
        # in that case. As model training is faster now, this should not be an issue, just something to flag out or
        # automate, or print a warning if the subject_id is found in normativeModel.subject for example.
        # This case-scenario doesn't affect the evaluate_singleSubject_allSes() function (this current function here) as
        # matching rows are selected here and subject can be discarded.

        normModel_matching_cols = np.array([i for i in range(len(self.normativeModel.covariate_names)) if
                                            self.normativeModel.covariate_names[i] in matching_covariates])

        # Tricky thing to do here: would be good to check if subject_id is in normativeModel.subject, as training model
        # might have been fit with a dataset that contains the subject (eg the subject came back for a scan and the user
        # would like to evaluate the new scan against the fit). However, depending on how things are done, a model might
        # be trained on a dataset with too generic IDs like sub-001... For now, printing a large warning message should
        # be enough.
        if subject_id in self.normativeModel.subject.keys():
            print("\n\n###########\n# WARNING #\n###########\n Subject %s being evaluated was found in"
                  " SOM.normativeModel.subject.keys(). Could be a coincidence for generic subject IDs (e.g. `sub-001`),"
                  " or you are evaluating a new scan for a subject that was included in the normative set. If this was"
                  " not on purpose, consider retraining the normativeModel on a dataset that does not contain the"
                  " subject.\n" % subject_id)

        # Extract current subject measured metrics and covariate values (all sessions as rows)
        subj_rows = self.get_subj_rows(subject_id)
        subj_measured_metrics = {}
        for k in self.measured_metrics.keys():
            subj_measured_metrics[k] = self.measured_metrics[k][subj_rows, :]
            if len(subj_measured_metrics[k].shape) == 1:
                subj_measured_metrics[k] = subj_measured_metrics[k][None, :]
        subj_covariate_values = self.covariate_values[subj_rows, :]
        if len(subj_covariate_values.shape) == 1:
            subj_covariate_values = subj_covariate_values[None, :]

        # Set output folder to session id or 'longitudinal' if using several timepoints
        if len(subj_rows) == 1:
            out_session_id = list(self.subject[subject_id].keys())[0]
        else:
            out_session_id = 'longitudinal_analysis'

        output = {}
        for k in self.normativeModel.measured_metrics.keys():
            # Get column matching between self.measured_metrics and self.normativeModel.measured_metrics, using metric_names
            # Intended to be used with arrays coming from new measurement, to match with ordering in normativeModel
            # eg subj_measured_metrics[:, normModel_cols[m]] - self.normativeModel.uncertainty[m]
            # All arrays intended to analyse the current subject should have same number of columns as normativeModel
            col_idxs_normativeModel = []
            col_idxs_self = []
            metric_names_union = [metric_name for metric_name in self.metric_names[:self.measured_metrics[k].shape[1]]
                                  if metric_name in self.normativeModel.metric_names[:self.normativeModel.measured_metrics[k].shape[1]]]
            for metric_name in metric_names_union:
                col_idxs_normativeModel.append(self.normativeModel.metric_names.index(metric_name))
                col_idxs_self.append(self.metric_names.index(metric_name))
            # Define outputs, should be n_session x len(self.normativeModel.measured_metrics)
            subj_ages = subj_covariate_values[:, self.covariate_names.index('age')]
            output[k] = {}
            output[k]['residuals'] = np.full((len(subj_ages), len(col_idxs_self)), np.nan)
            output[k]['devtn_dfs'] = np.ones((len(subj_ages), len(col_idxs_self)))
            output[k]['devtn_stats'] = np.full((len(subj_ages), len(col_idxs_self)), np.nan)
            output[k]['devtn_pvals'] = np.full((len(subj_ages), len(col_idxs_self)), np.nan)
            output[k]['devtn_logps'] = np.full((len(subj_ages), len(col_idxs_self)), np.nan)
            output[k]['devtall_dfs'] = np.ones(len(col_idxs_self))
            output[k]['devtall_stats'] = np.full(len(col_idxs_self), np.nan)
            output[k]['devtall_pvals'] = np.full(len(col_idxs_self), np.nan)
            output[k]['devtall_logps'] = np.full(len(col_idxs_self), np.nan)
            output[k]['trend_coefs'] = np.full((len(col_idxs_self), 2), np.nan)
            output[k]['trend_dfs'] = np.ones(len(col_idxs_self))
            output[k]['trend_stats'] = np.full(len(col_idxs_self), np.nan)
            output[k]['trend_pvals'] = np.full(len(col_idxs_self), np.nan)
            output[k]['trend_logps'] = np.full(len(col_idxs_self), np.nan)
            output[k]['metric_names'] = [self.metric_names[i] for i in col_idxs_self]

            # Compute subject residuals and residuals of matching subjects in the normativeModel
            matching_rows = []  # use list with append to add matches foreach session, and reduce with unique after
            predicted_values = self.normativeModel.predict_values(subj_covariate_values)
            i_sesAcq = 0
            for session_id in self.subject[subject_id].keys():
                for acq_id in self.subject[subject_id][session_id].keys():
                    # Extract matching subjSesAcq row indexes from normativeModel covariate_values
                    session_matches = np.where((self.normativeModel.covariate_values[:, normModel_matching_cols] ==
                                                   subj_covariate_values[i_sesAcq, normModel_matching_cols]).all(axis=1))[0]
                    output[k]['residuals'][i_sesAcq, :] = subj_measured_metrics[k][i_sesAcq, col_idxs_self] - \
                                                          predicted_values[k][i_sesAcq, col_idxs_normativeModel]
                    matching_residuals = self.normativeModel.fit_outputs[k]['residuals'][session_matches, :][:, col_idxs_normativeModel].copy()
                    outliers = self.normativeModel.outliers[k][session_matches, :][:, col_idxs_normativeModel].copy()
                    matching_residuals[outliers == 1] = np.nan
                    if len(matching_residuals.shape) == 1:
                        matching_residuals = matching_residuals[None, :]  # Make the array 2D to average over 1st dimension later
                    # Compute devtns
                    z = output[k]['residuals'][i_sesAcq, :] - np.nanmean(matching_residuals, axis=0)
                    z /= (np.sqrt(self.normativeModel.uncertainty[k][col_idxs_normativeModel] ** 2 + np.nanvar(matching_residuals, ddof=1, axis=0)))  # z-value for single subject/session against variance of matching residuals combined with uncertainties
                    p = 1 - norm.cdf(np.abs(z)) + norm.cdf(-np.abs(z))
                    output[k]['devtn_stats'][i_sesAcq, :] = z.copy()
                    output[k]['devtn_pvals'][i_sesAcq, :] = p.copy()
                    output[k]['devtn_logps'][i_sesAcq, :] = -np.sign(z) * np.log10(p)
                    matching_rows.append(session_matches)
                    i_sesAcq += 1
            # Compute devtalls
            if len(self.subject[subject_id]) > 1:
                matching_rows = np.unique(np.hstack(matching_rows))  # Group all sess matches and remove duplicates
                matching_residuals = self.normativeModel.fit_outputs[k]['residuals'][matching_rows, :][:, col_idxs_normativeModel].copy()
                outliers = self.normativeModel.outliers[k][matching_rows, :][:, col_idxs_normativeModel].copy()
                matching_residuals[outliers] = np.nan  # sets outliers to nan, so matching rows still includes outlier subjects, but not taken into account when using nan_policy='omit'
                if len(matching_residuals.shape) == 1:
                    matching_residuals = matching_residuals[None, :]  # Make the array 2D to average over 1st dimension later
                for m in range(len(col_idxs_self)):
                    if (~np.isnan(matching_residuals[:, m])).sum() >= min_num_ctrl:
                        t, p = ttest_ind(output[k]['residuals'][:, m], matching_residuals[:, m], nan_policy='omit')
                        output[k]['devtall_stats'][m] = t
                        output[k]['devtall_pvals'][m] = p
                        output[k]['devtall_logps'][m] = -np.sign(t) * np.log10(p)
                        # Line fit to measured_metric
                        if len(np.unique(subj_ages)) > 1:
                            c_metric = poly_model.fit(subj_ages, subj_measured_metrics[k][:, col_idxs_self[m]], deg=1)
                            res = subj_measured_metrics[k][:, col_idxs_self[m]] - c_metric(subj_ages)
                            T_fit = (len(res)-1)*np.var(res, ddof=1)/self.normativeModel.uncertainty[k][col_idxs_normativeModel[m]]**2
                            p_fit = 1 - chi2.cdf(T_fit, len(res) - 1)
                            if p_fit > alpha_uncorr:
                                output[k]['trend_coefs'][m, :] = c_metric.convert().coef.copy()
                            # Line fit to residuals
                            c_res = poly_model.fit(subj_ages, output[k]['residuals'][:, m], deg=1)
                            # Residue over predicted residuals
                            res = output[k]['residuals'][:, m] - c_res(subj_ages)
                            T_fit = (len(res)-1)*np.var(res, ddof=1)/self.normativeModel.uncertainty[k][col_idxs_normativeModel[m]]**2
                            p_fit = 1-chi2.cdf(T_fit, len(res)-1)
                            if p_fit > alpha_uncorr:
                                slope_est = c_res.convert().coef[1] if len(c_res.convert().coef) > 1 else 0.0
                                SSE = res.dot(res.T)
                                slope_std = np.sqrt(SSE/(len(res)-2)) / np.std(subj_ages, ddof=1)/np.sqrt(len(res))
                                z = slope_est/slope_std
                                p = 1-norm.cdf(np.abs(z))+norm.cdf(-np.abs(z))
                                output[k]['trend_dfs'][m] = len(res)
                                output[k]['trend_stats'][m] = z
                                output[k]['trend_pvals'][m] = p
                                output[k]['trend_logps'][m] = -np.sign(z)*np.log10(p)
                                # Weighted odd ratios can be computed here, by filtering/grouping odd ratios weighted by devtn_logp across
                                # morphological metrics (eg volume, thickness, and surface area, etc...)
                                # ...

        all_pvals = np.hstack(([output[k]['devtn_logps'].reshape(-1) for k in output.keys()]+[output[k]['trend_logps'] for k in output.keys()]))
        pID, _ = fdr(10 ** (-np.abs(all_pvals)), q_fdr)

        if create_html_report:
            for k in self.normativeModel.measured_metrics.keys():
                col_idxs_normativeModel = []
                col_idxs_self = []
                metric_names_union = [metric_name for metric_name in self.metric_names[:self.measured_metrics[k].shape[1]]
                                      if metric_name in self.normativeModel.metric_names[:self.normativeModel.measured_metrics[k].shape[1]]]
                for metric_name in metric_names_union:
                    col_idxs_normativeModel.append(self.normativeModel.metric_names.index(metric_name))
                    col_idxs_self.append(self.metric_names.index(metric_name))
                subj_ages = subj_covariate_values[:, self.covariate_names.index('age')]
                matching_rows = []  # use list with append to add matches foreach session, and reduce with unique after
                i_sesAcq = 0
                for session_id in self.subject[subject_id].keys():
                    for acq_id in self.subject[subject_id][session_id].keys():

                        # Copy html template
                        if create_html_report:
                            html_output = os.path.join(self.bids_database, 'derivatives', 'scanometrics',
                                                       self.normativeModel.model_dataset_id, subject_id, session_id, acq_id,
                                                       'html')
                            if not os.path.exists(html_output):
                                copytree(os.path.join(os.path.dirname(__file__), 'resources', 'html_template'), html_output)

                        # Find scans in normative model that match session
                        session_matches = np.where((self.normativeModel.covariate_values[:, normModel_matching_cols] ==
                                                    subj_covariate_values[i_sesAcq, normModel_matching_cols]).all(axis=1))[0]
                        matching_rows.append(session_matches)
                        for i_metric in range(len(col_idxs_self)):
                            # NB: where to save plots ? Should be in bids/derivatives/scanometrics/<model_dataset_id>/
                            # <subject_id>/<session_id>/html/figs'. Eg for LOOCV loop, the model_id remains the same
                            # (eg Polynomial) and dataset_id is original dataset (eg CHFIRST) with '-LOOCV-001'
                            # (-> dataset_model_id = 'polynomial_CHFIRST-LOOCV-001

                            bad = np.argwhere(self.normativeModel.outliers[k][:, col_idxs_normativeModel[i_metric]] == 1)
                            good = np.argwhere(self.normativeModel.outliers[k][:, col_idxs_normativeModel[i_metric]] == 0)
                            match_bad = np.intersect1d(session_matches, bad)
                            match_good = np.intersect1d(session_matches, good)
                            nonmatch_bad = np.intersect1d(np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), session_matches), bad)
                            nonmatch_good = np.intersect1d(np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), session_matches), good)

                            normModel_age = self.normativeModel.covariate_values[:, self.normativeModel.covariate_names.index('age')].copy()
                            fig = plt.figure()
                            ax1 = fig.gca()
                            ax1.plot(self.normativeModel.age_vec,
                                     self.normativeModel.fit_outputs[k]['fit_ave'][:, col_idxs_normativeModel[i_metric]], 'k')
                            # Commented out as normative dataset is large enough not to need CI boundaries
                            # ax1.plot(self.normativeModel.age_vec,
                            #          self.normativeModel.fit_outputs[k]['fit_lrg'][:, col_idxs_normativeModel[i_metric]], 'k-', alpha=0.3)
                            # ax1.plot(self.normativeModel.age_vec,
                            #          self.normativeModel.fit_outputs[k]['fit_sml'][:, col_idxs_normativeModel[i_metric]], 'k-', alpha=0.3)
                            # Plot matching/nonmatching outliers and non-outliers points from normative dataset
                            ax1.plot(normModel_age[match_good], self.normativeModel.measured_metrics[k][match_good, col_idxs_normativeModel[i_metric]], 'ko',
                                     markersize=6, markerfacecolor='none')
                            ax1.plot(normModel_age[match_bad], self.normativeModel.measured_metrics[k][match_bad, col_idxs_normativeModel[i_metric]], 'kx',
                                     markersize=6)
                            ax1.plot(normModel_age[nonmatch_good], self.normativeModel.measured_metrics[k][nonmatch_good, col_idxs_normativeModel[i_metric]],
                                     'ko', markersize=3, markerfacecolor='none', alpha=0.3)
                            ax1.plot(normModel_age[nonmatch_bad], self.normativeModel.measured_metrics[k][nonmatch_bad, col_idxs_normativeModel[i_metric]], 'kx',
                                     markersize=3, alpha=0.3)
                            # Plot subject of interest
                            ax1.scatter(subj_ages[i_sesAcq], subj_measured_metrics[k][i_sesAcq, col_idxs_self[i_metric]], s=50, facecolors='b',
                                        edgecolors='b')
                            ax1.errorbar(subj_ages[i_sesAcq], subj_measured_metrics[k][i_sesAcq, col_idxs_self[i_metric]],
                                         yerr=self.normativeModel.uncertainty[k][col_idxs_normativeModel[i_metric]], ecolor='b', ls='None')
                            # Color background
                            if not np.isnan(output[k]['devtn_pvals'][i_sesAcq, i_metric]) and output[k]['devtn_pvals'][i_sesAcq, i_metric] < pID:
                                ax1.set_facecolor([1, 0.5, 0.5])
                                ax1.set_alpha(0.5)
                            elif not np.isnan(output[k]['devtn_pvals'][i_sesAcq, i_metric]) and output[k]['devtn_pvals'][i_sesAcq, i_metric] < alpha_uncorr:
                                ax1.set_facecolor('yellow')
                                ax1.set_alpha(0.5)
                            # Arrange display
                            ax1.grid(visible=True, linewidth=0.5, alpha=0.5)
                            ax1.set_xlabel('Age (years)')
                            ax1.set_ylabel(self.metric_names[col_idxs_self[i_metric]])
                            if 'symmetryIndex' in self.metric_names[col_idxs_self[i_metric]]:
                                ax1.set_ylim((-1, 1))
                            # Add text information
                            ylim = ax1.get_ylim()
                            yrange = ylim[1] - ylim[0]
                            xlim = ax1.get_xlim()
                            xtext = xlim[0] + 0.05 * (xlim[1] - xlim[0])
                            ax1.text(xtext, ylim[1] + 0.10 * yrange,
                                     r'%s=%1.2f$\pm$%1.2f%s' % (self.metric_plotting_info['type'][col_idxs_self[i_metric]],
                                                                subj_measured_metrics[k][i_sesAcq, col_idxs_self[i_metric]],
                                                                self.normativeModel.uncertainty[k][col_idxs_normativeModel[i_metric]],
                                                                self.metric_plotting_info['units'][col_idxs_self[i_metric]]))
                            closest_agevec = np.argmin(np.abs(self.normativeModel.age_vec - subj_ages[i_sesAcq]))
                            ax1.text(xtext, ylim[1] + 0.05 * yrange, r'Norm=%1.2f%s, CI=[%1.2f-%1.2f]' % (
                              self.normativeModel.fit_outputs[k]['fit_ave'][closest_agevec, col_idxs_normativeModel[i_metric]],
                              self.metric_plotting_info['units'][col_idxs_self[i_metric]],
                              self.normativeModel.fit_outputs[k]['fit_sml'][closest_agevec, col_idxs_normativeModel[i_metric]],
                              self.normativeModel.fit_outputs[k]['fit_lrg'][closest_agevec, col_idxs_normativeModel[i_metric]]))
                            ax1.text(xtext, ylim[0] - 0.01 * yrange,
                                     r'z$_{dev}$=%1.3f, p$_{dev}$=%1.3f' % (output[k]['devtn_stats'][i_sesAcq, i_metric],
                                                                            output[k]['devtn_pvals'][i_sesAcq, i_metric]))

                            ax1.text(xtext, ylim[0] - 0.10 * yrange,
                                     r'Odds=%1.3f, outlier_frac=%1.3f, P(art.)=%1.3f' % (self.normativeModel.stats[k]['odds'][col_idxs_normativeModel[i_metric]],
                                                                                         self.normativeModel.stats[k]['fout'][col_idxs_normativeModel[i_metric]],
                                                                                         self.normativeModel.stats[k]['part'][col_idxs_normativeModel[i_metric]]))

                            ax1.set_ylim(ylim[0] - 0.15 * yrange, ylim[1] + 0.15 * yrange)
                            output_folder = os.path.join(html_output, 'figs')
                            if not os.path.exists(output_folder):
                                os.makedirs(output_folder)
                            fig.savefig(os.path.join(output_folder, '%s_%s.jpeg' % (self.normativeModel.metric_names[col_idxs_normativeModel[i_metric]], k)),
                                        dpi=300)
                            plt.close()
                        i_sesAcq += 1
                if len(self.subject[subject_id]) > 1:
                    # Create longitudinal report
                    html_output = os.path.join(self.bids_database, 'derivatives', 'scanometrics',
                                               self.normativeModel.model_dataset_id, subject_id, 'longitudinal_analysis',
                                               'html')
                    if not os.path.exists(html_output):
                        copytree(os.path.join(os.path.dirname(__file__), 'resources', 'html_template'), html_output)
                    # Group all sess matches and remove duplicates
                    matching_rows = np.unique(np.hstack(matching_rows))
                    for i_metric in range(len(col_idxs_self)):
                        bad = np.argwhere(self.normativeModel.outliers[k][:, col_idxs_normativeModel[i_metric]] == 1)
                        good = np.argwhere(self.normativeModel.outliers[k][:, col_idxs_normativeModel[i_metric]] == 0)
                        match_bad = np.intersect1d(matching_rows, bad)
                        match_good = np.intersect1d(matching_rows, good)
                        nonmatch_bad = np.intersect1d(np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), matching_rows), bad)
                        nonmatch_good = np.intersect1d(np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), matching_rows), good)

                        normModel_age = self.normativeModel.covariate_values[:, self.normativeModel.covariate_names.index('age')].copy()
                        fig = plt.figure()
                        ax1 = fig.gca()
                        ax1.plot(self.normativeModel.age_vec,
                                 self.normativeModel.fit_outputs[k]['fit_ave'][:, col_idxs_normativeModel[i_metric]], 'k')
                        # Commented out as normative dataset is large enough not to need CI boundaries
                        # ax1.plot(self.normativeModel.age_vec,
                        #          self.normativeModel.fit_outputs[k]['fit_lrg'][:, col_idxs_normativeModel[i_metric]], 'k-', alpha=0.3)
                        # ax1.plot(self.normativeModel.age_vec,
                        #          self.normativeModel.fit_outputs[k]['fit_sml'][:, col_idxs_normativeModel[i_metric]], 'k-', alpha=0.3)
                        # Plot matching/nonmatching outliers and non-outliers points from normative dataset
                        ax1.plot(normModel_age[match_good], self.normativeModel.measured_metrics[k][match_good, col_idxs_normativeModel[i_metric]], 'ko',
                                 markersize=6, markerfacecolor='none', alpha=0.5)
                        ax1.plot(normModel_age[match_bad], self.normativeModel.measured_metrics[k][match_bad, col_idxs_normativeModel[i_metric]], 'kx',
                                 markersize=6, alpha=0.5)
                        ax1.plot(normModel_age[nonmatch_good], self.normativeModel.measured_metrics[k][nonmatch_good, col_idxs_normativeModel[i_metric]],
                                 'ko', markersize=3, markerfacecolor='none', alpha=0.2)
                        ax1.plot(normModel_age[nonmatch_bad], self.normativeModel.measured_metrics[k][nonmatch_bad, col_idxs_normativeModel[i_metric]],
                                 'kx',
                                 markersize=3, alpha=0.2)
                        # Plot subject of interest
                        ax1.scatter(subj_ages, subj_measured_metrics[k][:, col_idxs_self[i_metric]], s=50,
                                    facecolors='b',
                                    edgecolors='b')
                        ax1.errorbar(subj_ages, subj_measured_metrics[k][:, col_idxs_self[i_metric]],
                                     yerr=self.normativeModel.uncertainty[k][col_idxs_normativeModel[i_metric]], ecolor='b', ls='None')
                        # If there are more than 2 tp, plot estimated longitudinal line and color background accordingly
                        if len(self.subject[subject_id]) > 2 and ~np.isnan(output[k]['trend_coefs'][i_metric, :].sum()):
                            c = poly_model.basis(2)
                            c.coef = output[k]['trend_coefs'][i_metric, :].copy()
                            ax1.plot(self.normativeModel.age_vec, c(self.normativeModel.age_vec), 'b')
                        # Color background
                        if not np.isnan(output[k]['devtall_pvals'][i_metric]) and output[k]['devtall_pvals'][i_metric] < pID:
                            ax1.set_facecolor([1, 0.5, 0.5])
                            ax1.set_alpha(0.5)
                        elif not np.isnan(output[k]['devtall_pvals'][i_metric]) and output[k]['devtall_pvals'][i_metric] < alpha_uncorr:
                            ax1.set_facecolor('yellow')
                            ax1.set_alpha(0.5)
                        # Arrange display
                        ax1.grid()
                        ax1.set_xlabel('Age (years)')
                        ax1.set_ylabel(self.metric_names[col_idxs_self[i_metric]])
                        if 'symmetryIndex' in self.metric_names[col_idxs_self[i_metric]]:
                            ax1.set_ylim((-1, 1))
                        # Add text information
                        ylim = ax1.get_ylim()
                        yrange = ylim[1] - ylim[0]
                        xlim = ax1.get_xlim()
                        xtext = xlim[0] + 0.05 * (xlim[1] - xlim[0])
                        ax1.text(xtext, ylim[1] + 0.10 * yrange,
                                 r'mean %s=%1.2f$\pm$%1.2f%s' % (self.metric_plotting_info['type'][col_idxs_self[i_metric]],
                                                            np.nanmean(subj_measured_metrics[k][:, col_idxs_self[i_metric]]),
                                                            np.nanstd(subj_measured_metrics[k][:, col_idxs_self[i_metric]], ddof=1),
                                                            self.metric_plotting_info['units'][col_idxs_self[i_metric]]))
                        closest_agevec = np.argmin(np.abs(self.normativeModel.age_vec - np.nanmean(subj_ages)))
                        ax1.text(xtext, ylim[1] + 0.05 * yrange, r'Norm=%1.2f%s, CI=[%1.2f-%1.2f]' % (
                            self.normativeModel.fit_outputs[k]['fit_ave'][closest_agevec, col_idxs_normativeModel[i_metric]],
                            self.metric_plotting_info['units'][col_idxs_self[i_metric]],
                            self.normativeModel.fit_outputs[k]['fit_sml'][closest_agevec, col_idxs_normativeModel[i_metric]],
                            self.normativeModel.fit_outputs[k]['fit_lrg'][closest_agevec, col_idxs_normativeModel[i_metric]]))
                        ax1.text(xtext, ylim[0] - 0.01 * yrange,
                                 r'z$_{dev}$=%1.3f, p$_{dev}$=%1.3f' % (output[k]['devtall_stats'][i_metric],
                                                                        output[k]['devtall_pvals'][i_metric]))

                        if len(self.subject[subject_id]) > 1:
                            ax1.text(xtext, ylim[0] - 0.05 * yrange,
                                     r'z$_{trend}$=%1.3f, p$_{trend}$=%1.3f' % (output[k]['trend_stats'][i_metric],
                                                                                output[k]['trend_pvals'][i_metric]))
                        ax1.text(xtext, ylim[0] - 0.10 * yrange,
                                 r'Odds=%1.3f, outlier_frac=%1.3f, P(art.)=%1.3f' % (
                                 self.normativeModel.stats[k]['odds'][col_idxs_normativeModel[i_metric]],
                                 self.normativeModel.stats[k]['fout'][col_idxs_normativeModel[i_metric]],
                                 self.normativeModel.stats[k]['part'][col_idxs_normativeModel[i_metric]]))

                        ax1.set_ylim(ylim[0] - 0.15 * yrange, ylim[1] + 0.15 * yrange)
                        output_folder = os.path.join(html_output, 'figs')
                        if not os.path.exists(output_folder):
                            os.makedirs(output_folder)
                        fig.savefig(os.path.join(output_folder, '%s_%s.jpeg' % (self.normativeModel.metric_names[col_idxs_normativeModel[i_metric]], k)),
                                    dpi=300)
                        plt.close()

        return output, pID

    def test_group_differences(self, matching_covariates, metric_names=None, normalizations=None, group_label=None, group_covariate=None):
        """
        Perform statistical test (Student t-test) between group residuals and residuals of matching normative scans.
        Consider different group options:
        
            * Loaded subjects against normative dataset
            * Groups inside loaded subjects, labeled by a variable, intended to be a single label to mask out subjects
              before testing. This function allows loading a complete dataset, and test a certain group against a normative
              model, one group at a time (for several groups, test_group_differences() should be called once per group,
              changing the group_label value).

        :param matching_covariates: list of covariate names to use to filter matching scans in the normative dataset
        :type matching_covariates: list
        :param metric_names: list of metric names to test for differences between groups. Defaults to None, which results
                             in testing all metrics in self.normativeModel.metric_names.
        :type metric_names: list of strings
        :param normalizations: list of normalization types to analyse. Defaults to None, which results into analysing
                               both 'orig' and 'norm' datasets in self.measured_metrics.
        :type normalizations: list of strings
        :param group_label: label of the group to test agains the normative dataset. Should correspond to a value in
                            group_covariate.
        :type group_label: int
        :param group_covariate: individual scan group label. Used with the group_label parameter to filter out subjects
                                to be tested against the normative dataset.
        :type group_covariate: array of ints
        """
        normModel_matching_cols = np.array([i for i in range(len(self.normativeModel.covariate_names)) if
                                            self.normativeModel.covariate_names[i] in matching_covariates])
        if metric_names is None:
            metric_names = self.metric_names.copy()
        output = {}
        if normalizations is None:
            normalizations = list(self.measured_metrics.keys())
        # Determine rows from self.covariate_values and self.measured_values[k] to use
        if group_label is None and group_covariate is None:
            subj_rows = np.arange(self.covariate_values.shape[0])
        else:
            subj_rows = np.argwhere(group_covariate == group_label)[:, 0]
        # Extract covariates and reorder according to normativeModel covariates
        subj_covariate_values = np.zeros((len(subj_rows), len(self.normativeModel.covariate_names)))
        for i, covariate_name in enumerate(self.normativeModel.covariate_names):
            subj_covariate_values[:, i] = self.covariate_values[subj_rows, :][:, self.covariate_names.index(covariate_name)]
        # Predict values for current subjects
        predicted_values = self.normativeModel.predict_values(subj_covariate_values)
        for k in normalizations:
            # Get column matching between self.measured_metrics and self.normativeModel.measured_metrics, using metric_names
            # Intended to be used with arrays coming from new measurement, to match with ordering in normativeModel
            # eg subj_measured_metrics[:, normModel_cols[m]] - self.normativeModel.uncertainty[m]
            # All arrays intended to analyse the current subject should have same number of columns as normativeModel
            normModel_cols = np.full(self.normativeModel.measured_metrics[k].shape[1], np.nan)
            for m in range(self.normativeModel.measured_metrics[k].shape[1]):
                if self.normativeModel.metric_names[m] in metric_names:
                    normModel_cols[m] = metric_names.index(self.normativeModel.metric_names[m])
            # Allocate output
            output[k] = {}
            output[k]['gp_residuals'] = np.full((len(subj_rows), len(normModel_cols)), np.nan)
            output[k]['gp-vs-gp_dfs'] = np.full((len(normModel_cols)), np.nan)
            output[k]['gp-vs-gp_ts'] = np.full((len(normModel_cols)), np.nan)
            output[k]['gp-vs-gp_pvals'] = np.full((len(normModel_cols)), np.nan)
            output[k]['gp-vs-gp_logps'] = np.full((len(normModel_cols)), np.nan)
            output[k]['gp-vs-gp_cohen-d'] = np.full((len(normModel_cols)), np.nan)
            output[k]['metric_names'] = [metric_names[int(m)] for m in normModel_cols if not np.isnan(m)]
            # Extract scan metrics
            subj_measured_metrics = self.measured_metrics[k][subj_rows, :]
            if len(subj_measured_metrics.shape) == 1:
                subj_measured_metrics = subj_measured_metrics[None, :]
            # Find matching controls in normative dataset
            matching_rows = []
            # Residuals are computed separately using precomputed predicted_value - scan value
            for i_scan in range(len(subj_rows)):
                # Extract matching subjSesAcq row indexes from normativeModel covariate_values
                scan_matches = np.where((self.normativeModel.covariate_values[:, normModel_matching_cols] ==
                                               subj_covariate_values[i_scan, normModel_matching_cols]).all(axis=1))[0]
                matching_rows.append(scan_matches)
            matching_rows = np.unique(np.hstack(matching_rows))  # Group all scan matches and remove duplicates
            matching_residuals = self.normativeModel.fit_outputs[k]['residuals'][matching_rows, :].copy()
            outliers = self.normativeModel.outliers[k][matching_rows, :].copy()
            matching_residuals[outliers] = np.nan  # sets outliers to nan, so matching rows still includes outlier subjects, but not taken into account when using nan_policy='omit'
            if len(matching_residuals.shape) == 1:
                matching_residuals = matching_residuals[None, :]  # Make the array 2D to average over 1st dimension later
            # Test for statistical difference between group residuals and matching residuals
            for m in np.argwhere(~np.isnan(normModel_cols))[:, 0]:
                output[k]['gp_residuals'][:, m] = subj_measured_metrics[:, int(normModel_cols[m])] - predicted_values[k][:, m]
                subj_noNans = output[k]['gp_residuals'][~np.isnan(output[k]['gp_residuals'][:, m]), m]
                matching_noNans = matching_residuals[~np.isnan(matching_residuals[:, m]), m]
                t, p = ttest_ind(subj_noNans, matching_noNans)
                output[k]['gp-vs-gp_dfs'][m] = len(subj_noNans)+len(matching_noNans)-1
                output[k]['gp-vs-gp_ts'][m] = t
                output[k]['gp-vs-gp_pvals'][m] = p
                output[k]['gp-vs-gp_logps'][m] = -np.sign(t) * np.log10(p)
                n1 = len(matching_noNans)
                n2 = len(subj_noNans)
                pooled_std = np.sqrt(((n1-1)*matching_noNans.std(ddof=1)**2+(n2-1)*subj_noNans.std(ddof=1)**2)/(n1+n2-2))
                output[k]['gp-vs-gp_cohen-d'][m] = (np.nanmean(output[k]['gp_residuals'][:, m])-np.nanmean(matching_residuals[:, m]))/pooled_std  # d<0 == test_group < norm, d>0 == test_group > norm
        return output
