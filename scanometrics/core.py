"""
Core ScanOMetrics classes and methods
"""
import io
import pickle
import numpy as np
import os
from scanometrics import normative_models, processing, __version__
from scanometrics.utils import logging
from scanometrics.utils.zenodo_api import list_normative_models, download_file
from scanometrics.utils.stats import fdr
from csv import DictReader
from scipy.stats import norm, ttest_ind, chi2
from numpy.polynomial import Polynomial as poly_model
from PyQt5.QtWidgets import QMessageBox
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
        'F': 1, 'f': 1}}, acq_pattern='*T1w', ses_delimiter="_", acq_delimiter="_", n_threads=-1, atlas='DesikanKilliany'):
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
        self.set_proc_pipeline(proc_pipeline, atlas)
        # Reset normative model
        self.normativeModel = None
        # Set dataset_id
        if dataset_id is None:
            dataset_id = '_'.join([self.metric_proc_pipeline.proc_pipeline_name, os.path.basename(os.path.normpath(self.bids_database))])
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
                    if row[k] in self.cov2float[k].keys():
                        row[k] = self.cov2float[k][row[k]]
                    else:
                        logging.ERROR("%s not in cov2float[%s] for ID=%s, ses_id=%s" % (row[k], k, ID, ses_id))
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
                        logging.PRINT("Removing %s from 'covariate_names'" % k)
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
                    row = {k: row[k] for k in self.covariate_names if k in row.keys()}
                    # Try converting everything to float, ValueErrors might rise from categorical variables, try cov2float and error exit if doesn't work
                    for k in self.covariate_names:
                        try:
                            row[k] = np.float64(row[k])
                        except ValueError as ve:
                            if k in self.cov2float.keys() and row[k] in self.cov2float[k].keys():
                                row[k] = self.cov2float[k][row[k]]
                            else:
                                logging.ERROR('ValueError for key %s, value "%s" and scan %s_%s_%s (error message was "' % (k, row[k], subj_id, ses_id, acq_label) + str(
                                    ve) + '") when running load_subjects(). Try adding categorical to'
                                          'numerical mapping in the cov2float dictionary when creating the ScanOMetrics_project'
                                          'and run load_subjects() again. Check that input covariate files do not have "NA", '
                                          '"na", "N/A" etc... entries instead of "NaN", "nan", or "NAN" and replace them '
                                          'accordingly')
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
                # Fill covariate values according to ordering of scans inside self.subject
                for subj_id in self.subject.keys():
                    for ses_id in self.subject[subj_id].keys():
                        for acq_label in self.subject[subj_id][ses_id].keys():
                            self.covariate_values.append(np.array([self.subject[subj_id][ses_id][acq_label][k] for k in self.covariate_names], dtype='float')[None, :])
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
                                ses_row = {k: ses_row[k] for k in self.covariate_names if k in ses_row.keys()}
                                logging.PRINT('ses_id=%s, sub_ses_exclude=%s, sub_ses_include=%s' % (ses_id, sub_ses_exclude, sub_ses_include))
                                if ((sub_ses_exclude is not None) and ((ID+'_'+ses_id) in sub_ses_exclude))\
                                or ((sub_ses_include is not None) and ((ID+'_'+ses_id) not in sub_ses_include)):
                                    continue
                                logging.PRINT('Adding subject %s (session %s)' % (ID, ses_id))
                                sub_acqs_matches = [os.path.basename(f).removeprefix('%s_%s_' % (ID, ses_id)).removesuffix('.gz').removesuffix('.nii') for f in glob(os.path.join(self.bids_database, ID, ses_id, 'anat', self.acq_pattern + '.nii*'))]
                                sub_acqs = []
                                if sub_ses_acq_include is not None:
                                    for acq_id in sub_acqs_matches:
                                        if self.get_subjSesAcq_id(ID, ses_id, acq_id) in sub_ses_acq_include:
                                            sub_acqs.append(acq_id)
                                else:
                                    sub_acqs = sub_acqs_matches.copy()
                                    if sub_ses_acq_exclude is not None:
                                        for acq_id in sub_acqs_matches:
                                            if self.get_subjSesAcq_id(ID, ses_id, acq_id) in sub_ses_acq_exclude:
                                                sub_acqs.remove(acq_id)
                                if len(sub_acqs) > 0:
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
        return len(self.subject)

    ############################
    #   DATA (PRE)PROCESSING   #
    ############################

    def set_proc_pipeline(self, metric_proc_pipeline, atlas):
        """
        Sets the pipeline used to process MRI scans and compute morphometric values for each participant. Should match
        `scanometrics/processing/<metric_proc_pipeline>.py`.
        """
        if hasattr(processing, metric_proc_pipeline):
            self.metric_proc_pipeline = getattr(processing, metric_proc_pipeline).proc_pipeline(self.bids_database, ses_delimiter=self.ses_delimiter, acq_delimiter=self.acq_delimiter, atlas=atlas)
        else:
            logging.ERROR("""Processing pipeline %s not found in scanometrics.processing module. Make sure the module file
exists, that it has been added to the __init__.py file, and that scanometrics is up-to-date"
with 'pip install -U .'""" % metric_proc_pipeline)

    def run_proc_pipeline(self, subject_id=None, n_threads=None):
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
        if n_threads is None:
            n_threads = -1
        logging.PRINT("scanometrics.core: n_threads set to %d" % n_threads)
        self.metric_proc_pipeline.run_pipeline(self.subject, n_threads, subject_id)

    def proc2table(self, n_threads=-1):
        """Method to convert pipeline specific outputs to a common table format. Has to be added to be able to import
        data processed externally but following the folder naming and organization implemented in the method
        self.metric_proc_pipeline.run_pipeline(). Externally processed data with different folder organization should be
        gathered into tables saved in <bids_database>/derivatives/<proc_pipeline>/<subject_id>_<session_id> folders
        by the user. Tables should match naming expected by self.metric_proc_pipeline.proc2metric(). In the case of
        freesurfer, this assumes that data to import is in a set of tables for each subject, located in
        <bids_directory>/derivatives/freesurfer_vX-X-X/<subject_id>_<session_id> and that subject IDs in the generated tables
        will be <subject_id>_<session_id>. If user processed freesurfer externally without using this format, they
        should change participants.tsv and folder structure to match it. Freesurfer's version should be specified by
        replacing vX-X-X with the appropriate value (done automatically when running processing from ScanOMetrics).
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
        Can be saved with save_proc_metrics() (eg a 75x1258 matrix requires 760 kB). Covariates is a list of variable
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
                                  self.metric_proc_pipeline.proc_pipeline_name, self.metric_proc_pipeline.version, self.cov2float, self.subject)
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
        if not os.path.exists(os.path.dirname(output_filename)):
            os.makedirs(os.path.dirname(output_filename))
        with open(output_filename, 'wb') as fout:
            pickle.dump(self.normativeModel, fout, pickle.HIGHEST_PROTOCOL)

    def load_normative_model(self, model_filename):
        """
        Loads normative model from pkl file into SOM.normativeModel structure, and overwrites SOM processing pipeline
        and cov2float to match those in the normative model.

        :param model_filename: filename of model to load. Can be one of the outputs of `list_normative_models()`.
                               Can also be the path to a pkl file (should contain the .pkl extension).
        :type model_filename: string
        """
        # if file exists locally, load it
        if os.path.exists(model_filename):
            logging.PRINT("%s found, loading it as normativeModel" % model_filename)
            with open(model_filename, 'rb') as fin:
                normative_model = pickle.load(fin)
        else:  # Check if model_filename is in SOM normative models folder, or on Zenodo
            # Retrieve Zenodo's and local normative models
            normative_models = list_normative_models()
            if model_filename in normative_models:
                if os.path.exists(normative_models[model_filename]['local_file']):
                    logging.PRINT("%s found in SOM normative models folder, loading it as normativeModel" % model_filename)
                else:
                    logging.PRINT("%s not found in SOM normative models folder, downloading it from Zenodo" % model_filename)
                    download_file(normative_models[model_filename]['download_url'], normative_models[model_filename]['local_file'])
                with open(normative_models[model_filename]['local_file'], 'rb') as fin:
                    normative_model = pickle.load(fin)
                if normative_model.proc_pipeline_version != normative_models[model_filename]['proc_pipeline_version']:
                    logging.ERROR('Mismatch in processing pipeline versions between downloaded file and the one queries from Zenodo, aborting')
            else:
                logging.ERROR("%s not found locally, nor in SOM normative models folders, nor on Zenodo, aborting..." % model_filename)
        self.normativeModel = normative_model
        proc_pipeline_id, _, atlas = self.normativeModel.proc_pipeline_id.split("_")
        logging.PRINT("\n\nAtlas set to %s" % atlas)
        # Set current processing pipeline to the same used to train the normative dataset
        self.set_proc_pipeline(proc_pipeline_id, atlas)  # atlas retrieved from last piece of self.normativeModel.proc_pipeline_id
        if self.metric_proc_pipeline.version != self.normativeModel.proc_pipeline_version:
            logging.WARNING("Processing pipeline in %s is %s but version installed on the system is %s. We recommend to process data with the same software version, please install version %s on your system before processing any MRIs on this system." % (model_filename, self.normativeModel.proc_pipeline_version, self.metric_proc_pipeline.version, self.normativeModel.proc_pipeline_version))
            self.metric_proc_pipeline.update_version(self.normativeModel.proc_pipeline_version, atlas)
            logging.WARNING("Processing version set to %s to ensure compatibility. Make sure the derivatives folder has a subfolder named '%s' before loading metrics." % (self.metric_proc_pipeline.version, self.metric_proc_pipeline.subjects_dir))

        self.cov2float = self.normativeModel.cov2float.copy()
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

    def evaluate_singleSubject_allSes(self, subject_id, matching_covariates, min_num_ctrl=5, alpha_uncorr=0.01, q_fdr=0.05):
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

        # set to order both subject and normativeModel matching covariate columns to match
        normModel_matching_cols = np.array([self.normativeModel.covariate_names.index(matching_name) for matching_name in matching_covariates])
        self_matching_cols = np.array([self.covariate_names.index(matching_name) for matching_name in matching_covariates])

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
            predicted_values = self.normativeModel.predict_values(subj_covariate_values[:, self.covariate_names.index('age')])
            i_sesAcq = 0
            for session_id in self.subject[subject_id].keys():
                for acq_id in self.subject[subject_id][session_id].keys():
                    # Extract matching subjSesAcq row indexes from normativeModel covariate_values
                    session_matches = np.where((self.normativeModel.covariate_values[:, normModel_matching_cols] ==
                                                   subj_covariate_values[i_sesAcq, self_matching_cols]).all(axis=1))[0]
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

        return subj_covariate_values, normModel_matching_cols, subj_measured_metrics, output, pID

    def plot_single_metric(self, subject_id, selected_metric, selected_session, selected_acquisition, output,
                           matching_covariates, pid, alpha_uncorr=0.01):
        """
        Plot the selected metric for the given subject, session, and acquisition, ensuring that both hemispheres are plotted together if available, and symmetry index in the middle.

        Parameters:
        subject_id (str): The ID of the subject.
        selected_metric (str): The metric to be plotted (eg aparc_lh_entorhinal_thickness).
        selected_session (str): The session of the subject.
        selected_acquisition (str): The acquisition label.
        output (dict): The output data containing metric names and deviation statistics.
        subj_covariate_values (ndarray): Covariate values for the subject.
        normModel_matching_cols (list): Columns for matching the normative model.
        subj_measured_metrics (ndarray): Measured metrics for the subject.
        pid (int): FDR threshold.
        alpha_uncorr (float, optional): Uncorrected alpha level. Default is 0.01.

        Returns:
        Figure: A matplotlib figure containing the plot.
        """
        plt.close("all")
        plt.style.use('classic')

        # Determine which data set to use based on the metric type
        k = 'norm'

        normModel_matching_cols = np.array([self.normativeModel.covariate_names.index(matching_name) for matching_name in matching_covariates])
        self_matching_cols = np.array([self.covariate_names.index(matching_name) for matching_name in matching_covariates])

        if selected_metric not in self.metric_names:
            return

        opposite_metric = None
        symmetry_metric = None

        # Determine the opposite metric for the other hemisphere
        if 'lh' in selected_metric:
            opposite_metric = selected_metric.replace('lh', 'rh')
            symmetry_metric = selected_metric.replace('lh', 'symmetryIndex')
        elif 'rh' in selected_metric:
            opposite_metric = selected_metric.replace('rh', 'lh').replace('entolhinal', 'entorhinal')
            symmetry_metric = selected_metric.replace('rh', 'symmetryIndex').replace('entolhinal', 'entorhinal')
        elif 'Left' in selected_metric:
            opposite_metric = selected_metric.replace('Left', 'Right')
            symmetry_metric = selected_metric.replace('Left', 'symmetryIndex')
        elif 'Right' in selected_metric:
            opposite_metric = selected_metric.replace('Right', 'Left')
            symmetry_metric = selected_metric.replace('Right', 'symmetryIndex')

        if (opposite_metric and opposite_metric not in self.metric_names) or (symmetry_metric and symmetry_metric not in self.metric_names):
            opposite_metric = None
            symmetry_metric = None

        # find index of metric in normativeModel metric_names
        col_idx_normativeModel = self.normativeModel.metric_names.index(selected_metric)
        # Find index of metric in self.metric_names
        col_idx_self = self.metric_names.index(selected_metric)
        # Find index of metric in output
        col_idx_output = output[k]["metric_names"].index(selected_metric)
        # Extract age of subject of interest
        subj_covariate_values = self.covariate_values[self.get_subjSesAcq_row(subject_id, selected_session, selected_acquisition), :]
        subj_ages = subj_covariate_values[self.covariate_names.index('age')]

        col_idx_normativeModel_opposite = None
        col_idx_self_opposite = None
        col_idx_output_opposite = None
        sym_idx_normativeModel = None
        sym_idx_self = None
        sym_idx_output = None

        if opposite_metric:
            col_idx_normativeModel_opposite = self.normativeModel.metric_names.index(opposite_metric)
            col_idx_self_opposite = self.metric_names.index(opposite_metric)
            col_idx_output_opposite = output[k]["metric_names"].index(opposite_metric)
            sym_idx_normativeModel = self.normativeModel.metric_names.index(symmetry_metric)
            sym_idx_self = self.metric_names.index(symmetry_metric)
            sym_idx_output = output['orig']["metric_names"].index(symmetry_metric)

        if selected_session not in self.subject[subject_id] or selected_acquisition not in self.subject[subject_id][
            selected_session]:
            logging.PRINT(
                f"Session {selected_session} or acquisition {selected_acquisition} not found for subject {subject_id}.")
            return

        # Extract subjSesAcq value
        subj_measured_metric = self.measured_metrics[k][self.get_subjSesAcq_row(subject_id, selected_session, selected_acquisition), col_idx_self]

        # Check if the metric has a side indication
        has_side = any(side in selected_metric for side in ['lh_', 'rh_', 'Left-', 'Right-']) or 'CortexVol' in selected_metric

        fig, ax = plt.subplots(nrows=1, ncols=3 if opposite_metric and has_side else 1,
                               figsize=(30 if opposite_metric and has_side else 10, 5))

        # Ensure 'Left' or 'lh' is on the left and 'Right' or 'rh' is on the right
        if has_side and opposite_metric and ('Right-' in selected_metric or 'rh' in selected_metric):
            ax = ax[::-1]  # Reverse the order of the axes

        run_idx = np.argwhere(self.get_subj_rows(subject_id) == self.get_subjSesAcq_row(subject_id, selected_session, selected_acquisition))[:, 0]
        logging.PRINT("Subject %s has %d scans (selected scan is scan number %d" % (subject_id, len(self.get_subj_rows(subject_id)), run_idx))
        for scan_id in [self.get_subjSesAcq_array()[i] for i in self.get_subj_rows(subject_id)]:
            logging.PRINT("    - %s" % scan_id)

        try:
            # Find matching and non-matching sessions
            session_matches = np.where((self.normativeModel.covariate_values[:, normModel_matching_cols] == subj_covariate_values[self_matching_cols]).all(axis=1))[0]
            bad = np.argwhere(self.normativeModel.outliers[k][:, col_idx_normativeModel] == 1)
            good = np.argwhere(self.normativeModel.outliers[k][:, col_idx_normativeModel] == 0)
            match_bad = np.intersect1d(session_matches, bad)
            match_good = np.intersect1d(session_matches, good)
            nonmatch_bad = np.intersect1d(
                np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), session_matches), bad)
            nonmatch_good = np.intersect1d(
                np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), session_matches), good)

            normModel_age = self.normativeModel.covariate_values[:,
                            self.normativeModel.covariate_names.index('age')].copy()

            # Plot for the selected metric
            ax1 = ax if not (opposite_metric and has_side) else ax[0]
            ax1.plot(self.normativeModel.age_vec,
                     self.normativeModel.fit_outputs[k]['fit_ave'][:, col_idx_normativeModel], 'k')
            ax1.plot(normModel_age[match_good],
                     self.normativeModel.measured_metrics[k][match_good, col_idx_normativeModel], 'ko',
                     markersize=6, markerfacecolor='none')
            ax1.plot(normModel_age[match_bad],
                     self.normativeModel.measured_metrics[k][match_bad, col_idx_normativeModel], 'kx',
                     markersize=6)
            ax1.plot(normModel_age[nonmatch_good],
                     self.normativeModel.measured_metrics[k][nonmatch_good, col_idx_normativeModel], 'ko',
                     markersize=3, markerfacecolor='none', alpha=0.3)
            ax1.plot(normModel_age[nonmatch_bad],
                     self.normativeModel.measured_metrics[k][nonmatch_bad, col_idx_normativeModel], 'kx',
                     markersize=3, alpha=0.3)
            ax1.scatter(subj_ages, subj_measured_metric, s=50,
                        facecolors='b', edgecolors='b')
            ax1.errorbar(subj_ages, subj_measured_metric,
                         yerr=self.normativeModel.uncertainty[k][col_idx_normativeModel], ecolor='b',
                         ls='None')

            pval = output[k]['devtn_pvals'][run_idx, col_idx_output]
            logging.PRINT(f"p-value: {pval}, pid: {pid}, alpha_uncorr: {alpha_uncorr}")
            if pval < pid:
                ax1.set_facecolor([1, 0.5, 0.5])
                ax1.set_alpha(0.5)
            elif pval < alpha_uncorr:
                ax1.set_facecolor('yellow')
                ax1.set_alpha(0.5)

            ax1.grid(visible=True, linewidth=0.5, alpha=0.5)
            ax1.set_xlabel('Age (years)')
            ax1.set_ylabel(self.metric_names[col_idx_self])
            ax1.set_title(f'Session: {selected_session}, Acquisition: {selected_acquisition}')
            if 'symmetryIndex' in self.metric_names[col_idx_self]:
                ax1.set_ylim((-1, 1))

            ylim = ax1.get_ylim()
            yrange = ylim[1] - ylim[0]
            xlim = ax1.get_xlim()
            xtext = xlim[0] + 0.05 * (xlim[1] - xlim[0])
            ax1.text(xtext, ylim[1] + 0.10 * yrange,
                     r'%s=%1.2f$\pm$%1.2f%s' % (self.metric_plotting_info['type'][col_idx_self],
                                                subj_measured_metric,
                                                self.normativeModel.uncertainty[k][col_idx_normativeModel],
                                                self.metric_plotting_info['units'][col_idx_self]))
            closest_agevec = np.argmin(np.abs(self.normativeModel.age_vec - subj_ages))
            ax1.text(xtext, ylim[1] + 0.05 * yrange, r'Norm=%1.2f%s, CI=[%1.2f-%1.2f]' % (
                self.normativeModel.fit_outputs[k]['fit_ave'][closest_agevec, col_idx_normativeModel],
                self.metric_plotting_info['units'][col_idx_self],
                self.normativeModel.fit_outputs[k]['fit_sml'][closest_agevec, col_idx_normativeModel],
                self.normativeModel.fit_outputs[k]['fit_lrg'][closest_agevec, col_idx_normativeModel]))
            ax1.text(xtext, ylim[0] - 0.01 * yrange,
                     r'z$_{dev}$=%1.3f, p$_{dev}$=%1.3f' % (
                         output[k]['devtn_stats'][run_idx, col_idx_output],
                         output[k]['devtn_pvals'][run_idx, col_idx_output]))
            ax1.text(xtext, ylim[0] - 0.10 * yrange,
                     r'Odds=%1.3f, outlier_frac=%1.3f, P(art.)=%1.3f' % (
                         self.normativeModel.stats[k]['odds'][col_idx_normativeModel],
                         self.normativeModel.stats[k]['fout'][col_idx_normativeModel],
                         self.normativeModel.stats[k]['part'][col_idx_normativeModel]))
            ax1.set_ylim(ylim[0] - 0.15 * yrange, ylim[1] + 0.15 * yrange)

            # Plot for the opposite hemisphere and symmetry index if available
            if opposite_metric and has_side:
                subj_opposite_metric = self.measured_metrics[k][self.get_subjSesAcq_row(subject_id, selected_session, selected_acquisition), col_idx_self_opposite]
                subj_sym_metric = self.measured_metrics['orig'][self.get_subjSesAcq_row(subject_id, selected_session, selected_acquisition), sym_idx_self]

                ax2 = ax[1]
                ax3 = ax[2]

                bad_sym = np.argwhere(self.normativeModel.outliers['orig'][:, sym_idx_normativeModel] == 1)
                good_sym = np.argwhere(self.normativeModel.outliers['orig'][:, sym_idx_normativeModel] == 0)
                match_bad_sym = np.intersect1d(session_matches, bad_sym)
                match_good_sym = np.intersect1d(session_matches, good_sym)
                nonmatch_bad_sym = np.intersect1d(
                    np.delete(np.arange(self.normativeModel.outliers['orig'].shape[0]), session_matches),
                    bad_sym)
                nonmatch_good_sym = np.intersect1d(
                    np.delete(np.arange(self.normativeModel.outliers['orig'].shape[0]), session_matches),
                    good_sym)

                # Plot symmetry index
                ax2.plot(self.normativeModel.age_vec,
                         self.normativeModel.fit_outputs['orig']['fit_ave'][:, sym_idx_normativeModel], 'k')
                ax2.plot(normModel_age[match_good_sym],
                         self.normativeModel.measured_metrics['orig'][
                             match_good_sym, sym_idx_normativeModel], 'ko',
                         markersize=6, markerfacecolor='none')
                ax2.plot(normModel_age[match_bad_sym],
                         self.normativeModel.measured_metrics['orig'][
                             match_bad_sym, sym_idx_normativeModel], 'kx',
                         markersize=6)
                ax2.plot(normModel_age[nonmatch_good_sym],
                         self.normativeModel.measured_metrics['orig'][
                             nonmatch_good_sym, sym_idx_normativeModel], 'ko',
                         markersize=3, markerfacecolor='none', alpha=0.3)
                ax2.plot(normModel_age[nonmatch_bad_sym],
                         self.normativeModel.measured_metrics['orig'][
                             nonmatch_bad_sym, sym_idx_normativeModel], 'kx',
                         markersize=3, alpha=0.3)

                ax2.scatter(subj_ages, subj_sym_metric, s=50,
                            facecolors='b', edgecolors='b')
                ax2.errorbar(subj_ages, subj_sym_metric,
                             yerr=self.normativeModel.uncertainty['orig'][sym_idx_normativeModel],
                             ecolor='b',
                             ls='None')

                pval = output['orig']['devtn_pvals'][run_idx, sym_idx_output]
                logging.PRINT(f"Symmetry index - p-value: {pval}, pid: {pid}, alpha_uncorr: {alpha_uncorr}")
                if pval < pid:
                    ax2.set_facecolor([1, 0.5, 0.5])
                    ax2.set_alpha(0.5)
                elif pval < alpha_uncorr:
                    ax2.set_facecolor('yellow')
                    ax2.set_alpha(0.5)

                ax2.grid(visible=True, linewidth=0.5, alpha=0.5)
                ax2.set_xlabel('Age (years)')
                ax2.set_ylabel(self.metric_names[sym_idx_self])

                ylim = [-1, 1]
                yrange = ylim[1] - ylim[0]
                xlim = ax2.get_xlim()
                xtext = xlim[0] + 0.05 * (xlim[1] - xlim[0])

                ax2.text(xtext, ylim[1] + 0.10 * yrange,
                         r'%s=%1.2f$\pm$%1.2f%s' % (self.metric_plotting_info['type'][sym_idx_self],
                                                    subj_sym_metric,
                                                    self.normativeModel.uncertainty['orig'][sym_idx_normativeModel],
                                                    self.metric_plotting_info['units'][sym_idx_self]))

                closest_agevec = np.argmin(np.abs(self.normativeModel.age_vec - subj_ages))
                ax2.text(xtext, ylim[1] + 0.05 * yrange, r'Norm=%1.2f%s, CI=[%1.2f-%1.2f]' % (
                    self.normativeModel.fit_outputs['orig']['fit_ave'][closest_agevec, sym_idx_normativeModel],
                    self.metric_plotting_info['units'][sym_idx_self],
                    self.normativeModel.fit_outputs['orig']['fit_sml'][closest_agevec, sym_idx_normativeModel],
                    self.normativeModel.fit_outputs['orig']['fit_lrg'][closest_agevec, sym_idx_normativeModel]))

                ax2.text(xtext, ylim[0] - 0.01 * yrange,
                         r'z$_{dev}$=%1.3f, p$_{dev}$=%1.3f' % (
                             output['orig']['devtn_stats'][run_idx, sym_idx_output],
                             output['orig']['devtn_pvals'][run_idx, sym_idx_output]))

                ax2.text(xtext, ylim[0] - 0.10 * yrange,
                         r'Odds=%1.3f, outlier_frac=%1.3f, P(art.)=%1.3f' % (
                             self.normativeModel.stats['orig']['odds'][sym_idx_normativeModel],
                             self.normativeModel.stats['orig']['fout'][sym_idx_normativeModel],
                             self.normativeModel.stats['orig']['part'][sym_idx_normativeModel]))

                ax2.set_ylim(ylim[0] - 0.15 * yrange, ylim[1] + 0.15 * yrange)

                # Plot contralateral pointcloud
                bad_opposite = np.argwhere(self.normativeModel.outliers[k][:, col_idx_normativeModel_opposite] == 1)
                good_opposite = np.argwhere(self.normativeModel.outliers[k][:, col_idx_normativeModel_opposite] == 0)
                match_bad_opposite = np.intersect1d(session_matches, bad_opposite)
                match_good_opposite = np.intersect1d(session_matches, good_opposite)
                nonmatch_bad_opposite = np.intersect1d(
                    np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), session_matches),
                    bad_opposite)
                nonmatch_good_opposite = np.intersect1d(
                    np.delete(np.arange(self.normativeModel.outliers[k].shape[0]), session_matches),
                    good_opposite)

                ax3.plot(self.normativeModel.age_vec,
                         self.normativeModel.fit_outputs[k]['fit_ave'][:, col_idx_normativeModel_opposite], 'k')
                ax3.plot(normModel_age[match_good_opposite],
                         self.normativeModel.measured_metrics[k][
                             match_good_opposite, col_idx_normativeModel_opposite], 'ko',
                         markersize=6, markerfacecolor='none')
                ax3.plot(normModel_age[match_bad_opposite],
                         self.normativeModel.measured_metrics[k][
                             match_bad_opposite, col_idx_normativeModel_opposite], 'kx',
                         markersize=6)
                ax3.plot(normModel_age[nonmatch_good_opposite],
                         self.normativeModel.measured_metrics[k][
                             nonmatch_good_opposite, col_idx_normativeModel_opposite], 'ko',
                         markersize=3, markerfacecolor='none', alpha=0.3)
                ax3.plot(normModel_age[nonmatch_bad_opposite],
                         self.normativeModel.measured_metrics[k][
                             nonmatch_bad_opposite, col_idx_normativeModel_opposite], 'kx',
                         markersize=3, alpha=0.3)

                ax3.scatter(subj_ages, subj_opposite_metric, s=50,
                            facecolors='b', edgecolors='b')
                ax3.errorbar(subj_ages, subj_opposite_metric,
                             yerr=self.normativeModel.uncertainty[k][col_idx_normativeModel_opposite],
                             ecolor='b',
                             ls='None')

                pval_opposite = output[k]['devtn_pvals'][run_idx, col_idx_output_opposite]
                logging.PRINT(f"Opposite hemisphere - p-value: {pval_opposite}, pid: {pid}, alpha_uncorr: {alpha_uncorr}")
                if pval_opposite < pid:
                    ax3.set_facecolor([1, 0.5, 0.5])
                    ax3.set_alpha(0.5)
                elif pval_opposite < alpha_uncorr:
                    ax3.set_facecolor('yellow')
                    ax3.set_alpha(0.5)

                ax3.grid(visible=True, linewidth=0.5, alpha=0.5)
                ax3.set_xlabel('Age (years)')
                ax3.set_ylabel(self.metric_names[col_idx_self_opposite])
                ax3.set_title(f'Session: {selected_session}, Acquisition: {selected_acquisition} (Right Hemisphere)')

                ylim = ax3.get_ylim()
                yrange = ylim[1] - ylim[0]
                xlim = ax3.get_xlim()
                xtext = xlim[0] + 0.05 * (xlim[1] - xlim[0])

                ax3.text(xtext, ylim[1] + 0.10 * yrange,
                         r'%s=%1.2f$\pm$%1.2f%s' % (self.metric_plotting_info['type'][col_idx_self_opposite],
                                                    subj_opposite_metric,
                                                    self.normativeModel.uncertainty[k][
                                                        col_idx_normativeModel_opposite],
                                                    self.metric_plotting_info['units'][col_idx_self_opposite]))

                closest_agevec = np.argmin(np.abs(self.normativeModel.age_vec - subj_ages))
                ax3.text(xtext, ylim[1] + 0.05 * yrange, r'Norm=%1.2f%s, CI=[%1.2f-%1.2f]' % (
                    self.normativeModel.fit_outputs[k]['fit_ave'][closest_agevec, col_idx_normativeModel_opposite],
                    self.metric_plotting_info['units'][col_idx_self_opposite],
                    self.normativeModel.fit_outputs[k]['fit_sml'][closest_agevec, col_idx_normativeModel_opposite],
                    self.normativeModel.fit_outputs[k]['fit_lrg'][closest_agevec, col_idx_normativeModel_opposite]))

                ax3.text(xtext, ylim[0] - 0.01 * yrange,
                         r'z$_{dev}$=%1.3f, p$_{dev}$=%1.3f' % (
                             output[k]['devtn_stats'][run_idx, col_idx_output_opposite],
                             output[k]['devtn_pvals'][run_idx, col_idx_output_opposite]))

                ax3.text(xtext, ylim[0] - 0.10 * yrange,
                         r'Odds=%1.3f, outlier_frac=%1.3f, P(art.)=%1.3f' % (
                             self.normativeModel.stats[k]['odds'][col_idx_normativeModel_opposite],
                             self.normativeModel.stats[k]['fout'][col_idx_normativeModel_opposite],
                             self.normativeModel.stats[k]['part'][col_idx_normativeModel_opposite]))

                ax3.set_ylim(ylim[0] - 0.15 * yrange, ylim[1] + 0.15 * yrange)

        except IndexError as e:
            error_message = f"IndexError for subject_id {subject_id}, session_id {selected_session}, acq_id {selected_acquisition}, metric {selected_metric}: {e}"
            print(error_message)
            QMessageBox.critical(None, "IndexError", error_message)

        except Exception as e:
            error_message = f"Unexpected error for subject_id {subject_id}, session_id {selected_session}, acq_id {selected_acquisition}, metric {selected_metric}: {e}"
            print(error_message)
            QMessageBox.critical(None, "Unexpected Error", error_message)

        return fig

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
                             in testing all metrics intersecting self.metric_names and self.normativeModel.metric_names.
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
        normModel_matching_cols = np.array( [self.normativeModel.covariate_names.index(matching_name) for matching_name in matching_covariates])
        self_matching_cols = np.array( [self.covariate_names.index(matching_name) for matching_name in matching_covariates])
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
        predicted_values = self.normativeModel.predict_values(self.covariate_values[:, self.covariate_names.index('age')])
        for k in normalizations:
            # Get column matching between self.measured_metrics and self.normativeModel.measured_metrics, using metric_names
            col_idxs_normativeModel = []
            col_idxs_self = []
            metric_names_union = [metric_name for metric_name in metric_names[:self.measured_metrics[k].shape[1]] if
                                  metric_name in self.normativeModel.metric_names[
                                                 :self.normativeModel.measured_metrics[k].shape[1]]]
            for metric_name in metric_names_union:
                col_idxs_normativeModel.append(self.normativeModel.metric_names.index(metric_name))
                col_idxs_self.append(metric_names.index(metric_name))
            # Allocate output
            output[k] = {}
            output[k]['gp_residuals'] = np.full((len(subj_rows), len(col_idxs_self)), np.nan)
            output[k]['gp-vs-gp_dfs'] = np.full((len(col_idxs_self)), np.nan)
            output[k]['gp-vs-gp_ts'] = np.full((len(col_idxs_self)), np.nan)
            output[k]['gp-vs-gp_pvals'] = np.full((len(col_idxs_self)), np.nan)
            output[k]['gp-vs-gp_logps'] = np.full((len(col_idxs_self)), np.nan)
            output[k]['gp-vs-gp_cohen-d'] = np.full((len(col_idxs_self)), np.nan)
            output[k]['metric_names'] = [metric_names[i] for i in col_idxs_self]
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
                                               subj_covariate_values[i_scan, self_matching_cols]).all(axis=1))[0]
                matching_rows.append(scan_matches)
            matching_rows = np.unique(np.hstack(matching_rows))  # Group all scan matches and remove duplicates
            matching_residuals = self.normativeModel.fit_outputs[k]['residuals'][matching_rows, :][:, col_idxs_normativeModel].copy()
            outliers = self.normativeModel.outliers[k][matching_rows, :][:, col_idxs_normativeModel].copy()
            matching_residuals[outliers] = np.nan  # sets outliers to nan, so matching rows still includes outlier subjects, but not taken into account when using nan_policy='omit'
            if len(matching_residuals.shape) == 1:
                matching_residuals = matching_residuals[None, :]  # Make the array 2D to average over 1st dimension later
            # Test for statistical difference between group residuals and matching residuals
            output[k]['gp_residuals'] = subj_measured_metrics[:, col_idxs_self] - predicted_values[k][:, col_idxs_normativeModel]
            for i_metric in range(len(col_idxs_self)):
                subj_noNans = output[k]['gp_residuals'][~np.isnan(output[k]['gp_residuals'][:, i_metric]), i_metric]
                matching_noNans = matching_residuals[~np.isnan(matching_residuals[:, i_metric]), i_metric]
                t, p = ttest_ind(subj_noNans, matching_noNans)
                output[k]['gp-vs-gp_dfs'][i_metric] = len(subj_noNans)+len(matching_noNans)-1
                output[k]['gp-vs-gp_ts'][i_metric] = t
                output[k]['gp-vs-gp_pvals'][i_metric] = p
                output[k]['gp-vs-gp_logps'][i_metric] = -np.sign(t) * np.log10(p)
                n1 = len(matching_noNans)
                n2 = len(subj_noNans)
                pooled_std = np.sqrt(((n1-1)*matching_noNans.std(ddof=1)**2+(n2-1)*subj_noNans.std(ddof=1)**2)/(n1+n2-2))
                output[k]['gp-vs-gp_cohen-d'][i_metric] = (np.nanmean(output[k]['gp_residuals'][:, i_metric])-np.nanmean(matching_residuals[:, i_metric]))/pooled_std  # d<0 == test_group < norm, d>0 == test_group > norm
        return output
