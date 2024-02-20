"""
Wrapper for DL+DiReCT scripts to be run on <subjects_dir>/<subjid>. Includes pipeline template parameters as
proc_pipeline_name, and methods as run(), and proc2metric().
"""


from shutil import rmtree, copyfile, which
import subprocess
import numpy as np
import os
from scanometrics.utils import logging
from csv import DictReader
from numpy.polynomial import Polynomial as poly_model
from nibabel.freesurfer.io import read_annot, write_annot
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull
from multiprocessing import Process, cpu_count
import re


class proc_pipeline:

    def __init__(self, bids_database, dldirect_version='1.0.0', load_from_ID=True, subjSesAcq_delimiter="_"):
        """
        DL+DiReCT processing pipeline constructor. Sets the bids path where derivatives/dldirect is or will be stored.
        Allows specifying the version of dldirect to use. The `load_from_ID` parameter can be set to False when loading
        of metrics is intended to be done from a single stats2table file with all subjects.

        :param bids_database: path to bids database folder.
        :type bids_database: string
        :param dldirect_version: version of dldirect to be used when processing, or loading atlas region names.
        :type dldirect_version: string
        :param load_from_ID: set to True to load data from each subject's folder. Set to False if data should be loaded
                             from files containing all subjects, and specify 'stats2table_folder' when calling the
                             proc2metric() method.
        :type load_from_ID: bool
        """
        self.bids_database = bids_database
        self.proc_pipeline_name = 'dldirect'
        self.dldirect_version = dldirect_version
        self.subjects_dir = os.path.join(bids_database, 'derivatives', self.proc_pipeline_name)
        self.load_from_ID = load_from_ID
        self.subjSesAcq_delimiter = subjSesAcq_delimiter
        self.setting = {'DL_DiReCT_condaenv': 'DL_DiReCT', 'model': 'v0'}
        self.seg_metrics = ['volume']
        self.parc35_metrics = ['volume', 'thickness', 'thicknessstd']
        self.parc75_metrics = ['volume', 'thickness', 'thicknessstd']

        # global variables
        self.atlases = {'1.0.0': {'DesikanKilliany': {'ROIs': ['bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'cuneus', 'entorhinal', 'fusiform', 'inferiorparietal', 'inferiortemporal', 'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal', 'lingual', 'medialorbitofrontal', 'middletemporal', 'parahippocampal', 'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral', 'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate', 'rostralmiddlefrontal', 'superiorfrontal', 'superiorparietal', 'superiortemporal', 'supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal', 'insula'],
                                              'views': ['lateral',                  'medial',             'lateral', 'medial',     'medial',   'medial',          'lateral',          'lateral',       'medial',              'lateral',              'lateral',  'medial',              'medial',        'lateral',          'medial',      'medial',         'lateral',       'lateral',          'lateral',        'medial',     'lateral',             'medial',    'lateral',    'medial',                   'medial',              'lateral',         'lateral',          'lateral',          'lateral',       'lateral',      'medial',       'medial',            'lateral', 'lateral']},
                                   'Destrieux': {'ROIs': ['G&S_frontomargin', 'G&S_occipital_inf', 'G&S_paracentral', 'G&S_subcentral', 'G&S_transv_frontopol', 'G&S_cingul-Ant', 'G&S_cingul-Mid-Ant', 'G&S_cingul-Mid-Post', 'G_cingul-Post-dorsal', 'G_cingul-Post-ventral', 'G_cuneus', 'G_front_inf-Opercular', 'G_front_inf-Orbital', 'G_front_inf-Triangul', 'G_front_middle', 'G_front_sup', 'G_Ins_lg&S_cent_ins', 'G_insular_short', 'G_occipital_middle', 'G_occipital_sup', 'G_oc-temp_lat-fusifor', 'G_oc-temp_med-Lingual', 'G_oc-temp_med-Parahip', 'G_orbital', 'G_pariet_inf-Angular', 'G_pariet_inf-Supramar', 'G_parietal_sup', 'G_postcentral', 'G_precentral', 'G_precuneus', 'G_rectus', 'G_subcallosal', 'G_temp_sup-G_T_transv', 'G_temp_sup-Lateral', 'G_temp_sup-Plan_polar', 'G_temp_sup-Plan_tempo', 'G_temporal_inf', 'G_temporal_middle', 'Lat_Fis-ant-Horizont', 'Lat_Fis-ant-Vertical', 'Lat_Fis-post', 'Pole_occipital', 'Pole_temporal', 'Medial_wall', 'S_calcarine', 'S_central', 'S_cingul-Marginalis', 'S_circular_insula_ant', 'S_circular_insula_inf', 'S_circular_insula_sup', 'S_collat_transv_ant', 'S_collat_transv_post', 'S_front_inf', 'S_front_middle', 'S_front_sup', 'S_interm_prim-Jensen', 'S_intrapariet&P_trans', 'S_oc_middle&Lunatus', 'S_oc_sup&transversal', 'S_occipital_ant', 'S_oc-temp_lat', 'S_oc-temp_med&Lingual', 'S_orbital_lateral', 'S_orbital_med-olfact', 'S_orbital-H_Shaped', 'S_parieto_occipital', 'S_pericallosal', 'S_postcentral', 'S_precentral-inf-part', 'S_precentral-sup-part', 'S_suborbital', 'S_subparietal', 'S_temporal_inf', 'S_temporal_sup', 'S_temporal_transverse'],
                                              'views': [        'lateral',           'lateral',          'medial',        'lateral',              'lateral',         'medial',             'medial',              'medial',               'medial',                'medial',   'medial',               'lateral',             'lateral',              'lateral',        'lateral',      'medial',             'lateral',         'lateral',            'lateral',          'medial',                'medial',                'medial',                'medial',   'lateral',              'lateral',               'lateral',        'lateral',       'lateral',      'lateral',      'medial',   'medial',        'medial',               'lateral',            'lateral',                'medial',               'lateral',        'lateral',           'lateral',              'lateral',              'lateral',      'lateral',        'lateral',       'lateral',     'medial' ,      'medial',   'lateral',              'medial',               'lateral',               'lateral',               'lateral',              'medial',               'medial',     'lateral',        'lateral',     'lateral',              'lateral',               'lateral',             'lateral',              'lateral',         'lateral',       'lateral',                'medial',           'lateral',              'lateral',            'lateral',              'medial',         'medial',       'lateral',               'lateral',               'lateral',       'medial',        'medial',        'lateral',        'lateral',               'lateral']}},
                        }

        self.lobe_ROIs = [['frontalLobe', ['caudalmiddlefrontal', 'lateralorbitofrontal', 'medialorbitofrontal',
                                           'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis',
                                           'precentral', 'rostralmiddlefrontal', 'superiorfrontal', 'frontalpole']],
                          ['parietalLobe', ['inferiorparietal', 'postcentral', 'precuneus', 'superiorparietal', 'supramarginal']],
                          ['occipitalLobe', ['cuneus', 'lateraloccipital', 'lingual', 'pericalcarine']],
                          ['temporalLobe', ['bankssts', 'entorhinal', 'fusiform', 'inferiortemporal', 'middletemporal',
                                            'parahippocampal', 'superiortemporal', 'temporalpole', 'transversetemporal']],
                          ['cingulateLobe', ['caudalanteriorcingulate', 'isthmuscingulate', 'posteriorcingulate',
                                             'rostralanteriorcingulate']],
                          ['insulaLobe', ['insula']]]

        self.TIVproxy_ROIs = ["Left-Cerebral-White-Matter","Right-Cerebral-White-Matter",
                              "Left-Ventricle-all","Right-Ventricle-all",
                              "Left-Thalamus-Proper","Right-Thalamus-Proper",
                              "Left-Caudate","Right-Caudate",
                              "Left-Putamen","Right-Putamen",
                              "Left-Pallidum","Right-Pallidum",
                              "Left-Hippocampus","Right-Hippocampus",
                              "Left-Amygdala","Right-Amygdala",
                              "Left-Accumbens-area","Right-Accumbens-area",
                              "3rd-Ventricle",
                              "4th-Ventricle",
                              "lhCortexVol","rhCortexVol"]

    def set(self, setting_key, setting_value):
        self.setting[setting_key] = setting_value

    def run_pipeline(self, subjects, n_threads):
        """Pipeline template method overided here. Calls DL+DiReCT processing scripts. Single or multi-subject
        processing controlled through n_threads (default of -1 allocates all cpus available through a call to
        multiprocessing.cpu_count(). Processes scans with batch-dl+direct after renaming all scans to T1w and saving them
        in a flattened folder structure in <bids_database>/derivatives/tmp_dldirectRawdata. The temporary directory is
        deleted afterwards, and the outputs of DL+DiReCT are saved into individual folders in <bids_database>/derivatives/
        dldirect/<sub-label_ses-label_acq-label>

        :param subjects: multi-layered dictionary, with first level being subject IDs, then session IDs, followed by the
                         acq labels. Used to identify T1w scans to process.
        :type subjects: dictionary
        :param n_threads: number of threads to use (default calls set to -1 to use all possible cpus). GPU calls all
                          available ressources, as available from torch.cuda.device_count().
        :type n_threads: int
        """
        import torch
        # os.system('eval "$(conda shell.bash hook)"')
        # os.system('conda activate DL_DiReCT')
        if which('dl+direct') is None:
            logging.ERROR("""dl+direct command not found. Make sure you have a conda environment that has DL+DiReCT
            properly installed, and that the environment can be activated through 'conda activate DL_DiReCT. You can
            change the name of the environment with SOM.processing_pipeline.set('DL_DiReCT_condaenv', <your_env_name>).
            If your environment is not saved in Conda's envs folder, you can specify the path to the environment instead
            of its name.""")
        # List subjects with acq_label in their BIDS anat folder, and path to file
        subjSesAcq_list, T1_list = self.get_subjSesAcq_T1s(subjects)
        # Copy T1s from BIDS anat folder to <bids_directory>/derivatives/tmp_dldirectRawdata
        for subjSesAcq_id, T1_file in zip(subjSesAcq_list, T1_list):
            subjSesAcq_folder = os.path.join(self.bids_database, 'derivatives', 'tmp_dldirectRawdata', subjSesAcq_id)
            if not os.path.exists(subjSesAcq_folder):
                os.makedirs(subjSesAcq_folder)
            copyfile(T1_file, os.path.join(subjSesAcq_folder, 'T1w.nii.gz'))  # rename all files to same name to get dldirect to work as intended (aqc label already in folder name)
        if n_threads == -1:
            n_threads = cpu_count()
        subprocess.run(['batch-dl+direct', '--cpu', str(n_threads),
                        '--gpu', str(torch.cuda.device_count()), '--bet',
                        '--t1', 'T1w.nii.gz', '--model', self.setting['model'],
                        os.path.join(self.bids_database, 'derivatives', 'tmp_dldirectRawdata'), self.subjects_dir])
        # subprocess.run(['conda', 'deactivate'])
        rmtree(os.path.join(self.bids_database, 'derivatives', 'tmp_dldirectRawdata'))

    def proc2table(self, subjects, n_threads):
        """
        n_thread is kept for compatibility with core.py functions, but DL+DiReCT stat2table script currently gathers all
        subject statistics into grouped files, so n_thread is ignored.

        :param subjects: dictionary with all subject IDs, session IDs, and acquisition labels used to identify scans.
        :type subjects: dictionary
        :param n_threads: ignored
        :type n_threads: int
        """
        # The following doesn't work, try using something similar to subprocess.run('conda init bash;source activate DL_DiReCT;which python', shell=True, executable='/bin/bash')
        # subprocess.run(['conda', 'activate', self.setting['DL_DiReCT_condaenv']])
        subjSesAcq_list = self.get_subjSesAcq_array(subjects)
        p = Process(target=self.run_proc2table, args=(subjSesAcq_list,))
        p.start()
        p.join()
        # subprocess.run(['conda', 'deactivate'])

    def run_proc2table(self, subjSesAcq_list):
        """Wrapper for DL+DiReCT stats2table. Assumes subjects were processed with DL+DiReCT.
        Creates new directories <subjects_dir>/<subjSesAcq_id>/stats2table (deletes old ones if present).

        :param subjSesAcq_list: list of subjects and sessions to include (lines with subjects and sessions not in the
                                list are ignored).
        :type subjSesAcq_list: list of strings
        """
        if which('stats2table') is None:
            logging.ERROR("""stats2table command not found. Make sure you have a conda environment that has DL+DiReCT
            properly installed, and that the environment can be activated through 'conda activate DL_DiReCT'. You can
            change the name of the environment with SOM.processing_pipeline.set('DL_DiReCT_condaenv', <your_env_name>).
            If your environment is not saved in Conda's envs folder, you can specify the path to the environment instead
            of its name.""")
        subprocess.run(['stats2table', self.subjects_dir, os.path.join(self.bids_database, 'derivatives', 'tmp_stats2tableAll')])
        split_subj_ses = re.compile('sub-(.*)_ses-(.*)')
        # Split stats2table file into several individual files
        for stats2table_in_file in ['aseg_stats_volume.txt'] + ['%s.aparc_stats_%s.txt' % (hemi, metric) for hemi in
                                                                ['lh', 'rh'] for metric in
                                                                ['volume', 'thickness', 'thicknessstd']]:
            with open(os.path.join(self.bids_database, 'derivatives', 'tmp_stats2tableAll', stats2table_in_file), 'r') as f_in:
                header = next(f_in)
                for line in f_in:
                    orig_id = line.split("\t")[0]
                    res = split_subj_ses.search(orig_id)
                    subj_id = res.group(1)
                    ses_id = res.group(2)
                    if ('sub-%s_ses-%s' % (subj_id, ses_id)) in subjSesAcq_list:
                        stats2table_folder = os.path.join(self.subjects_dir, 'sub-%s_ses-%s' % (subj_id, ses_id), 'stats2table')
                        if not os.path.exists(stats2table_folder):
                            os.makedirs(stats2table_folder)
                        with open(os.path.join(stats2table_folder, stats2table_in_file), 'w') as f_out:
                            f_out.write(header)
                            f_out.write('sub-%s_ses-%s\t' % (subj_id, ses_id))
                            f_out.write("\t".join(line.split("\t")[1:]))
        rmtree(os.path.join(self.bids_database, 'derivatives', 'tmp_stats2tableAll'))

    def proc2metric(self, subjects, covariate_values, covariate_names, stats2table_folder=None, ref_rows=None,
                    ref_metric_values=None, ref_metric_names=None, ref_covariate_values=None, ref_covariate_names=None,
                    metric_include=None, metric_exclude=None):
        """
        Loads all measurements from stats2table files. The `load_from_ID` option allows to load from individual files
        in the bids structure when set to True. If set to False, then `stats2table_folder` must be specified for SOM to
        load all files in there, which should contain all metrics for all subjects. The different ref_* parameters allow
        specifying the metric_values, metric_names, covariate_values and covariate_names of the reference dataset to use
        (typically a copy of the SOM.normativeModel.* variables).

        :param subjects: dictionary with subject IDs, session IDs, and acquisition labels.
        :type subjects: dictionary
        :param covariate_values: covariate values for subjects in the subjects variable.
        :type covariate_values: numpy matrix of size len(subjSesAqc_array)*n_covariates
        :param covariate_names: list of covariate names corresponding to columns in covariate_values.
        :type covariate_names: list
        :param stats2table_folder: path to folder containing files with all subjects measurements, to be loaded if the
                                   `load_from_ID` parameter in SOM.proc_pipeline is set to True.
        :type stats2table_folder: string
        :param ref_rows: array of row indexes in ref_metric_values and ref_covariate_values to use as reference when
                         normalizing metrics.
        :type ref_rows: 1D array
        :param ref_metric_values: numpy matrix with normative values to be used as reference when normalizing
        :type ref_metric_values: numpy matrix
        :param ref_metric_names: list of metric names corresponding to columns in ref_metric_values
        :type ref_metric_names: list of strings
        :param ref_covariate_values: numpy matrix with covariate_values used to select matching normative subjects in terms of
                                     age, sequence, scanner, etc... when normalizing
        :type ref_covariate_values: numpy matrix
        :param ref_covariate_names: list of covariate names corresponding to colums in ref_covariate_values
        :type ref_covariate_names: list of strings
        """
        # Error checks
        if 'sequence' not in covariate_names or 'scanner' not in covariate_names:
            logging.ERROR("The dldirect pipeline requires 'scanner' and 'sequence' as covariates, but were not"
                          " found in 'covariate_names' list. Please run the 'proc2metric' command again after including"
                          " them. You can specify them in the 'covariates_include' list when running SOM.load_subjects(),"
                          " or avoid specifying them in the 'covariates_exclude' list, or avoid specifying both inclusion"
                          " and exclusion lists to load all available covariates. You might want to check if the"
                          " covariate is included in the BIDS 'participants.tsv' file too.\n")
        # Metric filtering: add components of symmetry indexes and lobe indexes if they are in the include option
        # Dirty way is to add both lh/Left and rh/Right labels as it seems it's not breaking anything and would be
        # more tedious to figure out which label to add for which metric
        if metric_include is not None:
            metric_include_bckp = metric_include.copy()
            for metric_name in metric_include:
                if 'symmetryIndex' in metric_name:
                    for hemi in ['rh', 'Right', 'lh', 'Left']:
                        if metric_name.replace('symmetryIndex', hemi) not in metric_include:
                            metric_include.append(metric_name.replace('symmetryIndex', hemi))
                if 'Lobe' in metric_name:
                    hemi, lobe, metric = metric_name.split("_")
                    if metric == 'volume':  # Lobe thickness requires ROI area to do a weighted sum, we only consider the volume for now until DL+DiReCT outputs each ROI's area.
                        for ctx_label in self.lobe_ROIs[[self.lobe_ROIs[i][0] for i in range(len(self.lobe_ROIs))].index(lobe)][1]:
                            metric_include.append('aparc_%s_%s_%s' % (hemi, ctx_label, metric))
            # Add TIVproxyROI volumes for normalization if they are not included
            for TIVproxy_roi in self.TIVproxy_ROIs:
                metric_name = 'asegvolume_%s' % TIVproxy_roi
                if metric_name not in metric_include:
                    metric_include.append(metric_name)
        subject_list = self.get_subjSesAcq_array(subjects)
        if self.load_from_ID:
            metric_values, metric_names = self.load_subject_metrics_from_ID(subject_list[0], metric_include=metric_include, metric_exclude=metric_exclude)
            for ID in subject_list[1:]:
                subj_metric_values, subj_metric_names = self.load_subject_metrics_from_ID(ID, metric_include=metric_include, metric_exclude=metric_exclude)
                # Add potential new metric names at the end, and fill old metric_values with nans
                new_metric_names = [nm for nm in subj_metric_names if nm not in metric_names]
                metric_names += new_metric_names
                metric_values = np.hstack((metric_values, np.full((metric_values.shape[0], len(new_metric_names)), np.nan)))
                # Add nans to new subject values if metric does not exists, reorder everything to match metric_names
                missing_metric_names = [nm for nm in metric_names if nm not in subj_metric_names]
                subj_metric_names += missing_metric_names
                subj_metric_values = np.hstack((subj_metric_values, np.full((1, len(missing_metric_names)), np.nan)))
                metric_values = np.vstack((metric_values, subj_metric_values[:, np.array([subj_metric_names.index(m) for m in metric_names], dtype='int')]))
        else:
            if stats2table_folder is None:
                logging.ERROR('load_from_stats2tableFolder requires a path to the tables folder, which is currently set'
                              ' to None. Provide the path when calling the proc2metric() function.')
            metric_values, metric_names = self.load_subject_metrics_from_stats2tableFolder(stats2table_folder, subject_list, metric_include=metric_include, metric_exclude=metric_exclude)

        # Compute lobe metrics
        for hemi in ['lh', 'rh']:
            for lobe, ctx_labels in self.lobe_ROIs:
                for metric in self.parc35_metrics:
                    # Skip if not in metric_include or in metric_exclude
                    if metric_include is not None and '%s_%s_%s' % (hemi, lobe, metric) not in metric_include:
                        continue
                    if metric_exclude is not None and '%s_%s_%s' % (hemi, lobe, metric) in metric_exclude:
                        continue
                    # Lobe thickness requires ROI area to do a weighted sum, we only consider the volume for now until
                    # DL+DiReCT outputs each ROI's area.
                    if metric != 'volume':
                        continue
                    col_idx = np.array([metric_names.index('aparc_%s_%s_%s' % (hemi, ctx_label, metric)) for ctx_label in ctx_labels])
                    metric_names += ['%s_%s_%s' % (hemi, lobe, metric)]
                    metric_values = np.hstack((metric_values, np.nansum(metric_values[:, col_idx], axis=1, keepdims=True)))

        # Compute symmetric index
        sym_metric_names = []
        sym_metric_values = []
        for i, m in enumerate(metric_names):
            if 'lh' in m:  # lh_ label means there's a column with rh_ label that corresponds (look for lh instead of rh because rh matches entorhinal and other structures with rh in it)
                num = (metric_values[:, metric_names.index(m.replace('lh', 'rh'))] - metric_values[:, i])
                denom = (np.abs(metric_values[:, metric_names.index(m.replace('lh', 'rh'))]) + np.abs(metric_values[:, i]))
                if metric_include is not None and m.replace('lh', 'symmetryIndex') not in metric_include:
                    continue
                if metric_exclude is not None and m.replace('lh', 'symmetryIndex') in metric_exclude:
                    continue
                sym_metric_values.append(num / denom)  # Sym index is (rh-lh)/(rh+lh)
                sym_metric_names.append(m.replace('lh', 'symmetryIndex'))
            elif 'Left' in m:  # Left label means there's a column with Right label that corresponds
                num = (metric_values[:, metric_names.index(m.replace('Left', 'Right'))] - metric_values[:, i])
                denom = (np.abs(metric_values[:, metric_names.index(m.replace('Left', 'Right'))]) + np.abs(metric_values[:, i]))
                if metric_include is not None and m.replace('Left', 'symmetryIndex') not in metric_include:
                    continue
                if metric_exclude is not None and m.replace('Left', 'symmetryIndex') in metric_exclude:
                    continue
                sym_metric_values.append(num / denom)  # Sym index is (rh-lh)/(rh+lh)
                sym_metric_names.append(m.replace('Left', 'symmetryIndex'))

        # Normalize original data
        if ref_metric_values is None and ref_metric_names is None:
            ref_metric_values = metric_values.copy()
            ref_metric_names = metric_names.copy()
        # Set ref_covariate_values and ref_covariate_names
        if ref_covariate_values is None and ref_covariate_names is None:
            ref_covariate_values = covariate_values.copy()
            ref_covariate_names = covariate_names.copy()
        # Extract ref_age vector
        ref_age = ref_covariate_values[:, ref_covariate_names.index('age')].copy()
        if ref_rows is None:  # Set ref_rows to all subjSesAcq if not specified
            ref_rows = np.arange(ref_metric_values.shape[0])
        ref_metric_values = ref_metric_values[ref_rows, :]
        ref_covariate_values = ref_covariate_values[ref_rows, :]
        ref_eTIV = np.nanmean(ref_metric_values[:, ref_metric_names.index('asegvolume_EstimatedTotalIntraCranialVol')])
        # Future work: split into a different submodule with a normalization dict argument.
        norm_coef = np.full(len(metric_names), np.nan)
        metric_unit = []
        metric_type = []
        for j, m in enumerate(metric_names+sym_metric_names):
            if 'symmetryIndex' in m:
                metric_unit.append(r'')
                metric_type.append('Symmetry index')
                continue
            if m in ['asegvolume_EstimatedTotalIntraCranialVol', 'aparc_eTIV']:
                metric_unit.append(r'')
                metric_type.append('Estimated total ICV')
                continue
            if 'volume' in m:
                norm_coef[j] = 1.0
                metric_unit.append(r'mm$^3$')
                metric_type.append('Volume')
            elif 'area' in m:
                norm_coef[j] = 0.6660
                metric_unit.append(r'mm$^2$')
                metric_type.append('Area')
            elif 'thicknessstd' in m:
                norm_coef[j] = 0.3330
                metric_unit.append(r'mm')
                metric_type.append('Thickness std')
            elif 'thickness' in m:
                norm_coef[j] = 0.3330
                metric_unit.append(r'mm')
                metric_type.append('Thickness')
            elif 'meancurv' in m:
                norm_coef[j] = -0.3330
                metric_unit.append(r'mm$^{-1}$')
                metric_type.append('Mean curvature')
            elif 'gauscurv' in m:
                norm_coef[j] = -0.6660
                metric_unit.append(r'mm$^{-2}$')
                metric_type.append('Gaus. curvature')
            elif 'foldind' in m:
                norm_coef[j] = -0.6660
                metric_unit.append(r'mm$^{-2}$')
                metric_type.append('Folding index')
            elif 'curvind' in m:
                norm_coef[j] = -0.6660
                metric_unit.append(r'mm$^{-2}$')
                metric_type.append('Curve index')
            elif 'pctmean' in m:
                metric_unit.append(r'%')
                metric_type.append('mean PCT')
            elif 'pctstd' in m:
                metric_unit.append(r'%')
                metric_type.append('std PCT')
            elif 'pctsnr' in m:
                metric_unit.append(r'%')
                metric_type.append('PCT snr')
            else:
                metric_unit.append(r'unknown')
                metric_type.append('unknown')
                logging.PRINT('%s does not match any metric type in Freesurfer.py metric list' % m)

        # Get indexes of unique combinations of scanner and sequence in current subjects
        matching_cols = np.array([covariate_names.index('sequence'), covariate_names.index('scanner')])
        ref_matching_cols = np.array([ref_covariate_names.index('sequence'), ref_covariate_names.index('scanner')])
        scanseq_combinations = np.unique(ref_covariate_values[:, ref_matching_cols], axis=0)

        # Normalize metrics
        normalized_metrics = np.full(metric_values.shape, np.nan)
        TIV = metric_values[:, metric_names.index('asegvolume_EstimatedTotalIntraCranialVol')].copy()  # has length # subjects
        TIV /= ref_eTIV
        vol_scaling_factor = (np.tile(TIV[:, None], (1, len(metric_names))))**(-np.tile(norm_coef[None, :], (metric_values.shape[0], 1)))
        volnorm_col_idxs = np.argwhere(~np.isnan(norm_coef))
        normalized_metrics[:, volnorm_col_idxs] = metric_values[:, volnorm_col_idxs] * vol_scaling_factor[:, volnorm_col_idxs]
        for j in np.argwhere(np.isnan(norm_coef))[:, 0]:
            slopes = []
            intercepts = []
            weights = []
            # Loop through combinations of scanseq.
            # NB: differs from octave implementation as handles cases when only one subject has the combination (then single
            #     value is taken as the intercept
            for scanseq in scanseq_combinations:
                ref_matching_rows = (np.where((ref_covariate_values[:, ref_matching_cols] == scanseq).all(axis=1))[0]).astype('int')
                x = ref_age[ref_matching_rows]
                y = ref_metric_values[ref_matching_rows, ref_metric_names.index(metric_names[j])]
                if len(ref_matching_rows) > 1 and ~np.any(np.isnan(y)) and len(np.unique(x)) > 1:
                    c = poly_model.fit(x, y, deg=1).convert()
                    intercepts.append(c.coef[0])
                    if len(c.coef) > 1:
                        slopes.append(c.coef[1])  # skips slope estimation when y has the same value for all matching subj
                    else:
                        slopes.append(np.nan)
                    weights.append(len(ref_matching_rows) / ref_metric_values.shape[0])
                elif len(ref_matching_rows) == 1:
                    intercepts.append(y[0])
                    slopes.append(0.0)
                    weights.append(1/ref_metric_values.shape[0])
            slopes = np.array(slopes)
            intercepts = np.array(intercepts)
            weights = np.array(weights)
            mean_slope = np.nansum(slopes * weights) / np.nansum(weights)
            mean_intercept = np.nansum(intercepts * weights) / np.nansum(weights)
            c = poly_model([mean_intercept, mean_slope])
            # Normalize by substracting shift for matching groups
            # NB: this leads to negative values if shift is greater than subject predicted value
            for scanseq in scanseq_combinations:
                matching_rows = np.where((covariate_values[:, matching_cols] == scanseq).all(axis=1))[0].astype('int')
                ref_matching_rows = (np.where((ref_covariate_values[:, ref_matching_cols] == scanseq).all(axis=1))[0]).astype('int')
                if len(matching_rows) > 0:
                    shift = np.nanmean(ref_metric_values[ref_matching_rows, ref_metric_names.index(metric_names[j])] - c(ref_age[ref_matching_rows]))
                    normalized_metrics[matching_rows, j] = metric_values[matching_rows, j] - shift
        # Filter out metrics (remove those metrics that were temporarily added to compute other ones, and reshuffle to match order from inclusion input)
        sym_metric_values = np.array(sym_metric_values, dtype='float').T
        if metric_include is not None:
            metric_include = metric_include_bckp.copy()
            idxs = np.array([metric_names.index(metric_name) for idx, metric_name in enumerate(metric_include) if metric_name in metric_names])
            if len(idxs) > 0:
                metric_names = [metric_names[i] for i in idxs]
                metric_values = metric_values[:, idxs]
                normalized_metrics = normalized_metrics[:, idxs]
            idxs = np.array([sym_metric_names.index(metric_name) for idx, metric_name in enumerate(metric_include) if metric_name in sym_metric_names])
            if len(idxs) > 0:
                sym_metric_values = sym_metric_values[:, idxs]
                sym_metric_names = [sym_metric_names[i] for i in idxs]
        if metric_exclude is not None:
            idxs = np.array([idx for idx, metric_name in enumerate(metric_names) if metric_name not in metric_exclude])
            if len(idxs) > 0:
                metric_names = [metric_names[i] for i in idxs]
                metric_values = metric_values[:, idxs]
                normalized_metrics = normalized_metrics[:, idxs]
            idxs = np.array([idx for idx, metric_name in enumerate(sym_metric_names) if metric_name not in metric_exclude])
            if len(idxs) > 0:
                sym_metric_values = sym_metric_values[:, idxs]
                sym_metric_names = [sym_metric_names[i] for i in idxs]
        # Combined orig and normalized matrices into a single dictionary
        if len(metric_values.shape) == 2 and len(sym_metric_values.shape) == 2:
            orig_metrics = np.hstack((metric_values, sym_metric_values))
        if len(metric_values.shape) == 1 and len(sym_metric_values).shape == 2:
            orig_metrics = sym_metric_values
        if len(metric_values.shape) == 2 and len(sym_metric_values.shape) == 1:
            orig_metrics = metric_values
        measured_metrics = {'orig': orig_metrics,
                            'norm': normalized_metrics}
        metric_plotting_info = {'units': metric_unit, 'type': metric_type}
        return metric_names+sym_metric_names, measured_metrics, metric_plotting_info

    def load_subject_metrics_from_ID(self, ID, metric_include=None, metric_exclude=None):
        """
        Load metrics from stats2table folder in <subjects_dir>/ID/stats2table. Here, ID is the combination of subject ID,
        session ID, and acquisition label, usually obtained when looping through SOM.subjects dictionary. All metrics are
        loaded into a single row which is returned along with metric names.

        :param ID: subject ID corresponding to the folder in <subjects_dir> that contains all metric files. The metric
                   files are supposed to contain a single line corresponding to ID, which should in turn correspond to
                   the combination of subject ID, session ID and acquisition label of the scan used to compute the metrics.
        :type ID: string
        """
        metric_names = []
        metric_values = []
        stats2table_folder = os.path.join(self.subjects_dir, ID, 'stats2table')
        logging.PRINT('stats2table_folder=%s' % stats2table_folder)
        for metric in self.seg_metrics:
            aseg_file = os.path.join(stats2table_folder, 'aseg_stats_%s.txt' % metric)
            logging.PRINT('aseg_file=%s' % aseg_file)
            if os.path.exists(aseg_file):
                logging.PRINT('aseg_file exists, loading it...')
                with open(aseg_file, 'r') as f:
                    dict_reader = DictReader(f, delimiter='\t')
                    tmp_metric_names = ['aseg%s_%s' % (metric, m) for m in dict_reader.fieldnames[1:]]  # Read metric names, skipping 1st as corresponds to participant ID
                    filtering_idxs = np.arange(len(tmp_metric_names))
                    if metric_include is not None:
                        filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name in metric_include]
                        filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name in metric_include])
                        tmp_metric_names = filtered_metric_names
                    if metric_exclude is not None:
                        filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name not in metric_exclude]
                        filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name not in metric_exclude])
                        tmp_metric_names = filtered_metric_names
                    for row in dict_reader:
                        logging.PRINT('Loaded metrics are %s' % (', '.join(tmp_metric_names)))
                        new_subj = list(row.values())[0]
                        if new_subj != ID:
                            logging.PRINT('Entered 1st if statement')
                            logging.WARNING("Subject %s is labeled as %s in DL+DiReCT output (aseg_stats_%s.txt)."
                                            " The subject was discarded." % (ID, new_subj, metric))
                        else:
                            logging.PRINT('Entered 2nd if statement')
                            metric_values.append(np.array(list(row.values())[1:], dtype='float')[filtering_idxs])  # skip participant ID and filter out metrics
                            metric_names += tmp_metric_names
                        logging.PRINT('Exiting if statement')
                logging.PRINT('Loaded metric_names are %s' % (', '.join(metric_names)))
                logging.PRINT('len(metric_names)=%d' % (len(metric_names)))
        metric_values = np.array(metric_values, dtype='float')
        TIVproxy_columns = []
        for TIVproxy_roi in self.TIVproxy_ROIs:
            TIVproxy_columns.append(metric_names.index('asegvolume_%s' % TIVproxy_roi))
        metric_values = np.hstack((metric_values, np.nansum(metric_values[:, TIVproxy_columns], axis=1)[:, None]))
        metric_names.append('asegvolume_EstimatedTotalIntraCranialVol')
        logging.PRINT('Index of asegvolume_EstimatedTotalIntraCranialVol is %d' % (metric_names.index('asegvolume_EstimatedTotalIntraCranialVol')))
        # Make it 2D
        if len(metric_values.shape) == 1:
            metric_values = metric_values[None, :]
        for hemi in ['rh', 'lh']:
            for metric in self.parc35_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.aparc_stats_%s.txt' % (hemi, metric))
                if os.path.exists(stats2table_file):
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['aparc_%s' % m for m in dict_reader.fieldnames[1:]]  # skip <hemi>.aparc.<metric> entry at start
                        filtering_idxs = np.arange(len(tmp_metric_names))
                        if metric_include is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name in metric_include]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name in metric_include])
                            tmp_metric_names = filtered_metric_names
                        if metric_exclude is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name not in metric_exclude]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name not in metric_exclude])
                            tmp_metric_names = filtered_metric_names
                        if len(tmp_metric_names) > 0:
                            for row in dict_reader:
                                new_subj = list(row.values())[0]
                                if new_subj != ID:
                                    logging.WARNING("Subject %s is labeled as %s in freesurfer output (%s.aparc_stats_%s.txt)."
                                                    " The subject was discarded." % (ID, new_subj, hemi, metric))
                                else:
                                    metric_values = np.hstack((metric_values, np.array(list(row.values())[1:], dtype='float')[filtering_idxs][None, :]))  # skip participant ID
                                    metric_names += tmp_metric_names
            for metric in self.parc75_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.aparca2009s_stats_%s.txt' % (hemi, metric))
                if os.path.exists(stats2table_file):
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['aparc.a2009s_%s' % m for m in dict_reader.fieldnames[1:]]  # skip <hemi>.aparc.a2009s.<metric> entry
                        filtering_idxs = np.arange(len(tmp_metric_names))
                        if metric_include is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name in metric_include]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name in metric_include])
                            tmp_metric_names = filtered_metric_names
                        if metric_exclude is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name not in metric_exclude]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name not in metric_exclude])
                            tmp_metric_names = filtered_metric_names
                        if len(tmp_metric_names) > 0:
                            for row in dict_reader:
                                new_subj = list(row.values())[0]
                                if new_subj != ID:
                                    logging.WARNING(
                                        "Subject %s is labeled as %s in freesurfer output (%s.aparca2009s_stats_%s.txt)."
                                        " The subject was discarded." % (ID, new_subj, hemi, metric))
                                else:
                                    metric_values = np.hstack((metric_values, np.array(list(row.values())[1:], dtype='float')[filtering_idxs][None, :]))  # skip participant ID
                                    metric_names += tmp_metric_names
        if len(metric_values) == 0:
            logging.WARNING('No metric values were found for subject %s, will likely have only NaNs' % ID)
        return metric_values, metric_names

    def load_subject_metrics_from_stats2tableFolder(self, stats2table_folder, subjSesAcq_list, metric_include=None, metric_exclude=None):
        """
        Loads tables with multiple subjects. Requires SOM.proc_pipeline.load_from_ID to be set to False for this to be
        called in proc2metric(). The function loops through files defined in seg_metrics, parc35_metrics and parc75_metrics,
        and loads all subjects found in the files, reorganizing order according to subjSesAcq_list. Metric_include and
        metric_exclude lists can be used to filter metrics by name.

        :param stats2table_folder: path to folder containing all the stats2metric files, with all subjects in each file.
        :type stats2table_folder: string
        :param subjSesAcq_list: list of IDs, used to reorder loaded metrics according to the desired order. Should also
                                contain all subjects in the metric files, as loading a subject that is not in the list
                                will probably trigger a ValueError when trying to find the index in this list.
        :type subjSesAcq_list: list
        :param metric_include: list of metric names to keep
        :type metric_include: list
        :param metric_exclude: list of metric names to remove
        :type metric_exclude: list
        """

        metric_names = []
        metric_values = np.zeros((len(subjSesAcq_list), 1))  # Initialize array, to be trimmed out after 1st set of metrics
        # Loop over first table to get ID order, to be used to reorder tables without IDs according to subjSesAcq_list
        loaded_IDs = []
        i = 0
        aseg_file = os.path.join(stats2table_folder, 'aseg_stats_%s.txt' % self.seg_metrics[i])
        while not os.path.exists(aseg_file):
            i += 1
            aseg_file = os.path.join(stats2table_folder, 'aseg_stats_%s.txt' % self.seg_metrics[i])
        with open(aseg_file, 'r') as f:
            for row in DictReader(f, delimiter='\t'):
                loaded_IDs.append(list(row.values())[0])
        logging.PRINT('stats2table_folder=%s' % stats2table_folder)
        for metric in self.seg_metrics:  # Loop through tables
            aseg_file = os.path.join(stats2table_folder, 'aseg_stats_%s.txt' % metric)
            logging.PRINT('aseg_file=%s' % aseg_file)
            if os.path.exists(aseg_file):
                logging.PRINT('aseg_file exists, loading it...')
                with open(aseg_file, 'r') as f:
                    dict_reader = DictReader(f, delimiter='\t')
                    tmp_metric_names = ['aseg%s_%s' % (metric, m) for m in dict_reader.fieldnames[1:]]  # skip first field
                    filtering_idxs = np.arange(len(tmp_metric_names))
                    if metric_include is not None:
                        filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name in metric_include]
                        filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name in metric_include])
                        tmp_metric_names = filtered_metric_names
                    if metric_exclude is not None:
                        filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name not in metric_exclude]
                        filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name not in metric_exclude])
                        tmp_metric_names = filtered_metric_names
                    if len(tmp_metric_names) > 0:
                        table_to_fill = np.full((len(subjSesAcq_list), len(tmp_metric_names)), np.nan)
                        for i, row in enumerate(dict_reader):
                            logging.PRINT('Loaded metrics are %s' % (', '.join(list(row.keys()))))
                            ID = list(row.values())[0]  # First value is subject ID, no common key between tables
                            if ID != loaded_IDs[i]:
                                logging.WARNING("%s in aseg_stats_%s.txt (line %d) doesn't match previous order !" % (ID, metric, i))
                            if ID not in subjSesAcq_list:
                                logging.PRINT('Entered 1st if statement')
                                logging.WARNING("Subject %s in freesurfer output (aseg_stats_%s.txt) is not in subjSesAcq_list."
                                                " The subject was discarded." % (ID, metric))
                            else:
                                logging.PRINT('Entered 2nd if statement')
                                table_to_fill[subjSesAcq_list.index(ID), :] = np.array(list(row.values())[1:])[filtering_idxs]  # skip ID value, fill row in table_to_fill that corresponds to subjSesAcq_list
                            logging.PRINT('Exiting if statement')
                        metric_names += tmp_metric_names
                        metric_values = np.hstack((metric_values, table_to_fill))  # append all values in current table file to metric_values
        logging.PRINT('Loaded metric_names are %s' % (', '.join(metric_names)))
        logging.PRINT('len(metric_names)=%d' % (len(metric_names)))
        TIVproxy_columns = []
        for TIVproxy_roi in self.TIVproxy_ROIs:
            TIVproxy_columns.append(metric_names.index('asegvolume_%s' % TIVproxy_roi))
        metric_values = np.hstack((metric_values, np.nansum(metric_values[:, TIVproxy_columns], axis=1)[:, None]))
        metric_names.append('asegvolume_EstimatedTotalIntraCranialVol')
        logging.PRINT('Index of asegvolume_EstimatedTotalIntraCranialVol is %d' % (metric_names.index('asegvolume_EstimatedTotalIntraCranialVol')))
        metric_values = metric_values[:, 1:]  # Trim out 1st column used to initialize array
        for hemi in ['rh', 'lh']:
            for metric in self.parc35_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.aparc_stats_%s.txt' % (hemi, metric))
                if os.path.exists(stats2table_file):
                    logging.PRINT("File %s exists, opening it" % stats2table_file)
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['aparc_%s' % m for m in dict_reader.fieldnames[1:]]
                        logging.PRINT("tmp_metric_names has length %d" % (len(tmp_metric_names)))
                        filtering_idxs = np.arange(len(tmp_metric_names))
                        if metric_include is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name in metric_include]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name in metric_include])
                            tmp_metric_names = filtered_metric_names
                        if metric_exclude is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name not in metric_exclude]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name not in metric_exclude])
                            tmp_metric_names = filtered_metric_names
                        if len(tmp_metric_names) > 0:
                            table_to_fill = np.full((len(subjSesAcq_list), len(tmp_metric_names)), np.nan)
                            logging.PRINT("table_to_fill has shape %dx%d" % (table_to_fill.shape[0], table_to_fill.shape[1]))
                            for i, row in enumerate(dict_reader):
                                ID = list(row.values())[0]  # First value is subject ID, no common key between tables
                                logging.PRINT("Subject ID is %s" % ID)
                                if ID != loaded_IDs[i]:
                                    logging.WARNING("%s in %s.aparc_stats_%s.txt (line %d) doesn't match previous order !" % (ID, hemi, metric, i))
                                if ID not in subjSesAcq_list:
                                    logging.WARNING("Subject %s in freesurfer output (%s.aparc_stats_%s.txt) is not in"
                                                    " subjSesAcq_list. The subject was discarded." % (ID, hemi, metric))
                                else:
                                    logging.PRINT("list(row.values())[1:] has %d values" % (len(list(row.values())[1:])))
                                    logging.PRINT("filtering_idxs has length %d" % len(filtering_idxs))
                                    logging.PRINT("Subject idx is %d" % (subjSesAcq_list.index(ID)))
                                    table_to_fill[subjSesAcq_list.index(ID), :] = np.array(list(row.values())[1:])[filtering_idxs]  # skip ID value, fill row in table_to_fill that corresponds to subjSesAcq_list
                            metric_names += tmp_metric_names
                            metric_values = np.hstack((metric_values, table_to_fill))
            for metric in self.parc75_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.aparca2009s_stats_%s.txt' % (hemi, metric))
                if os.path.exists(stats2table_file):
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['aparc.a2009s_%s' % m for m in dict_reader.fieldnames[1:]]
                        filtering_idxs = np.arange(len(tmp_metric_names))
                        if metric_include is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name in metric_include]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name in metric_include])
                            tmp_metric_names = filtered_metric_names
                        if metric_exclude is not None:
                            filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name not in metric_exclude]
                            filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name not in metric_exclude])
                            tmp_metric_names = filtered_metric_names
                        if len(tmp_metric_names) > 0:
                            table_to_fill = np.full((len(subjSesAcq_list), len(tmp_metric_names)), np.nan)
                            for i, row in enumerate(dict_reader):
                                ID = list(row.values())[0]  # First value is subject ID, no common key between tables
                                if ID != loaded_IDs[i]:
                                    logging.WARNING("%s in %s.aparca2009s_stats_%s.txt (line %d) doesn't match previous order !" % (ID, hemi, metric, i))
                                if ID not in subjSesAcq_list:
                                    logging.WARNING(
                                        "Subject %s in freesurfer output (%s.aparca2009s_stats_%s.txt) is not in"
                                        " subjSesAcq_list. The subject was discarded." % (ID, hemi, metric))
                                else:
                                    table_to_fill[subjSesAcq_list.index(ID), :] = np.array(list(row.values())[1:])[filtering_idxs]  # skip ID value, fill row in table_to_fill that corresponds to subjSesAcq_list
                            metric_names += tmp_metric_names
                            metric_values = np.hstack((metric_values, table_to_fill))
        return metric_values, metric_names

    def generate_location_plots(self, output_folder):
        my_env = os.environ.copy()
        if 'FREESURFER_HOME' not in my_env.keys():
            logging.ERROR('$FREESURFER_HOME environment variable not found. Set $FREESURFER_HOME to your FREESURFER '
                          'installation path before running ScanOMetrics.')
        fs_home = my_env['FREESURFER_HOME']
        """dk_ROIs = [['bankssts', 'lateral'],
                   ['caudalanteriorcingulate', 'medial'],
                   ['caudalmiddlefrontal', 'lateral'],
                   ['cuneus', 'medial'],
                   ['entorhinal', 'medial'],
                   ['fusiform', 'medial'],
                   ['inferiorparietal', 'lateral'],
                   ['inferiortemporal', 'lateral'],
                   ['isthmuscingulate', 'medial'],
                   ['lateraloccipital', 'lateral'],
                   ['lateralorbitofrontal', 'lateral'],
                   ['lingual', 'medial'],
                   ['medialorbitofrontal', 'medial'],
                   ['middletemporal', 'lateral'],
                   ['parahippocampal', 'medial'],
                   ['paracentral', 'medial'],
                   ['parsopercularis', 'lateral'],
                   ['parsorbitalis', 'lateral'],
                   ['parstriangularis', 'lateral'],
                   ['pericalcarine', 'medial'],
                   ['postcentral', 'lateral'],
                   ['posteriorcingulate', 'medial'],
                   ['precentral', 'lateral'],
                   ['precuneus', 'medial'],
                   ['rostralanteriorcingulate', 'medial'],
                   ['rostralmiddlefrontal', 'lateral'],
                   ['superiorfrontal', 'lateral'],
                   ['superiorparietal', 'lateral'],
                   ['superiortemporal', 'lateral'],
                   ['supramarginal', 'lateral'],
                   ['frontalpole', 'medial'],
                   ['temporalpole', 'medial'],
                   ['transversetemporal', 'lateral'],
                   ['insula', 'lateral']]"""
        hemis = ['lh', 'rh']
        for roi, roi_view in [self.atlases[self.fs_version]['DesikanKilliany']['ROIs'], self.atlases[self.fs_version]['DesikanKilliany']['views']]:
            for hemi in hemis:
                orig_annot = read_annot(os.path.join(fs_home, 'subjects', 'fsaverage', 'label', '%s.aparc.annot' % hemi))
                new_vertices = np.ones(orig_annot[0].shape, dtype='int32') * -1
                new_colors = [[0, 0, 0, 0, 23234234]]
                new_names = [b'unknown']
                roi_idx = orig_annot[2].index(roi.encode('ascii'))
                new_vertices[orig_annot[0] == roi_idx] = 1
                new_colors.append([0, 0, 255, 0, 2647065])
                new_names.append(orig_annot[2][roi_idx])
                roi_annot_file = '/tmp/tmp_roi.annot'
                write_annot(roi_annot_file, new_vertices, np.array(new_colors, dtype='int32'), new_names)
                if hemi == 'lh' and roi_view == 'lateral':
                    rotation = '0'
                elif hemi == 'lh' and roi_view == 'medial':
                    rotation = '180'
                elif hemi == 'rh' and roi_view == 'lateral':
                    rotation = '180'
                elif hemi == 'rh' and roi_view == 'medial':
                    rotation = '0'
                cmd = ['freeview',
                       '-f', os.path.join(fs_home, 'subjects', 'fsaverage', 'surf', '%s.pial' % hemi) +
                           ':edgethickness=0'+
                           ':curvature_method=off' +
                           ':annot=%s' % roi_annot_file,
                       # '-view', view,
                       # '-zoom', '1.5',
                       '-layout', '1', '-viewport', '3d',
                       '-cam', 'Dolly', '1.5', 'azimuth', rotation,
                       '-ss', os.path.join(output_folder, 'location_%s_%s.jpeg' % (roi, hemi)), '0.5']
                # logging.PRINT(' '.join(cmd))
                subprocess.call(cmd)

        # Destrieux ROIs
        """destrieux_ROIs = [['G_and_S_frontomargin', 'lateral'],
                          ['G_and_S_occipital_inf', 'lateral'],
                          ['G_and_S_paracentral', 'medial'],
                          ['G_and_S_subcentral', 'lateral'],
                          ['G_and_S_transv_frontopol', 'lateral'],
                          ['G_and_S_cingul-Ant', 'medial'],
                          ['G_and_S_cingul-Mid-Ant', 'medial'],
                          ['G_and_S_cingul-Mid-Post', 'medial'],
                          ['G_cingul-Post-dorsal', 'medial'],
                          ['G_cingul-Post-ventral', 'medial'],
                          ['G_cuneus', 'medial'],
                          ['G_front_inf-Opercular', 'lateral'],
                          ['G_front_inf-Orbital', 'lateral'],
                          ['G_front_inf-Triangul', 'lateral'],
                          ['G_front_middle', 'lateral'],
                          ['G_front_sup', 'medial'],
                          ['G_Ins_lg_and_S_cent_ins', 'lateral'],
                          ['G_insular_short', 'lateral'],
                          ['G_occipital_middle', 'lateral'],
                          ['G_occipital_sup', 'medial'],
                          ['G_oc-temp_lat-fusifor', 'medial'],
                          ['G_oc-temp_med-Lingual', 'medial'],
                          ['G_oc-temp_med-Parahip', 'medial'],
                          ['G_orbital', 'lateral'],
                          ['G_pariet_inf-Angular', 'lateral'],
                          ['G_pariet_inf-Supramar', 'lateral'],
                          ['G_parietal_sup', 'lateral'],
                          ['G_postcentral', 'lateral'],
                          ['G_precentral', 'lateral'],
                          ['G_precuneus', 'medial'],
                          ['G_rectus', 'medial'],
                          ['G_subcallosal', 'medial'],
                          ['G_temp_sup-G_T_transv', 'lateral'],
                          ['G_temp_sup-Lateral', 'lateral'],
                          ['G_temp_sup-Plan_polar', 'medial'],
                          ['G_temp_sup-Plan_tempo', 'lateral'],
                          ['G_temporal_inf', 'lateral'],
                          ['G_temporal_middle', 'lateral'],
                          ['Lat_Fis-ant-Horizont', 'lateral'],
                          ['Lat_Fis-ant-Vertical', 'lateral'],
                          ['Lat_Fis-post', 'lateral'],
                          # ['Medial_wall', 'lateral'],
                          ['Pole_occipital', 'lateral'],
                          ['Pole_temporal', 'lateral'],
                          ['S_calcarine', 'medial'],
                          ['S_central', 'lateral'],
                          ['S_cingul-Marginalis', 'medial'],
                          ['S_circular_insula_ant', 'lateral'],
                          ['S_circular_insula_inf', 'lateral'],
                          ['S_circular_insula_sup', 'lateral'],
                          ['S_collat_transv_ant', 'medial'],
                          ['S_collat_transv_post', 'medial'],
                          ['S_front_inf', 'lateral'],
                          ['S_front_middle', 'lateral'],
                          ['S_front_sup', 'lateral'],
                          ['S_interm_prim-Jensen', 'lateral'],
                          ['S_intrapariet_and_P_trans', 'lateral'],
                          ['S_oc_middle_and_Lunatus', 'lateral'],
                          ['S_oc_sup_and_transversal', 'lateral'],
                          ['S_occipital_ant', 'lateral'],
                          ['S_oc-temp_lat', 'lateral'],
                          ['S_oc-temp_med_and_Lingual', 'medial'],
                          ['S_orbital_lateral', 'lateral'],
                          ['S_orbital_med-olfact', 'lateral'],
                          ['S_orbital-H_Shaped', 'lateral'],
                          ['S_parieto_occipital', 'medial'],
                          ['S_pericallosal', 'medial'],
                          ['S_postcentral', 'lateral'],
                          ['S_precentral-inf-part', 'lateral'],
                          ['S_precentral-sup-part', 'lateral'],
                          ['S_suborbital', 'medial'],
                          ['S_subparietal', 'medial'],
                          ['S_temporal_inf', 'lateral'],
                          ['S_temporal_sup', 'lateral'],
                          ['S_temporal_transverse', 'lateral']]"""
        for roi, roi_view in [self.atlases[self.fs_version]['Destrieux']['ROIs'], self.atlases[self.fs_version]['Destrieux']['views']]:
            for hemi in hemis:
                orig_annot = read_annot(os.path.join(fs_home, 'subjects', 'fsaverage', 'label', '%s.aparc.a2009s.annot' % hemi))
                new_vertices = np.ones(orig_annot[0].shape, dtype='int32') * -1
                new_colors = [[0, 0, 0, 0, 23234234]]
                new_names = [b'unknown']
                roi_idx = orig_annot[2].index(roi.encode('ascii'))
                new_vertices[orig_annot[0] == roi_idx] = 1
                new_colors.append([0, 0, 255, 0, 2647065])
                new_names.append(orig_annot[2][roi_idx])
                roi_annot_file = '/tmp/tmp_roi.annot'
                write_annot(roi_annot_file, new_vertices, np.array(new_colors, dtype='int32'), new_names)
                if hemi == 'lh' and roi_view == 'lateral':
                    rotation = '0'
                elif hemi == 'lh' and roi_view == 'medial':
                    rotation = '180'
                elif hemi == 'rh' and roi_view == 'lateral':
                    rotation = '180'
                elif hemi == 'rh' and roi_view == 'medial':
                    rotation = '0'
                if roi[:2] == 'S_':
                    surf = os.path.join(fs_home, 'subjects', 'fsaverage', 'surf', '%s.inflated' % hemi)+':curvature_method=binary'
                else:
                    surf = os.path.join(fs_home, 'subjects', 'fsaverage', 'surf', '%s.pial' % hemi)+':curvature_method=off'
                cmd = ['freeview',
                       '-f', surf +
                       ':edgethickness=0' +
                       ':annot=%s' % roi_annot_file,
                       # '-view', view,
                       # '-zoom', '1.5',
                       '-layout', '1', '-viewport', '3d',
                       '-cam', 'Dolly', '1.5', 'azimuth', rotation,
                       '-ss', os.path.join(output_folder, 'location_%s_%s.jpeg' % (roi, hemi)), '0.5']
                # logging.PRINT(' '.join(cmd))
                subprocess.call(cmd)

        for lobe, rois in self.lobe_ROIs:
            for hemi in hemis:
                orig_annot = read_annot(os.path.join(fs_home, 'subjects', 'fsaverage', 'label', '%s.aparc.annot' % hemi))
                new_vertices = np.ones(orig_annot[0].shape, dtype='int32') * -1
                new_colors = [[0, 0, 0, 0, 23234234]]
                new_names = [b'unknown']
                for roi in rois:
                    roi_idx = orig_annot[2].index(roi.encode('ascii'))
                    new_vertices[orig_annot[0] == roi_idx] = 1
                new_colors.append([0, 0, 255, 0, 2647065])
                new_names.append(lobe.encode('ascii'))
                roi_annot_file = '/tmp/tmp_roi.annot'
                write_annot(roi_annot_file, new_vertices, np.array(new_colors, dtype='int32'), new_names)
                if hemi == 'lh':
                    for rotation, roi_view in [['0', 'lateral'], ['180', 'medial']]:
                        cmd = ['freeview',
                               '-f', os.path.join(fs_home, 'subjects', 'fsaverage', 'surf', '%s.pial' % hemi) +
                               ':edgethickness=0' +
                               ':curvature_method=off' +
                               ':annot=%s' % roi_annot_file,
                               '-layout', '1', '-viewport', '3d',
                               '-cam', 'Dolly', '1.5', 'azimuth', rotation,
                               '-ss', os.path.join(output_folder, 'location_%s_%s.jpeg' % (lobe, hemi)), '0.5']
                        # logging.PRINT(' '.join(cmd))
                        subprocess.call(cmd)
                elif hemi == 'rh':
                    for rotation, roi_view in [['0', 'medial'], ['180', 'lateral']]:
                        cmd = ['freeview',
                               '-f', os.path.join(fs_home, 'subjects', 'fsaverage', 'surf', '%s.pial' % hemi) +
                               ':edgethickness=0' +
                               ':curvature_method=off' +
                               ':annot=%s' % roi_annot_file,
                               '-layout', '1', '-viewport', '3d',
                               '-cam', 'Dolly', '1.5', 'azimuth', rotation,
                               '-ss', os.path.join(output_folder, 'location_%s_%s.jpeg' % (lobe, hemi)), '0.5']
                        # logging.PRINT(' '.join(cmd))
                        subprocess.call(cmd)

        # Prepare baseline data for aseg ROIs
        import nibabel as nib
        aseg = nib.load(os.path.join(fs_home, 'subjects', 'bert', 'mri', 'aparc+aseg.mgz')).get_fdata()
        affine = nib.load(os.path.join(fs_home, 'subjects', 'bert', 'mri', 'aparc+aseg.mgz')).affine
        background = np.zeros(aseg.shape)
        background[aseg > 0] = 0.5
        for label in [2, 41, 7, 46]:
            background[aseg == label] = 1
        for label in [4, 43, 5, 44, 14, 15]:
            background[aseg == label] = 0
        nib.save(nib.Nifti1Image(background, affine), '/tmp/background.mgz')
        subcort_ROIs = [['Lateral-Ventricle', [4, 43]],
                        ['Inf-Lat-Vent', [5, 44]],
                        ['3rd-Ventricle', [14]],
                        ['4th-Ventricle', [15]],
                        # ['5th-Ventricle',[]], <- not mentioned in FS LUT...
                        ['Thalamus-Proper', [10, 49]],
                        ['Accumbens-area', [26, 58]],
                        ['Caudate', [11, 50]],
                        ['Putamen', [12, 51]],
                        ['Pallidum', [13, 52]],
                        ['Amygdala', [18, 54]],
                        ['Hippocampus', [17, 53]],
                        ['Brain-Stem', [16]],
                        ['VentralDC', [28, 60]],
                        ['choroid-plexus', [31, 63]]]
        for roi_name, roi_labels in subcort_ROIs:
            for i_hemi, roi_label in enumerate(roi_labels):
                roi = np.zeros(aseg.shape)
                roi[aseg == roi_label] = 1
                nib.save(nib.Nifti1Image(roi, affine), '/tmp/roi.mgz')
                for viewport, slice in [['x', 'sagittal'], ['y', 'coronal'], ['z', 'axial']]:
                    if len(roi_labels) == 1:
                        cmd = ['freeview', '-layout', '1', '-v', '/tmp/background.mgz',
                               '-v', os.path.join(fs_home,'subjects', 'bert', 'mri', 'aparc+aseg.mgz')+':outline=true'
                                    + ':smoothed=true:colormap=nih:colorscale=0,1',
                               '-v', '/tmp/roi.mgz:colormap=jet:colorscale=0,1:structure=1',
                               '-zoom', '1.2',
                               '-viewport', viewport, '-ss', os.path.join(output_folder, 'location_%s_%s.jpeg' % (roi_name, slice)),
                               ]
                        subprocess(cmd)
                    else:
                        cmd = ['freeview', '-layout', '1', '-v', '/tmp/background.mgz',
                               '-v', os.path.join(fs_home, 'subjects', 'bert', 'mri', 'aparc+aseg.mgz') + ':outline=true'
                               + ':smoothed=true:colormap=nih:colorscale=0,1',
                               '-v', '/tmp/roi.mgz:colormap=jet:colorscale=0,1:structure=1',
                               '-zoom', '1.2',
                               '-viewport', viewport, '-ss', os.path.join(output_folder, 'location_%s_%s_sagittal.jpeg' % (roi_name, hemis[i_hemi], slice)),
                               ]
                        subprocess(cmd)

    def pial2outer(self, subjects_dir, subject, atlas, hemi):
        """
        Labels outer surface with closest pial labels.
        :param subjects_dir: Freesurfer subjet directory
        :param subject: subject ID (folder name in subjects_dir)
        :param atlas: atlas to use (DesikanKilliany or Destrieux)
        :param hemi: hemispheres to process ('lh' or 'rh')
        """
        # Compute neighbouring distances for all vertices in hemisphere outer surface
        with open(os.path.join(subjects_dir, subject, 'surf', '%s.pial.asc' % hemi), 'r') as f_in:
            next(f_in)
            N_vtx_pial = int(next(f_in).split(' ')[0])
        coord_pial = np.loadtxt(os.path.join(subjects_dir, subject, 'surf', '%s.pial.asc' % hemi), skiprows=2)
        coord_pial = coord_pial[:N_vtx_pial, :3]
        with open(os.path.join(subjects_dir, subject, 'surf', '%s.pial_outer_smoothed.asc' % hemi), 'r') as f_in:
            next(f_in)
            N_vtx_outer = int(next(f_in).split(' ')[0])
        coord_outer = np.loadtxt(os.path.join(subjects_dir, subject, 'surf', '%s.pial_outer_smoothed.asc' % hemi), skiprows=2)
        coord_outer = coord_outer[:N_vtx_outer, :3]
        dist_neighbour = np.zeros(N_vtx_outer)
        vtx_neighbour = np.zeros(N_vtx_outer, dtype='int')  # indexes of outer vertices' closest pial point
        for n in range(N_vtx_outer):
            # all_dist = np.sqrt(np.sum((np.tile(coord_outer[n, :], (coord_pial.shape[0],1)) - coord_pial)**2, axis=1))  # about 6 times slower than cdist
            all_dist = cdist(coord_outer[n, :][None, :], coord_pial)
            dist_neighbour[n] = all_dist.min()
            vtx_neighbour[n] = np.argmin(all_dist)
        # Loop through specified ROIs
        for roi in self.atlases[self.fs_version][atlas]['ROIs']:
            label_file = os.path.join(subjects_dir, subject, 'label_%s' % atlas, 'pial', '%s.%s.label' % (hemi, roi))
            if os.path.exists(label_file):
                vtx_chosen_ROI = np.loadtxt(label_file, skiprows=2)[:, 0].astype(int)  # list of pial indexes corresponding to ROI
                idx = np.argwhere(np.isin(vtx_neighbour, vtx_chosen_ROI))[:, 0]  # indexes of outer surface vertices corresponding to ROI
                N = len(idx)
                if N > 0:
                    with open(os.path.join(subjects_dir, subject, 'label_%s' % atlas, 'pial_outer_smoothed', '%s.%s.label' % (hemi, roi)), 'w') as f_out:
                        f_out.write('#!ascii label  , from subject %s vox2ras=TkReg\n' % subject)
                        f_out.write('%d\n' % N)
                        for n in range(N):
                            f_out.write('%d\t%f\t%f\t%f\t%f\n' % (idx[n], coord_outer[idx[n], 0], coord_outer[idx[n], 1], coord_outer[idx[n], 2], 0))

    def calc_area_gauscurv(self, subject, atlas, hemi):
        surfnames = ['pial', 'pial_outer_smoothed']
        for surfname in surfnames:
            surf_file = os.path.join(self.subjects_dir, subject, 'surf', '%s.%s.asc' % (hemi, surfname))
            with open(surf_file, 'r') as f_in:
                next(f_in)
                line = next(f_in).rstrip()
            N_vtx, N_tri = [int(n) for n in line.split(' ')]
            coord_triangles = np.loadtxt(surf_file, skiprows=2)
            coord = coord_triangles[:N_vtx, :3]
            triangles = coord_triangles[N_vtx:N_vtx+N_tri, :3].astype('int')
            vtx_convHull = ConvexHull(coord, qhull_options='Qt').simplices
            gauscurv_file = os.path.join(self.subjects_dir, subject, 'surf', '%s.%s_gauscurv.asc' % (hemi, surfname))
            vtx_coord_curv = np.loadtxt(gauscurv_file)
            curv = vtx_coord_curv[:, 4]
            ROIs = self.atlases[self.fs_version][atlas]['ROIs']
            surf_area = np.full(len(ROIs), np.nan)
            gauscurv = np.full(len(ROIs), np.nan)
            for i, roi in enumerate(ROIs):
                label_file = os.path.join(self.subjects_dir, subject, 'label_%s' % atlas, surfname, '%s.%s.label' % (hemi, roi))
                if os.path.exists(label_file):
                    vtx_chosen_ROI = np.loadtxt(label_file, skiprows=2)
                    if len(vtx_chosen_ROI.shape) == 1:
                        vtx_chosen_ROI = vtx_chosen_ROI[:, None]  # Can happen that only one pial point is part of the roi..
                    vtx_chosen_ROI = vtx_chosen_ROI[:, 0].astype(int)  # list of pial indexes corresponding to ROI
                    idx1 = np.argwhere(np.isin(triangles[:, 0], vtx_chosen_ROI))
                    idx2 = np.argwhere(np.isin(triangles[:, 1], vtx_chosen_ROI))
                    idx3 = np.argwhere(np.isin(triangles[:, 2], vtx_chosen_ROI))
                    idx_sel = np.vstack((idx1, idx2, idx3))[:, 0]
                    N_sel = len(idx_sel)
                    triangles_sel = triangles[idx_sel, :]
                    areas = np.zeros(N_sel)
                    gausdens = np.zeros(N_sel)
                    for n in range(N_sel):
                        A = np.sqrt(((coord[triangles_sel[n, 0], :] - coord[triangles_sel[n, 1], :])**2).sum())
                        B = np.sqrt(((coord[triangles_sel[n, 0], :] - coord[triangles_sel[n, 2], :])**2).sum())
                        C = np.sqrt(((coord[triangles_sel[n, 1], :] - coord[triangles_sel[n, 2], :])**2).sum())
                        S = (A+B+C)/2.0
                        areas[n] = np.sqrt(S*(S-A)*(S-B)*(S-C))
                        if triangles_sel[n, 0] in vtx_convHull:
                            gausdens[n] += curv[triangles_sel[n, 0]]
                        if triangles_sel[n, 1] in vtx_convHull:
                            gausdens[n] += curv[triangles_sel[n, 1]]
                        if triangles_sel[n, 2] in vtx_convHull:
                            gausdens[n] += curv[triangles_sel[n, 2]]
                        gausdens[n] /= 3
                    surf_area[i] = areas.sum()/3
                    gauscurv[i] = np.sum(areas*gausdens)/3

            if atlas == 'DesikanKilliany':
                out_file = os.path.join(self.subjects_dir, subject, 'stats2table', '%s.aparc_stats_area_%s.txt' % (hemi, surfname))
            elif atlas == 'Destrieux':
                out_file = os.path.join(self.subjects_dir, subject, 'stats2table', '%s.aparca2009s_stats_area_%s.txt' % (hemi, surfname))
            with open(out_file, 'w') as f_out:
                if atlas == 'DesikanKilliany':
                    f_out.write('%s.aparc.area' % (hemi))
                elif atlas == 'Destrieux':
                    f_out.write('%s.aparc.a2009s.area' % (hemi))
                for i, roi in enumerate(ROIs):
                    f_out.write('\t%s_%s_area_%s' % (hemi, roi, surfname))
                f_out.write('\t%s_TotalSurfArea_%s\n' % (hemi, surfname))
                f_out.write(subject)
                for i, roi in enumerate(ROIs):
                    f_out.write('\t%1.6f' % surf_area[i])
                f_out.write('\t%1.6f\n' % np.nansum(surf_area))

            if atlas == 'DesikanKilliany':
                out_file = os.path.join(self.subjects_dir, subject, 'stats2table', '%s.aparc_stats_gauscurv_%s_FS.txt' % (hemi, surfname))
            elif atlas == 'Destrieux':
                out_file = os.path.join(self.subjects_dir, subject, 'stats2table', '%s.aparca2009s_stats_gauscurv_%s_FS.txt' % (hemi, surfname))
            with open(out_file, 'w') as f_out:
                if atlas == 'DesikanKilliany':
                    f_out.write('%s.aparc.gauscurv' % hemi)
                elif atlas == 'Destrieux':
                    f_out.write('%s.aparc.a2009s.gauscurv' % hemi)
                for i, roi in enumerate(ROIs):
                    f_out.write('\t%s_%s_gauscurv_%s_FS' % (hemi, roi, surfname))
                f_out.write('\t%s_TotalGausCurv_%s_FS\n' % (hemi, surfname))
                f_out.write(subject)
                for i, roi in enumerate(ROIs):
                    f_out.write('\t%1.6f' % gauscurv[i])
                f_out.write('\t%1.6f\n' % np.nansum(gauscurv))

    def get_subjSesAcq_array(self, subjects):
        """
        Loops through a subjects dictionary and returns array of combinations of subject ID, session ID, as well as the
        acquisition label. Usually used to retrieve linear IDs used for processing individual scans and recover the
        corresponding metrics in the same order as covariate_value arrays.

        :param subjects: dictionary with subject IDs, session IDs, and acquisition labels.
        :type subjects: dictionary
        """
        subject_IDmergedSes = []
        for subj_id in subjects.keys():
            for ses_id in subjects[subj_id].keys():
                for acq_id in subjects[subj_id][ses_id].keys():
                    """if ses_id == '':
                        subject_IDmergedSes.append('%s_%s' % (subj_id, acq_id))
                    else:
                        subject_IDmergedSes.append('%s_%s_%s' % (subj_id, ses_id, acq_id))"""
                    subjSesAcq_id = [subj_id, ses_id, acq_id]
                    while("" in subjSesAcq_id):
                        subjSesAcq_id.remove("")
                    subject_IDmergedSes.append(self.subjSesAcq_delimiter.join(subjSesAcq_id))
        return subject_IDmergedSes

    def get_subjSesAcq_row(self, subjects, subject_id, session_id, acq_label):
        """Quick and dirty way of getting row index for measured_metrics and covariate_values"""
        subject_IDmergedSes = self.get_subjSesAcq_array(subjects)
        """if session_id == '':
            subjSesAcq_id = '%s_%s' % (subject_id, acq_label)
        else:
            subjSesAcq_id = '%s_%s_%s' % (subject_id, session_id, acq_label)"""
        subjSesAcq_id = [subject_id, session_id, acq_label]
        while("" in subjSesAcq_id):
            subjSesAcq_id.remove("")
        subjSesAcq_id = self.subjSesAcq_delimiter.join(subjSesAcq_id)
        return subject_IDmergedSes.index(subjSesAcq_id)

    def get_subjSesAcq_T1s(self, subjects):
        subjSesAcq_list = []
        subjSesAcqs_T1s = []
        for subj_id in subjects.keys():
            for ses_id in subjects[subj_id].keys():
                for acq_label in subjects[subj_id][ses_id].keys():
                    subjSesAcq_id = [subj_id, ses_id, acq_label]
                    while("" in subjSesAcq_id):
                        subjSesAcq_id.remove("")
                    subjSesAcq_id = self.subjSesAcq_delimiter.join(subjSesAcq_id)
                    T1_file = os.path.join(self.bids_database, subj_id, ses_id, 'anat', '%s.nii.gz' % subjSesAcq_id)
                    if os.path.exists(T1_file):
                        subjSesAcq_list.append(subjSesAcq_id)
                        subjSesAcqs_T1s.append(T1_file)
                    else:
                        logging.WARNING("""File %s not found: batch-dl+direct will skip subject %s, you might want
                        to process this subject manually.""" % (T1_file, subjSesAcq_id))
        return subjSesAcq_list, subjSesAcqs_T1s
