"""
Wrapper for freesurfer scripts to be run on <subjects_dir>/<subjid>. Includes pipeline template parameters as
proc_pipeline_name, and methods as run(), and proc2metric().
"""


from shutil import rmtree, which
from glob import glob
import subprocess
import numpy as np
import os
import re
from scanometrics.utils import logging
from csv import DictReader
from numpy.polynomial import Polynomial as poly_model
from nibabel.freesurfer.io import read_annot, write_annot
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull
from multiprocessing import Process, cpu_count
from packaging.version import Version


aparc_code = {'DesikanKilliany': 'aparc', 'Destrieux': 'aparca2009s'}

class proc_pipeline:

    def __init__(self, bids_database, compute_lgi=False, load_from_ID=True, ses_delimiter="_", acq_delimiter="_", atlas="DesikanKilliany"):
        self.bids_database = bids_database
        try:
            regex = r'(\d+\.\d+\.\d+)'
            self.version = re.search(regex, subprocess.check_output(["recon-all", "-version"]).decode("utf-8")).group(1)
        except FileNotFoundError:
            logging.ERROR("Path to `recon-all` not found, make sure to set FREESURFER_HOME to the appropriate path, and source $FREESURFER_HOME/SetUpFreeSurfer.sh")
        self.proc_pipeline_name = 'freesurfer_v%s_%s' % (self.version.replace('.', '-'), atlas)
        self.subjects_dir = os.path.join(bids_database, 'derivatives', self.proc_pipeline_name)
        self.load_from_ID = load_from_ID
        self.ses_delimiter = ses_delimiter
        self.acq_delimiter = acq_delimiter
        # global variables
        self.atlas = atlas
        self.atlases = {'6.0.0': {'DesikanKilliany': {'ROIs': ['bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'cuneus', 'entorhinal', 'fusiform', 'inferiorparietal', 'inferiortemporal', 'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal', 'lingual', 'medialorbitofrontal', 'middletemporal', 'parahippocampal', 'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral', 'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate', 'rostralmiddlefrontal', 'superiorfrontal', 'superiorparietal', 'superiortemporal', 'supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal', 'insula'],
                                              'views': ['lateral',                  'medial',             'lateral', 'medial',     'medial',   'medial',          'lateral',          'lateral',       'medial',              'lateral',              'lateral',  'medial',              'medial',        'lateral',          'medial',      'medial',         'lateral',       'lateral',          'lateral',        'medial',     'lateral',             'medial',    'lateral',    'medial',                   'medial',              'lateral',         'lateral',          'lateral',          'lateral',       'lateral',      'medial',       'medial',            'lateral', 'lateral']},
                                   'Destrieux': {'ROIs': ['G&S_frontomargin', 'G&S_occipital_inf', 'G&S_paracentral', 'G&S_subcentral', 'G&S_transv_frontopol', 'G&S_cingul-Ant', 'G&S_cingul-Mid-Ant', 'G&S_cingul-Mid-Post', 'G_cingul-Post-dorsal', 'G_cingul-Post-ventral', 'G_cuneus', 'G_front_inf-Opercular', 'G_front_inf-Orbital', 'G_front_inf-Triangul', 'G_front_middle', 'G_front_sup', 'G_Ins_lg&S_cent_ins', 'G_insular_short', 'G_occipital_middle', 'G_occipital_sup', 'G_oc-temp_lat-fusifor', 'G_oc-temp_med-Lingual', 'G_oc-temp_med-Parahip', 'G_orbital', 'G_pariet_inf-Angular', 'G_pariet_inf-Supramar', 'G_parietal_sup', 'G_postcentral', 'G_precentral', 'G_precuneus', 'G_rectus', 'G_subcallosal', 'G_temp_sup-G_T_transv', 'G_temp_sup-Lateral', 'G_temp_sup-Plan_polar', 'G_temp_sup-Plan_tempo', 'G_temporal_inf', 'G_temporal_middle', 'Lat_Fis-ant-Horizont', 'Lat_Fis-ant-Vertical', 'Lat_Fis-post', 'Pole_occipital', 'Pole_temporal', 'Medial_wall', 'S_calcarine', 'S_central', 'S_cingul-Marginalis', 'S_circular_insula_ant', 'S_circular_insula_inf', 'S_circular_insula_sup', 'S_collat_transv_ant', 'S_collat_transv_post', 'S_front_inf', 'S_front_middle', 'S_front_sup', 'S_interm_prim-Jensen', 'S_intrapariet&P_trans', 'S_oc_middle&Lunatus', 'S_oc_sup&transversal', 'S_occipital_ant', 'S_oc-temp_lat', 'S_oc-temp_med&Lingual', 'S_orbital_lateral', 'S_orbital_med-olfact', 'S_orbital-H_Shaped', 'S_parieto_occipital', 'S_pericallosal', 'S_postcentral', 'S_precentral-inf-part', 'S_precentral-sup-part', 'S_suborbital', 'S_subparietal', 'S_temporal_inf', 'S_temporal_sup', 'S_temporal_transverse'],
                                              'views': [        'lateral',           'lateral',          'medial',        'lateral',              'lateral',         'medial',             'medial',              'medial',               'medial',                'medial',   'medial',               'lateral',             'lateral',              'lateral',        'lateral',      'medial',             'lateral',         'lateral',            'lateral',          'medial',                'medial',                'medial',                'medial',   'lateral',              'lateral',               'lateral',        'lateral',       'lateral',      'lateral',      'medial',   'medial',        'medial',               'lateral',            'lateral',                'medial',               'lateral',        'lateral',           'lateral',              'lateral',              'lateral',      'lateral',        'lateral',       'lateral',     'medial' ,      'medial',   'lateral',              'medial',               'lateral',               'lateral',               'lateral',              'medial',               'medial',     'lateral',        'lateral',     'lateral',              'lateral',               'lateral',             'lateral',              'lateral',         'lateral',       'lateral',                'medial',           'lateral',              'lateral',            'lateral',              'medial',         'medial',       'lateral',               'lateral',               'lateral',       'medial',        'medial',        'lateral',        'lateral',               'lateral']}},
                        '6.0.1': {'DesikanKilliany': {'ROIs': ['bankssts', 'caudalanteriorcingulate', 'caudalmiddlefrontal', 'cuneus', 'entorhinal','fusiform', 'inferiorparietal', 'inferiortemporal', 'isthmuscingulate', 'lateraloccipital','lateralorbitofrontal', 'lingual', 'medialorbitofrontal', 'middletemporal', 'parahippocampal','paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine','postcentral', 'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate','rostralmiddlefrontal', 'superiorfrontal', 'superiorparietal', 'superiortemporal','supramarginal', 'frontalpole', 'temporalpole', 'transversetemporal', 'insula'],
                                              'views': ['lateral', 'medial', 'lateral', 'medial', 'medial', 'medial', 'lateral', 'lateral', 'medial','lateral', 'lateral', 'medial', 'medial', 'lateral', 'medial', 'medial', 'lateral', 'lateral','lateral', 'medial', 'lateral', 'medial', 'lateral', 'medial', 'medial', 'lateral', 'lateral','lateral', 'lateral', 'lateral', 'medial', 'medial', 'lateral', 'lateral']},
                                   'Destrieux': {'ROIs': ['G_and_S_frontomargin', 'G_and_S_occipital_inf', 'G_and_S_paracentral', 'G_and_S_subcentral', 'G_and_S_transv_frontopol', 'G_and_S_cingul-Ant', 'G_and_S_cingul-Mid-Ant', 'G_and_S_cingul-Mid-Post', 'G_cingul-Post-dorsal', 'G_cingul-Post-ventral', 'G_cuneus', 'G_front_inf-Opercular', 'G_front_inf-Orbital', 'G_front_inf-Triangul', 'G_front_middle', 'G_front_sup', 'G_Ins_lg_and_S_cent_ins', 'G_insular_short', 'G_occipital_middle', 'G_occipital_sup', 'G_oc-temp_lat-fusifor', 'G_oc-temp_med-Lingual', 'G_oc-temp_med-Parahip', 'G_orbital', 'G_pariet_inf-Angular', 'G_pariet_inf-Supramar', 'G_parietal_sup', 'G_postcentral', 'G_precentral', 'G_precuneus', 'G_rectus', 'G_subcallosal', 'G_temp_sup-G_T_transv', 'G_temp_sup-Lateral', 'G_temp_sup-Plan_polar', 'G_temp_sup-Plan_tempo', 'G_temporal_inf', 'G_temporal_middle', 'Lat_Fis-ant-Horizont', 'Lat_Fis-ant-Vertical', 'Lat_Fis-post', 'Pole_occipital', 'Pole_temporal', 'S_calcarine', 'S_central', 'S_cingul-Marginalis', 'S_circular_insula_ant', 'S_circular_insula_inf', 'S_circular_insula_sup', 'S_collat_transv_ant', 'S_collat_transv_post', 'S_front_inf', 'S_front_middle', 'S_front_sup', 'S_interm_prim-Jensen', 'S_intrapariet_and_P_trans', 'S_oc_middle_and_Lunatus', 'S_oc_sup_and_transversal', 'S_occipital_ant', 'S_oc-temp_lat', 'S_oc-temp_med_and_Lingual', 'S_orbital_lateral', 'S_orbital_med-olfact', 'S_orbital-H_Shaped', 'S_parieto_occipital', 'S_pericallosal', 'S_postcentral', 'S_precentral-inf-part', 'S_precentral-sup-part', 'S_suborbital', 'S_subparietal', 'S_temporal_inf', 'S_temporal_sup', 'S_temporal_transverse'],
                                             'views': ['lateral', 'lateral', 'medial', 'lateral', 'lateral', 'medial', 'medial', 'medial', 'medial', 'medial', 'medial', 'lateral', 'lateral', 'lateral', 'lateral', 'medial', 'lateral', 'lateral', 'lateral', 'medial', 'medial', 'medial', 'medial', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'medial', 'medial', 'medial', 'lateral', 'lateral', 'medial', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'medial', 'lateral', 'medial', 'lateral', 'lateral', 'lateral', 'medial', 'medial', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'lateral', 'medial', 'lateral', 'lateral', 'lateral', 'medial', 'medial', 'lateral', 'lateral', 'lateral', 'medial', 'medial', 'lateral', 'lateral', 'lateral']}},

                        }
        self.atlases['7.4.1'] = self.atlases['6.0.1']

        self.lobe_ROIs = {'DesikanKilliany':
                              [['frontalLobe', ['caudalmiddlefrontal', 'lateralorbitofrontal', 'medialorbitofrontal', 'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis', 'precentral', 'rostralmiddlefrontal', 'superiorfrontal', 'frontalpole']],
                               ['parietalLobe', ['inferiorparietal', 'postcentral', 'precuneus', 'superiorparietal', 'supramarginal']],
                               ['occipitalLobe', ['cuneus', 'lateraloccipital', 'lingual', 'pericalcarine']],
                               ['temporalLobe', ['bankssts', 'entorhinal', 'fusiform', 'inferiortemporal', 'middletemporal', 'parahippocampal', 'superiortemporal', 'temporalpole', 'transversetemporal']],
                               ['cingulateLobe', ['caudalanteriorcingulate', 'isthmuscingulate', 'posteriorcingulate', 'rostralanteriorcingulate']],
                               ['insulaLobe', ['insula']]],
                          'Destrieux':
                              [['frontalLobe', ['G&S_frontomargin', 'G&S_frontomargin', 'G&S_paracentral', 'G&S_subcentral', 'G&S_transv_frontopol', # 'G&S_cingul-Ant' and 'G&S_cingul-Mid-Post'moved to cingulate lobe
                                                'G_front_inf-Opercular', 'G_front_inf-Orbital', 'G_front_inf-Triangul', 'G_front_middle', 'G_front_sup', 'G_orbital', 'G_precentral', 'G_rectus', 'G_subcallosal', 'Lat_Fis-ant-Horizont', 'Lat_Fis-ant-Vertical', 'S_central', 'S_front_inf', 'S_front_middle', 'S_front_sup', 'S_orbital_lateral', 'S_orbital_med-olfact', # 'S_circular_insula_ant' moved to insula lobe
                                                'S_orbital-H_Shaped', 'S_precentral-inf-part', 'S_precentral-sup-part', 'S_suborbital']],
                               ['parietalLobe', ['G_occipital_sup', 'G_pariet_inf-Angular', 'G_pariet_inf-Supramar', 'G_parietal_sup', 'G_postcentral', 'G_precuneus', 'G_temp_sup-Plan_tempo', 'S_calcarine', # 'S_cingul-Marginalis' moved to cingulate lobe
                                                 'S_interm_prim-Jensen', 'S_intrapariet&P_trans', 'S_oc_sup&transversal', 'S_parieto_occipital', 'S_postcentral', 'S_subparietal']],
                               ['occipitalLobe', ['G&S_occipital_inf', 'G_cuneus', 'G_occipital_middle', 'G_oc-temp_med-Lingual', 'Pole_occipital', 'S_oc_middle&Lunatus', 'S_occipital_ant', 'S_oc-temp_med&Lingual']],
                               ['temporalLobe', ['G_oc-temp_lat-fusifor', 'G_oc-temp_med-Parahip', 'G_temp_sup-G_T_transv', 'G_temp_sup-Lateral', 'G_temp_sup-Plan_polar', 'G_temporal_inf', 'G_temporal_middle', 'Pole_temporal', 'S_collat_transv_ant', 'S_collat_transv_post', 'S_oc-temp_lat', 'S_temporal_inf', 'S_temporal_sup', 'S_temporal_transverse']], # 'S_circular_insula_inf' moved to insula lobe
                               ['cingulateLobe', ['G&S_cingul-Ant', 'G&S_cingul-Mid-Ant', 'G&S_cingul-Mid-Post', 'G_cingul-Post-dorsal', 'G_cingul-Post-ventral', 'S_cingul-Marginalis']],
                               ['insula', ['G_Ins_lg&S_cent_ins', 'G_insular_short', 'S_circular_insula_ant', 'S_circular_insula_inf', 'S_circular_insula_sup']]]
                          }

        self.set_settings(compute_lgi)

    def update_version(self, version, atlas):
        """Added function to update version, used mainly when loading a model that was trained with another version, to update software specific variables as subjects dir for dldirect"""
        self.version = version
        self.proc_pipeline_name = 'freesurfer_v%s_%s' % (self.version.replace('.', '-'), atlas)
        self.subjects_dir = os.path.join(self.bids_database, 'derivatives', self.proc_pipeline_name)

    def set_settings(self, compute_lgi):
        self.recon_all_settings = {'compute_lgi': compute_lgi}
        self.seg_metrics = ['volume']
        self.parc35_metrics = ['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind']
        self.parc75_metrics = ['volume', 'area', 'thickness', 'thicknessstd', 'meancurv', 'gauscurv', 'foldind', 'curvind']
        self.parcSTD_metrics = ['mean', 'std', 'snr']
        if self.recon_all_settings['compute_lgi']:
            self.parc35_metrics += ['area_pial', 'area_pial_outer_smoothed', 'gauscurv_pial_FS', 'gauscurv_pial_outer_smoothed_FS']
            self.parc75_metrics += ['area_pial', 'area_pial_outer_smoothed', 'gauscurv_pial_FS', 'gauscurv_pial_outer_smoothed_FS']

    def run_pipeline(self, subjects, n_threads, subject_id=None):
        """Pipeline template method overrided here. Calls freesurfer processing scripts. Single or multi-subject
        processing controlled through n_threads (default of -1 allocates all cpus available through a call to
        multiprocessing.cpu_count(). If T1_file is specified, checks that a .mgz file is present in
        <subject_dir>/<subj_id>/mri/orig, and stops with an error otherwise. Subject directory intended
        to be the actual subject folder, and the subject id is taken as the session id to follow bids structure.

        :param subj_id: participant code for subject to be analysed.
        :type subj_id: string
        """
        regex = r'(\d+\.\d+\.\d+)'
        current_version = re.search(regex, subprocess.check_output(["recon-all", "-version"]).decode("utf-8")).group(1)
        if current_version != self.version:
            logging.ERROR("Version installed on this system is %s, although processing version was set to %s (probably from loaded normative dataset). Please install required version to avoid discrepancies between normative model and newly processed data." % (current_version, self.version))
        if subject_id is None:
            logging.PRINT("Evaluating all subjects (should not happen from GUI)")
            subjSesAcq_list, T1_files = self.get_subjSesAcq_T1s(subjects)
        else:
            subjSesAcq_list, T1_files = self.get_subjSesAcq_T1s(subjects)
        if n_threads == -1:
            n_threads = cpu_count()
        for thread in [range(k, min(k + n_threads, len(subjSesAcq_list))) for k in range(0, len(subjSesAcq_list), n_threads)]:
            p_thread = []
            for i in thread:
                subjSesAcq_id = subjSesAcq_list[i]
                T1_file = T1_files[i]
                logging.PRINT('Running recon all for scan %s' % subjSesAcq_id)
                # Runs recon-all with <subjects_dir>/<participant_id>_<session_id>_<acq_label> structure. Same structure
                # is assumed for run_stats2tables(). Should stats2table folder be in derivatives/freesurfer or in
                # derivatives/scanometrics ? As most are freesurfer metrics, we choose to put them in
                # derivatives/freesurfer for now. Also makes sense to have one table per subject, to be able to load
                # different sets at runtime without having to parse a full table. Also makes it easier to add a couple
                # of subjects without having to regenerate the whole table again.
                # FR: what to do when user has processed fs data without bids like structure ? The best would be to have a
                # participants.tsv file with participant_id that matches the name of freesurfer subjects. To be handled by users
                # in that case.
                p = Process(target=self.run_recon_all, args=(subjSesAcq_id, T1_file))
                p.start()
                p_thread.append(p)
            for p in p_thread:
                p.join()

    def proc2table(self, subjects, n_threads):
        """
        Calls run_proc2table to convert Freesurfer segmentation and parcellation statistics into text files

        :param subjects: dictionary with subjects considered in the study
        :type subjects: dictionary
        :param n_threads: number of threads to use
        :type n_threads: integer
        """
        subjSesAcq_list = self.get_subjSesAcq_array(subjects)
        if n_threads == -1:
            n_threads = cpu_count()
        for thread in [range(k, min(k + n_threads, len(subjSesAcq_list))) for k in range(0, len(subjSesAcq_list), n_threads)]:
            p_thread = []
            for i in thread:
                subjSesAcq_id = subjSesAcq_list[i]
                p = Process(target=self.run_proc2table, args=(subjSesAcq_id,))
                p.start()
                p_thread.append(p)
            for p in p_thread:
                p.join()

    def run_recon_all(self, subjSesAcq_id, T1_file):
        """Wrapper for freesurfer recon-all. Takes subj_id, ses_id and acq_label as input to create a flattened folder
        structure in self.subjects_dir (i.e. bids/derivatives/freesurfer) and overwrites previous recon-all outputs.
        Computation of LGI requires matlab to be installed.

        :param subj_id: subject ID (usually taken from participants.tsv).
        :type subj_id: string
        :param ses_id: session ID (usually taken from <subj_id>_session.tsv).
        :type ses_id: string
        :param acq_label: label of scan to process (usually taken from glob() output on acq_pattern matches
        :type acq_label: string
        """
        os_env = os.environ.copy()
        os_env['SUBJECTS_DIR'] = self.subjects_dir
        if which('recon-all') is None:
            logging.ERROR("""'recon-all command not found. Please set $FREESURFER_HOME to your Freesurfer installation
            directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
        # Check MCR availability
        if not os.path.exists(os_env['FREESURFER_HOME'] + '/MCRv80'):
            logging.WARNING('Matlab Compiler Runtime (MCR) 2012b not found in $FREESURFER_HOME/MCRv80. '
                          'Make sure $FREESURFER_HOME is set to the correct location. To install MCR 2012b, follow the'
                          ' instructions in https://surfer.nmr.mgh.harvard.edu/fswiki/MatlabRuntime')
        # Check matlab availability if lgi computation is activated
        if self.recon_all_settings['compute_lgi'] and which('matlab') is None:
            logging.ERROR("'matlab' command not found, but required to compute LGI. Please add matlab location to $PATH, "
                          "or set 'compute_lgi' to False in scanometrics.processing.freesurfer.recon_all_settings")
        # Clean up output dir if not empty
        if os.path.exists(os.path.join(self.subjects_dir, subjSesAcq_id)):
            rmtree(os.path.join(self.subjects_dir, subjSesAcq_id))
        os.makedirs(os.path.join(self.subjects_dir, subjSesAcq_id, 'mri', 'orig'))
        cmd = ["mri_convert", T1_file, os.path.join(self.subjects_dir, subjSesAcq_id, 'mri', 'orig', '001.mgz')]
        p = subprocess.run(cmd, env=os_env)
        if p.returncode:
            logging.ERROR('Command failed : %s' % (' '.join(cmd)))
        if self.recon_all_settings['compute_lgi']:
            os.makedirs(os.path.join(self.subjects_dir, subjSesAcq_id, 'label_%s' % self.atlas, 'pial'))
            os.makedirs(os.path.join(self.subjects_dir, subjSesAcq_id, 'label_%s' % self.atlas, 'pial_outer_smoothed'))
        os.makedirs(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts'))
        with open(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts', 'recon-all.log'), 'w') as fs_log_file:
            # Run recon-all
            cmd = ['recon-all', '-subjid', subjSesAcq_id]
            if self.atlas == 'Destrieux':
                cmd += ['-noparcstats', '-parcstats2']
            cmd.append('-all')
            p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
            if p.returncode:
                logging.ERROR('Command recon-all failed with code %d (see %s for details).' % (p.returncode, fs_log_file))

        # Run LGI computation if selected
        if self.recon_all_settings['compute_lgi']:
            with open(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts', 'recon-all_lgi.log'), 'w') as fs_log_file:
                # Run compute_lgi script
                p = subprocess.run(['recon-all',
                                    '-subjid', subjSesAcq_id,
                                    '-localGI'], env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                if p.returncode:
                    logging.ERROR('Command recon-all localGI computation failed with code %d (see %s for details).' % (
                                   p.returncode, fs_log_file))

        # Run qcache
        with open(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts', 'recon-all_qcache.log'), 'w') as fs_log_file:
            cmd = ["recon-all", "-subject", subjSesAcq_id, '-qcache']
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
            if p.returncode:
                logging.ERROR('Command recon-all qcache failed with code %d (see %s for details).' % (p.returncode, fs_log_file))

        # Process brain-stem structures
        with open(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts', 'recon-all_brainstem.log'), 'w') as fs_log_file:
            if Version(self.version) >= Version("7.3.0"):
                cmd = ["segment_subregions", "--cross", subjSesAcq_id, "brainstem"]
            else:
                cmd = ['recon-all', '-subject', subjSesAcq_id, '-brainstem-structures']
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
            if p.returncode:
                logging.ERROR('Command recon-all brainstem-structures failed with code %d (see %s for details).' % (p.returncode, fs_log_file))

        # Process hippocampus
        with open(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts', 'recon-all_hippocampus.log'), 'w') as fs_log_file:
            if Version(self.version) >= Version("7.3.0"):
                cmd = ["segment_subregions", "--cross", subjSesAcq_id, "hippo-amygdala"]
            else:
                cmd = ["recon-all", '-subject', subjSesAcq_id, '-hippocampal-subfields-T1']
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
            if p.returncode:
                logging.ERROR('Command recon-all hippocampal-subfields failed with code %d (see %s for details).' % (p.returncode, fs_log_file))

    def run_proc2table(self, subjSesAcq_id):
        """Wrapper for asegstats2table and aparcstats2table. Assumes subject was processed with recon-all -parcstats2.
         Creates new directory <subjects_dir>/<subjSesAcq_id>/stats2table (deletes old one if present)."""
        if which('asegstats2table') is None:
            logging.ERROR("""asegstats2table command not found. Please set $FREESURFER_HOME to your Freesurfer installation
            directory and run 'source $FREESURFER_HOME/SetUpFreeSurfer.sh'.""")
        stats_folder = os.path.join(self.subjects_dir, subjSesAcq_id, 'stats2table')
        if os.path.exists(stats_folder):
            rmtree(stats_folder)
        os.makedirs(stats_folder)
        if not os.path.exists(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts')):
            os.makedirs(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts'))
        os_env = os.environ.copy()
        os_env['SUBJECTS_DIR'] = self.subjects_dir
        if Version(self.version).major == 6:  # v6 requires to run with python2.7, which is not necessarily in the conda environment
            asegstats_cmd = ['python2.7', os.path.join(os_env["FREESURFER_HOME"], 'bin', 'asegstats2table')]
            aparcstats_cmd = ['python2.7', os.path.join(os_env["FREESURFER_HOME"], 'bin', 'aparcstats2table')]
        else:
            asegstats_cmd = ['asegstats2table']
            aparcstats_cmd = ['aparcstats2table']

        # Post processing on each hemisphere
        with open(os.path.join(self.subjects_dir, subjSesAcq_id, 'scripts', 'stats2table.log'), 'w') as fs_log_file:
            for hemi in ['lh', 'rh']:
                if self.recon_all_settings['compute_lgi']:
                    # Convert surface files to ASCII for further manual processing
                    p = subprocess.run(['mris_convert',
                                        os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial' % hemi),
                                        os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial.asc' % hemi)
                                        ], env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('Conversion of %s.pial to ascii failed with code %d (see %s for details).' % (hemi, p.returncode, fs_log_file))

                    p = subprocess.run(['mris_convert',
                                        os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial-outer-smoothed' % hemi),
                                        os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial-outer-smoothed.asc' % hemi)
                                        ], env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('Conversion of %s.pial-outer-smoothed to ascii failed with code %d (see %s for details).' % (hemi, p.returncode, fs_log_file))

                    # Convert to labels
                    p = subprocess.run(['mri_annotation2label',
                                        '--subject', subjSesAcq_id,
                                        '--hemi', hemi,
                                        '--annotation', aparc_code[self.atlas],
                                        '--surface', 'pial',
                                        '--outdir', os.path.join(self.subjects_dir, subjSesAcq_id, 'label_DesikanKilliany' if self.atlas=='DesikanKilliany' else 'label_Destrieux', 'pial')
                                        ], env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('Annotation2label for %s pial failed with code %d (see %s for details).' % (aparc_code[self.atlas], p.returncode, fs_log_file))

                    # Label and annotate outer surface
                    self.pial2outer(self.subjects_dir, subjSesAcq_id, self.atlas, hemi)
                    annot_file = os.path.join(self.subjects_dir, subjSesAcq_id, 'label', '%s.aparc-outer-smoothed.annot' % hemi)
                    if os.path.exists(annot_file):
                        os.remove(annot_file)
                    cmd = ['mris_label2annot',
                           '--s', subjSesAcq_id,
                           '--h', hemi,
                           '--ctab', os.path.join(self.subjects_dir, subjSesAcq_id, 'label', 'aparc.annot.ctab'),
                           '--surf', 'pial-outer-smoothed',
                           '--ldir', os.path.join(self.subjects_dir, subjSesAcq_id, 'label_DesikanKilliany', 'pial_outer_smoothed'),
                           '--no-unknown',
                           '--a', 'aparc-outer-smoothed']
                    p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('label2annot for DesikanKilliany pial-outer-smoothed failed with code %d (see %s for'
                                      ' details).' % (p.returncode, fs_log_file))
                    annot_file = os.path.join(self.subjects_dir, subjSesAcq_id, 'label', '%s.aparc-outer-smoothed.a2009s.annot' % hemi)
                    if os.path.exists(annot_file):
                        os.remove(annot_file)
                    cmd = ['mris_label2annot',
                           '--s', subjSesAcq_id,
                           '--h', hemi,
                           '--ctab', os.path.join(self.subjects_dir, subjSesAcq_id, 'label', 'aparc.annot.a2009s.ctab'),
                           '--surf', 'pial-outer-smoothed',
                           '--ldir', os.path.join(self.subjects_dir, subjSesAcq_id, 'label_Destrieux', 'pial_outer_smoothed'),
                           '--no-unknown',
                           '--a', 'aparc-outer-smoothed.a2009s']
                    p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('label2annot for Destrieux pial-outer-smoothed failed with code %d (see %s for'
                                      ' details).' % (p.returncode, fs_log_file))

                    # Compute curvature
                    p = subprocess.run(['mris_curvature',
                                        '-w', '-a', '10', os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial' % hemi)],
                                       env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('mris_curvature for %s.pial failed with code %d (see %s for details).' % (hemi,
                                       p.returncode, fs_log_file))
                    p = subprocess.run(['mris_curvature',
                                        '-w', '-a', '10',
                                        os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial-outer-smoothed' % hemi)],
                                       env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('mris_curvature for %s.pial-outer-smoothed failed with code %d (see %s for details).'
                                      % (hemi, p.returncode, fs_log_file))

                    # Convert %s.pial_gauscurv.asc
                    subprocess.run(['mris_convert',
                                    '-c', os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial.K' % hemi),
                                    os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial' % hemi),
                                    os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial_gauscurv.asc' % hemi)],
                                   env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('Conversion of %s.pial_gauscurv.asc failed with code %d (see %s for details).'
                                      % (hemi, p.returncode, fs_log_file))
                    subprocess.run(['mris_convert',
                                    '-c', os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial-outer-smoothed.K' % hemi),
                                    os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial-outer-smoothed' % hemi),
                                    os.path.join(self.subjects_dir, subjSesAcq_id, 'surf',
                                                 '%s.pial_outer_smoothed_gauscurv.asc' % hemi)],
                                   env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR(
                            'Conversion of %s.pial_outer_smoothed_gauscurv.asc failed with code %d (see %s for details).'
                            % (hemi, p.returncode, fs_log_file))

                    # Run internal area and gauscurv computation
                    self.calc_area_gauscurv(subjSesAcq_id, atlas, hemi)

                # WM-GM PCT
                if self.atlas == 'DesikanKilliany':
                    stats_file = '%s.w-g.pct.stats' % hemi
                elif self.atlas == 'Destrieux':
                    stats_file = '%s.w-g.pct.a2009s.stats' % hemi
                cmd = ['mri_segstats', '--in', os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.w-g.pct.mgh' % hemi),
                       '--annot', subjSesAcq_id, hemi, 'aparc' if self.atlas=='DesikanKilliany' else 'aparc.a2009s',
                       '--sum', os.path.join(self.subjects_dir, subjSesAcq_id, 'stats', stats_file), '--snr']
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                if p.returncode:
                    logging.ERROR('mri_segstats failed with code %d (see %s for details).' % (p.returncode, fs_log_file))

                # LGI segstats if selected
                if self.recon_all_settings['compute_lgi']:
                    if self.atlas == 'DesikanKilliany':
                        stats_file = '%s.pial_lgi.stats' % hemi
                    elif self.atlas == 'Destrieux':
                        stats_file = '%s.pial_lgi.a2009s.stats' % hemi
                    cmd = ['mri_segstats', '--in', os.path.join(self.subjects_dir, subjSesAcq_id, 'surf', '%s.pial_lgi' % hemi),
                           '--annot', subjSesAcq_id, hemi, 'aparc' if self.atlas=='DesikanKilliany' else 'aparc.a2009s',
                           '--sum', os.path.join(self.subjects_dir, subjSesAcq_id, 'stats', stats_file)]
                    logging.PRINT('Running command: %s' % (' '.join(cmd)))
                    p = subprocess.run(cmd, env=os_env, stdout=fs_log_file, stderr=fs_log_file)
                    if p.returncode:
                        logging.ERROR('mri_segstats failed with code %d (see %s for details).' % (p.returncode, fs_log_file))

        for metric in self.seg_metrics:
            cmd = asegstats_cmd + [
                   '--subjects', subjSesAcq_id,
                   '--meas', metric,
                   '--all-segs',
                   '--tablefile', os.path.join(stats_folder, 'aseg_stats_%s.txt' % metric)
                   ]
            logging.PRINT('Running command: %s' % (' '.join(cmd)))
            p = subprocess.run(cmd, env=os_env)
            if p.returncode:
                logging.ERROR('asegstats2table failed with code %d.' % p.returncode)

        # Export brainstem structures if available (using tab as delimiter to match aseg_stats_volume for example)
        brainstem_volumes_file = glob(os.path.join(self.subjects_dir, subjSesAcq_id, 'mri', 'brainstemSsVolumes.v*.txt'))[0]
        if os.path.exists(brainstem_volumes_file):
            file_content = []
            with open(brainstem_volumes_file, 'r') as f_in:
                for line in f_in:
                    file_content.append(line.strip().split(' '))
            with open(os.path.join(stats_folder, 'brainstem_volume.txt'), 'w') as f_out:
                f_out.write('Subject\t%s\n' % ('\t'.join([tmp[0] for tmp in file_content])))
                f_out.write('%s\t%s\n' % (subjSesAcq_id, '\t'.join([tmp[1] for tmp in file_content])))

        # Export hippocampal subfields stats if file exists
        for hemi in ['rh', 'lh']:
            hipp_volume_file = glob(os.path.join(self.subjects_dir, subjSesAcq_id, 'mri', '%s.hippoSfVolumes*.txt' % hemi))[0]
            if os.path.exists(hipp_volume_file):
                file_content = []
                with open(hipp_volume_file, 'r') as f_in:
                    for line in f_in:
                        file_content.append(line.strip().split(' '))
                with open(os.path.join(stats_folder, '%s.hippoSf_volume.txt' % hemi), 'w') as f_out:
                    f_out.write('Subject\t%s\n' % ('\t'.join([hemi+'.'+tmp[0] for tmp in file_content])))
                    f_out.write('%s\t%s\n' % (subjSesAcq_id, '\t'.join([tmp[1] for tmp in file_content])))

        # Export aparc stats
        for hemi in ['rh', 'lh']:
            if self.atlas == 'DesikanKilliany':
                parc_metrics = self.parc35_metrics
            elif self.atlas == 'Destrieux':
                parc_metrics = self.parc75_metrics
            for metric in parc_metrics:
                cmd = aparcstats_cmd + [
                                "--subjects", subjSesAcq_id,
                                "--parc", 'aparc' if self.atlas=='DesikanKilliany' else 'aparc.a2009s',
                                "--hemi", hemi,
                                "--measure", metric,
                                "--parcs-from-file", os.path.join(os.path.dirname(__file__), 'resources', '%s_ROIs_select.txt' % self.atlas),
                                "--tablefile", os.path.join(stats_folder, '%s.%s_stats_%s.txt' % (hemi, aparc_code[self.atlas], metric))]
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env)
                if p.returncode:
                    logging.ERROR("aparcstats2table failed with code %d." % p.returncode)

            # GM-WM contrast needs asegstats (table format a little different)
            for metric in self.parcSTD_metrics:
                if self.atlas == 'DesikanKilliany':
                    stats_file = '%s.w-g.pct.stats' % hemi
                elif self.atlas == 'Destrieux':
                    stats_file = '%s.w-g.pct.a2009s.stats' % hemi
                cmd = asegstats_cmd + [
                    "--inputs", os.path.join(self.subjects_dir, subjSesAcq_id, 'stats', stats_file),
                    "--meas", metric,
                    "--tablefile", os.path.join(stats_folder, '%s.%s_stats_pct%s.txt' % (hemi, aparc_code[self.atlas], metric))]
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env)

            # Check if lgi values where computed
            if self.recon_all_settings['compute_lgi']:
                if self.atlas == 'DesikanKilliany':
                    stats_file = '%s.pial_lgi.stats' % hemi
                elif self.atlas == 'Destrieux':
                    stats_file = '%s.pial_lgi.a2009s.stats' % hemi
                if not os.path.exists(os.path.join(self.subjects_dir, subjSesAcq_id, 'stats', stats_file)):
                    logging.ERROR(f'{stats_file} file not found for subject {subjSesAcq_id}. Run recon-all with SOM.recon_all_settings["compute_lgi"]=True or disable lgi in scanometrics.processing.freesurfer.run_stats2table')
                cmd = asegstats_cmd + [
                    "--inputs", os.path.join(self.subjects_dir, subjSesAcq_id, 'stats', stats_file),
                    "--meas", "mean",
                    "--tablefile", os.path.join(stats_folder, '%s.%s_stats_lgimean.txt' % (hemi, aparc_code[self.atlas]))]
                logging.PRINT('Running command: %s' % (' '.join(cmd)))
                p = subprocess.run(cmd, env=os_env)

    def proc2metric(self, subjects, covariate_values, covariate_names, stats2table_folder=None, ref_rows=None,
                    ref_metric_values=None, ref_metric_names=None, ref_covariate_values=None, ref_covariate_names=None,
                    metric_include=None, metric_exclude=None):
        """Function defined for each preprocessing pipeline, in order to get all the variables from the preprocessing
        in a list of metrics, usually a table saved as text file that can then be read by load_proc_metrics(). Subjects
        with missing values get assigned a np.nan value. Scan duplicates should be checked in the future. Computes
        lobe metrics and asymmetric indexes, and adds results at end of the measured_metric matri. Normalized values
        are computed on with respect to averages in the normative dataset. Some variables might not
        exist for all subjects, such scans get assigned a np.nan value.

        :param subjects: multilevel dictionary of subjects (e.g.: SOM.subject) starting with subject IDs, then session
                         IDs, and finally acq labels. Used to generate subjSesAcq_array to analyse.
        :type subjects: dictionary
        :param covariate_values: covariate matrix for subjects being processed (e.g.: SOM.covariate_values.copy())
        :type covariate_values: numpy array
        :param covariate_names: names of covariates for the scans being processed (e.g.: SOM.covariate_names.copy())
        :type covariate_names: list of strings
        :param stats2table_folder: path to folder containing stat files with all subjects inside. Used when the flag
                                   self.load_from_ID is set to False.
        :type stats2table_folder: string
        :param ref_rows: indexes of rows in ref_metric_values and ref_covariate_values to keep before finding matches.
        :type ref_rows: numpy array
        :param ref_metric_values: matrix of metric_values to use as reference when normalizing
                                  (e.g.: SOM.measured_metrics.copy())
        :type ref_metric_values: numpy matrix
        :param ref_metric_names: list of metric names in ref_metric_values
        :type ref_metric_names: list of strings
        :param ref_covariate_values: array with covariates of reference dataset, used to find matches
        :type ref_covariate_values: numpy matrix
        :param ref_covariate_names: name of covariates in ref_covariate_values
        :type ref_covariate_names: list of strings
        """
        # Error checks
        if 'sequence' not in covariate_names or 'scanner' not in covariate_names:
            logging.ERROR("The freesurfer pipeline requires 'scanner' and 'sequence' as covariates, but were not"
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
                    if metric in ['area', 'area_pial', 'area_pial_outer_smoothed', 'volume', 'gauscurv']:
                        for ctx_label in self.lobe_ROIs[self.atlas][[self.lobe_ROIs[self.atlas][i][0] for i in range(len(self.lobe_ROIs[self.atlas]))].index(lobe)][1]:
                            metric_include.append('%s_%s_%s_%s' % (aparc_code[self.atlas], hemi, ctx_label, metric))
                    else:
                        metric_include.append('%s_%s_%s_area' % (aparc_code[self.atlas], hemi, ctx_label))
                        if metric == 'pctmean':
                            metric_include.append('%s_%s_pctmean' % (hemi, ctx_label))
                        else:
                            metric_include.append('%s_%s_%s_%s' % (aparc_code[self.atlas], hemi, ctx_label, metric))

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
            for lobe, ctx_labels in self.lobe_ROIs[self.atlas]:
                for metric in self.parc35_metrics:
                    if metric_include is not None and '%s_%s_%s' % (hemi, lobe, metric) not in metric_include:
                        continue
                    if metric_exclude is not None and '%s_%s_%s' % (hemi, lobe, metric) in metric_exclude:
                        continue
                    if metric in ['area', 'area_pial', 'area_pial_outer_smoothed', 'volume', 'gauscurv']:  # Might add 'gauscurv_pial', 'gauscurv_pial_outer_smoothed' if needed, but check if we should normalise, and by what
                        metric_names += ['%s_%s_%s' % (hemi, lobe, metric)]
                        col_idx = np.array([metric_names.index("_".join([aparc_code[self.atlas], hemi, ctx_label, metric])) for ctx_label in ctx_labels])
                        metric_values = np.hstack((metric_values, np.nansum(metric_values[:, col_idx], axis=1, keepdims=True)))
                    else:
                        col_idx = []
                        for ctx_label in ctx_labels:
                            key_name = '%s_%s_%s_area' % (aparc_code[self.atlas], hemi, ctx_label)
                            if key_name in metric_names:
                                col_idx.append(metric_names.index(key_name))
                        col_idx = np.array(col_idx)
                        if len(col_idx) > 0:
                            norm_area = metric_values[:, col_idx].copy()
                            norm_area /= norm_area.sum(axis=1, keepdims=True)
                        if metric == 'pctmean':
                            col_idx = np.array([metric_names.index('%s_%s_pctmean' % (hemi, ctx_label)) for ctx_label in ctx_labels])
                        else:
                            col_idx = np.array([metric_names.index("_".join([aparc_code[self.atlas], hemi, ctx_label, metric])) for ctx_label in ctx_labels])
                        metric_names += ['%s_%s_%s' % (hemi, lobe, metric)]
                        metric_values = np.hstack((metric_values, np.nansum(metric_values[:, col_idx]*norm_area, axis=1, keepdims=True)))

        # Compute symmetric index
        sym_metric_names = []
        sym_metric_values = []
        for i, m in enumerate(metric_names):
            if metric_include is not None and ((m.replace('lh', 'symmetryIndex') not in metric_include) or (m.replace('Left', 'symmetryIndex') not in metric_include)):
                continue
            if metric_exclude is not None and ((m.replace('lh', 'symmetryIndex') in metric_exclude) or (m.replace('Left', 'symmetryIndex') in metric_exclude)):
                continue
            if 'lh' in m:  # lh_ label means there's a column with rh_ label that corresponds (look for lh instead of rh because rh matches entorhinal and other structures with rh in it)
                if m.replace('lh', 'rh') in metric_names:
                    num = (metric_values[:, metric_names.index(m.replace('lh', 'rh'))] - metric_values[:, i])
                    denom = (np.abs(metric_values[:, metric_names.index(m.replace('lh', 'rh'))]) + np.abs(metric_values[:, i]))
                    sym_metric_values.append(num / denom)  # Sym index is (rh-lh)/(rh+lh)
                    sym_metric_names.append(m.replace('lh', 'symmetryIndex'))
            elif 'Left' in m:  # Left label means there's a column with Right label that corresponds
                if m.replace('Left', 'Right') in metric_names:
                    num = (metric_values[:, metric_names.index(m.replace('Left', 'Right'))] - metric_values[:, i])
                    denom = (np.abs(metric_values[:, metric_names.index(m.replace('Left', 'Right'))]) + np.abs(metric_values[:, i]))
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
            if m in ['asegvolume_EstimatedTotalIntraCranialVol', '%s_eTIV' % aparc_code[self.atlas]]:
                metric_unit.append(r'')
                metric_type.append('Estimated total ICV')
                continue
            if 'volume' in m or 'BrainSegVolNotVent' in m:
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
                print('%s does not match any metric type in Freesurfer.py metric list' % m)

        # Get indexes of unique combinations of scanner and sequence in reference subjects
        matching_cols = np.array([covariate_names.index('sequence'), covariate_names.index('scanner')])
        ref_matching_cols = np.array([ref_covariate_names.index('sequence'), ref_covariate_names.index('scanner')])
        scanseq_combinations = np.unique(ref_covariate_values[:, ref_matching_cols], axis=0)

        # Normalize metrics
        normalized_metrics = np.full(metric_values.shape, np.nan)
        TIV = metric_values[:, metric_names.index('asegvolume_EstimatedTotalIntraCranialVol')].copy()  # has length # scans
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
                    filtering_idxs = np.arange(len(dict_reader.fieldnames)-1)
                    if metric_include is not None:
                        filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name in metric_include]
                        filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name in metric_include])
                        tmp_metric_names = filtered_metric_names
                    if metric_exclude is not None:
                        filtered_metric_names = [metric_name for metric_name in tmp_metric_names if metric_name not in metric_exclude]
                        filtering_idxs = np.array([idx for idx, metric_name in enumerate(tmp_metric_names) if metric_name not in metric_exclude])
                        tmp_metric_names = filtered_metric_names
                    if len(tmp_metric_names) > 0:
                        logging.PRINT('Loaded metrics are %s' % (', '.join(list(tmp_metric_names))))
                        for row in dict_reader:
                            new_subj = list(row.values())[0]
                            if new_subj != ID:
                                logging.PRINT('Entered 1st if statement')
                                logging.WARNING("Subject %s is labeled as %s in freesurfer output (aseg_stats_%s.txt)."
                                                " The subject was discarded." % (ID, new_subj, metric))
                            else:
                                logging.PRINT('Entered 2nd if statement')
                                row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                                metric_values.append(np.array(list(row.values())[1:], dtype='float')[filtering_idxs])  # skip participant ID and filter metrics
                                metric_names += tmp_metric_names
                        logging.PRINT('Exiting if statement')
                logging.PRINT('Loaded metric_names are %s' % (', '.join(metric_names)))
                logging.PRINT('len(metric_names)=%d' % (len(metric_names)))
                logging.PRINT('Index of asegvolume_EstimatedTotalIntraCranialVol is %d' % (metric_names.index('asegvolume_EstimatedTotalIntraCranialVol')))
        metric_values = np.array(metric_values, dtype='float')
        # Make it 2D
        if len(metric_values.shape) == 1:
            metric_values = metric_values[None, :]
        for hemi in ['rh', 'lh']:
            if self.atlas == 'DesikanKilliany':
                parc_metrics = self.parc35_metrics
            elif self.atlas == 'Destrieux':
                parc_metrics = self.parc75_metrics
            for metric in parc_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.%s_stats_%s.txt' % (hemi, aparc_code[self.atlas], metric))
                if os.path.exists(stats2table_file):
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['%s_%s' % (aparc_code[self.atlas], m) for m in dict_reader.fieldnames[1:]]  # skip <hemi>.aparc.<metric> entry at start
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
                                    logging.WARNING("Subject %s is labeled as %s in freesurfer output (%s.%s_stats_%s.txt)."
                                                    " The subject was discarded." % (ID, new_subj, hemi, aparc_code[self.atlas], metric))
                                else:
                                    row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                                    metric_values = np.hstack((metric_values, np.array(list(row.values())[1:], dtype='float')[filtering_idxs][None, :]))  # skip participant ID and filter metrics
                                    metric_names += tmp_metric_names

            # GM-WM contrast needs asegstats (table format a little different, no ID at each line in CH-first for example)
            for metric in self.parcSTD_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.%s_stats_pct%s.txt' % (hemi, aparc_code[self.atlas], metric))
                if os.path.exists(stats2table_file):
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['%s_%s_pct%s' % (hemi, m, metric) for m in dict_reader.fieldnames[1:]]  # skip <hemi>.aparc.<metric> entry
                        filtering_idxs = np.arange(len(dict_reader.fieldnames)-1)
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
                                row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                                metric_values = np.hstack((metric_values, np.array(list(row.values())[1:], dtype='float')[filtering_idxs][None, :]))  # skip subject index as doesn't provide any info, and filter metrics
                                # loaded_subjects.append(list(row.values())[0])  <- commented as no ID in pct files (used --inputs option in asegstats2table), using same index as other FS output tables
                                metric_names += tmp_metric_names

            # Optional lgi files if processed with matlab script
            lgi_file = os.path.join(stats2table_folder, '%s.%s_stats_lgimean.txt' % (hemi, aparc_code[self.atlas]))
            if os.path.exists(lgi_file):
                with open(lgi_file, 'r') as f:
                    dict_reader = DictReader(f, delimiter='\t')
                    tmp_metric_names = ['%s_%s_lgimean' % (hemi, m) for m in dict_reader.fieldnames[1:]]  # skip <hemi>.aparc.<metric> entry
                    filtering_idxs = np.arange(len(dict_reader.fieldnames)-1)
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
                            row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                            metric_values = np.hstack((metric_values, np.array(list(row.values())[1:], dtype='float')[filtering_idxs][None, :]))  # skip participant ID and filter metrics
                            metric_names += tmp_metric_names

        if len(metric_values) == 0:
            logging.WARNING('No metric values were found for subject %s, will likely have only NaNs' % ID)
        else:
            # Remove repeats
            repeats = [i for i, metric_name in enumerate(metric_names) if metric_name in metric_names[i+1:]]
            while(len(repeats) > 0):
                metric_names.pop(repeats[0])
                metric_values = np.delete(metric_values, repeats[0], 1)
                repeats = [i for i, metric_name in enumerate(metric_names) if metric_name in metric_names[i + 1:]]
        return metric_values, metric_names

    def load_subject_metrics_from_stats2tableFolder(self, stats2table_folder, subjSesAcq_list, metric_include=None,
                                                    metric_exclude=None):
        """
        Loads tables with multiple subjects. Numpy arrays have to be hstacked from table to table.
        Test: allocate numpy array of NaNs of shape N_subj (known from subjSesAcq_list) x N_metrics (known from dictreader)
        and fill it when reading metric table line by line using subjSesAcq_list.index(ID) or something similar

        :param stats2table_folder: path to folder containing files with all subject metrics
        :type stats2table_folder: string
        :param subjSesAcq_list: list of <subj_id>/<ses_id>/<aqc_label> combinations to load
        :type subjSesAcq_list: list of strings
        :param metric_include: list of metric names to keep
        :type metric_include: list
        :param metric_exclude: list of metric names to remove
        :type metric_exclude: list
        :return:
        """

        metric_names = ["Dummy_col"]
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
                    filtering_idxs = np.arange(len(dict_reader.fieldnames)-1)
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
                            logging.PRINT('Loaded metrics are %s' % (', '.join(tmp_metric_names)))
                            ID = list(row.values())[0]  # First value is subject ID, no common key between tables
                            if ID != loaded_IDs[i]:
                                logging.WARNING("%s in aseg_stats_%s.txt (line %d) doesn't match previous order !" % (ID, metric, i))
                            if ID not in subjSesAcq_list:
                                logging.PRINT('Entered 1st if statement')
                                logging.WARNING("Subject %s in freesurfer output (aseg_stats_%s.txt) is not in subjSesAcq_list."
                                                " The subject was discarded." % (ID, metric))
                            else:
                                logging.PRINT('Entered 2nd if statement')
                                row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                                table_to_fill[subjSesAcq_list.index(ID), :] = np.array(list(row.values())[1:])[filtering_idxs]  # skip ID and filter metrics
                            logging.PRINT('Exiting if statement')
                        metric_names += tmp_metric_names
                        metric_values = np.hstack((metric_values, table_to_fill))
        logging.PRINT('Loaded metric_names are %s' % (', '.join(metric_names)))
        logging.PRINT('len(metric_names)=%d' % (len(metric_names)))
        logging.PRINT('Index of asegvolume_EstimatedTotalIntraCranialVol is %d' % (metric_names.index('asegvolume_EstimatedTotalIntraCranialVol')))
        metric_values = metric_values[:, 1:]  # Trim out 1st column used to initialize array
        metric_names = metric_names[1:]
        for hemi in ['rh', 'lh']:
            if self.atlas == 'DesikanKilliany':
                parc_metrics = self.parc35_metrics
            elif self.atlas == 'Destrieux':
                parc_metrics = self.parc75_metrics
            for metric in parc_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.%s_stats_%s.txt' % (hemi, aparc_code[self.atlas], metric))
                if os.path.exists(stats2table_file):
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['%s_%s' % (aparc_code[self.atlas], m) for m in dict_reader.fieldnames[1:]]
                        filtering_idxs = np.arange(len(dict_reader.fieldnames)-1)
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
                                    logging.WARNING("%s in %s.%s_stats_%s.txt (line %d) doesn't match previous order !" % (ID, hemi, aparc_code[self.atlas], metric, i))
                                if ID not in subjSesAcq_list:
                                    logging.WARNING("Subject %s in freesurfer output (%s.%s_stats_%s.txt) is not in"
                                                    " subjSesAcq_list. The subject was discarded." % (ID, hemi, aparc_code[self.atlas], metric))
                                else:
                                    row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                                    table_to_fill[subjSesAcq_list.index(ID), :] = np.array(list(row.values())[1:])[filtering_idxs]
                            metric_names += tmp_metric_names
                            metric_values = np.hstack((metric_values, table_to_fill))

            # GM-WM contrast needs asegstats (table format a little different, no ID at each line in CH-first for example)
            for metric in self.parcSTD_metrics:
                stats2table_file = os.path.join(stats2table_folder, '%s.%s_stats_pct%s.txt' % (hemi, aparc_code[self.atlas], metric))
                if os.path.exists(stats2table_file):
                    with open(stats2table_file, 'r') as f:
                        dict_reader = DictReader(f, delimiter='\t')
                        tmp_metric_names = ['%s_%s_pct%s' % (hemi, m, metric) for m in dict_reader.fieldnames[1:]]
                        filtering_idxs = np.arange(len(dict_reader.fieldnames)-1)
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
                                if loaded_IDs[i] not in subjSesAcq_list:
                                    logging.WARNING(
                                        "Line %d in %s.%s_stats_pct%s.txt (%s in tables with IDs) is not in"
                                        " subjSesAcq_list. The subject was discarded." % (i, hemi, aparc_code[self.atlas], metric, loaded_IDs[i]))
                                else:
                                    row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                                    table_to_fill[subjSesAcq_list.index(loaded_IDs[i]), :] = np.array(list(row.values())[1:])[filtering_idxs]
                            metric_names += tmp_metric_names
                            metric_values = np.hstack((metric_values, table_to_fill))

            # Optional lgi files if processed with matlab script
            lgi_file = os.path.join(stats2table_folder, '%s.%s_stats_lgimean.txt' % (hemi, aparc_code[self.atlas]))
            if os.path.exists(lgi_file):
                with open(lgi_file, 'r') as f:
                    dict_reader = DictReader(f, delimiter='\t')
                    tmp_metric_names = ['%s_%s_lgimean' % (hemi, m) for m in dict_reader.fieldnames[1:]]
                    filtering_idxs = np.arange(len(dict_reader.fieldnames)-1)
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
                            if loaded_IDs[i] not in subjSesAcq_list:
                                logging.WARNING(
                                    "Line %d in %s.%s_stats_lgimean.txt (%s in tables with IDs) is not in"
                                    " subjSesAcq_list. The subject was discarded." % (i, hemi, aparc_code[self.atlas], loaded_IDs[i]))
                            else:
                                row = {key: np.nan if value == "" else value for key, value in row.items()}  # Convert empty strings to nan
                                table_to_fill[subjSesAcq_list.index(loaded_IDs[i]), :] = np.array(list(row.values())[1:])[filtering_idxs]
                        metric_names += tmp_metric_names
                        metric_values = np.hstack((metric_values, table_to_fill))

        return metric_values, metric_names

    def generate_location_plots(self, output_folder):
        my_env = os.environ.copy()
        if 'FREESURFER_HOME' not in my_env.keys():
            logging.ERROR('$FREESURFER_HOME environment variable not found. Set $FREESURFER_HOME to your FREESURFER '
                          'installation path before running ScanOMetrics.')
        fs_home = my_env['FREESURFER_HOME']
        hemis = ['lh', 'rh']
        for atlas in self.atlases[self.version].keys():
            for roi, roi_view in [self.atlases[self.version][atlas]['ROIs'], self.atlases[self.version][atlas]['views']]:
                for hemi in hemis:
                    orig_annot = read_annot(os.path.join(fs_home, 'subjects', 'fsaverage', 'label', '%s.%s.annot' % (hemi, aparc_code[atlas])))
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
                    # print(' '.join(cmd))
                    subprocess.call(cmd)

            for lobe, rois in self.lobe_ROIs[atlas]:
                for hemi in hemis:
                    orig_annot = read_annot(os.path.join(fs_home, 'subjects', 'fsaverage', 'label', '%s.%s.annot' % (hemi, aparc_code[atlas])))
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
                            # print(' '.join(cmd))
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
                            # print(' '.join(cmd))
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

        :param subjects_dir: path to Freesurfer subjet directory
        :type subjects_dir: string
        :param subject: subject ID (folder name in subjects_dir)
        :type subject: string
        :param atlas: atlas to use (DesikanKilliany or Destrieux)
        :type atlas: string
        :param hemi: hemispheres to process ('lh' or 'rh')
        :type hemi: string
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
        for roi in self.atlases[self.version][atlas]['ROIs']:
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
            ROIs = self.atlases[self.version][atlas]['ROIs']
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

            out_file = os.path.join(self.subjects_dir, subject, 'stats2table', '%s.%s_stats_area_%s.txt' % (hemi, aparc_code[atlas], surfname))
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

            out_file = os.path.join(self.subjects_dir, subject, 'stats2table', '%s.%s_stats_gauscurv_%s_FS.txt' % (hemi, aparc_code[atlas], surfname))
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
        subject_IDmergedSes = []
        for subj_id in subjects.keys():
            for ses_id in subjects[subj_id].keys():
                for acq_id in subjects[subj_id][ses_id].keys():
                    subjSesAcq_id = self.get_subjSesAcq_id(subj_id, ses_id, acq_id)
                    subject_IDmergedSes.append(subjSesAcq_id)
        return subject_IDmergedSes

    def get_subjSesAcq_row(self, subjects, subject_id, session_id, acq_id):
        """Quick and dirty way of getting row index for measured_metrics and covariate_values"""
        subject_IDmergedSes = self.get_subjSesAcq_array(subjects)
        subjSesAcq_id = self.get_subjSesAcq_id(subject_id, session_id, acq_id)
        return subject_IDmergedSes.index(subjSesAcq_id)

    def get_subjSesAcq_id(self, subject_id, session_id, acq_id):
        subjSesAcq_id = subject_id
        if session_id != "":
            subjSesAcq_id += self.ses_delimiter+session_id
        if acq_id != "":
            subjSesAcq_id += self.acq_delimiter+acq_id
        return subjSesAcq_id

    def get_subjSesAcq_T1s(self, subjects):
        subjSesAcq_list = []
        subjSesAcqs_T1s = []
        for subj_id in subjects.keys():
            for ses_id in subjects[subj_id].keys():
                for acq_label in subjects[subj_id][ses_id].keys():
                    subjSesAcq_id = self.get_subjSesAcq_id(subj_id, ses_id, acq_label)
                    T1_file = os.path.join(self.bids_database, subj_id, ses_id, 'anat', '%s.nii.gz' % subjSesAcq_id)
                    if os.path.exists(T1_file):
                        subjSesAcq_list.append(subjSesAcq_id)
                        subjSesAcqs_T1s.append(T1_file)
                    else:
                        logging.WARNING("""File %s not found: batch-dl+direct will skip subject %s, you might want
                        to process this subject manually.""" % (T1_file, subjSesAcq_id))
        return subjSesAcq_list, subjSesAcqs_T1s
