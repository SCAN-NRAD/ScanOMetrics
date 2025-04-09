#!/usr/bin/env python
"""
Module to convert DICOMs organised in the sourcefile directory of a BIDS database folder.
WIP: hs_pipeline functions need to be copied here, or reimplemented for users outside of SCAN to use this script
"""

import os
import argparse
import tempfile
from shutil import which, copyfile
import subprocess
from csv import DictReader
from pydicom import dcmread as pydicom_read_file
from datetime import datetime
import json
import hashlib

seriesmapping = {# 't1_cs_mp2rage_1mm_INV1': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': ''},
                 # 't1_cs_mp2rage_1mm_INV2': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w'},
                 't1_cs_mp2rage_1mm_UNI_Images': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MP2RAGE'},
                 't1_mprage_sag_1mm': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRstd'},
                 't1_mprage_tra': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRstd'},
                 't1_mprage_sag_1mm_ms': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRstd'},
                 't1_mpr_adni_sag_iso_1mm_ns': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 'MPR_MS_fischl_freesurfer': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRfischl'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TR 2.3_IR 0.9': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TI 0.9_TR 2': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TR 1.7_IR 0.9': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TR 2.6_IR 1.1': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TR 1.84_IR 1.1': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TI 0.8_TR 1.7': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TI 1.1_TR 2.3': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TI 1.1_TR 2': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni3_sag_iso_1mm_ns_TR 2.3_IR 0.9': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mprage_sag_1mm_ns': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRstd'},
                 't1_mprage_sag_1mm_we': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRstd'},
                 't1_mprage_sag_1mm_ns_fischl_freesurfer': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRfischl'},
                 # 'flair_space_sag_1mm': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'FLAIR', 'sequence': 'MPRstd'},
                 't1_mprage_sag_p2_iso': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRstd'},
                 't1_mprage_sag_1mm_adni_ep': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_sag_iso_1mm_we': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mpr_adni_sag_iso_we_1mm': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 't1_mprage_sag_1mm_fischl_dm': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRfischl'},
                 't1_mprage_sag_1mm_ce_fs_dm': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'ce-Gd_T1w', 'sequence': 'MPRfischl'},
                 't1_cs_mp2rage_1mm_UNI_Images': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRadni'},
                 'Ax 3D T1 Gado': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'Ax3DGE'},
                 'cs_T1W_3D_TFE_km': {'modality': 'anat', 'file_in': 'T1', 'file_suffix': 'T1w', 'sequence': 'MPRtfe'}
                 }

def _generate_safeReading_random_code(length=10, seed=None):
    """Generate random alphanumerical code from a set without commonly misread characters
    (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3541865), starting with a letter for proper sorting in file browsers.

    :param length: specifies length of random string, defaults to 10
    :type length: int
    :param seed: seed for random string generator, defaults to None
    :type seed: int, optional
    :return: random alphanumerical string with given length
    :rtype: string
    """

    safe_alpha_set = 'ACGHJKLMNPQRUVWXY'
    safe_alphanumerical_set = 'ACGHJKLMNPQRUVWXY469'
    random.seed(seed)
    return ''.join(random.choices(safe_alpha_set)+random.choices(safe_alphanumerical_set, k=length-1))

def accession2hash(patient_id, accession_number, bids_directory, hash_length=8):
    """Function to retrieve hash code from local file, create and add random hash if not found. Automatically adds a
    .gitignore file to avoid unintended addition to a git repo."""
    # If local file exists, check if input already has a random key
    key_file = os.path.join(bids_directory, 'mapping.keys')
    sub_key = ''
    ses_key = ''
    if os.path.exists(key_file):
        with open(key_file, 'r') as key_in:
            for row in DictReader(key_in, delimiter='\t'):
                if row['PatientID'] == patient_id:
                    sub_key = row['SubKey']
                if row['AccessionNumber'] == accession_number:
                    if row['PatientID'] != patient_id:
                        print('ERROR: session key found but PatientID does not match')
                    else:
                        ses_key = row['SesKey']
                        print('INFO: existing key pair found: ("%s", "%s") => ("%s", "%s")' % (patient_id, accession_number, sub_key, ses_key))
    # If file does't exist, create it along with a .gitignore file to avoid unintended upload
    else:
        with open(key_file, 'w') as key_out:
            key_out.write('PatientID\tAccessionNumber\tSubKey\tSesKey\n')
        gitignore_file = os.path.join(bids_directory, '.gitignore')
        ignored_file = False
        if os.path.exists(gitignore_file):
            with open(gitignore_file, 'r') as git_in:
                for row in git_in:
                    if 'mapping.keys' in row:
                        ignored_file = True
            if not ignored_file:
                with open(gitignore_file, 'a') as git_out:
                    git_out.write('mapping.keys')
        else:
            with open(gitignore_file, 'w') as git_out:
                git_out.write('mapping.keys')
    # If ses_key=='' or sub_key=='', create missing ones and add them to the local file
    if ses_key == '' or sub_key == '':
        if sub_key == '':
            sub_key = 'sub-'+_generate_safeReading_random_code(length=hash_length)
        if ses_key == '':
            ses_key = 'ses-'+_generate_safeReading_random_code(length=hash_length)
        print('INFO: new key pair generated: ("%s", "%s") => ("%s","%s")' % (patient_id, accession_number, sub_key, ses_key))
        with open(key_file, 'a') as key_out:
            key_out.write('%s\t%s\t%s\t%s\n' % (patient_id, accession_number, sub_key, ses_key))
    return sub_key, ses_key

def scannercode(institution, device):
    if device == '27272' or device == '169589':
        return 2
    elif device == '41437':
        return 3
    elif device == '40186' or device == '170029':
        return 4
    elif device == '40270':
        return 5
    elif device == '175711':
        return 51
    elif device == '20608' or device == '35028':
        return 6
    elif device == '67081':
        return 61
    elif device == '166117':
        return 8
    else:
        return -1  # Sets scanner to -1 when it is not preset in scannercode (was set to "${INST: -1}" in dicom_extract_demographics.sh)

def main():
    # Process inputs
    parser = argparse.ArgumentParser(description='Organize data from DICOM folder to bids database')
    parser.add_argument('bidsDatabase', type=str, help='Path to bids database')
    parser.add_argument('dicomFolder', type=str, help='Folder with patient DICOM data')
    args = parser.parse_args()
    bids_directory = args.bidsDatabase
    
    # Create working directory
    entry_dir = os.getcwd()
    pacs_directory = os.path.join(bids_directory, 'sourcedata', 'pacs_imports')
    if not os.path.exists(pacs_directory):
        os.makedirs(pacs_directory)
    # working_directory = tempfile.TemporaryDirectory(dir=pacs_directory, prefix='tmp-')
    working_directory = os.path.join(pacs_directory, 'tmp-'+_generate_safeReading_random_code(5))
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
    os.chdir(working_directory)
    print("Changed to %s directory" % (os.getcwd()))
    
    # Create mapping.csv
    morpho_package_folder = os.path.dirname(which("hs_pipeline.sh"))
    CASE_ID = os.path.basename(args.dicomFolder)
    cmd = ['python', f'{morpho_package_folder}/../data_curation/dcmsort.py', '--skip_derived', '1', '--symlink', '1', '--mapping_csv', 'mapping.csv', '--pattern', f'{CASE_ID}_%MagneticFieldStrength%T%IS_CE:-ce%_%StudyDate%_%ProtocolName|SeriesDescription%', args.dicomFolder, 'rawdata-export']
    subprocess.run(cmd)

    # Convert to nifti
    raw_in_folder = os.path.join(working_directory, "rawdata-export")
    raw_out_folder = os.path.join(working_directory, "rawdata")
    subprocess.run([os.path.join(morpho_package_folder, '..', 'data_curation', 'dcm2nifti-rawdata.sh'), raw_in_folder, raw_out_folder])

    # Copy and rename images in mapping.csv to bids format
    participants_file = os.path.join(bids_directory, 'participants.tsv')
    if not os.path.exists(participants_file):
        with open(participants_file, 'w') as f2:
            f2.write("participant_id\tsex\n")
    mapping_file = os.path.join(working_directory, 'mapping.csv')  # Currently loops through selected T1 folders, not sure how mapping.csv behaves when other modalities are selected
    with open(mapping_file, 'r') as f_in:
        dict_reader = DictReader(f_in, delimiter=';')
        for scan_info in dict_reader:
            dicom_folder = os.path.join(raw_in_folder, scan_info['SUBJECT_ID'])
            dicom_file = os.listdir(dicom_folder)[0]
            dicom = pydicom_read_file(os.path.join(dicom_folder, dicom_file), stop_before_pixels=True)
            sex = dicom['PatientSex'].value
            if dicom['SeriesDescription'].value == '':
                series_desc = dicom['ProtocolName'].value
            else:
                series_desc = dicom['SeriesDescription'].value
            sequence = ''
            try:
                modality = seriesmapping[series_desc]['modality']
                sequence = seriesmapping[series_desc]['sequence']
            except KeyError:
                print('ERROR: %s is not in seriesmapping dict. Please update pacs2bids.py accordingly and reinstall.' % (series_desc))
                continue
            age_days = (datetime.strptime(dicom['StudyDate'].value, '%Y%m%d') - datetime.strptime(dicom['PatientBirthDate'].value, '%Y%m%d')).days
            institution = dicom['InstitutionName'].value
            device_serial = dicom['DeviceSerialNumber'].value
            scanner = scannercode(institution, device_serial)
            subj_id, ses_id = accession2hash(scan_info['PatientID'], scan_info['AccessionNumber'], bids_directory)  # PatientID can be the same btwn different accession numbers
            subj_folder = os.path.join(bids_directory, subj_id)
            sessions_file = os.path.join(subj_folder, '%s_sessions.tsv' % subj_id)
            ses_folder = os.path.join(subj_folder, ses_id, modality)
            acq_id = 'acq-%s' % sequence
            nifti_input = os.path.join(raw_out_folder, scan_info['SUBJECT_ID'], '%s.nii.gz' % (seriesmapping[series_desc]['file_in']))
            json_input = os.path.join(raw_out_folder, scan_info['SUBJECT_ID'], '%s.json' % (seriesmapping[series_desc]['file_in']))
            with open(json_input, 'r') as json_in:
                json_content = json.load(json_in)
                if "ContrastBolusAgent" in json_content.keys():
                    print("\n\n\nInfo: %s_%s_%s has %s as ContrastBolusAgent value (%s)\n\n\n" % (subj_id, ses_id, acq_id, json_content["ContrastBolusAgent"], json_input))
                    acq_id += "_ce-Gd"
            run = 1
            # Copy T1 file, check how mapping.csv behaves for other modalities to copy other files
            nifti_output = os.path.join(ses_folder, '%s_%s_%s_run-%d_%s.nii.gz' % (subj_id, ses_id, acq_id, run, seriesmapping[series_desc]['file_suffix']))
            json_output = os.path.join(ses_folder, '%s_%s_%s_run-%d_%s.json' % (subj_id, ses_id, acq_id, run, seriesmapping[series_desc]['file_suffix']))
            skip_scan = False
            while(os.path.exists(nifti_output)):
                old_hash = hashlib.md5(open(nifti_output, 'rb').read()).hexdigest()
                new_hash = hashlib.md5(open(nifti_input, 'rb').read()).hexdigest()
                if old_hash == new_hash:
                    print('WARNING: %s already exists, skipping...' % nifti_output)
                    skip_scan = True
                    break
                else:
                    print('INFO: %s already exists, but %s has a different hash, updating run...' % (nifti_output, nifti_input))
                    run += 1
                    nifti_output = os.path.join(ses_folder, '%s_%s_%s_run-%d_T1w.nii.gz' % (subj_id, ses_id, acq_id, run))
                    json_output = os.path.join(ses_folder, '%s_%s_%s_run-%d_T1w.json' % (subj_id, ses_id, acq_id, run))
            if skip_scan:
                continue
            if not os.path.exists(subj_folder):
                subj_added = False
                with open(participants_file, 'r') as f2:
                    participant_reader = DictReader(f2, delimiter='\t')
                    for line in participant_reader:
                        if line['participant_id'] == subj_id:
                            subj_added = True
                os.makedirs(subj_folder)
                if not subj_added:
                    with open(participants_file, 'a') as f2:
                        f2.write('%s\t%s\n' % (subj_id, sex))
            if not os.path.exists(ses_folder):
                ses_added = False
                if not os.path.exists(sessions_file):
                    with open(sessions_file, 'w') as f3:
                        f3.write('session_id\tage\tscanner\tsequence\n')
                with open(sessions_file, 'r') as f3:
                    session_reader = DictReader(f3, delimiter='\t')
                    for line in session_reader:
                        if line['session_id'] == ses_id:
                            ses_added = True
                os.makedirs(ses_folder)
                if not ses_added:
                    with open(sessions_file, 'a') as f3:
                        f3.write('%s\t%1.2f\t%s\t%s\n' % (ses_id, age_days/365.25, scanner, sequence))
            copyfile(nifti_input, nifti_output)
            copyfile(json_input, json_output)

    os.chdir(entry_dir)
    print("Working dir is %s" % working_directory)
    print("Finished")

if __name__ == "__main__":
    main()
