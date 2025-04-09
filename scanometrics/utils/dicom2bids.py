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
from .pacs2bids import seriesmapping, _generate_safeReading_random_code, accession2hash, scannercode

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
