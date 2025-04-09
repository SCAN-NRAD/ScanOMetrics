dicom2bids (utility tool)
=========================

We provide a simple script to convert a dicom folder into nifti files, and export into a bids-compatible folder
structure to be read by ScanOMetrics.

The script uses the dicom `PatientID` and `AccessionNumber` to create anonymous participant and session IDs. The tool
keeps track of the conversion through a local `mapping.keys' file, which can simply be deleted when sharing the BIDS
database in an anonymized way.

ScanOMetrics evaluates new scans by comparing participant metrics to matching scans in the normative dataset. The
matching currently involves finding scans with similar contrast. We assign sequence labels according to the dicom
`SeriesDescription` value. The labeling is done based on the `seriesmapping` dictionary inside the `pacs2bids.py` file.
Sequence labels should match the ones used in the normative model used for scan evaluation. You should make sure that
`seriesmapping` maps `SeriesDescription` values to the same sequence labels used in the normative model. If you were to
use a different mapping dictionary, make sure to recompile the ScanOMetrics package and re-run the dicom2bids export.

Usage
-----

To add new scans from a dicom folder to your BIDS database, you can run the following::

    dicom2bids /path/to/bids/database /path/to/dicom/folder

The script assumes the dicom folder contains scans from a single session. Repeated scans will lead to nifti files with
increasing `run-` labels.

BIDS folder structure
---------------------

Assuming there are two participants `sub-001` and `sub-002`, each with 2 sessions called `ses-1` and `ses-2`, the output
bids folder will be populated as follows:::

    BIDS_database/
    |
    ├ code/
    |
    ├ derivatives/
    |
    ├ participants.tsv: tsv file with the list of participant IDs and sex
    |
    ├ sub-001/
        |
        ├ sub-001_sessions.tsv: tsv file with the list of sub-001 session IDs, age, and scan sequences
        |
        ├ sub-001_ses-1/
            └ anat/
                ├ sub-001_ses-1_acq-<seqLabel>_run-1_T1w.json
                └ sub-001_ses-1_acq-<seqLabel>_run-1_T1w.nii.gz
        ¦
        └ sub-001_ses-2/
            └ anat/
                ├ sub-001_ses-2_acq-<seqLabel>_run-1_T1w.json
                └ sub-001_ses-2_acq-<seqLabel>_run-1_T1w.nii.gz
    ¦
    └ sub-002
        |
        ├ sub-002_sessions.tsv: tsv file with the list of sub-002 session IDs, age, and scan sequences
        |
        ├ sub-002_ses-1/
            └ anat/
                ├ sub-002_ses-1_acq-<seqLabel>_run-1_T1w.json
                └ sub-002_ses-1_acq-<seqLabel>_run-1_T1w.nii.gz
        ¦
        └ sub-002_ses-2/
            └ anat/
                ├ sub-002_ses-2_acq-<seqLabel>_run-1_T1w.json
                └ sub-002_ses-2_acq-<seqLabel>_run-1_T1w.nii.gz

Such a folder structure contains the minimal information required by ScanOMetrics to process a participant (age, sex,
sequence, contrast or native T1 scan, etc...).