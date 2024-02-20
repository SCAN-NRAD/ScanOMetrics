"""
Module to run FSL processing scripts.
"""

import os
from scanometrics.utils import logging
import subprocess
import nibabel as nib
from shutil import which


def fsl_pve_vols(fsl_dir, subj_id, T1_file):
    """Run FSL's BET and FAST scripts on subject subj_id, and extract Partial Volume Estimates (PVEs) for CSF, GM and WM.
    PVEs are given in float, as opposed to the integer result in the original octave implementation.

    :param fsl_dir: path to FSL derivatives folder (eg bids/derivatives/fsl)
    :type fsl_dir: string
    :param subj_id: ID of subject to process
    :type subj_id: string
    :param T1_file: path for T1.nii.gz file used as input
    :type T1_file: string
    """
    if which('fslreorient2std') is None:
        logging.ERROR('fslreorient2std command not found. Please make sure FSL is installed, and that $FSLDIR is properly'
                      ' set up. You may want to source the $FSLDIR/etc/fslconf/fsl.sh script too.')
    # Change to absolute path if not already specified as such
    T1_file = os.path.abspath(T1_file)
    subj_folder = os.path.join(fsl_dir, subj_id)
    if not os.path.exists(subj_folder):
        os.makedirs(subj_folder)
    os.chdir(subj_folder)
    # Run fslreorient2std
    cmd = ["fslreorient2std", T1_file, 'T1_reorient2std.nii.gz']
    logging.PRINT('Running command %s' % (' '.join(cmd)))
    p = subprocess.run(cmd)
    if p.returncode:
        logging.ERROR('%s command exited with code %d' % (cmd[0], p.returncode))
    # Run brain extraction with BET
    cmd = ['bet', 'T1_reorient2std.nii.gz', 'T1_reorient2std_brain.nii.gz', '-R', '-S', '-B', '-f', '0.3', '-v']
    logging.PRINT('Running command %s' % (' '.join(cmd)))
    p = subprocess.run(cmd)
    if p.returncode:
        logging.ERROR('%s command exited with code %d' % (cmd[0], p.returncode))
    # Run FSL's FAST
    cmd = ['fast', '-v', '-t', '1', '-n', '3', '-H', '0.1', '-I', '4', '-l', '20.0', '-o', 'T1_reorient2std_brain', 'T1_reorient2std_brain']
    logging.PRINT('Running command %s' % (' '.join(cmd)))
    p = subprocess.run(cmd)
    if p.returncode:
        logging.ERROR('%s command exited with code %d' % (cmd[0], p.returncode))
    # Recover partial volume estimate values
    stats_folder = 'stats2table'  # we are already in subj_folder
    if not os.path.exists(stats_folder):
        os.makedirs(stats_folder)
    with open(os.path.join(stats_folder, 'fsl_pve_vols.txt'), 'w') as f1:
        f1.write('subject\tcsf\tgm\twm\n')
        pves = [subj_id]
        with open('T1_reorient2std_brain_pve.txt', 'w') as f2:
            for i in range(3):
                pve = nib.load('T1_reorient2std_brain_pve_%d.nii.gz' % i).get_fdata()
                vox_size = nib.load('T1_reorient2std_brain_pve_%d.nii.gz' % i).header.get_zooms()
                vol_pve = pve[pve > 0].sum()*vox_size[0]*vox_size[1]*vox_size[2]
                pves.append(str(vol_pve))
                f2.write('%d    %1.9f\n' % (i, vol_pve))
        f1.write('\t'.join(pves)+'\n')
