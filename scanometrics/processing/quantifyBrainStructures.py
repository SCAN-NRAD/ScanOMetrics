"""
Methods to gather results from different subjects
"""


import os
from scanometrics.utils import logging

def quantifybrainstem(subjects_dir, subj_id, output_filename='brainstem_volume.txt'):
    """
    Collects relevant brainstem information from morphometry output

    :param subjects_dir: absolute or relative path to directory with all freesurfer outputs, by subject
    :type subjects_dir: string
    :param subj_id: name of subject being processed
    :type subj_id: string
    :param output_filename: filename of output (defaults to brainstem_volume.txt, written to <subjects_dir>/<subj_id>/stats2table/<output_filename>)
    :type output_filename: string
    :return:
    """
    input_filename = os.path.join(subjects_dir, subj_id, 'mri', 'brainstemSsVolumes.v10.txt')
    if os.path.exists(input_filename):
        with open(os.path.join(subjects_dir, subj_id, 'stats2table', output_filename), 'w') as f_out:
            logging.PRINT("Collecting data for subject: %s" % subj_id)
            # gather names of structures and volume values
            hdr = ["Subject"]
            values = [subj_id]
            with open(input_filename,'r') as f:
                for line in f:
                    tmp = line.rstrip().split(' ')
                    hdr.append(tmp[0])
                    values.append(tmp[1])
            # Write to output as two lines
            f_out.write(' '.join(hdr)+'\n')
            f_out.write(' '.join(values)+'\n')
    else:
        logging.PRINT('No brainstemSsVolumes.v10.txt file found for subject %s (skipped).' % subj_id)

def quantifyhippocampalsubfields(subjects_dir, subj_id, hemi, output_filename='hippoSf_volume.txt',suffix='T1'):
    """
    Collects relevant information for hippocampal subfields from morphometry output

    :param subjects_dir: relative or absolute path to freesurfer subjects output
    :type subjects_dir: string
    :param subj_id: subject name
    :type subj_id: string
    :param hemi: hemisphere to recover (left=lh, right=rh)
    :type hemi: string
    :param output_filename: output name, without hemisphere prefix to avoid call errors.
    :type output_filename: string
    :return: no return value
    """

    output_filename = hemi + '.' + output_filename
    input_filename = os.path.join(subjects_dir, subj_id, 'mri', '%s.hippoSfVolumes-%s.v10.txt' % (hemi, suffix))
    if os.path.exists(input_filename):
        with open(os.path.join(subjects_dir, subj_id, 'stats2table', output_filename), 'w') as f_out:
            logging.PRINT("Collecting data for subject: %s" % subj_id)
            # gather names of structures and volume values
            hdr = ["Subject"]
            values = [subj_id]
            with open(input_filename, 'r') as f:
                for line in f:
                    tmp = line.rstrip().split(' ')
                    hdr.append(hemi + '.' + tmp[0])
                    values.append(tmp[1])
            # Write to output as two lines
            f_out.write(' '.join(hdr) + '\n')
            f_out.write(' '.join(values) + '\n')
    else:
        logging.PRINT('No brainstemSsVolumes.v10.txt file found for subject %s (skipped).' % subj_id)
