"""
Template for processing pipeline. Should contain a run() generic function to run the processing, as well as a proc2metric
function in order to gather the generated outputs and combine them into a matrix. The template provides the metrics as
a table to bypass processing.
Remember to add new processing modules to the __init__.py for import and recognition as submodule by scanometrics.
"""


import os
import numpy as np
from csv import DictReader
import sys
from scanometrics.utils import logging


def proc2metric(subjects_dir, subj_ids, metric_tables=['metrics.txt'], sep=','):
    """
    Template function to convert tables with processing outputs to a list of variable names and corresponding numpy
    array, intended to be stacked with other subjects in a single matrix for a scanometrics study. Functions build from
    this template should be callable through proc2metric(subjects_dir, subj_ids).
    :param subjects_dir: path to folder containing the different subjects
    :param subj_ids: ids of subjects to process in parallel
    :param metric_tables: list of table files. Should have a header with variable names, and a single row with values.
    Path should be specified relative to <subjects_dir>/<subj_id>
    :param sep: character used as separator in metric_tables
    :return: [metric_names, metric_values] as a list of variable names, and numpy array
    """
    metric_names = []
    metric_values = []
    for metric_table in metric_tables:
        for subj_id in subj_ids:
            with open(os.path.join(subjects_dir, subj_id, metric_table), 'r') as f:
                for row in DictReader(f, delimiter=sep):
                    new_metrics = [new_metric for new_metric in list(row.keys()) if new_metric not in metric_names]
                    if len(new_metrics):
                        logging.PRINT('Adding metrics %s to metric_names' % new_metrics)
                        metric_names += [new_metric for new_metric in list(row.keys()) if new_metric not in metric_names]
                    subject_values = []
                    for metric in metric_names:
                        if metric not in row.keys():
                            subject_values.append('nan')
                        else:
                            subject_values.append(row[metric])
                    metric_values.append(subject_values)
    try:
        metric_values = np.array(metric_values, dtype='float')
    except ValueError as e:
        print('Metric values cannot be casted to float, make sure metric tables contain only numerical data.')
        sys.exit(e)
    return metric_names, np.array(metric_values, dtype='float')
