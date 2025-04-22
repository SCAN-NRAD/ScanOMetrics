"""
Methods to query zenodo's REST API.
Separate datasets are uploaded for each normative model, with versions corresponding to the software used to process the data used for training.
This allows keeping track of the different models trained with different processing software versions, without having files all over the place.
Eg: 'Polynomial_dldirect_cleanOASIS3.pkl' has record id 14712271, whose version 1.0.0 has been fitted with dldirect v1.0.0
"""

import requests
import os
import scanometrics.processing
from . import logging
from tqdm import tqdm
from typing import Optional
from glob import glob
import pickle

__ZENODO_URL__ = 'https://zenodo.org/api/records/15178429'  # zenodos base url
proc_software_list = dir(scanometrics.processing)


def get_scannrad_records():
    local_path = os.path.join(os.path.dirname(__file__), '..', 'resources', 'normative_models')
    response = requests.get(__ZENODO_URL__)
    record = response.json()
    models = {}
    for fileset in record['files']:
        file_name = fileset['key']
        if os.path.splitext(file_name)[1] != '.pkl':
            # skip file in case we upload non pkl files at some point
            continue
        models[file_name] = decode_model_filename(file_name)
        models[file_name]['download_url'] = f"https://zenodo.org/records/{record['id']}/files/{file_name}?download=1"
        models[file_name]['local_file'] = os.path.join(local_path, file_name)
    return models

def decode_model_filename(model_filename):
    # splits model filename into normative_model_id, proc_pipeline_id, proc_pipeline_version
    norm_model_id, proc_pipeline_name, proc_pipeline_version, atlas, dataset_id, tmp = model_filename.split("_")  # eg norm_model_id='Polynomial',
    #                                                                                                                  proc_pipeline_name='dldirect
    #                                                                                                                  proc_pipeline_version='v1-0-3'
    #                                                                                                                  atlas='DesikanKilliany'
    #                                                                                                                  dataset_id='OASIS3'
    #                                                                                                                  tmp='som-v1-0-3'
    #
    proc_pipeline_version = proc_pipeline_version[1:].replace("-", ".")  # eg proc_pipeline_version=1.0.3
    som_version = tmp[5:10].replace('-','.')  # eg som_version="1.0.0"
    return {'norm_model_id': norm_model_id,
            'proc_pipeline_id': proc_pipeline_name,
            'proc_pipeline_version': proc_pipeline_version,
            'atlas': atlas,
            'dataset_id': dataset_id,
            'som_version': som_version}

def list_normative_models():
    """
    Lists available models on default ScanOMetrics normative model folder and on scan-nrad Zenodo repository.
    This will list model_ids like 'Polynomial_dldirect-v1.0.3_NormColl_som-v1.0.0.pkl'
    Output: list of models with ['filename', 'proc_pipeline_id', 'version']
    """
    normative_models = get_scannrad_records()
    local_path = os.path.join(os.path.dirname(__file__), '..', 'resources', 'normative_models')
    local_files = glob(os.path.join(local_path, '*.pkl'))  # expected to match Polynomial_dldirect_v1-0-3_DesikanKilliany_NormColl_som-v0-1-0.pkl
    for f in local_files:
        normative_model_id = os.path.basename(f)
        if normative_model_id not in normative_models.keys():
            normative_models[normative_model_id] = decode_model_filename(normative_model_id)
            normative_models[normative_model_id]['local_file'] = f
    logging.PRINT('Available models:')
    for f in normative_models.keys():
        logging.PRINT('\t- %s' % f)
    logging.PRINT("Eg: to load the 1st model, run `SOM.load_normative_model('%s')`" % (list(normative_models.keys())[0]))
    return normative_models

def download_file(url: str, local_filename: str, chunk_size: Optional[int] = 8192 * 16) -> str:
    # Inspired from https://github.com/MIC-DKFZ/HD-BET/blob/678e44d546a84de0f2a7fc245f176b82b7d912fd/HD_BET/checkpoint_download.py
    # which itself borrowed from https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
    # NOTE the stream=True parameter below
    logging.PRINT(f"Downloading {url} to {local_filename}")
    with requests.get(url, stream=True, timeout=100) as r:
        r.raise_for_status()
        with tqdm.wrapattr(open(local_filename, 'wb'), "write", total=int(r.headers.get("Content-Length"))) as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                f.write(chunk)
    return local_filename