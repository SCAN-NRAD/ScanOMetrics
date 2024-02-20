"""
Module for reading, writing data methods
"""


def _read_fs_table(bids_path, subject_code):
    fs_table = {'fs_version': np.loadtxt(os.path.join(bids_path,'derivatives','freesurfer','sub-%s' % subject_code,'ses-%d','freesurfer_version'))}
    if fs_table['version'] == '6.0.0':
        for hemi in ['lh', 'rh']:
            for ROI in ['bankssts','caudalanteriorcingulate','caudalmiddlefrontal','corpuscallosum','cuneus',
                        'entorhinal','fusiform','inferiorparietal','inferiortemporal','isthmuscingulate',
                        'lateraloccipital','lateralorbitofrontal','lingual','medialorbitofrontal','middletemporal',
                        'parahippocampal','paracentral','parsopercularis','parsorbitalis','parstriangularis',
                        'pericalcarine','postcentral','posteriorcingulate','precentral','precuneus',
                        'rostralanteriorcingulate','rostralmiddlefrontal','superiorfrontal','superiorparietal',
                        'superiortemporal','supramarginal','frontalpole','temporalpole','transversetemporal','insula']:
                fs_table[ROI]

    return fs_table


def save_project(som_project, project_file=None):
    """Method to save a ScanOMetrics_project instance

    :param som_project: ScanOMetrics instance to save
    :type som_project: ScanOMetrics_project instance
    :param project_file: path of file to save information to (defaults to None to save to SOM_project.bids_database)
    :type project_file: string
    """
    if project_file is None:
        project_file = os.path.join(som_project.bids_database,'ScanOMetrics','ScanOMetrics_project.pkl')
    with open(project_file, 'wb') as f:
        pickle.dump(som_project, f)


def load_project(project_file):
    """Method to load a ScanOMetrics_project instance

    :param project_file: path to file with ScanOMetrics_project info
    :type project_file: string
    :returns: ScanOMetrics_project instance
    """
    with open(project_file, 'rb') as f:
        return pickle.load(f)


