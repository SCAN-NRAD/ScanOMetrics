Processing of MRI scans (tutorial)
==================================

Scan-O-Metrics can process MRI scans using wrappers around Freesurfer and DL+DiReCT. If you want to specify your own
processing pipeline, you can follow the corresponding tutorial `here <add_processing_module.html>`_.

To process MRI scans, you just have to configure a `ScanOMetrics_project()` with the processing pipeline you want
to use, load the subjects you want to process, and start the processing pipeline.

Configure the processing pipeline
---------------------------------

Scan-O-Metrics has two built-in processing pipelines wrappers, `freesurfer` and `dldirect`, from which you can choose
when starting a `ScanOMetrics_project`. Both processing pipelines support cortical segmentations following either the
`DesikanKilliany` or `Destrieux` atlases::

    from scanometrics.core import ScanOMetrics_project
    proc_pipeline = 'dldirect'  # <- can be either 'freesurfer' or 'dldirect'
    bids_directory = "~/Downloads/OASIS3"  # <-- example to be edited
    atlas = 'DesikanKilliany'
    SOM = ScanOMetrics_project(bids_directory, proc_pipeline=proc_pipeline, atlas=atlas)

By default, all scans in the BIDS directory are loaded. You can exclude and include sets of subjects, and also include or
exclude sets of scans, based on the respective IDs and acqusition labels. For more details, check the documentation of
:meth:`scanometrics.core.ScanOMetrics_project.__init__`. Processing can be multi-threaded by setting `n_threads` to the
desired number of threads::

    n_threads = 10
    SOM.load_subjects()
    SOM.run_proc_pipeline(n_threads)
    SOM.proc2table()

Processing outputs are saved in a "linearized" way, meaning each scan has its output saved in a folder combining
subject name, session, and acquisition label inside the `<bids_directory>/derivatives/<processing_pipeline>`.


