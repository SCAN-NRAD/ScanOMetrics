Setup
=====

.. _installation:

Create virtual environment (optional)
-------------------------------------
Download and install `Miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/linux.html>`_ and
create a new conda environment:

.. code-block:: console

    conda create -y -n scanometrics python=3.10
    source activate scanometrics

Installation
------------

You can install ScanOMetrics by cloning the code from  and using `pip`. First `cd` to the
folder you want to download ScanOMetrics to, and run the following (a `scanometrics` folder
in your current directory will be created):

.. code-block:: console

    (scanometrics)$ git clone https://github.com/SCAN-NRAD/scanometrics.git
    (scanometrics)$ cd scanometrics
    (scanometrics)$ pip install .

Data formatting
---------------

ScanOMetrics heavily relies on the dataset structure to follow the format recommended by
the `BIDS initiative on Brain Imaging Data Structure <https://bids.neuroimaging.io/>`_. Make
sure your data is properly organized before running ScanOMetrics.

Processing MRI data
-------------------

ScanOMetrics relies on either `Freesurfer <https://surfer.nmr.mgh.harvard.edu/>`_ or
`DL+DiReCT <https://github.com/SCAN-NRAD/DL-DiReCT>`_ to process MRI scans available in BIDS directories.
Check out the tutorial on `MRI processing <./tutorials/process_MRI_scans.html>`_ for more details.

Single subject evaluation
-------------------------

If you have processed MRI scans from a subject, you can compare them to a normative dataset by
following the tutorial on `Subject evaluation <./tutorials/evaluate_single_subject.html>`_.

Group evaluation
----------------

ScanOMetrics is mainly intended for single subject evaluation. However, we also provide a way
to perform group comparisons, by following `this tutorial <./tutorials/evaluate_group.html>`_.

Fit a normative model
---------------------

We provide a selection of normative models for subject and group evaluation. However, if you
want to create your own normative model for ScanOMetrics, you can do so by following the
tutorial on `Fitting normative data <./tutorials/fit_normative_data.html>`_.

Adding a processing module
--------------------------

We provide default processing modules, as wrappers around Freesurfer and DL+DiReCT. If you'd
like to implement your own processing pipeline, you can follow some basic ideas `here <./tutorials/add_processing_module.html>`_.
