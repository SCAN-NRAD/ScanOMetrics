Setup
=====

.. _installation:

Create virtual environment (optional)
-------------------------------------
Download and install `Miniconda <https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions>`_ and
create a new conda environment:

.. code-block:: console

    conda create -y -n scanometrics python=3.10
    source activate scanometrics

Installation
------------

You can install ScanOMetrics by cloning the code from our `GitHub repo <https://github.com/SCAN-NRAD/scanometrics.git>`_ and using `pip`. First `cd` to the
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

MRI processing tools
--------------------

ScanOMetrics relies on either `Freesurfer <https://surfer.nmr.mgh.harvard.edu/>`_ or
`DL+DiReCT <https://github.com/SCAN-NRAD/DL-DiReCT>`_ to process MRI scans available in BIDS directories.
Installation instructions are found on each tool's websites.

Graphical User Interface
------------------------
ScanOMetrics main functionalities can be accessed through the dedicated `GUI <./tutorials/gui.html>`_.

Command line usage
------------------
Advanced features are available through our `Command line implementation <./tutorials/command_line.html>`_.