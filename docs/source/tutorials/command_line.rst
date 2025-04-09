Command line implementation
===========================

This page summarizes advanced features of ScanOMetrics, accessible through command line.

Processing MRI data
-------------------

To process new MRI scans, check the tutorial on `MRI processing <./process_MRI_scans.html>`_ for more details.

Single subject evaluation
-------------------------

If you have processed MRI scans from a subject, you can compare them to a normative dataset by
following the tutorial on `Subject evaluation <./evaluate_single_subject.html>`_.

Group evaluation
----------------

ScanOMetrics is mainly intended for single subject evaluation. However, we also provide a way
to perform group comparisons, by following `this tutorial <./evaluate_group.html>`_.

Fit a normative model
---------------------

We provide a selection of normative models for subject and group evaluation. However, if you
want to create your own normative model for ScanOMetrics, you can do so by following the
tutorial on `Fitting normative data <./fit_normative_data.html>`_.

Adding a processing module
--------------------------

We provide default processing modules, as wrappers around Freesurfer and DL+DiReCT. If you'd
like to implement your own processing pipeline, you can follow some basic ideas `here <./add_processing_module.html>`_.
