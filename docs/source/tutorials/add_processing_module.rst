Adding a processing module to Scan-O-Metrics
============================================

This tutorial summarizes how to add a new processing module to Scan-O-Metrics.

Add processing module to scanometrics.processing.<new_module>.py
****************************************************************

You can copy the pipeline_template.py file and rename as your <new_module>.py
Edit the file to overide default functions.

Add processing pipeline to scanometrics.processing.__init__.py file
*******************************************************************

Add the new module to the imports done in __init__.py

