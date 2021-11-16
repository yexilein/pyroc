This repository contains the code for the analyses and figures for our preprint ["Defining the extent of gene function using ROC curvature"](https://doi.org/10.1101/2021.09.03.458825).


Requirements
------------

The code contains a mixture of shell scripts (downloading data), Python (analysis of ROC curve) and R (data handling, analysis of results and figure generation). The code was tested and run on a Linux server with CentOS7, Python 3.7.3 and R 4.0.5.

To run scripts on a command line, use `sh <shell_script>.sh`, `python <python_script>.py`, `Rscript <shell_script>.R`.

Required packages for Python are: numpy, pandas, seaborn, matplotlib.pyplot, scipy.stats, scipy.io, scipy.sparse, sklean.metrics, h5py.

Required packages for R are: tidyverse, Matrix, org.Mm.eg.db, rhdf5, rjson, ontologyIndex, ggupset, MetaMarkers (https://github.com/gillislab/MetaMarkers), SingleCellExperiment.


Description of files and directories
------------------------------------

Here is an overall description of the elements found in the root directory:
 - 01_data: this directory contains all the scripts needed to download and pre-process the data for subsequent analyses.
 - pyroc: this directory contains the Python module used to analyze ROC curves (e.g., extraction of linear segments). The module is used by other scripts in the project, there is nothing to run in this directory.