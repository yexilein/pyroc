This directory contains the scripts generating ROC curves and extracting FECs for single-cell markers based on the Tabula Muris dataset. These scripts assume that you have downloaded the Tabula Muris dataset and converted it to SingleCellExperiment format. They also assume the existence of the following (empty) directories: de_lists, fec, figs.

Run the scripts in following order:
 - 01_compute_markers.R: compute markers for all cell types. Results are stored in "de_lists".
 - 02_plot_fecs.R: extract and plot FECs for B cells (Fig. 5). Results are stored in "fec" and "figs".