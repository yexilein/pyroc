
This directory contains the scripts generating ROC curves and extracting FECs from Protein-Protein Interaction networks and co-expression networks. The scripts assume the existence of the following (empty) directories: figs, panels, curves, curves/coexp, curves/ppi, curves/ppi_binary, stats, stats/coexp, stats/ppi, stats/ppi_binary.

Run in the following order:
 - 01_make_go_curves.py: Generate ROC curves for all GO terms from the PPI and co-expression networks. Results are stored in "curves".
 - 02_summarize_go_curves.py: Extract FECs and other summary statistics from ROC curves. Results are stored in "stats".
 - 03_plot_panels.R: Generate panels for figures 4 and 5. Results are stored in "panels".
 - 04_plot_permutations.R: Generate panels for Figure 3 (robustness of FECs to label permutations). Results are stored in "figs".