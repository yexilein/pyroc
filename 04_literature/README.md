
This directory contains the scripts extracting FECs from ROC curves from published articles. The scripts assume the existence of the following (empty) directories: figs, stats. The "curves" directory contains the ROC curves extracted from the literature.

Run in the following order:
 - 01_summarize_curves.py: Extract FECs and other summary statistics from ROC curves. Results are stored in "stats".
 - 02_plot_panels.R: Generate panels for figures 1. Results are stored in "figs".
