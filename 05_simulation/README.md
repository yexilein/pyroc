
This directory contains all the scripts needed to reproduce all analyses related to FECs as seen in simulation models (gaussian mixture model, network model). The scripts assume the existence of the following (empty) directories: figs, curves, simulated_networks.

Run in the following order:
 - 01_simulate_functionality.py: Generate ROC curves from gaussian mixture panel, then generate panels for figure 1.
 - 02_simulate_networks.R: Generate weighted networks of genes following pre-defined community structures and functional states. Results are stored "simulated_networks".
 - 03_make_roc_curves.py: Generate ROC curves from the simulated networks using EGAD. Results are stored in "curves".
 - 04_plot_curves.R: Generate panels for supplementary figure 1. Results are stored in "figs".
