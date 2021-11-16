
Files used to download and pre-process data.

Run in the following order:
 - 01_download_data.sh: download data from GO, MGI, BIOGRID and CoCoCoNET.
 - 02_go_annotations.R: filter and transform GO data in sparse matrix format.
 - 03_biogrid.R: filter and transform BIOGRID data in sparse matrix format.
 - 04_make_id_map.R: generate mapping from EnsemblID to Gene symbols for CoCoCoNET (which uses EnsemblID).
 - 05_make_common_genes.py: identify the set of genes that is common to the PPI data, GO annotations and co-expression network.
 
When running the scripts, additional files are generated:
 - 4 files containing the raw data (GO, GO annotations, BIOGRID, CoCoCoNET).
 - 2 files containing the parsed biogrid matrix (one for the matrix, one for gene names).
 - 3 files containing the parsed go annotations (matrix, column labels, row labels).
 - 1 file containing the EnsemblID to Symbol mapping.
 - 1 file containing the set of common genes.