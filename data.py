
import os.path
import numpy as np
import scipy.io
import pandas as pd
import h5py


DATA_DIR = "/home/fischer/projects/ROC/01_data/"

def load_coexp(file_prefix):
    network = np.loadtxt(file_prefix + ".txt.gz")
    with open(file_prefix + "_genes.txt") as fp:
        genes = fp.read().splitlines()
    return network, genes

def mouse_ppi():
    return load_ppi(os.path.join(DATA_DIR, "biogrid_mouse"))

def load_ppi(file_prefix):
    network = scipy.io.mmread(file_prefix + "_ppi.txt").A
    with open(file_prefix + "_genes.txt") as fp:
        genes = fp.read().splitlines()
    return network, genes

def load_go_mouse():
    result = load_go(os.path.join(DATA_DIR, "go_mouse"))
    return result[0].T, result[2], result[1]

def load_go(file_prefix):
    network = scipy.io.mmread(file_prefix + ".txt").A
    with open(file_prefix + "_row_labels.txt") as fp:
        genes = fp.read().splitlines()
    with open(file_prefix + "_col_labels.txt") as fp:
        terms = fp.read().splitlines()
    return network, genes, terms

def mouse_coexp():
    f = h5py.File(os.path.join(DATA_DIR, "mouse_cococonet_210420.hdf5"), "r")
    network = np.array(f["agg"])
    genes = list(g.decode("utf-8") for g in f["row"])
    mouse_ids = mouse_symbols()
    known_ids = set(mouse_ids.keys())
    keep_gene = [g in known_ids for g in genes]
    subnetwork = network[np.ix_(keep_gene,keep_gene)]
    symbols = list(mouse_ids.get(g, None) for g in genes)
    symbols = list(s for s in symbols if s is not None)
    return subnetwork, symbols

def mouse_symbols():
    id_map = pd.read_csv(os.path.join(DATA_DIR, "mouse_ids_210421.csv")).dropna()
    result = dict(zip(list(id_map["ensembl"]), list(id_map["symbol"])))
    return result

def mouse_common_genes():
    with open(os.path.join(DATA_DIR, "mouse_common_genes.txt")) as fp:
        genes = fp.read().splitlines()
    return genes
