
import numpy as np
import scipy.io

from pyroc import compute_neighbor_voting, Roc


def main():
    curves = np.loadtxt("new_data/ppi_curves_legacy.txt")
    export_stats(curves, "new_data/ppi_legacy_")
    curves = np.loadtxt("new_data/coexp_curves_legacy.txt")
    export_stats(curves, "new_data/coexp_legacy_")
    

def generate_ppi_curves():
    ppi_used = scipy.io.mmread("legacy/new_data/ppi_used.mtx").A
    go_used = scipy.io.mmread("legacy/new_data/GO_used.mtx").A
    stats = compute_neighbor_voting(ppi_used, go_used)
    rocs, aurocs = zip(*stats)
    np.savetxt("new_data/ppi_curves_legacy.txt", rocs)
    np.savetxt("new_data/ppi_aurocs_legacy.txt", aurocs)


def generate_coexp_curves():
    coexp_used = np.loadtxt("legacy/new_data/coexp_used.txt")
    go_used = scipy.io.mmread("legacy/new_data/GO_used.mtx").A
    stats = compute_neighbor_voting(coexp_used, go_used)
    aurocs, rocs = zip(*stats)
    np.savetxt("new_data/coexp_curves_legacy.txt", rocs)
    np.savetxt("new_data/coexp_aurocs_legacy.txt", aurocs)
    

def export_stats(curves, fileprefix):
    stats = [Roc(c) for c in curves]
    np.savetxt(fileprefix + "aurocs.txt", [s.roc() for s in stats])
    np.savetxt(fileprefix + "all_null_segments.txt", np.vstack(s.null_segments for s in stats))
    np.savetxt(fileprefix + "longest_null_segment.txt", [s.longest_null_segment for s in stats])
    np.savetxt(fileprefix + "auroc_longest_only.txt", [s.roc_longest_segment() for s in stats])
    np.savetxt(fileprefix + "auroc_excluding_longest.txt", [s.roc_excluding_longest() for s in stats])
    np.savetxt(fileprefix + "auroc_best_flip.txt", [s.roc_flip for s in stats])
    np.savetxt(fileprefix + "all_quasilinear_segments.txt", np.vstack(s.quasilinear_segments for s in stats))
    np.savetxt(fileprefix + "longest_quasilinear_segment.txt", [s.quasilinear_segments[0] for s in stats])
    

if __name__ == "__main__":
    main()