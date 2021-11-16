
import numpy as np
import sys
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append("..")
import pyroc
import data


def main():
    np.random.seed(17)
    
    genes = data.mouse_common_genes()
    go_sets = data.load_go_mouse()
    go_used = filter_go(go_sets, genes)
    go_prior = make_go_prior(go_used)

    network = data.mouse_coexp()
    generate_curves(network, genes, go_used, go_prior, "curves/coexp/")
    
    network = data.mouse_ppi()
    generate_curves(network, genes, go_used, go_prior, "curves/ppi_binary/")

    network = (propagate_ppi(network[0]), network[1])
    generate_curves(network, genes, go_used, go_prior, "curves/ppi/")
    
    
def filter_go(go_sets, common_genes, min_genes=20):
    go_used, go_genes, go_terms = go_sets
    
    indices_go = find_indices(go_genes, common_genes)
    go_used = go_used[indices_go,]
    
    n_genes = np.sum(go_used, axis = 0)
    keep_term = n_genes >= min_genes
    go_used = go_used[:, keep_term]
    go_terms = [t for keep, t in zip(keep_term, go_terms) if keep]
    n_genes = n_genes[keep_term]
    with open("curves/go_terms.txt", "w") as fp:
        fp.write("Curve index\tGO term\tNumber of genes\n")
        for i, t, g in zip(range(len(n_genes)), go_terms, n_genes):
            if g > 2:
                fp.write("\t".join((str(i), t, str(g))) + "\n")
    with open("curves/genes.txt", "w") as fp:
        fp.write("\n".join(common_genes) + "\n")
    return go_used


def find_indices(genes, ref_genes):
    return [genes.index(g) for g in ref_genes]


def make_go_prior(go_net):
    n_positives = np.sum(go_net, axis=0)
    n_negatives = go_net.shape[0]-n_positives
    term_score = 1/(n_positives*n_negatives)
    go_prior = np.sum(go_net * term_score, axis=1)
    np.savetxt("curves/go_prior.txt", go_prior)
    return(go_prior)


def generate_curves(network, genes, go_used, go_prior, output_prefix):
    net_adj, net_genes = network
    
    indices_net = find_indices(net_genes, genes)
    net_used = net_adj[np.ix_(indices_net, indices_net)]
    
    rocs, aurocs, predictors = pyroc.compute_neighbor_voting(net_used, go_used)
    np.savetxt(output_prefix + "curves.txt", rocs)
    np.savetxt(output_prefix + "aurocs.txt", aurocs)
    np.savetxt(output_prefix + "predictors.txt", predictors)
    
    rocs, aurocs, predictors = pyroc.compute_node_degree_performance(net_used, go_used)
    np.savetxt(output_prefix + "node_degree_curves.txt", rocs)
    np.savetxt(output_prefix + "node_degree_aurocs.txt", aurocs)
    np.savetxt(output_prefix + "node_degree_predictor.txt", predictors)

    rocs, aurocs = pyroc.compute_prior_performance(go_prior, go_used)
    np.savetxt(output_prefix + "go_prior_curves.txt", rocs)
    np.savetxt(output_prefix + "go_prior_aurocs.txt", aurocs)


def propagate_ppi(network):
    graph = csr_matrix(network)
    dist_matrix = shortest_path(csgraph=graph, directed=False, unweighted=True)
    np.fill_diagonal(dist_matrix, 1)
    return 1/dist_matrix


if __name__ == "__main__":
    main()
