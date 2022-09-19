
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
    coexp = data.mouse_coexp()
    ppi_binary = data.mouse_ppi()
    ppi = (propagate_ppi(ppi_binary[0]), ppi_binary[1])
    codomain_sparse = data.mouse_codomain()
    codomain = (propagate_codomain(codomain_sparse[0]), codomain_sparse[1])
    
    # protein function prediction
    genes = data.mouse_common_genes()
    go_sets = data.load_go_mouse()
    go_used = filter_gene_sets(go_sets, genes, prefix="curves/go")
    go_prior = make_go_prior(go_used, prefix="curves/go")
    np.random.seed(17)
    generate_curves(coexp, genes, go_used, go_prior, "curves/coexp/")
    generate_curves(ppi_binary, genes, go_used, go_prior, "curves/ppi_binary/")
    generate_curves(ppi, genes, go_used, go_prior, "curves/ppi/")
    
    # protein function prediction based on domain structure
    genes = data.mouse_common_genes()
    genes = set(codomain[1]).intersection(genes)
    go_used = filter_gene_sets(go_sets, genes, prefix="curves/go_cod")
    go_prior = make_go_prior(go_used, prefix="curves/go_cod")
    np.random.seed(17)
    generate_curves(codomain_sparse, genes, go_used, go_prior, "curves/codomain_sparse/")
    generate_curves(codomain, genes, go_used, go_prior, "curves/codomain/")
    
    # drug-target interaction prediction
    net_genes = data.mouse_common_genes()
    dti_sets = data.load_dti_mouse()
    genes = set(dti_sets[1]).intersection(net_genes)
    dti_used = filter_gene_sets(dti_sets, genes, prefix="curves/dti")
    dti_prior = make_go_prior(dti_used, prefix="curves/dti")
    np.random.seed(17)
    generate_curves(coexp, genes, dti_used, dti_prior, "curves/coexp_dti/")    
    generate_curves(ppi, genes, dti_used, dti_prior, "curves/ppi_dti/")
    
    
def filter_gene_sets(gene_sets, common_genes, min_genes=20, prefix="curves/go"):
    gs_used, gs_genes, gs_terms = gene_sets
    
    indices_gs = find_indices(gs_genes, common_genes)
    gs_used = gs_used[indices_gs,]
    
    n_genes = np.sum(gs_used, axis = 0)
    keep_term = n_genes >= min_genes
    gs_used = gs_used[:, keep_term]
    gs_terms = [t for keep, t in zip(keep_term, gs_terms) if keep]
    n_genes = n_genes[keep_term]
    with open(prefix + "_terms.txt", "w") as fp:
        fp.write("Curve index\tGene set\tNumber of genes\n")
        for i, t, g in zip(range(len(n_genes)), gs_terms, n_genes):
            if g > 2:
                fp.write("\t".join((str(i), t, str(g))) + "\n")
    with open(prefix + "_genes.txt", "w") as fp:
        fp.write("\n".join(common_genes) + "\n")
    return gs_used


def find_indices(genes, ref_genes):
    return [genes.index(g) for g in ref_genes]


def make_go_prior(go_net, prefix="curves/go"):
    n_positives = np.sum(go_net, axis=0)
    n_negatives = go_net.shape[0]-n_positives
    term_score = 1/(n_positives*n_negatives)
    go_prior = np.sum(go_net * term_score, axis=1)
    np.savetxt(prefix + "_prior.txt", go_prior)
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


def propagate_codomain(network):
    graph = csr_matrix(1/network)
    dist_matrix = shortest_path(csgraph=graph, directed=False)
    np.fill_diagonal(dist_matrix, 1)
    return 1/dist_matrix


if __name__ == "__main__":
    main()
