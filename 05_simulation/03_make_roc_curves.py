
import os.path
import numpy as np
import sys

sys.path.append("..")
import pyroc


def main():
    np.random.seed(17)
    for net_name in ["net1","net2","net3","net4"]:
        network, genes = load_network(net_name)
        labels = make_labels(genes, label_parameters())
        generate_curves(network, labels, os.path.join("curves", net_name + "_"))

def load_network(net_name, input_dir = "simulated_networks"):
    net = np.loadtxt(os.path.join(input_dir, net_name + ".txt"))
    with open(os.path.join(input_dir, net_name + "_genes.txt")) as fp:
        genes = fp.read().splitlines()
    return net, genes

def make_labels(genes, params):
    genes = np.array(genes)
    result = []
    for my_params in params:
        block_index = np.array([int(g.rsplit("|")[1]) for g in genes])
        annotation_strength = np.random.uniform(0,1,len(block_index))
        p_annotated = np.array(my_params["p_annotated"])[block_index-1]
        result.append(annotation_strength <= p_annotated)
    result = np.array(result)*1
    return result.T

def label_parameters():
    return [
        {"p_annotated": [0.5, 0.5, 0.5, 0.5]},
        {"p_annotated": [0.8, 0.1, 0.1, 0.1]},
        {"p_annotated": [0.8, 0.5, 0.1, 0.1]},
        {"p_annotated": [0.4, 0.3, 0.2, 0.1]}
    ]

def generate_curves(network, labels, output_prefix):    
    rocs, aurocs, predictors = pyroc.compute_neighbor_voting(network, labels)
    np.savetxt(output_prefix + "curves.txt", rocs)
    np.savetxt(output_prefix + "aurocs.txt", aurocs)
    np.savetxt(output_prefix + "predictors.txt", predictors)


if __name__ == "__main__":
    main()
