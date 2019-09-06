import numpy as np
import scipy.stats
from collections import Counter


def compute_neighbor_voting(network, positive_matrix, n_folds = 3):
    node_degrees = np.sum(network, 0)
    result = [compute_neighbor_voting_one_positive(network, node_degrees, p, n_folds)
              for p in positive_matrix.T]
    return result

def compute_neighbor_voting_one_positive(network, node_degrees, positives, n_folds = 3):
    train_positives, test_positives = make_cv_vectors(positives, n_folds)
    votes = train_positives.dot(network) / node_degrees
    roc = [compute_roc(v[is_train == 0], p[is_train == 0])
           for v, p, is_train in zip(votes, test_positives, train_positives)]
    roc = np.mean(roc, axis = 0)
    auroc = [compute_auroc(v[is_train == 0], p[is_train == 0])
             for v, p, is_train in zip(votes, test_positives, train_positives)]
    auroc = np.mean(auroc)
    return roc, auroc
    
def make_cv_vectors(positives, n_folds = 3):
    nonzero_indices = np.flatnonzero(positives)
    np.random.shuffle(nonzero_indices)
    test_nonzeros = np.array_split(nonzero_indices, n_folds)

    test_positives = np.zeros(n_folds*len(positives)).reshape(n_folds, len(positives))
    for i, indices in enumerate(test_nonzeros):
        test_positives[i, indices] = 1

    train_positives = positives - test_positives
    return train_positives, test_positives
    
def compute_roc(predictor, positives, n_points = 200):
    positive_structure = order_positives(predictor, positives)
    result = convert_ordered_positives_to_roc(positive_structure)
    return(standardize_roc(result, n_points))
    
def order_positives(predictor, positives):
    ranked_predictor = scipy.stats.rankdata(predictor)
    n_positives = np.sum(positives)
    n_negatives = len(positives) - n_positives
    
    result = np.zeros(int(n_positives + n_negatives))
    positive_vals = ranked_predictor[positives != 0]
    number_ties = Counter(ranked_predictor)
    for v in positive_vals:
        n_ties = number_ties[v]
        tie_min = int(np.floor(v - (n_ties-1)/2))
        tie_max = int(np.floor(v + (n_ties-1)/2))
        tie_range = range(tie_min-1, tie_max)
        result[tie_range] = result[tie_range] + 1/n_ties;
    result = np.flip(result)
    return result

def convert_ordered_positives_to_roc(ordered_positives):
    tpr = np.cumsum(ordered_positives) / np.sum(ordered_positives)
    fpr = np.cumsum(1-ordered_positives) / np.sum(1-ordered_positives)
    return np.array([fpr, tpr])

def standardize_roc(roc, resolution = 200):
    x = np.concatenate(([0], roc[0], [1]), axis = None)
    y = np.concatenate(([0], roc[1], [1]), axis = None)
    result = np.interp(np.linspace(0, 1, resolution+1), x, y)
    return result

def compute_auroc(predictor, positives):
    ranked_predictor = scipy.stats.rankdata(predictor)
    n_positives = np.sum(positives)
    n_negatives = len(positives) - n_positives
    return (np.sum(ranked_predictor[positives != 0]) / n_positives - (n_positives+1)/2) / n_negatives