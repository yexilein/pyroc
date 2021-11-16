import numpy as np
from collections import Counter

from pyroc.roc import compute_roc, compute_auroc


def compute_neighbor_voting(network, positive_matrix, n_folds = 3):
    node_degrees = np.sum(network, 0)
    # avoid division by zero later
    node_degrees[node_degrees == 0] = 1
    result = [compute_neighbor_voting_one_positive(network, node_degrees, p, n_folds)
              for p in positive_matrix.T]
    roc, auroc, predictors = zip(*result)
    return roc, auroc, predictors

def compute_neighbor_voting_one_positive(network, node_degrees, positives, n_folds = 3):
    train_positives, test_positives = make_random_cv_vectors(positives, n_folds)
    votes = train_positives.dot(network) / node_degrees
    roc = [compute_roc(v[is_train == 0], p[is_train == 0])
           for v, p, is_train in zip(votes, test_positives, train_positives)]
    roc = np.mean(roc, axis = 0)
    auroc = [compute_auroc(v[is_train == 0], p[is_train == 0])
             for v, p, is_train in zip(votes, test_positives, train_positives)]
    auroc = np.mean(auroc)
    predictor = votes
    for i in range(predictor.shape[0]):
        predictor[i,train_positives[i,:]==1] = np.nan
    predictor = np.nanmean(predictor, axis=0)
    return roc, auroc, predictor
    
def make_random_cv_vectors(positives, n_folds = 3):
    nonzero_indices = np.flatnonzero(positives)
    np.random.shuffle(nonzero_indices)
    test_nonzeros = np.array_split(nonzero_indices, n_folds)

    test_positives = np.zeros(n_folds*len(positives)).reshape(n_folds, len(positives))
    for i, indices in enumerate(test_nonzeros):
        test_positives[i, indices] = 1

    train_positives = positives - test_positives
    return train_positives, test_positives

def make_deterministic_cv_vectors(positives, n_folds = 3):
    test_positives = np.zeros(n_folds*len(positives)).reshape(n_folds, len(positives))
    nonzero_indices = np.flatnonzero(positives)
    for i in range(n_folds):
        test_positives[i, nonzero_indices[range(i, len(nonzero_indices), n_folds)]] =  1
    train_positives = positives - test_positives
    return train_positives, test_positives
    
def compute_node_degree_performance(network, positive_matrix):
    node_degrees = np.sum(network, 0)
    # avoid division by zero later
    node_degrees[node_degrees == 0] = 1
    roc, auroc = compute_prior_performance(node_degrees, positive_matrix)
    return roc, auroc, node_degrees

def compute_prior_performance(prior, positive_matrix):
    roc = [compute_roc(prior, p) for p in positive_matrix.T]
    auroc = [compute_auroc(prior, p) for p in positive_matrix.T]
    return roc, auroc
