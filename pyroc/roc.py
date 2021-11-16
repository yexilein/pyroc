
import numpy as np
import scipy.stats
from sklearn.metrics import roc_curve
import pandas as pd
import seaborn as sns


class Roc(object):
    
    def __init__(self, curve,
                 n_positives = 100,
                 local_roc_threshold = 0.51,
                 min_segment_length = 0.05,
                 quasilinear_threshold = 1,
                 min_quasilinear_length = 0.05):
        self._curve = curve
        self.n_positives = n_positives
        self.local_roc_threshold = local_roc_threshold
        self.min_segment_length = min_segment_length
        self.min_quasilinear_length = min_quasilinear_length
        self.quasilinear_threshold = quasilinear_threshold
        self.recompute_fecs()
        
    @classmethod
    def from_curve(cls, curve, *args, **kwargs):
        return cls(curve, *args, **kwargs)

    @classmethod
    def from_predictor(cls, predictor, labels, *args, **kwargs):
        curve = compute_roc(predictor, labels)
        return cls(curve, n_positives = np.sum(labels), *args, **kwargs)
        
    def recompute_fecs(self):
        self._total_length = len(self._curve)
        self.local_roc = self._compute_local_auroc()
        self._linearity_deviation, self._normalized_linearity = self._compute_linearity_deviation()
        self.null_segments = self._extract_segments(self._is_null(), self.min_segment_length)
        self.linear_segments = self._extract_segments(np.logical_and(self._is_quasilinear(), self._is_null()), self.min_quasilinear_length)
        self.flip_interval, self.roc_flip = self._find_best_flip()
            
    def _compute_local_auroc(self):
        partial_auc = np.cumsum(np.concatenate(([0], (self._curve[:-1] + self._curve[1:])/2)))
        partial_auc = partial_auc - partial_auc[np.newaxis].T
        indices = np.arange(self._total_length)[np.newaxis]
        delta_x = indices - indices.T
        delta_y = self._curve - self._curve[np.newaxis].T
        with np.errstate(divide="ignore", invalid="ignore"):
            # take partial AUC, normalize x to [0,1], remove offset, normalize y to [0,1]
            result = (partial_auc / delta_x - self._curve[np.newaxis].T) / delta_y
        result[np.logical_or(delta_x == 0, delta_y == 0)] = 0.5
        result[np.tril_indices(result.shape[0])] = result.T[np.tril_indices(result.shape[0])]
        return result
    
    def _compute_linearity_deviation(self):
        deviation = np.zeros((self._total_length, self._total_length))
        norm_deviation = np.zeros((self._total_length, self._total_length))
        for start in range(self._total_length-1):
            delta_y = self._curve[start+1:] - self._curve[start]
            delta_x = np.arange(len(delta_y)) + 1
            slope = delta_y / delta_x
            # distance to interpolated curve: we only care about the [0, delta_x] interval (lower triangle)
            local_deviation = (delta_y - slope[np.newaxis].T*delta_x) * np.tri(len(delta_y))
            # KS statistic
            with np.errstate(divide="ignore", invalid="ignore"):
                ks = np.max(abs(local_deviation / delta_y[np.newaxis].T), axis=1)
                ks[delta_y == 0] = 1
                deviation[start, start+1:] = ks
                # normalized KS statistic (depends on the number of positives)
                norm_deviation[start, start+1:] = ks * np.sqrt(abs(delta_y)*self.n_positives)
            # Mean absolute deviation
            #deviation[start, start+1:] = np.sum(abs(local_deviation), axis=1) / delta_x
        deviation[np.tril_indices(deviation.shape[0])] = deviation.T[np.tril_indices(deviation.shape[0])]
        norm_deviation[np.tril_indices(deviation.shape[0])] = norm_deviation.T[np.tril_indices(deviation.shape[0])]
        return deviation, norm_deviation
                
    def _extract_segments(self, is_valid_segment, min_length):
        indices = np.arange(self._total_length)[np.newaxis]
        segment_length = indices - indices.T
        segment_length[np.logical_not(is_valid_segment)] = 0
        intervals = []
        stats = []
        while np.max(segment_length) > 0:
            i, j = np.unravel_index(np.argmax(segment_length), segment_length.shape)
            if (j-i)/self._total_length >= min_length and self._curve[j] - self._curve[i] >= min_length:
                intervals.append((i,j))
            segment_length[i:(j+1),:] = 0
            segment_length[:,i:(j+1)] = 0

        return intervals
    
    def _is_null(self, abs_deviations=True):
        if abs_deviations:
            return abs(self.local_roc-0.5) + 0.5 <= self.local_roc_threshold
        else:
            return self.local_roc <= self.local_roc_threshold

    def _is_quasilinear(self):
        return self._normalized_linearity <= self.quasilinear_threshold
    
    def _find_best_flip(self):
        """Find flip that maximizes AUROC improvement (if ties, longest such flip)"""
        indices = np.arange(self._total_length)[np.newaxis]
        delta_x = (indices - indices.T)  / (self._total_length-1)
        delta_y = self._curve - self._curve[np.newaxis].T
        rocs = self.roc() + (1-2*self.local_roc)*delta_y*delta_x
        
        new_roc = np.max(rocs)
        max_indices = np.argwhere(rocs == new_roc)
        longest_segment = np.argmax(np.diff(max_indices))
        interval = max_indices[longest_segment, :]
        return interval, new_roc
                
    def flip_segment(self, start, end):
        result = self._curve.copy()
        result[start:(end+1)] = self._curve[start] + self._curve[end] - np.flip(self._curve[start:(end+1)])
        return result
    
    def roc(self):
        return np.trapz(self._curve, dx = 1/(self._total_length-1))
        
    def roc_longest_segment(self):
        """Reduce ROC to 4 points [0,i,j,1] where i,j are start and end of the longest straight line (first interv)"""
        if not self.longest_null_segment():
            return None
        start, end = self.longest_null_segment()
        simplified_x = np.array([0, start, end, self._total_length-1]) / (self._total_length-1)
        simplified_y = np.array([0, self._curve[start], self._curve[end], 1])
        return np.trapz(simplified_y, simplified_x)

    def roc_excluding_longest(self):
        if not self.longest_null_segment():
            return None
        start, end = self.longest_null_segment()
        curve_excluding_segment = self._curve.copy()
        np.delete(curve_excluding_segment, range(start, end+1))
        curve_excluding_segment[end:] += self._curve[start] - self._curve[end]
        curve_excluding_segment /= np.max(curve_excluding_segment)
        return np.trapz(curve_excluding_segment, dx = 1 / (len(curve_excluding_segment)-1))
    
    def flipped_curve(self):
        start, end = self.flip_interval
        return self.flip_segment(*self.flip_interval)
    
    def fraction_null(self):
        if self.null_segments:
            return sum(e - s + 1 for s, e in self.null_segments) / self._total_length
        else:
            return 0
    
    def longest_null_segment(self):
        if self.null_segments:
            return self.null_segments[0]
        else:
            return None

    def fraction_linear(self):
        if self.linear_segments:
            return sum(e - s + 1 for s, e in self.linear_segments) / self._total_length
        else:
            return 0
    
    def longest_linear_segment(self):
        if self.linear_segments:
            return self.linear_segments[0]
        else:
            return None

    def recompute_null_segments(self):
        self.null_segments = self._extract_segments(self._is_null(), self.min_segment_length)
        
    def recompute_linear_segments(self):
        self.linear_segments = self._extract_segments(np.logical_and(self._is_quasilinear(), self._is_null()), self.min_quasilinear_length)
        
    def plot(self):
        result = sns.lineplot(x=np.linspace(0,1,len(self._curve)), y=self._curve)
        return result
        
    def plot_linearity(self, ref_fpr = 0, n_positives = None):
        if n_positives is None:
            n_positives = self.n_positives
        ref_point = round(ref_fpr*len(self._curve))
        ref_tpr = self._curve[ref_point]
        to_plot = pd.DataFrame(dict(
            fpr = np.linspace(0,1,len(self._curve)),
            ks_stat = self._linearity_deviation[ref_point,:]
        ))
        threshold = 1.358 / np.sqrt(abs(ref_tpr - self._curve)*n_positives)
        result = sns.lineplot(data=to_plot, x="fpr", y="ks_stat")
        result.plot(to_plot["fpr"], threshold, color="k", dashes=(2,1))
        return result

    def plot_normalized_linearity(self, ref_fpr = 0):
        ref_point = round(ref_fpr*len(self._curve))
        to_plot = pd.DataFrame(dict(
            fpr = np.linspace(0,1,len(self._curve)),
            ks_stat = self._normalized_linearity[ref_point,:]
        ))
        result = sns.lineplot(data=to_plot, x="fpr", y="ks_stat")
        result.axhline(1.358, color="k", dashes=(2,1))
        return result

    def plot_linear_segments(self):
        result = self.plot()
        L = len(self._curve)
        for s in self.linear_segments:
            result.plot([s[0]/L, s[1]/L], [self._curve[s[0]], self._curve[s[1]]], color="k", dashes=(2,1))
        return result

    def plot_linear_segment(self, i):
        s = self.linear_segments[i]
        result = sns.lineplot(x=np.linspace(0,1,len(self._curve))[s[0]:(s[1]+1)], y=self._curve[s[0]:(s[1]+1)])
        L = len(self._curve)
        result.plot([s[0]/L, s[1]/L], [self._curve[s[0]], self._curve[s[1]]], color="k", dashes=(2,1))
        return result

    def plot_local_auroc(self, ref_fpr = 0):
        ref_point = round(ref_fpr*len(self._curve))
        to_plot = self.local_roc[ref_point,:]
        result = sns.lineplot(x=np.linspace(0,1,len(to_plot)), y=to_plot)
        result.axhline(0.5, color="k", dashes = (2,1))
        return result
    

def compute_roc(predictor, positives, n_points = 200):
    result = roc_curve(positives, predictor)
    return(standardize_roc(result, n_points))
    
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

def compute_roc_old(predictor, positives, n_points = 200):
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

