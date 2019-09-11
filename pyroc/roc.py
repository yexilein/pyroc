
import numpy as np


class Roc(object):
    local_roc_threshold = 0.55
    min_segment_length = 0.05
    
    def __init__(self, curve):
        self._curve = curve
        self._total_length = len(self._curve)
        self.local_roc = self._create_local_auroc_matrix()
        self.null_segments = self._extract_null_segments()[0]
        if self.null_segments:
            self.longest_null_segment = self.null_segments[0]
        else:
            self.longest_null_segment = None
        self.flip_interval, self.roc_flip = self._find_best_flip()
            
    def _create_local_auroc_matrix(self):
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
    
    def _extract_null_segments(self, abs_deviations=True):
        """Extract longest non-overlapping straight lines in curve"""
        if abs_deviations:
            local_roc = abs(self.local_roc-0.5) + 0.5
        else:
            local_roc = self.local_roc

        return self._extract_segments(local_roc, self.local_roc_threshold, self.min_segment_length)
        
    def _extract_segments(self, stat, threshold, min_length):
        indices = np.arange(self._total_length)[np.newaxis]
        segment_length = indices - indices.T
        segment_length[stat > threshold] = 0
        intervals = []
        stats = []
        while np.max(segment_length) > 0:
            i, j = np.unravel_index(np.argmax(segment_length), segment_length.shape)
            if (j-i)/self._total_length + self._curve[j] - self._curve[i] >= min_length:
                intervals.append((i,j))
                stats.append(stat[i,j])
            segment_length[i:(j+1),:] = 0
            segment_length[:,i:(j+1)] = 0

        return intervals, stats

    def extract_quasilinear_segments(self, quasilinear_threshold=0.005, min_quasilinear_length=0.1, absolute_deviation=True):
        deviation = np.zeros((self._total_length, self._total_length))
        for start in range(self._total_length-1):
            for end in range(start+1, self._total_length):
                interpolated_curve = np.interp(np.arange(start, end+1), [start, end], self._curve[[start, end]])
                if absolute_deviation:
                    deviation[start, end] = max(abs(self._curve[start:(end+1)] - interpolated_curve))
                else:
                    deviation[start, end] = max(self._curve[start:(end+1)] - interpolated_curve)
                deviation[end, start] = deviation[start, end]
        return self._extract_segments(deviation, quasilinear_threshold, min_quasilinear_length)[0]
    
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
        if not self.longest_null_segment:
            return None
        start, end = self.longest_null_segment
        simplified_x = np.array([0, start, end, self._total_length-1]) / (self._total_length-1)
        simplified_y = np.array([0, self._curve[start], self._curve[end], 1])
        return np.trapz(simplified_y, simplified_x)

    def roc_excluding_longest(self):
        if not self.longest_null_segment:
            return None
        start, end = self.longest_null_segment
        curve_excluding_segment = self._curve.copy()
        np.delete(curve_excluding_segment, range(start, end+1))
        curve_excluding_segment[end:] += self._curve[start] - self._curve[end]
        curve_excluding_segment /= np.max(curve_excluding_segment)
        return np.trapz(curve_excluding_segment, dx = 1 / (len(curve_excluding_segment)-1))
    
    def flipped_curve(self):
        start, end = self.flip_interval
        return self.flip_segment(*self.flip_interval)

