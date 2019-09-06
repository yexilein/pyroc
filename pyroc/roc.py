
import numpy as np


class Roc(object):
    local_roc_threshold = 0.55
    min_segment_length = 0.05
    quasilinear_threshold = 0.005
    min_quasilinear_length = 0.1
    
    def __init__(self, curve):
        self._curve = curve
        self._total_length = len(self._curve)
        self.null_segments = self._extract_null_segments()[0]
        if self.null_segments:
            self.longest_null_segment = self.null_segments[0]
        else:
            self.longest_null_segment = None
        self.quasilinear_segments = self._extract_quasilinear_segments()[0]
        self.flip_interval, self.roc_flip = self._find_best_flip()
        
    def _extract_null_segments(self, abs_deviations = True):
        """Extract longest non-overlapping straight lines in curve"""
        roc = np.zeros((self._total_length, self._total_length))
        # TODO: replace that part of the code with (approximately)
        #  contribution[i] = (y[i] + y[i+1]) / 2
        #  result = cumsum(contribution)
        #  result = result - contribution.T
        #
        for start in range(self._total_length-1):
            for end in range(start+1, self._total_length):
                segment_length = end - start
                if (self._curve[end] - self._curve[start] > 0):
                    local_roc = np.trapz((self._curve[start:(end+1)] - self._curve[start]) /
                                         (self._curve[end] - self._curve[start]), dx = 1/segment_length)
                else:
                    local_roc = 0.5
                roc[start, end] = local_roc
                roc[end, start] = local_roc
        np.fill_diagonal(roc, 0.5)

        if abs_deviations:
            roc = abs(roc-0.5) + 0.5

        return self._extract_segments(roc, self.local_roc_threshold, self.min_segment_length)
        
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

    def _extract_quasilinear_segments(self, absolute_deviation = True):
        deviation = np.zeros((self._total_length, self._total_length))
        for start in range(self._total_length-1):
            for end in range(start+1, self._total_length):
                interpolated_curve = np.interp(np.arange(start, end+1), [start, end], self._curve[[start, end]])
                if absolute_deviation:
                    deviation[start, end] = max(abs(self._curve[start:(end+1)] - interpolated_curve))
                else:
                    deviation[start, end] = max(self._curve[start:(end+1)] - interpolated_curve)
                deviation[end, start] = deviation[start, end]
        return  self._extract_segments(deviation, self.quasilinear_threshold, self.min_quasilinear_length)
    
    def _find_best_flip(self):
        """Find flip that maximizes AUROC improvement (if ties, longest such flip)"""
        rocs = np.zeros((self._total_length, self._total_length))
        for start in range(self._total_length-1):
            for end in range(start+1, self._total_length):
                flipped_curve = self.flip_segment(start, end)
                rocs[start, end] = np.trapz(flipped_curve, dx = 1/(len(flipped_curve)-1))
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

