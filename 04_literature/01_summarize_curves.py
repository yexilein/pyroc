
import sys
import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append("..")
import pyroc
import data


def main():
    export_stats("curves", "")

    
def export_stats(input_dir, file_prefix, output_dir="stats"):
    output_prefix = os.path.join(output_dir, file_prefix)
    curves = np.loadtxt(os.path.join(input_dir, file_prefix + "curves.txt"))
    n = 100
    stats = [pyroc.Roc.from_curve(c, n) for c in curves]
    np.savetxt(output_prefix + "all_null_segments.txt",
               np.vstack([np.c_[i*np.ones(len(s.null_segments)), np.array(s.null_segments)] for i, s in enumerate(stats) if s.null_segments]))
    np.savetxt(output_prefix + "all_linear_segments.txt",
               np.vstack([np.c_[i*np.ones(len(s.linear_segments)), np.array(s.linear_segments)] for i, s in enumerate(stats) if s.linear_segments]))
       
    df = pd.DataFrame({
        "auroc": [s.roc() for s in stats],
        "fraction_null": [s.fraction_null() for s in stats],
        "number_null_segments": [len(s.null_segments) for s in stats],
        "longest_null_segment": [int(np.diff(s.longest_null_segment())) if s.longest_null_segment() else 0 for s in stats],
        "fraction_linear": [s.fraction_linear() for s in stats],
        "number_linear_segments": [len(s.linear_segments) for s in stats],
        "longest_linear_segment": [int(np.diff(s.longest_linear_segment())) if s.longest_linear_segment() else 0 for s in stats],
        "reduced_curve_auroc": [s.roc_longest_segment() for s in stats],
        "roc_excluding_longest": [s.roc_excluding_longest() for s in stats],
        "flip_size": [int(np.diff(s.flip_interval)) for s in stats],
        "flip_start": [int(s.flip_interval[0]) for s in stats],
        "flip_end": [int(s.flip_interval[1]) for s in stats],
        "flip_performance_increase": [s.roc_flip - s.roc() for s in stats]
    })
    df.to_csv(output_prefix + "stats.csv")        

    
if __name__ == "__main__":
    main()