import h5py
import pandas as pd
import anndata as ad

def get_df(path,name):
    with h5py.File(path, "r") as f:
        r = ad._io.h5ad.read_elem(f[name])
    return r

import numpy as np


import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def compute_overlap_matrices(
    dict_a: dict[str, list | set],
    dict_b: dict[str, list | set],
) -> tuple[np.ndarray, np.ndarray, list[str], list[str], pd.DataFrame]:
    """
    Compute pairwise overlap and Jaccard index matrices between two dicts of lists/sets,
    and find the optimal one-to-one matching using the Hungarian algorithm.

    Args:
        dict_a: First dict mapping names to lists or sets (rows)
        dict_b: Second dict mapping names to lists or sets (columns)

    Returns:
        overlap_matrix:  (len_a x len_b) int array — size of intersection for each pair
        jaccard_matrix:  (len_a x len_b) float array — Jaccard index for each pair
        keys_a:          row labels (keys of dict_a)
        keys_b:          column labels (keys of dict_b)
        matches_df:      DataFrame with best matches and statistics
    """
    keys_a = list(dict_a.keys())
    keys_b = list(dict_b.keys())

    sets_a = [set(dict_a[k]) for k in keys_a]
    sets_b = [set(dict_b[k]) for k in keys_b]

    n, m = len(keys_a), len(keys_b)
    overlap_matrix = np.zeros((n, m), dtype=int)
    jaccard_matrix = np.zeros((n, m), dtype=float)

    for i, sa in enumerate(sets_a):
        for j, sb in enumerate(sets_b):
            intersection = len(sa & sb)
            union = len(sa | sb)
            overlap_matrix[i, j] = intersection
            jaccard_matrix[i, j] = intersection / union if union > 0 else 0.0

    # Hungarian algorithm — maximise overlap (negate for minimisation)
    row_ind, col_ind = linear_sum_assignment(-overlap_matrix)

    matches_df = pd.DataFrame({
        "key_a":    [keys_a[i] for i in row_ind],
        "key_b":    [keys_b[j] for j in col_ind],
        "size_a":   [len(sets_a[i]) for i in row_ind],
        "size_b":   [len(sets_b[j]) for j in col_ind],
        "overlap":  [overlap_matrix[i, j] for i, j in zip(row_ind, col_ind)],
        "jaccard":  [jaccard_matrix[i, j] for i, j in zip(row_ind, col_ind)],
    }).sort_values("overlap", ascending=False).reset_index(drop=True)
    matches_df["szymkiewicz–simpson"] = matches_df["overlap"]/matches_df[["size_a", "size_b"]].min(axis=1)

    return overlap_matrix, jaccard_matrix, keys_a, keys_b, matches_df