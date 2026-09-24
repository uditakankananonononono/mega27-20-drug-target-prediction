"""DTI benchmark metrics: MSE and the DeepDTA concordance index."""
from __future__ import annotations

import numpy as np


def mse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.mean((np.asarray(y_true) - np.asarray(y_pred)) ** 2))


def concordance_index(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """CI over comparable pairs (DeepDTA definition): pair (i,j) comparable
    when y_i != y_j; concordant if sign(y_i - y_j) == sign(p_i - p_j)."""
    y = np.asarray(y_true, dtype=float)
    p = np.asarray(y_pred, dtype=float)
    n = len(y)
    total = 0.0
    score = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dy = y[i] - y[j]
            if dy == 0:
                continue
            dp = p[i] - p[j]
            total += 1
            if dy * dp > 0:
                score += 1
            elif dp == 0:
                score += 0.5
    return score / total if total else float("nan")


def concordance_index_fast(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """O(n log n) CI via pairwise-difference sign counting on sorted values."""
    y = np.asarray(y_true, dtype=float)
    p = np.asarray(y_pred, dtype=float)
    total = score = 0.0
    # vectorized per-i against suffix after sorting by y
    order = np.argsort(y, kind="stable")
    ys, ps = y[order], p[order]
    for i in range(len(ys)):
        gt = ys[i + 1:] > ys[i]  # strictly greater true affinity
        if not gt.any():
            continue
        pj = ps[i + 1:][gt]
        total += gt.sum()
        score += (pj > ps[i]).sum() + 0.5 * (pj == ps[i]).sum()
    return score / total if total else float("nan")
