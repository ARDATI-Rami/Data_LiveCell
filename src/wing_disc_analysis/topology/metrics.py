# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""
Topology metrics and statistical analysis.

This module provides functions for computing topological metrics on cell shape
distributions, including Jensen-Shannon divergence, Wasserstein distance, and
topological charge.
"""

import numpy as np
from scipy.stats import entropy
from typing import Dict, List, Tuple, Optional


def infer_support(all_data: List[np.ndarray],
                 min_global_pct: float = 0.0) -> np.ndarray:
    """
    Infer polygon classes present globally and filter rare classes.

    :param all_data: List of 1D arrays (per frame) with polygon numbers (e.g., 3..10).
    :type all_data: List[np.ndarray]
    :param min_global_pct: Minimum global percentage to keep a class (0..100).
    :type min_global_pct: float
    :return: Sorted unique polygon classes retained.
    :rtype: np.ndarray
    """
    if not all_data:
        return np.array([], dtype=int)
    
    g = np.concatenate(all_data)
    if g.size == 0:
        return np.array([], dtype=int)
    
    classes, counts = np.unique(g, return_counts=True)
    pct = 100 * counts / counts.sum()
    keep = pct >= min_global_pct
    
    return classes[keep]


def tally_roi_counts(roi_data_list: List[Dict],
                    support: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert ROI polygon lists to count matrices over a fixed support.

    :param roi_data_list: List of ROI dicts, each with 'polygon_numbers' as List[np.ndarray].
    :type roi_data_list: List[Dict]
    :param support: Sorted polygon classes to count (e.g., np.array([3,4,5,6,7])).
    :type support: np.ndarray
    :return: Tuple (C, N), where C[i, k] = counts of support[k] in ROI i,
             N[i] = total counts in ROI i.
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    num_rois = len(roi_data_list)
    K = len(support)
    C = np.zeros((num_rois, K), dtype=np.int64)
    
    for i, roi in enumerate(roi_data_list):
        # Concatenate all frames for this ROI
        arr = np.concatenate(roi.get('polygon_numbers', [])) if roi.get('polygon_numbers') else np.array([], dtype=int)
        if arr.size:
            # Vectorized counting per class
            C[i, :] = np.array([(arr == cls).sum() for cls in support], dtype=np.int64)
    
    N = C.sum(axis=1)
    return C, N


def dirichlet_posterior_mean(counts: np.ndarray, alpha: float = 0.5) -> np.ndarray:
    """
    Compute Dirichlet posterior mean probabilities with symmetric prior.

    Uses Bayesian smoothing to estimate probabilities, avoiding zero probabilities.
    Jeffreys prior (α=0.5) is the default; α=1.0 gives Laplace smoothing.

    :param counts: Count vector or matrix [..., K] for K classes.
    :type counts: np.ndarray
    :param alpha: Symmetric Dirichlet prior concentration per class (default 0.5).
    :type alpha: float
    :return: Posterior mean probabilities with shape matching counts.
    :rtype: np.ndarray
    """
    counts = np.asarray(counts, dtype=np.float64)
    K = counts.shape[-1]
    num = counts + alpha
    den = counts.sum(axis=-1, keepdims=True) + alpha * K
    return num / np.clip(den, 1e-12, None)


def js_divergence(P: np.ndarray, Q: np.ndarray, base: float = 2.0) -> np.ndarray:
    """
    Jensen–Shannon divergence between rows of P and a single Q.

    JSD is a symmetric, bounded measure of distribution similarity.
    Returns 0 for identical distributions, up to log₂(2)=1 bit for completely different.

    :param P: Array [M, K] of probability vectors (one per ROI).
    :type P: np.ndarray
    :param Q: Array [K] global probability vector.
    :type Q: np.ndarray
    :param base: Logarithm base (2 for bits, e for nats).
    :type base: float
    :return: JSD values per row of P.
    :rtype: np.ndarray
    """
    Q = np.asarray(Q, dtype=np.float64)
    P = np.asarray(P, dtype=np.float64)
    M = 0.5 * (P + Q[None, :])
    
    # Compute KL(P||M) + KL(Q||M) for each row
    jsd = np.empty(P.shape[0], dtype=np.float64)
    for i in range(P.shape[0]):
        jsd[i] = 0.5 * entropy(P[i], M[i], base=base) + 0.5 * entropy(Q, M[i], base=base)
    
    return jsd


def wasserstein1_discrete(P: np.ndarray, Q: np.ndarray, support: np.ndarray) -> np.ndarray:
    """
    1D Wasserstein-1 distance (Earth Mover's Distance) between distributions.

    For discrete measures on ordered points x₁ < ... < xₖ:
    W₁ = Σ |CDF_P(k) - CDF_Q(k)| * (x_{k+1} - x_k)

    :param P: Array [M, K] of probability vectors (ROIs).
    :type P: np.ndarray
    :param Q: Array [K] global probability vector.
    :type Q: np.ndarray
    :param support: Sorted support values (polygon classes, e.g., [3,4,5,6,7]).
    :type support: np.ndarray
    :return: Wasserstein-1 distance per row of P (units = support units, i.e., "sides").
    :rtype: np.ndarray
    """
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    
    cdf_P = np.cumsum(P, axis=1)
    cdf_Q = np.cumsum(Q)
    gaps = np.diff(support.astype(np.float64))
    
    # |CDF difference| at internal cutpoints times gap sizes
    return np.sum(np.abs(cdf_P[:, :-1] - cdf_Q[:-1][None, :]) * gaps[None, :], axis=1)


def multinomial_residuals(counts: np.ndarray, p0: np.ndarray) -> np.ndarray:
    """
    Per-class standardized residuals under multinomial(n, p0).

    z_k = (x_k - n·p0_k) / sqrt(n·p0_k·(1 - p0_k))

    These are not independent across classes but useful for diagnostics.

    :param counts: Array [M, K] of counts per ROI.
    :type counts: np.ndarray
    :param p0: Array [K] of baseline probabilities.
    :type p0: np.ndarray
    :return: Array [M, K] of z-scores; NaN for n=0 or p0_k≈0.
    :rtype: np.ndarray
    """
    counts = counts.astype(np.float64)
    n = counts.sum(axis=1, keepdims=True)
    exp = n * p0[None, :]
    var = n * p0[None, :] * (1.0 - p0[None, :])
    
    z = (counts - exp) / np.sqrt(np.maximum(var, 1e-12))
    z[np.where(n == 0)[0], :] = np.nan
    z[:, p0 < 1e-9] = np.nan
    
    return z


def topological_charge(counts: np.ndarray, support: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Mean topological charge E[n-6] and E[|n-6|] per ROI.

    In hexagonal packing, cells have 6 neighbors (n=6). Charge measures
    deviation from this ideal.

    :param counts: Array [M, K] of counts per ROI.
    :type counts: np.ndarray
    :param support: Sorted polygon classes (e.g., [3,4,5,6,7]).
    :type support: np.ndarray
    :return: Tuple (mean_q, mean_abs_q) each [M].
             mean_q = E[n-6], mean_abs_q = E[|n-6|].
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    n = counts.sum(axis=1, keepdims=True)
    
    with np.errstate(invalid='ignore', divide='ignore'):
        P = counts / np.clip(n, 1e-12, None)
    
    q = support.astype(np.float64) - 6.0
    mean_q = np.nansum(P * q[None, :], axis=1)
    mean_abs_q = np.nansum(P * np.abs(q)[None, :], axis=1)
    
    mean_q[n.ravel() == 0] = np.nan
    mean_abs_q[n.ravel() == 0] = np.nan
    
    return mean_q, mean_abs_q


def compute_topology_deviation_metrics(
    roi_data_list: List[Dict],
    all_polygon_numbers: List[np.ndarray],
    polygons: Optional[np.ndarray] = None,
    alpha: float = 0.5,
    min_global_pct: float = 0.0
) -> Dict[str, np.ndarray]:
    """
    Compute comprehensive topology deviation metrics per ROI.

    This function:
      1) Pools counts across frames
      2) Applies Dirichlet(α) smoothing (Jeffreys α≈0.5 by default)
      3) Compares ROI posterior means to global posterior mean with:
         - Jensen–Shannon divergence (bits)
         - 1D Wasserstein distance over polygon classes (sides)
         - Per-class standardized residuals (z-scores)
         - Mean topological charge E[n-6] and E[|n-6|]

    :param roi_data_list: Each ROI dict must have 'polygon_numbers': List[np.ndarray].
    :type roi_data_list: List[Dict]
    :param all_polygon_numbers: Per-frame arrays of polygon numbers (whole tissue).
    :type all_polygon_numbers: List[np.ndarray]
    :param polygons: Optional sorted array of polygon classes; if None, inferred globally.
    :type polygons: Optional[np.ndarray]
    :param alpha: Symmetric Dirichlet prior concentration per class.
    :type alpha: float
    :param min_global_pct: Drop classes with global percentage below this threshold.
    :type min_global_pct: float
    :return: Dictionary with keys:
             'support', 'C', 'N', 'p_roi', 'p_global', 'jsd', 'w1', 'z',
             'mean_q', 'mean_abs_q'.
    :rtype: Dict[str, np.ndarray]
    """
    # 1) Determine support
    if polygons is None:
        support = infer_support(all_polygon_numbers, min_global_pct=min_global_pct)
    else:
        support = np.asarray(polygons, dtype=int)
    
    if support.size == 0:
        raise ValueError("No polygon classes found/kept. Check inputs or min_global_pct.")
    
    # 2) Tally counts per ROI and globally
    C, N = tally_roi_counts(roi_data_list, support)
    G_counts = C.sum(axis=0)
    
    # 3) Posterior means with Dirichlet smoothing
    p_roi = dirichlet_posterior_mean(C, alpha=alpha)
    p_global = dirichlet_posterior_mean(G_counts, alpha=alpha)
    
    # 4) Compute metrics
    jsd = js_divergence(p_roi, p_global, base=2.0)
    w1 = wasserstein1_discrete(p_roi, p_global, support)
    z = multinomial_residuals(C, p_global)
    mean_q, mean_abs_q = topological_charge(C, support)
    
    return dict(
        support=support,
        C=C,
        N=N,
        p_roi=p_roi,
        p_global=p_global,
        jsd=jsd,
        w1=w1,
        z=z,
        mean_q=mean_q,
        mean_abs_q=mean_abs_q
    )

