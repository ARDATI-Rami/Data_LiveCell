"""
Spatial operations and utilities.

This module provides functions for spatial queries, centroid computation,
and nearest-neighbor searches using KDTree.
"""

import numpy as np
from scipy.spatial import KDTree
from typing import Tuple, Optional


def compute_centroids(x: np.ndarray, y: np.ndarray,
                     roi_indices: np.ndarray,
                     num_rois: int) -> np.ndarray:
    """
    Compute centroids for each ROI.

    :param x: X coordinates of all points.
    :type x: np.ndarray
    :param y: Y coordinates of all points.
    :type y: np.ndarray
    :param roi_indices: ROI assignment for each point.
    :type roi_indices: np.ndarray
    :param num_rois: Total number of ROIs.
    :type num_rois: int
    :return: (num_rois, 2) array of centroids. NaN if ROI is empty.
    :rtype: np.ndarray
    """
    centroids = np.full((num_rois, 2), np.nan, dtype=float)

    for roi_id in range(num_rois):
        mask = roi_indices == roi_id
        if np.any(mask):
            centroids[roi_id, 0] = np.mean(x[mask])
            centroids[roi_id, 1] = np.mean(y[mask])

    return centroids


def build_kdtree(points: np.ndarray) -> KDTree:
    """
    Build a KDTree from point coordinates for fast spatial queries.

    :param points: (N, 2) array of point coordinates.
    :type points: np.ndarray
    :return: KDTree object.
    :rtype: KDTree
    """
    return KDTree(points)


def nearest_neighbors(tree: KDTree, query_points: np.ndarray,
                     k: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find k nearest neighbors for query points.

    :param tree: KDTree built from reference points.
    :type tree: KDTree
    :param query_points: (M, 2) array of query coordinates.
    :type query_points: np.ndarray
    :param k: Number of nearest neighbors to find.
    :type k: int
    :return: Tuple of (distances, indices) arrays.
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    distances, indices = tree.query(query_points, k=k)
    return distances, indices


def farthest_point_sampling(points: np.ndarray, k: int,
                            rng: Optional[np.random.Generator] = None) -> list:
    """
    Select k seed points using farthest-point sampling.

    This greedy algorithm selects points that are maximally separated,
    useful for initializing spatially distributed ROIs.

    :param points: (N, 2) array of point coordinates.
    :type points: np.ndarray
    :param k: Number of seed points to select.
    :type k: int
    :param rng: NumPy random generator for reproducibility.
    :type rng: Optional[np.random.Generator]
    :return: List of selected point indices.
    :rtype: list
    """
    if rng is None:
        rng = np.random.default_rng()

    n = points.shape[0]
    if k >= n:
        return list(range(n))

    # Start with a random point
    seeds = [rng.integers(0, n)]
    d2 = np.full(n, np.inf)

    for _ in range(1, k):
        last = seeds[-1]
        diff = points - points[last]
        d2 = np.minimum(d2, np.einsum('ij,ij->i', diff, diff))
        d2[seeds] = -1.0
        nxt = int(np.argmax(d2))
        seeds.append(nxt)

    return seeds


def radial_distance(x: np.ndarray, y: np.ndarray,
                   cx: float, cy: float) -> np.ndarray:
    """
    Compute radial distance from points to a center.

    :param x: X coordinates of points.
    :type x: np.ndarray
    :param y: Y coordinates of points.
    :type y: np.ndarray
    :param cx: X coordinate of center.
    :type cx: float
    :param cy: Y coordinate of center.
    :type cy: float
    :return: Array of distances.
    :rtype: np.ndarray
    """
    return np.hypot(x - cx, y - cy)


def radial_cv(x: np.ndarray, y: np.ndarray,
             cx: float, cy: float) -> float:
    """
    Compute coefficient of variation of radial distances.

    This measures circularity: lower CV = more circular.

    :param x: X coordinates of points.
    :type x: np.ndarray
    :param y: Y coordinates of points.
    :type y: np.ndarray
    :param cx: X coordinate of center.
    :type cx: float
    :param cy: Y coordinate of center.
    :type cy: float
    :return: Coefficient of variation (std/mean).
    :rtype: float
    """
    if len(x) < 3:
        return 1.0

    r = radial_distance(x, y, cx, cy)
    mu = np.mean(r)
    sd = np.std(r)

    return float(sd / max(mu, 1e-12))


def robust_outlier_threshold(distances: np.ndarray,
                             k_mad: float = 2.5) -> float:
    """
    Compute robust outlier threshold using Median Absolute Deviation (MAD).

    :param distances: Array of distance values.
    :type distances: np.ndarray
    :param k_mad: Multiplier for MAD (default 2.5 for moderate outlier detection).
    :type k_mad: float
    :return: Threshold value.
    :rtype: float
    """
    med = np.median(distances)
    mad = 1.4826 * np.median(np.abs(distances - med))
    scale = mad if mad > 1e-12 else (np.std(distances) if len(distances) > 1 else 0.0)
    return med + k_mad * scale

