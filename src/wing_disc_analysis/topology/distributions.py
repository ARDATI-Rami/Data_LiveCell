"""
Distribution analysis and global statistics.

This module provides functions for computing polygon shape distributions
and global statistics across frames.
"""

import numpy as np
from typing import List, Tuple, Dict


# Standard shape names
SHAPE_NAMES = {
    3: 'Triangle',
    4: 'Quadrilateral',
    5: 'Pentagon',
    6: 'Hexagon',
    7: 'Heptagon',
    8: 'Octagon',
    9: 'Nonagon',
    10: 'Decagon'
}


def get_shape_name(n: int) -> str:
    """
    Get the name of a polygon with n sides.

    :param n: Number of sides.
    :type n: int
    :return: Shape name (e.g., "Hexagon" for 6).
    :rtype: str
    """
    return SHAPE_NAMES.get(n, f'{n}-gon')


def compute_global_distribution(all_polygon_numbers: List[np.ndarray],
                                threshold: float = 0.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute global mean and std of polygon distribution across frames.

    Filters out rare shapes below the threshold percentage.

    :param all_polygon_numbers: List of per-frame polygon arrays.
    :type all_polygon_numbers: List[np.ndarray]
    :param threshold: Minimum mean percentage to keep a polygon type.
    :type threshold: float
    :return: Tuple of (polygons, mean_pct, std_pct) filtered.
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    if not all_polygon_numbers:
        return np.array([]), np.array([]), np.array([])
    
    # Get all unique polygon types across frames
    all_data_flat = np.concatenate(all_polygon_numbers)
    unique = np.unique(all_data_flat)
    
    mean_pct = []
    std_pct = []
    
    for polygon in unique:
        percentages = []
        for data in all_polygon_numbers:
            if len(data) > 0:
                pct = (data == polygon).sum() / len(data) * 100
                percentages.append(pct)
            else:
                percentages.append(0.0)
        
        mean_pct.append(np.mean(percentages))
        std_pct.append(np.std(percentages))
    
    unique = np.array(unique, dtype=int)
    mean_pct = np.array(mean_pct)
    std_pct = np.array(std_pct)
    
    # Filter by threshold
    mask = mean_pct > threshold
    return unique[mask], mean_pct[mask], std_pct[mask]


def compute_frame_distribution(polygon_numbers: np.ndarray) -> Dict[int, float]:
    """
    Compute polygon distribution for a single frame.

    :param polygon_numbers: Array of polygon numbers for one frame.
    :type polygon_numbers: np.ndarray
    :return: Dictionary mapping polygon number to percentage.
    :rtype: Dict[int, float]
    """
    if len(polygon_numbers) == 0:
        return {}
    
    unique, counts = np.unique(polygon_numbers, return_counts=True)
    total = counts.sum()
    
    distribution = {}
    for poly, count in zip(unique, counts):
        distribution[int(poly)] = (count / total) * 100
    
    return distribution


def mean_sidedness_per_roi(metrics: dict,
                           posterior: bool = True,
                           empty_policy: str = "nan") -> np.ndarray:
    """
    Compute mean polygon sidedness E[n] per ROI.

    For ROIs with zero total count, handles according to empty_policy.

    :param metrics: Output dict from compute_topology_deviation_metrics.
    :type metrics: dict
    :param posterior: If True, use Dirichlet-smoothed p_roi; else use empirical.
    :type posterior: bool
    :param empty_policy: How to handle N==0 ROIs:
                         "nan" (default), "global", "uniform", or numeric string.
    :type empty_policy: str
    :return: Array [R] with mean sidedness per ROI.
    :rtype: np.ndarray
    """
    support = metrics['support'].astype(float)
    
    if posterior:
        mean_n = metrics['p_roi'] @ support
    else:
        mean_n = 6.0 + metrics['mean_q']  # NaN when N==0
    
    # Identify empty ROIs
    empties = (metrics['N'] == 0)
    
    if empty_policy == "nan":
        mean_n = mean_n.astype(float)
        mean_n[empties] = np.nan
    elif empty_policy == "global":
        g_mean = float(metrics['p_global'] @ support)
        mean_n = mean_n.astype(float)
        mean_n[empties] = g_mean
    elif empty_policy == "uniform":
        # Keep as-is (posterior uniform → mean(support))
        pass
    else:
        # Try numeric override
        try:
            fill = float(empty_policy)
            mean_n = mean_n.astype(float)
            mean_n[empties] = fill
        except Exception:
            # Fallback to NaN
            mean_n = mean_n.astype(float)
            mean_n[empties] = np.nan
    
    return mean_n


def compute_roi_deviations(roi_data_list: List[Dict],
                          filtered_polygons: np.ndarray,
                          global_mean: np.ndarray) -> np.ndarray:
    """
    Compute absolute deviations of ROI distributions from global.

    :param roi_data_list: Per-ROI collected polygons over time.
    :type roi_data_list: List[Dict]
    :param filtered_polygons: Polygon categories used.
    :type filtered_polygons: np.ndarray
    :param global_mean: Global mean percentages aligned with filtered_polygons.
    :type global_mean: np.ndarray
    :return: Deviations[roi, shape_index] as absolute percentage differences.
    :rtype: np.ndarray
    """
    n_roi = len(roi_data_list)
    n_shapes = len(filtered_polygons)
    deviations = np.zeros((n_roi, n_shapes), dtype=float)
    
    global_dict = {int(pg): float(mu) for pg, mu in zip(filtered_polygons, global_mean)}
    
    for r in range(n_roi):
        # Compute mean percentage within ROI over frames
        means = []
        for pg in filtered_polygons:
            percs = []
            for arr in roi_data_list[r]['polygon_numbers']:
                if arr.size > 0:
                    percs.append((arr == pg).sum() / len(arr) * 100)
            means.append(np.mean(percs) if percs else 0.0)
        
        roi_vec = np.array(means, dtype=float)
        glob_vec = np.array([global_dict[int(pg)] for pg in filtered_polygons], dtype=float)
        deviations[r, :] = np.abs(roi_vec - glob_vec)
    
    return deviations

