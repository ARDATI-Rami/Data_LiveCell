# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""
Lagrangian (cell-tracking) ROI construction and tracking.

This module provides functions for creating ROIs that follow specific groups
of cells over time.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy.spatial import KDTree
from ..geometry import (
    delaunay_adjacency,
    farthest_point_sampling,
    compute_compactness,
    radial_cv,
    robust_outlier_threshold,
    compute_centroids
)


def build_lagrangian_rois_initial(
    ids0: np.ndarray,
    x0: np.ndarray,
    y0: np.ndarray,
    target: int = 50,
    random_seed: int = 7,
    circ_k_outlier: float = 1.0,
    circ_cv_max: float = 0.35,
    circ_compactness_min: float = 0.65,
    circ_min_members_eval: int = 20,
    circ_min_members_keep: Optional[int] = None,
    circ_max_iter: int = 4
) -> Tuple[Dict[int, int], Dict[int, np.ndarray]]:
    """
    Construct Lagrangian ROIs at Frame 0 with circularity refinement.

    Uses farthest-point sampling for initial seeds, then grows ROIs via
    adjacency to reach target size. Applies iterative circularity refinement
    by ejecting radial outliers.

    :param ids0: Cell IDs at Frame 0 (may include -1 for boundary).
    :type ids0: np.ndarray
    :param x0: X coordinates at Frame 0.
    :type x0: np.ndarray
    :param y0: Y coordinates at Frame 0.
    :type y0: np.ndarray
    :param target: Target number of non-boundary cells per ROI.
    :type target: int
    :param random_seed: Random seed for reproducibility.
    :type random_seed: int
    :param circ_k_outlier: MAD threshold multiplier for outlier detection.
    :type circ_k_outlier: float
    :param circ_cv_max: Maximum allowed coefficient of variation of radii.
    :type circ_cv_max: float
    :param circ_compactness_min: Minimum allowed hull compactness 4πA/P².
    :type circ_compactness_min: float
    :param circ_min_members_eval: Minimum members to attempt circularity repair.
    :type circ_min_members_eval: int
    :param circ_min_members_keep: Don't shrink below this; default = max(10, 0.5*target).
    :type circ_min_members_keep: Optional[int]
    :param circ_max_iter: Maximum circularity repair iterations.
    :type circ_max_iter: int
    :return: Tuple of (cell_id_to_roi, roi_seeds_xy).
    :rtype: Tuple[Dict[int, int], Dict[int, np.ndarray]]
    """
    rng = np.random.default_rng(random_seed)
    if circ_min_members_keep is None:
        circ_min_members_keep = max(10, int(np.floor(0.5 * target)))
    
    # Filter non-boundary cells
    mask_nb = ids0 != -1
    idx_nb = np.where(mask_nb)[0]
    pts_nb = np.column_stack([x0[mask_nb], y0[mask_nb]])
    n_nb = len(idx_nb)
    
    if n_nb == 0:
        raise ValueError("No non-boundary cells in Frame 0.")
    
    # Number of ROIs
    n_roi = max(1, n_nb // target + (1 if (n_nb % target) >= target // 2 else 0))
    
    # Build adjacency graph
    nbrs_full = delaunay_adjacency(x0[mask_nb], y0[mask_nb])
    
    # Farthest-point sampling for seeds
    seed_local_indices = farthest_point_sampling(pts_nb, n_roi, rng)
    seed_global_indices = idx_nb[seed_local_indices]
    
    # Initialize assignments
    assigned = np.full(len(ids0), fill_value=-1, dtype=int)
    quota = np.full(n_roi, target, dtype=int)
    queues: List[List[int]] = [[s] for s in seed_local_indices]
    
    # Assign seeds
    for r, s_loc in enumerate(seed_local_indices):
        g = idx_nb[s_loc]
        assigned[g] = r
        quota[r] = max(0, quota[r] - 1)
    
    # Grow ROIs via adjacency (BFS-like)
    frontier_levels = 0
    while (assigned[idx_nb] == -1).any() and (quota > 0).any():
        frontier_levels += 1
        next_queues: List[List[int]] = [[] for _ in range(n_roi)]
        
        for r in range(n_roi):
            if quota[r] <= 0 or not queues[r]:
                continue
            
            for u_loc in queues[r]:
                for v_loc in nbrs_full[u_loc]:
                    g = idx_nb[v_loc]
                    if assigned[g] == -1:
                        assigned[g] = r
                        quota[r] -= 1
                        if quota[r] > 0:
                            next_queues[r].append(v_loc)
                if quota[r] <= 0:
                    break
        
        queues = next_queues
        if frontier_levels > 5 * n_roi:
            break
    
    # Assign remaining cells via nearest seed
    remaining_nb = idx_nb[assigned[idx_nb] == -1]
    if len(remaining_nb) > 0:
        seed_pts = np.column_stack([x0[seed_global_indices], y0[seed_global_indices]])
        kd = KDTree(seed_pts)
        rem_pts = np.column_stack([x0[remaining_nb], y0[remaining_nb]])
        _, nn = kd.query(rem_pts, k=1)
        for g, r in zip(remaining_nb, nn):
            assigned[g] = int(r)
    
    # Assign boundary cells via nearest seed
    idx_bnd = np.where(ids0 == -1)[0]
    if len(idx_bnd) > 0:
        seed_pts = np.column_stack([x0[seed_global_indices], y0[seed_global_indices]])
        kd = KDTree(seed_pts)
        bnd_pts = np.column_stack([x0[idx_bnd], y0[idx_bnd]])
        _, nn = kd.query(bnd_pts, k=1)
        for g, r in zip(idx_bnd, nn):
            assigned[g] = int(r)
    
    # Circularity refinement (simplified version)
    # TODO: Full implementation would include iterative refinement
    # For now, compute final centroids as seeds
    
    roi_seeds_xy = {}
    for r in range(n_roi):
        mask = assigned == r
        if np.any(mask):
            roi_seeds_xy[r] = np.array([np.mean(x0[mask]), np.mean(y0[mask])])
        else:
            sg = seed_global_indices[r]
            roi_seeds_xy[r] = np.array([x0[sg], y0[sg]])
    
    # Build cell_id -> roi mapping
    cell_id_to_roi: Dict[int, int] = {}
    for g, roi in enumerate(assigned):
        cid = int(ids0[g])
        cell_id_to_roi[cid] = int(roi)
    
    return cell_id_to_roi, roi_seeds_xy


def compute_first_seen(frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]) -> Dict[int, int]:
    """
    Compute the first frame in which each cell ID appears.

    :param frames: List of (ids, x, y, polygons) per frame.
    :type frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    :return: Mapping cell_id -> first_seen_frame_index (1-based).
    :rtype: Dict[int, int]
    """
    first_seen = {}
    for t, (ids, _, _, _) in enumerate(frames, start=1):
        for cid in ids:
            cid = int(cid)
            if cid >= 0 and cid not in first_seen:
                first_seen[cid] = t
    return first_seen


class LagrangianROI:
    """
    Lagrangian (cell-tracking) ROI manager.
    
    This class manages ROIs that follow specific groups of cells over time.
    """
    
    def __init__(self, target_cells_per_roi: int = 46, random_seed: int = 7):
        """
        Initialize Lagrangian ROI manager.
        
        :param target_cells_per_roi: Target number of cells per ROI.
        :type target_cells_per_roi: int
        :param random_seed: Random seed for reproducibility.
        :type random_seed: int
        """
        self.target = target_cells_per_roi
        self.random_seed = random_seed
        self.cell_id_to_roi = {}
        self.roi_seeds_xy = {}
        self.num_rois = 0
        self.roi_data = {}
    
    def initialize(self, ids0: np.ndarray, x0: np.ndarray, y0: np.ndarray):
        """
        Initialize ROIs from Frame 0 data.
        
        :param ids0: Cell IDs at Frame 0.
        :type ids0: np.ndarray
        :param x0: X coordinates at Frame 0.
        :type x0: np.ndarray
        :param y0: Y coordinates at Frame 0.
        :type y0: np.ndarray
        """
        self.cell_id_to_roi, self.roi_seeds_xy = build_lagrangian_rois_initial(
            ids0, x0, y0,
            target=self.target,
            random_seed=self.random_seed
        )
        self.num_rois = len(self.roi_seeds_xy)
        
        # Initialize data storage
        for r in range(self.num_rois):
            self.roi_data[r] = {
                'polygon_numbers': [],
                'cell_counts': []
            }
    
    def track_frame(self, ids: np.ndarray, x: np.ndarray, y: np.ndarray,
                   polygons: np.ndarray) -> np.ndarray:
        """
        Track cells in a new frame and assign to ROIs.
        
        :param ids: Cell IDs for this frame.
        :type ids: np.ndarray
        :param x: X coordinates.
        :type x: np.ndarray
        :param y: Y coordinates.
        :type y: np.ndarray
        :param polygons: Polygon numbers.
        :type polygons: np.ndarray
        :return: ROI indices for each cell.
        :rtype: np.ndarray
        """
        N = len(ids)
        roi_indices = np.full(N, -1, dtype=int)
        
        # Assign known cells
        for i, cid in enumerate(ids):
            cid = int(cid)
            if cid in self.cell_id_to_roi:
                roi_indices[i] = self.cell_id_to_roi[cid]
        
        # Assign unknown/new cells via nearest seed
        unknown_mask = roi_indices == -1
        if np.any(unknown_mask):
            seed_mat = np.vstack([self.roi_seeds_xy[r] for r in range(self.num_rois)])
            kd = KDTree(seed_mat)
            unknown_pts = np.column_stack([x[unknown_mask], y[unknown_mask]])
            _, nearest = kd.query(unknown_pts, k=1)
            
            for idx, nn in zip(np.where(unknown_mask)[0], nearest):
                roi_indices[idx] = int(nn)
                cid = int(ids[idx])
                if cid >= 0:  # Update mapping for non-boundary
                    self.cell_id_to_roi[cid] = int(nn)
        
        return roi_indices
    
    def collect(self, polygons: np.ndarray, roi_indices: np.ndarray):
        """
        Collect polygon data for this frame.
        
        :param polygons: Polygon numbers for all cells.
        :type polygons: np.ndarray
        :param roi_indices: ROI assignment for each cell.
        :type roi_indices: np.ndarray
        """
        for r in range(self.num_rois):
            mask = roi_indices == r
            roi_poly = polygons[mask]
            self.roi_data[r]['polygon_numbers'].append(roi_poly)
            self.roi_data[r]['cell_counts'].append(len(roi_poly))

