# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""
Eulerian (fixed-grid) ROI construction and analysis.

This module provides functions for creating and analyzing fixed spatial grid ROIs.
"""

import numpy as np
from typing import Tuple, List, Dict


def create_grid(x_min: float, x_max: float,
               y_min: float, y_max: float,
               num_grid_x: int, num_grid_y: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a fixed spatial grid for Eulerian ROIs.

    :param x_min: Minimum X coordinate.
    :type x_min: float
    :param x_max: Maximum X coordinate.
    :type x_max: float
    :param y_min: Minimum Y coordinate.
    :type y_min: float
    :param y_max: Maximum Y coordinate.
    :type y_max: float
    :param num_grid_x: Number of grid divisions along X axis.
    :type num_grid_x: int
    :param num_grid_y: Number of grid divisions along Y axis.
    :type num_grid_y: int
    :return: Tuple of (x_grid_edges, y_grid_edges).
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    x_grid_edges = np.linspace(x_min, x_max, num_grid_x + 1)
    y_grid_edges = np.linspace(y_min, y_max, num_grid_y + 1)
    return x_grid_edges, y_grid_edges


def compute_bounds_from_frames(all_x_coords: List[np.ndarray],
                               all_y_coords: List[np.ndarray]) -> Tuple[float, float, float, float]:
    """
    Compute overall min/max bounds from multiple frames.

    :param all_x_coords: List of X coordinate arrays (one per frame).
    :type all_x_coords: List[np.ndarray]
    :param all_y_coords: List of Y coordinate arrays (one per frame).
    :type all_y_coords: List[np.ndarray]
    :return: Tuple of (min_x, max_x, min_y, max_y).
    :rtype: Tuple[float, float, float, float]
    """
    min_x = min([x.min() for x in all_x_coords if len(x) > 0])
    max_x = max([x.max() for x in all_x_coords if len(x) > 0])
    min_y = min([y.min() for y in all_y_coords if len(y) > 0])
    max_y = max([y.max() for y in all_y_coords if len(y) > 0])
    return min_x, max_x, min_y, max_y


def assign_cells_to_grid(x_coords: np.ndarray, y_coords: np.ndarray,
                        x_grid_edges: np.ndarray, y_grid_edges: np.ndarray,
                        num_grid_x: int, num_grid_y: int) -> np.ndarray:
    """
    Assign cells to grid ROIs based on coordinates.

    :param x_coords: X coordinates of cells.
    :type x_coords: np.ndarray
    :param y_coords: Y coordinates of cells.
    :type y_coords: np.ndarray
    :param x_grid_edges: X edges of grid (length = num_grid_x + 1).
    :type x_grid_edges: np.ndarray
    :param y_grid_edges: Y edges of grid (length = num_grid_y + 1).
    :type y_grid_edges: np.ndarray
    :param num_grid_x: Number of grid divisions along X.
    :type num_grid_x: int
    :param num_grid_y: Number of grid divisions along Y.
    :type num_grid_y: int
    :return: Array of ROI indices for each cell.
    :rtype: np.ndarray
    """
    # Digitize coordinates to grid indices
    x_indices = np.digitize(x_coords, x_grid_edges) - 1
    y_indices = np.digitize(y_coords, y_grid_edges) - 1
    
    # Clip to valid range
    x_indices = np.clip(x_indices, 0, num_grid_x - 1)
    y_indices = np.clip(y_indices, 0, num_grid_y - 1)
    
    # Compute flat ROI index
    roi_indices = x_indices * num_grid_y + y_indices
    
    return roi_indices


def initialize_roi_data(num_rois: int) -> List[Dict]:
    """
    Initialize ROI data storage.

    :param num_rois: Number of ROIs.
    :type num_rois: int
    :return: List of dictionaries with 'polygon_numbers' and 'cell_counts' lists.
    :rtype: List[Dict]
    """
    roi_data_list = []
    for _ in range(num_rois):
        roi_data_list.append({
            'polygon_numbers': [],
            'cell_counts': []
        })
    return roi_data_list


def collect_roi_data(polygon_numbers: np.ndarray,
                    roi_indices: np.ndarray,
                    num_rois: int) -> Tuple[List[np.ndarray], List[int]]:
    """
    Collect polygon numbers and counts for each ROI.

    :param polygon_numbers: Polygon numbers for all cells.
    :type polygon_numbers: np.ndarray
    :param roi_indices: ROI assignment for each cell.
    :type roi_indices: np.ndarray
    :param num_rois: Total number of ROIs.
    :type num_rois: int
    :return: Tuple of (roi_polygons, roi_counts) lists.
    :rtype: Tuple[List[np.ndarray], List[int]]
    """
    roi_polygons = []
    roi_counts = []
    
    for roi_id in range(num_rois):
        mask = roi_indices == roi_id
        roi_poly = polygon_numbers[mask]
        roi_polygons.append(roi_poly)
        roi_counts.append(len(roi_poly))
    
    return roi_polygons, roi_counts


def generate_roi_labels(num_grid_x: int, num_grid_y: int) -> List[str]:
    """
    Generate ROI labels in format "ROI (i, j)".

    :param num_grid_x: Number of grid divisions along X.
    :type num_grid_x: int
    :param num_grid_y: Number of grid divisions along Y.
    :type num_grid_y: int
    :return: List of ROI labels.
    :rtype: List[str]
    """
    labels = []
    for i in range(num_grid_x):
        for j in range(num_grid_y):
            labels.append(f'ROI ({i}, {j})')
    return labels


class EulerianGrid:
    """
    Eulerian (fixed spatial grid) ROI manager.
    
    This class encapsulates the creation and management of a fixed spatial grid
    for analyzing cell topology at specific tissue locations.
    """
    
    def __init__(self, num_grid_x: int = 20, num_grid_y: int = 10):
        """
        Initialize Eulerian grid.
        
        :param num_grid_x: Number of grid divisions along X axis (default 20).
        :type num_grid_x: int
        :param num_grid_y: Number of grid divisions along Y axis (default 10).
        :type num_grid_y: int
        """
        self.num_grid_x = num_grid_x
        self.num_grid_y = num_grid_y
        self.num_rois = num_grid_x * num_grid_y
        
        self.x_grid_edges = None
        self.y_grid_edges = None
        self.bounds = None
        self.roi_data_list = initialize_roi_data(self.num_rois)
        self.roi_labels = generate_roi_labels(num_grid_x, num_grid_y)
    
    def fit(self, all_x_coords: List[np.ndarray], all_y_coords: List[np.ndarray]):
        """
        Fit the grid to data bounds.
        
        :param all_x_coords: List of X coordinate arrays (one per frame).
        :type all_x_coords: List[np.ndarray]
        :param all_y_coords: List of Y coordinate arrays (one per frame).
        :type all_y_coords: List[np.ndarray]
        """
        self.bounds = compute_bounds_from_frames(all_x_coords, all_y_coords)
        min_x, max_x, min_y, max_y = self.bounds
        self.x_grid_edges, self.y_grid_edges = create_grid(
            min_x, max_x, min_y, max_y, self.num_grid_x, self.num_grid_y
        )
    
    def assign(self, x_coords: np.ndarray, y_coords: np.ndarray) -> np.ndarray:
        """
        Assign cells to grid ROIs.
        
        :param x_coords: X coordinates of cells.
        :type x_coords: np.ndarray
        :param y_coords: Y coordinates of cells.
        :type y_coords: np.ndarray
        :return: ROI indices for each cell.
        :rtype: np.ndarray
        """
        if self.x_grid_edges is None:
            raise ValueError("Grid not fitted. Call fit() first.")
        
        return assign_cells_to_grid(
            x_coords, y_coords,
            self.x_grid_edges, self.y_grid_edges,
            self.num_grid_x, self.num_grid_y
        )
    
    def collect(self, polygon_numbers: np.ndarray, roi_indices: np.ndarray):
        """
        Collect polygon data for this frame.
        
        :param polygon_numbers: Polygon numbers for all cells.
        :type polygon_numbers: np.ndarray
        :param roi_indices: ROI assignment for each cell.
        :type roi_indices: np.ndarray
        """
        roi_polygons, roi_counts = collect_roi_data(
            polygon_numbers, roi_indices, self.num_rois
        )
        
        for roi_id in range(self.num_rois):
            self.roi_data_list[roi_id]['polygon_numbers'].append(roi_polygons[roi_id])
            self.roi_data_list[roi_id]['cell_counts'].append(roi_counts[roi_id])

