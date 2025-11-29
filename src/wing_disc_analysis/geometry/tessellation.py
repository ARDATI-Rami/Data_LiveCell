"""
Delaunay triangulation and adjacency graph utilities.

This module provides functions for computing Delaunay triangulations
and building cell adjacency graphs from point coordinates.
"""

import numpy as np
from scipy.spatial import Delaunay
from typing import List, Set, Tuple


def compute_delaunay(x: np.ndarray, y: np.ndarray) -> Delaunay:
    """
    Compute Delaunay triangulation from x and y coordinates.

    :param x: X coordinates of points.
    :type x: np.ndarray
    :param y: Y coordinates of points.
    :type y: np.ndarray
    :return: Delaunay triangulation object.
    :rtype: Delaunay
    """
    points = np.column_stack([x, y])
    return Delaunay(points)


def delaunay_adjacency(x: np.ndarray, y: np.ndarray) -> List[Set[int]]:
    """
    Build adjacency graph from Delaunay triangulation.

    Each simplex (triangle) connects three points; this function
    identifies all neighboring point pairs.

    :param x: X coordinates of points.
    :type x: np.ndarray
    :param y: Y coordinates of points.
    :type y: np.ndarray
    :return: List of neighbor sets, one per point (index-aligned).
    :rtype: List[Set[int]]
    """
    points = np.column_stack([x, y])
    tri = Delaunay(points)
    
    n = len(points)
    neighbors = [set() for _ in range(n)]
    
    # Each simplex is a triangle with 3 vertices
    for simplex in tri.simplices:
        a, b, c = simplex
        neighbors[a].update([b, c])
        neighbors[b].update([a, c])
        neighbors[c].update([a, b])
    
    return neighbors


def filter_triangles_by_circumradius(points: np.ndarray,
                                     tri: Delaunay,
                                     scale: float = 15.0) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Filter Delaunay triangles by circumradius threshold.

    Keeps triangles with circumradius R <= scale * median_edge_length / sqrt(3).
    This helps create concave boundaries by removing large "bridging" triangles.

    :param points: (N, 2) array of point coordinates.
    :type points: np.ndarray
    :param tri: Delaunay triangulation object.
    :type tri: Delaunay
    :param scale: Scaling factor for threshold (larger = more triangles kept).
    :type scale: float
    :return: List of triangle vertices as (A, B, C) coordinate tuples.
    :rtype: List[Tuple[np.ndarray, np.ndarray, np.ndarray]]
    """
    T = points[tri.simplices]  # (M, 3, 2) - M triangles, 3 vertices each, 2D coords
    
    # Compute edge lengths
    e1 = np.linalg.norm(T[:, 0] - T[:, 1], axis=1)
    e2 = np.linalg.norm(T[:, 1] - T[:, 2], axis=1)
    e3 = np.linalg.norm(T[:, 2] - T[:, 0], axis=1)
    edge_lengths = np.concatenate([e1, e2, e3])
    
    # Median edge length as characteristic scale
    L = np.median(edge_lengths)
    R_thresh = (scale * L) / np.sqrt(3)
    
    # Filter triangles by circumradius
    kept_triangles = []
    for A, B, C in T:
        a = np.linalg.norm(B - C)
        b = np.linalg.norm(A - C)
        c = np.linalg.norm(A - B)
        s = 0.5 * (a + b + c)
        area_sq = s * (s - a) * (s - b) * (s - c)
        
        if area_sq <= 0:
            continue
        
        area = np.sqrt(area_sq)
        R = (a * b * c) / (4.0 * area)
        
        if R <= R_thresh:
            kept_triangles.append((A, B, C))
    
    return kept_triangles


def estimate_edge_length_scale(tri: Delaunay, points: np.ndarray) -> float:
    """
    Estimate characteristic edge length from triangulation.

    :param tri: Delaunay triangulation.
    :type tri: Delaunay
    :param points: (N, 2) point coordinates.
    :type points: np.ndarray
    :return: Median edge length.
    :rtype: float
    """
    T = points[tri.simplices]
    e1 = np.linalg.norm(T[:, 0] - T[:, 1], axis=1)
    e2 = np.linalg.norm(T[:, 1] - T[:, 2], axis=1)
    e3 = np.linalg.norm(T[:, 2] - T[:, 0], axis=1)
    edge_lengths = np.concatenate([e1, e2, e3])
    return np.median(edge_lengths)

