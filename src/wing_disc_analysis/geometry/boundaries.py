# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""
Boundary detection and alpha-shape utilities.

This module provides functions for detecting tissue boundaries using
convex hull and alpha-shape (concave hull) algorithms.
"""

import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from shapely.geometry import Polygon, MultiLineString
from shapely.ops import polygonize, unary_union
from typing import Optional, List, Tuple


def compute_convex_hull(points: np.ndarray) -> Tuple[ConvexHull, Polygon]:
    """
    Compute convex hull from point coordinates.

    :param points: (N, 2) array of point coordinates.
    :type points: np.ndarray
    :return: Tuple of (ConvexHull object, Shapely Polygon).
    :rtype: Tuple[ConvexHull, Polygon]
    """
    hull = ConvexHull(points)
    hull_polygon = Polygon(points[hull.vertices])
    return hull, hull_polygon


def compute_alpha_shape(points: np.ndarray,
                       scale: float = 15.0) -> Optional[Polygon]:
    """
    Compute alpha-shape (concave hull) from points using Delaunay triangulation.

    Filters triangles by circumradius and polygonizes the edges to create
    a concave boundary that follows the tissue outline.

    :param points: (N, 2) array of point coordinates.
    :type points: np.ndarray
    :param scale: Scaling factor for circumradius threshold (default 15.0).
                  Larger values = smoother boundary; smaller = tighter fit.
    :type scale: float
    :return: Shapely Polygon representing the boundary, or None if failed.
    :rtype: Optional[Polygon]
    """
    try:
        # Delaunay triangulation
        tri = Delaunay(points)
        T = points[tri.simplices]  # (M, 3, 2)
        
        # Estimate length scale
        e1 = np.linalg.norm(T[:, 0] - T[:, 1], axis=1)
        e2 = np.linalg.norm(T[:, 1] - T[:, 2], axis=1)
        e3 = np.linalg.norm(T[:, 2] - T[:, 0], axis=1)
        edge_lengths = np.concatenate([e1, e2, e3])
        L = np.median(edge_lengths)
        
        # Threshold for circumradius
        R_thresh = (scale * L) / np.sqrt(3)
        
        # Collect edges from triangles with R <= threshold
        edges = []
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
                edges.append((tuple(A), tuple(B)))
                edges.append((tuple(B), tuple(C)))
                edges.append((tuple(C), tuple(A)))
        
        if not edges:
            # Fallback to convex hull
            hull = ConvexHull(points)
            return Polygon(points[hull.vertices]).buffer(0)
        
        # Polygonize edges
        mls = MultiLineString(edges)
        polys = list(polygonize(mls))
        
        if polys:
            merged = unary_union(polys)
            if merged.geom_type == "MultiPolygon":
                poly = max(merged.geoms, key=lambda g: g.area).buffer(0)
            else:
                poly = merged.buffer(0)
            return poly
        else:
            # Fallback to convex hull
            hull = ConvexHull(points)
            return Polygon(points[hull.vertices]).buffer(0)
    
    except Exception:
        # Fallback to convex hull on any error
        hull = ConvexHull(points)
        return Polygon(points[hull.vertices]).buffer(0)


def get_boundary_points(polygon: Polygon) -> np.ndarray:
    """
    Extract boundary points from a Shapely polygon.

    :param polygon: Shapely Polygon object.
    :type polygon: Polygon
    :return: (N, 2) array of boundary coordinates.
    :rtype: np.ndarray
    """
    coords = list(polygon.exterior.coords)
    return np.array(coords)


def compute_compactness(points: np.ndarray) -> float:
    """
    Compute compactness metric 4πA/P² for a set of points.

    Returns the compactness of the convex hull. Value of 1 is a perfect circle.

    :param points: (N, 2) array of point coordinates.
    :type points: np.ndarray
    :return: Compactness value in range (0, 1].
    :rtype: float
    """
    if len(points) < 3:
        return 0.0
    
    try:
        hull = ConvexHull(points)
        hull_vertices = points[hull.vertices]
        
        # Compute area and perimeter using shoelace formula
        x, y = hull_vertices[:, 0], hull_vertices[:, 1]
        area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
        perimeter = np.sum(np.hypot(np.diff(x, append=x[0]), np.diff(y, append=y[0])))
        
        if perimeter <= 1e-12:
            return 0.0
        
        return float(4.0 * np.pi * area / (perimeter * perimeter))
    except Exception:
        return 0.0

