# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""Geometry module for spatial operations and boundary detection."""

from .tessellation import (
    compute_delaunay,
    delaunay_adjacency,
    filter_triangles_by_circumradius,
    estimate_edge_length_scale
)

from .boundaries import (
    compute_convex_hull,
    compute_alpha_shape,
    get_boundary_points,
    compute_compactness
)

from .spatial import (
    compute_centroids,
    build_kdtree,
    nearest_neighbors,
    farthest_point_sampling,
    radial_distance,
    radial_cv,
    robust_outlier_threshold
)

__all__ = [
    # Tessellation
    'compute_delaunay',
    'delaunay_adjacency',
    'filter_triangles_by_circumradius',
    'estimate_edge_length_scale',
    # Boundaries
    'compute_convex_hull',
    'compute_alpha_shape',
    'get_boundary_points',
    'compute_compactness',
    # Spatial
    'compute_centroids',
    'build_kdtree',
    'nearest_neighbors',
    'farthest_point_sampling',
    'radial_distance',
    'radial_cv',
    'robust_outlier_threshold',
]

