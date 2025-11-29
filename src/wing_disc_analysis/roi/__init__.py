# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""ROI module for region-of-interest construction and tracking."""

from .eulerian import (
    create_grid,
    compute_bounds_from_frames,
    assign_cells_to_grid,
    initialize_roi_data,
    collect_roi_data,
    generate_roi_labels,
    EulerianGrid
)

from .lagrangian import (
    build_lagrangian_rois_initial,
    compute_first_seen,
    LagrangianROI
)

__all__ = [
    # Eulerian
    'create_grid',
    'compute_bounds_from_frames',
    'assign_cells_to_grid',
    'initialize_roi_data',
    'collect_roi_data',
    'generate_roi_labels',
    'EulerianGrid',
    # Lagrangian
    'build_lagrangian_rois_initial',
    'compute_first_seen',
    'LagrangianROI',
]

