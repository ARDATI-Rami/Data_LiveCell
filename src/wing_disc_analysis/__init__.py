# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""
Wing Disc Analysis Package

A modular Python package for analyzing epithelial cell topology in Drosophila
wing disc tissue. Provides both Eulerian (fixed-grid) and Lagrangian (cell-tracking)
analysis frameworks.

Modules:
--------
- io: XLS/ODS file reading utilities
- geometry: Spatial operations, tessellation, and boundary detection
- topology: Topology metrics and statistical analysis
- roi: Region-of-interest construction (Eulerian and Lagrangian)
- visualization: Plotting and figure generation
- utils: Helper functions and utilities

Example Usage:
--------------
    import wing_disc_analysis as wda

    # Load data
    wb = wda.io.load_workbook('data/all_frames_wing_discs.xls')
    frames = wda.io.read_all_frames(wb)

    # Eulerian analysis
    grid = wda.roi.EulerianGrid(num_grid_x=20, num_grid_y=10)
    # ... (see scripts/ directory for complete examples)

    # Lagrangian analysis
    lag_roi = wda.roi.LagrangianROI(target_cells_per_roi=46)
    # ... (see scripts/ directory for complete examples)
"""

__version__ = '0.1.0'
__author__ = 'Wing Disc Analysis Team'

# Import key classes and functions for convenience
from . import io
from . import geometry
from . import topology
from . import roi
from . import visualization
from . import utils

__all__ = [
    'io',
    'geometry',
    'topology',
    'roi',
    'visualization',
    'utils',
]

