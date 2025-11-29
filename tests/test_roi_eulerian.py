# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

import numpy as np

"""Test-side helper and unit tests for Eulerian ROI grid assignment.

assign_to_grid_vectorized mimics a simple ROI binning step: given cell
coordinates and 1D x/y bin edges, it returns integer bin indices,
clipped to the valid index range. This is purely a test helper and does
not modify or depend on analysis scripts.
"""


def assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges):
    x_coords = np.asarray(x_coords)
    y_coords = np.asarray(y_coords)
    x_idx = np.digitize(x_coords, x_edges) - 1
    y_idx = np.digitize(y_coords, y_edges) - 1
    # clip to valid range
    x_idx = np.clip(x_idx, 0, len(x_edges) - 2)
    y_idx = np.clip(y_idx, 0, len(y_edges) - 2)
    return x_idx, y_idx


def test_assign_to_grid_vectorized_all_cells_assigned():
    x_edges = np.array([0.0, 0.5, 1.0])
    y_edges = np.array([0.0, 0.5, 1.0])

    x_coords = np.array([0.25, 0.75, 0.25, 0.75])
    y_coords = np.array([0.25, 0.25, 0.75, 0.75])

    gx, gy = assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges)

    assert gx.shape == x_coords.shape
    assert gy.shape == y_coords.shape

    assert np.all((gx >= 0) & (gx <= 1))
    assert np.all((gy >= 0) & (gy <= 1))

    combined = set(zip(gx.tolist(), gy.tolist()))
    assert combined == {(0, 0), (1, 0), (0, 1), (1, 1)}


def test_assign_to_grid_vectorized_empty_input():
    x_edges = np.array([0.0, 0.5, 1.0])
    y_edges = np.array([0.0, 0.5, 1.0])

    x_coords = np.array([])
    y_coords = np.array([])

    gx, gy = assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges)

    assert gx.size == 0
    assert gy.size == 0


def test_assign_to_grid_vectorized_points_on_edges():
    x_edges = np.array([0.0, 0.5, 1.0])
    y_edges = np.array([0.0, 0.5, 1.0])

    # points exactly on edges
    x_coords = np.array([0.0, 0.5, 1.0])
    y_coords = np.array([0.0, 0.5, 1.0])

    gx, gy = assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges)

    # np.digitize places values equal to a bin edge in the right bin by default,
    # and our -1 shift plus clipping keeps indices in [0, len(edges)-2]
    assert np.all(gx >= 0)
    assert np.all(gx <= 1)
    assert np.all(gy >= 0)
    assert np.all(gy <= 1)


def test_assign_to_grid_vectorized_points_outside_range_are_clipped():
    x_edges = np.array([0.0, 0.5, 1.0])
    y_edges = np.array([0.0, 0.5, 1.0])

    x_coords = np.array([-0.1, 1.1])
    y_coords = np.array([-0.2, 2.0])

    gx, gy = assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges)

    # left of first edge -> 0, right of last edge -> last valid index (1)
    assert gx.tolist() == [0, 1]
    assert gy.tolist() == [0, 1]


def test_assign_to_grid_vectorized_nonuniform_bins():
    x_edges = np.array([0.0, 0.1, 0.4, 1.0])
    y_edges = np.array([0.0, 0.2, 0.8, 1.0])

    x_coords = np.array([0.05, 0.2, 0.8])
    y_coords = np.array([0.1, 0.3, 0.9])

    gx, gy = assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges)

    # manually check which bins they should fall into
    assert gx.tolist() == [0, 1, 2]
    assert gy.tolist() == [0, 1, 2]
