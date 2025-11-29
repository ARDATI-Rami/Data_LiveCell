# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

import numpy as np
import matplotlib.pyplot as plt
import pytest

from test_roi_eulerian import assign_to_grid_vectorized
from test_roi_lagrangian import simple_nearest_neighbor_track

"""Interactive visual tests for manual inspection of plots.

These tests are marked with @pytest.mark.interactive and are intended to
be run manually, not in automated CI. Each test shows a deterministic
plot and prompts the user to confirm that it looks correct.
Run them with: pytest -m interactive
"""


@pytest.mark.interactive
def test_visual_histogram_confirmation():
    # use deterministic synthetic data rather than randomness
    data = np.array([3] * 10 + [4] * 20 + [5] * 30 + [6] * 80 + [7] * 60)

    fig, ax = plt.subplots()
    ax.hist(data,
            bins=[2.5, 3.5, 4.5, 5.5, 6.5, 7.5],
            edgecolor="black")
    ax.set_xlabel("Polygon sides")
    ax.set_ylabel("Count")
    ax.set_title("Interactive check: polygon distribution histogram")

    plt.tight_layout()
    plt.show(block=True)

    answer = input("Does the histogram look correct and match the expected skew towards hexagons? [y/N]: ").strip().lower()
    assert answer == "y"


@pytest.mark.interactive
def test_interactive_eulerian_grid_overlay():
    # small synthetic set of cell coordinates in [0, 1] x [0, 1]
    x_coords = np.array([0.1, 0.3, 0.7, 0.9, 0.2, 0.8])
    y_coords = np.array([0.1, 0.8, 0.2, 0.7, 0.4, 0.6])

    # 4x4 grid
    x_edges = np.linspace(0.0, 1.0, 5)
    y_edges = np.linspace(0.0, 1.0, 5)

    gx, gy = assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges)

    fig, ax = plt.subplots()

    # draw grid lines
    for x in x_edges:
        ax.axvline(x, color="lightgray", linewidth=0.8)
    for y in y_edges:
        ax.axhline(y, color="lightgray", linewidth=0.8)

    # color points by their grid index
    cell_indices = gx + gy * (len(x_edges) - 1)
    sc = ax.scatter(x_coords, y_coords, c=cell_indices, cmap="viridis", s=60, edgecolor="black")
    plt.colorbar(sc, ax=ax, label="Grid cell index")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Interactive check: Eulerian grid assignment")

    plt.tight_layout()
    plt.show(block=True)

    answer = input(
        "Do the points appear in the expected grid cells with consistent colouring? [y/N]: "
    ).strip().lower()
    assert answer == "y"


@pytest.mark.interactive
def test_interactive_lagrangian_tracking_paths():
    # deterministic synthetic frames: frame1 is a small shift of frame0
    frame0 = np.array([
        [0.1, 0.1],
        [0.4, 0.2],
        [0.7, 0.3],
        [0.3, 0.7],
        [0.8, 0.8],
    ])
    frame1 = frame0 + np.array([0.1, 0.0])  # shift to the right

    tracks = simple_nearest_neighbor_track(frame0, frame1, max_distance=0.3)

    fig, ax = plt.subplots()
    ax.scatter(frame0[:, 0], frame0[:, 1], c="blue", label="Frame 0", s=40)
    ax.scatter(frame1[:, 0], frame1[:, 1], c="orange", label="Frame 1", s=40)

    # draw arrows/lines for tracked pairs
    for i0, j1 in tracks.items():
        x0, y0 = frame0[i0]
        x1, y1 = frame1[j1]
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="->", color="black", linewidth=1.0),
        )

    ax.set_xlim(0.0, 1.2)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal", adjustable="box")
    ax.legend()
    ax.set_title("Interactive check: Lagrangian nearest-neighbor trajectories")

    plt.tight_layout()
    plt.show(block=True)

    answer = input(
        "Do the arrows correctly connect each point in Frame 0 to its shifted partner in Frame 1? [y/N]: "
    ).strip().lower()
    assert answer == "y"
