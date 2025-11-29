import numpy as np

"""Test-side helper and unit tests for simple Lagrangian tracking.

simple_nearest_neighbor_track implements a tiny nearest-neighbor matcher
between two frames of 2D points. It is intentionally simple and only
used from tests to exercise Lagrangian-style logic without touching
analysis scripts.
"""


def simple_nearest_neighbor_track(frame0, frame1, max_distance=1.0):
    tracks = {}
    for i0, p0 in enumerate(frame0):
        if frame1.size == 0:
            continue
        dists = np.linalg.norm(frame1 - p0, axis=1)
        j = np.argmin(dists)
        if dists[j] <= max_distance:
            tracks[i0] = j
    return tracks


def test_track_cells_preserves_ids_for_small_displacements():
    frame0 = np.array([[0.0, 0.0], [1.0, 0.0]])
    frame1 = np.array([[0.1, 0.0], [0.9, 0.0]])

    tracks = simple_nearest_neighbor_track(frame0, frame1, max_distance=0.5)

    assert tracks[0] == 0
    assert tracks[1] == 1


def test_track_cells_handles_empty_initial_frame():
    frame0 = np.empty((0, 2))
    frame1 = np.array([[0.0, 0.0]])

    tracks = simple_nearest_neighbor_track(frame0, frame1, max_distance=1.0)

    assert tracks == {}


def test_track_cells_ignores_large_displacements():
    frame0 = np.array([[0.0, 0.0], [1.0, 1.0]])
    frame1 = np.array([[10.0, 10.0], [20.0, 20.0]])

    tracks = simple_nearest_neighbor_track(frame0, frame1, max_distance=0.5)

    assert tracks == {}


def test_track_cells_many_to_one_mapping_allowed():
    frame0 = np.array([[0.0, 0.0], [0.0, 0.1]])
    frame1 = np.array([[0.0, 0.05]])

    tracks = simple_nearest_neighbor_track(frame0, frame1, max_distance=0.2)

    # both initial points are closest to the single target point
    assert tracks == {0: 0, 1: 0}


def test_track_cells_tie_breaker_is_first_min_index():
    frame0 = np.array([[0.0, 0.0]])
    # two points at equal distance from the origin
    frame1 = np.array([[1.0, 0.0], [0.0, 1.0]])

    tracks = simple_nearest_neighbor_track(frame0, frame1, max_distance=2.0)

    # np.argmin will pick the first index in case of a tie
    assert tracks[0] == 0


def test_track_cells_handles_empty_second_frame():
    frame0 = np.array([[0.0, 0.0], [1.0, 1.0]])
    frame1 = np.empty((0, 2))

    tracks = simple_nearest_neighbor_track(frame0, frame1, max_distance=1.0)

    assert tracks == {}
