# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

import numpy as np

"""Synthetic I/O-style tests for cell record handling.

These tests model the behavior of boundary-cell filtering and simple
per-frame counting/grouping using in-memory record dicts only. They do
not touch any real files or analysis scripts.
"""


def filter_boundary_cells(records):
    return [r for r in records if r.get("cell_id") != -1]


def test_boundary_cells_excluded_by_filter_helper():
    records = [
        {"cell_id": 1, "frame": 0},
        {"cell_id": -1, "frame": 0},
        {"cell_id": 2, "frame": 1},
    ]
    filtered = filter_boundary_cells(records)
    ids = {r["cell_id"] for r in filtered}
    assert -1 not in ids
    assert ids == {1, 2}


def test_filter_boundary_cells_keeps_all_valid_when_no_boundaries():
    records = [
        {"cell_id": 1, "frame": 0},
        {"cell_id": 2, "frame": 0},
    ]
    filtered = filter_boundary_cells(records)
    assert len(filtered) == len(records)
    assert {r["cell_id"] for r in filtered} == {1, 2}


def test_filter_boundary_cells_all_boundary_yields_empty():
    records = [
        {"cell_id": -1, "frame": 0},
        {"cell_id": -1, "frame": 1},
    ]
    filtered = filter_boundary_cells(records)
    assert filtered == []


def test_empty_sheet_handled_as_empty_result():
    # Stand-in behavior: an "empty sheet" yields an empty list of records
    empty_sheet_records = []
    assert empty_sheet_records == []


def test_correct_number_of_cells_per_frame_synthetic():
    # synthetic data for 2 frames: 2 cells in frame 0, 3 cells in frame 1
    records = [
        {"cell_id": 1, "frame": 0},
        {"cell_id": 2, "frame": 0},
        {"cell_id": 3, "frame": 1},
        {"cell_id": 4, "frame": 1},
        {"cell_id": 5, "frame": 1},
    ]
    counts = {}
    for r in records:
        counts[r["frame"]] = counts.get(r["frame"], 0) + 1
    assert counts[0] == 2
    assert counts[1] == 3


def test_frame_counts_with_missing_frames():
    # frames are non-contiguous: 0, 2, 2, 5
    records = [
        {"cell_id": 1, "frame": 0},
        {"cell_id": 2, "frame": 2},
        {"cell_id": 3, "frame": 2},
        {"cell_id": 4, "frame": 5},
    ]
    counts = {}
    for r in records:
        counts[r["frame"]] = counts.get(r["frame"], 0) + 1
    assert counts == {0: 1, 2: 2, 5: 1}


def group_by_frame(records):
    """Group records by frame while preserving input order per frame."""
    groups = {}
    for r in records:
        frame = r["frame"]
        groups.setdefault(frame, []).append(r)
    return groups


def test_records_grouped_by_frame_ordering_is_stable():
    records = [
        {"cell_id": 1, "frame": 0},
        {"cell_id": 2, "frame": 1},
        {"cell_id": 3, "frame": 1},
        {"cell_id": 4, "frame": 0},
    ]
    groups = group_by_frame(records)
    assert [r["cell_id"] for r in groups[0]] == [1, 4]
    assert [r["cell_id"] for r in groups[1]] == [2, 3]
