# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""I/O module for wing disc analysis package."""

from .xls_reader import (
    load_workbook,
    extract_data,
    extract_frame,
    read_all_frames,
    get_frame_count
)

from .ods_reader import (
    load_ods,
    extract_column_data,
    extract_neighbor_distribution,
    compute_neighbor_statistics
)

__all__ = [
    'load_workbook',
    'extract_data',
    'extract_frame',
    'read_all_frames',
    'get_frame_count',
    'load_ods',
    'extract_column_data',
    'extract_neighbor_distribution',
    'compute_neighbor_statistics',
]

