# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""
Utility functions and helpers.

This module contains various utility functions that don't fit neatly into
other categories.
"""

import numpy as np
from typing import Tuple, Optional


# Standard page dimensions (from original scripts)
PAGE_WIDTH_PT = 455.24411  # LaTeX page width in points
PAGE_HEIGHT_PT = 702.78308  # LaTeX page height in points
POINTS_PER_INCH = 72.27

PAGE_WIDTH_IN = PAGE_WIDTH_PT / POINTS_PER_INCH  # ≈ 6.30 inches
PAGE_HEIGHT_IN = PAGE_HEIGHT_PT / POINTS_PER_INCH  # ≈ 9.73 inches


def get_page_dimensions() -> Tuple[float, float]:
    """
    Get standard page dimensions in inches.

    :return: Tuple of (width, height) in inches.
    :rtype: Tuple[float, float]
    """
    return PAGE_WIDTH_IN, PAGE_HEIGHT_IN


def compute_figure_size(width_fraction: float = 1.0,
                       height_fraction: float = 1.0) -> Tuple[float, float]:
    """
    Compute figure size as fraction of page dimensions.

    :param width_fraction: Fraction of page width (default 1.0).
    :type width_fraction: float
    :param height_fraction: Fraction of page height (default 1.0).
    :type height_fraction: float
    :return: Tuple of (width, height) in inches.
    :rtype: Tuple[float, float]
    """
    return PAGE_WIDTH_IN * width_fraction, PAGE_HEIGHT_IN * height_fraction


def ensure_dir_exists(path: str) -> None:
    """
    Ensure directory exists, create if necessary.

    :param path: Directory path.
    :type path: str
    """
    import os
    os.makedirs(path, exist_ok=True)


def generate_random_colors(n: int, seed: Optional[int] = None) -> list:
    """
    Generate n distinct random colors.

    :param n: Number of colors to generate.
    :type n: int
    :param seed: Random seed for reproducibility.
    :type seed: Optional[int]
    :return: List of RGBA tuples.
    :rtype: list
    """
    import matplotlib.pyplot as plt
    import random

    cmap = plt.get_cmap('gist_ncar', n)
    colors = [cmap(i / n) for i in range(n)]

    if seed is not None:
        random.seed(seed)
    random.shuffle(colors)

    return colors


def print_progress(current: int, total: int, prefix: str = 'Progress') -> None:
    """
    Print a simple progress indicator.

    :param current: Current iteration (0-based or 1-based).
    :type current: int
    :param total: Total number of iterations.
    :type total: int
    :param prefix: Prefix text for progress message.
    :type prefix: str
    """
    percentage = (current / total) * 100 if total > 0 else 0
    print(f"{prefix}: {current}/{total} ({percentage:.1f}%)")


