# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025  Rami Ardati

"""
Histogram and distribution plotting utilities.

This module provides functions for plotting cell shape distribution histograms.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional
import os


def plot_global_distribution(
    filtered_polygons: np.ndarray,
    mean_pct: np.ndarray,
    std_pct: np.ndarray,
    output_path: str,
    title: str = 'Cell shape distribution across all frames',
    page_w: float = 6.3,
    page_h: float = 3.0,
    dpi: int = 300
) -> None:
    """
    Plot global polygon distribution with error bars.

    :param filtered_polygons: Polygon categories kept.
    :type filtered_polygons: np.ndarray
    :param mean_pct: Mean percentages.
    :type mean_pct: np.ndarray
    :param std_pct: Standard deviations.
    :type std_pct: np.ndarray
    :param output_path: Output filepath (PNG).
    :type output_path: str
    :param title: Plot title.
    :type title: str
    :param page_w: Figure width in inches.
    :type page_w: float
    :param page_h: Figure height in inches.
    :type page_h: float
    :param dpi: Figure DPI.
    :type dpi: int
    """
    from ..topology import get_shape_name
    
    fig, ax = plt.subplots(figsize=(page_w, page_h))
    
    bars = ax.bar(filtered_polygons, mean_pct, yerr=std_pct,
                 color=plt.cm.tab20.colors, capsize=2)
    
    ax.set_title(title, weight="bold")
    ax.set_xlabel('Shape')
    ax.set_ylabel('Percentage of Cells')
    ax.set_ylim(0, max(mean_pct + std_pct) * 1.2)
    ax.set_xticks(filtered_polygons)
    ax.set_xticklabels([get_shape_name(i) for i in filtered_polygons])
    ax.grid(True)
    
    # Annotate bars with mean ± std
    for bar, m, s in zip(bars, mean_pct, std_pct):
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, yval + s + 0.5,
               f'{m:.2f}±{s:.2f}%', ha='center', va='bottom', fontsize=8)
    
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_frame_histogram(
    polygon_numbers: np.ndarray,
    frame_index: int,
    output_path: str,
    base_name: str = 'frame',
    page_w: float = 10.0,
    page_h: float = 6.0,
    dpi: int = 300
) -> None:
    """
    Plot histogram for a single frame.

    :param polygon_numbers: Polygon numbers for cells in this frame.
    :type polygon_numbers: np.ndarray
    :param frame_index: Frame number (for labeling).
    :type frame_index: int
    :param output_path: Output path.
    :type output_path: str
    :param base_name: Base name for file.
    :type base_name: str
    :param page_w: Figure width.
    :type page_w: float
    :param page_h: Figure height.
    :type page_h: float
    :param dpi: Figure DPI.
    :type dpi: int
    """
    from ..topology import get_shape_name
    
    # Calculate percentages
    unique, counts = np.unique(polygon_numbers, return_counts=True)
    percentages = counts / counts.sum() * 100
    
    fig, ax = plt.subplots(figsize=(page_w, page_h))
    
    bars = ax.bar(unique, percentages, color=plt.cm.tab20.colors)
    
    ax.set_title(f'Percentage of Cell Shapes (Polygon Count) - Frame {frame_index}')
    ax.set_xlabel('Shape')
    ax.set_ylabel('Percentage of Cells')
    ax.set_xticks(unique)
    ax.set_xticklabels([get_shape_name(i) for i in unique])
    ax.grid(True)
    
    # Show percentages on top of bars
    for bar, percentage in zip(bars, percentages):
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, yval + 0.5,
               f'{percentage:.2f}%', ha='center', va='bottom')
    
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_neighbor_distribution(
    neighbors: list,
    output_path: Optional[str] = None,
    max_neighbors: int = 14,
    page_w: float = 10.0,
    page_h: float = 6.0
) -> None:
    """
    Plot neighbor count distribution histogram.

    :param neighbors: List of neighbor counts.
    :type neighbors: list
    :param output_path: Optional output path (displays if None).
    :type output_path: Optional[str]
    :param max_neighbors: Maximum number of neighbors to show.
    :type max_neighbors: int
    :param page_w: Figure width.
    :type page_w: float
    :param page_h: Figure height.
    :type page_h: float
    """
    fig, ax = plt.subplots(figsize=(page_w, page_h))
    
    ax.hist(neighbors, bins=max_neighbors, edgecolor='black')
    ax.set_title('Distribution of Neighbors')
    ax.set_xlabel('Number of Neighbors')
    ax.set_ylabel('Frequency')
    ax.grid(True)
    
    if output_path:
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        plt.savefig(output_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()

