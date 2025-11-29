"""
Heatmap and spatial visualization utilities.

This module provides functions for creating heatmaps showing topology metrics
and cell distributions.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.spatial import ConvexHull
from typing import Optional, Tuple
import os


def plot_deviation_heatmap(
    deviation_grid: np.ndarray,
    x_grid_edges: np.ndarray,
    y_grid_edges: np.ndarray,
    title: str,
    output_path: str,
    colormap: str = 'viridis',
    levels: Optional[list] = None,
    hull_points: Optional[np.ndarray] = None,
    page_w: float = 6.3,
    page_h: float = 3.0,
    dpi: int = 300
) -> None:
    """
    Plot a deviation heatmap with optional convex hull boundary.

    :param deviation_grid: 2D array of deviation values (transposed for pcolormesh).
    :type deviation_grid: np.ndarray
    :param x_grid_edges: X edges of grid.
    :type x_grid_edges: np.ndarray
    :param y_grid_edges: Y edges of grid.
    :type y_grid_edges: np.ndarray
    :param title: Plot title.
    :type title: str
    :param output_path: Path to save the figure.
    :type output_path: str
    :param colormap: Matplotlib colormap name.
    :type colormap: str
    :param levels: Custom colormap levels for BoundaryNorm.
    :type levels: Optional[list]
    :param hull_points: Optional (N, 2) array for boundary contour.
    :type hull_points: Optional[np.ndarray]
    :param page_w: Figure width in inches.
    :type page_w: float
    :param page_h: Figure height in inches.
    :type page_h: float
    :param dpi: Figure DPI.
    :type dpi: int
    """
    X, Y = np.meshgrid(x_grid_edges, y_grid_edges)
    
    fig, ax = plt.subplots(figsize=(page_w, page_h))
    
    # Apply custom normalization if levels provided
    if levels is not None:
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256, clip=False)
        pcm = ax.pcolormesh(X, Y, deviation_grid, cmap=colormap, norm=norm, shading='auto')
        cbar = fig.colorbar(pcm, ax=ax, extend='max', orientation='horizontal', ticks=levels)
        tick_labels = [f'{levels[i]}-{levels[i+1]}' for i in range(len(levels)-1)] + [f'>{levels[-1]}']
        cbar.set_ticklabels(tick_labels)
    else:
        pcm = ax.pcolormesh(X, Y, deviation_grid, cmap=colormap, shading='auto')
        cbar = fig.colorbar(pcm, ax=ax, orientation='horizontal')
    
    # Draw hull boundary
    if hull_points is not None and len(hull_points) >= 3:
        try:
            hull = ConvexHull(hull_points)
            for simplex in hull.simplices:
                ax.plot(hull_points[simplex, 0], hull_points[simplex, 1], 'red', linewidth=3)
        except Exception:
            pass
    
    ax.set_title(title)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal')
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close(fig)


def plot_mean_sidedness_numbers(
    metrics: dict,
    x_grid_edges: np.ndarray,
    y_grid_edges: np.ndarray,
    num_grid_x: int,
    num_grid_y: int,
    page_w: float,
    page_h: float,
    save_path: str,
    hull_points: Optional[np.ndarray] = None,
    posterior: bool = True,
    empty_policy: str = "nan",
    text_fmt: str = "{:.2f}",
    na_label: str = "",
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    dpi: int = 300
) -> None:
    """
    Plot mean sidedness E[n] as numbers with colored background.

    :param metrics: Output dict from compute_topology_deviation_metrics.
    :type metrics: dict
    :param x_grid_edges: X edges of grid.
    :type x_grid_edges: np.ndarray
    :param y_grid_edges: Y edges of grid.
    :type y_grid_edges: np.ndarray
    :param num_grid_x: Number of grid divisions along X.
    :type num_grid_x: int
    :param num_grid_y: Number of grid divisions along Y.
    :type num_grid_y: int
    :param page_w: Figure width in inches.
    :type page_w: float
    :param page_h: Figure height in inches.
    :type page_h: float
    :param save_path: Output path.
    :type save_path: str
    :param hull_points: Optional boundary points.
    :type hull_points: Optional[np.ndarray]
    :param posterior: If True, use posterior-smoothed E[n].
    :type posterior: bool
    :param empty_policy: Policy for empty ROIs ("nan", "global", etc.).
    :type empty_policy: str
    :param text_fmt: Format string for values.
    :type text_fmt: str
    :param na_label: Label for NaN values.
    :type na_label: str
    :param cmap: Colormap name.
    :type cmap: str
    :param vmin: Min value for colormap.
    :type vmin: Optional[float]
    :param vmax: Max value for colormap.
    :type vmax: Optional[float]
    :param dpi: Figure DPI.
    :type dpi: int
    """
    from ..topology import mean_sidedness_per_roi
    
    mean_n = mean_sidedness_per_roi(metrics, posterior=posterior, empty_policy=empty_policy)
    mean_n_grid = mean_n.reshape((num_grid_x, num_grid_y))
    mean_n_grid_masked = np.ma.masked_invalid(mean_n_grid.T)
    
    fig, ax = plt.subplots(figsize=(page_w, page_h/3))
    
    X, Y = np.meshgrid(x_grid_edges, y_grid_edges)
    pcm = ax.pcolormesh(X, Y, mean_n_grid_masked, cmap=cmap,
                       shading='auto', vmin=vmin, vmax=vmax)
    
    cbar = fig.colorbar(pcm, ax=ax, label='Mean polygon sidedness E[n]', shrink=0.9)
    cbar.ax.tick_params(labelsize=8)
    
    # Add text labels
    dx = (x_grid_edges[1] - x_grid_edges[0]) / 2.0
    dy = (y_grid_edges[1] - y_grid_edges[0]) / 2.0
    
    for i in range(num_grid_x):
        for j in range(num_grid_y):
            val = mean_n_grid[i, j]
            cx = x_grid_edges[i] + dx
            cy = y_grid_edges[j] + dy
            
            if np.isfinite(val):
                # Determine text color based on background
                norm_val = (val - (vmin if vmin is not None else np.nanmin(mean_n_grid))) / \
                          ((vmax if vmax is not None else np.nanmax(mean_n_grid)) -
                           (vmin if vmin is not None else np.nanmin(mean_n_grid)))
                text_color = 'black' if norm_val > 0.5 else 'white'
                ax.text(cx, cy, text_fmt.format(val), ha='center', va='center',
                       fontsize=4.2, color=text_color, weight='bold')
            elif na_label:
                ax.text(cx, cy, na_label, ha='center', va='center',
                       fontsize=4.2, color='black')
    
    # Grid lines
    for x in x_grid_edges:
        ax.plot([x, x], [y_grid_edges[0], y_grid_edges[-1]], 'k--', lw=0.5, alpha=0.7)
    for y in y_grid_edges:
        ax.plot([x_grid_edges[0], x_grid_edges[-1]], [y, y], 'k--', lw=0.5, alpha=0.7)
    
    # Hull boundary
    if hull_points is not None and len(hull_points) >= 3:
        try:
            hull = ConvexHull(hull_points)
            for simplex in hull.simplices:
                ax.plot(hull_points[simplex, 0], hull_points[simplex, 1], 'k-', lw=2)
        except Exception:
            pass
    
    ax.set_title('Mean Polygon Sidedness per ROI')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal')
    
    os.makedirs(os.path.dirname(save_path) if os.path.dirname(save_path) else '.', exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=dpi)
    plt.close(fig)


def plot_metric_map(
    metric_values: np.ndarray,
    x_grid_edges: np.ndarray,
    y_grid_edges: np.ndarray,
    num_grid_x: int,
    num_grid_y: int,
    title: str,
    colorbar_label: str,
    output_path: str,
    page_w: float = 6.3,
    page_h: float = 3.0,
    cmap: str = 'viridis',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    dpi: int = 300
) -> None:
    """
    Generic metric map plotter (JSD, W1, charge, etc.).

    :param metric_values: 1D array of metric values per ROI.
    :type metric_values: np.ndarray
    :param x_grid_edges: X grid edges.
    :type x_grid_edges: np.ndarray
    :param y_grid_edges: Y grid edges.
    :type y_grid_edges: np.ndarray
    :param num_grid_x: Grid divisions along X.
    :type num_grid_x: int
    :param num_grid_y: Grid divisions along Y.
    :type num_grid_y: int
    :param title: Plot title.
    :type title: str
    :param colorbar_label: Colorbar label.
    :type colorbar_label: str
    :param output_path: Save path.
    :type output_path: str
    :param page_w: Figure width.
    :type page_w: float
    :param page_h: Figure height.
    :type page_h: float
    :param cmap: Colormap.
    :type cmap: str
    :param vmin: Min value for colormap.
    :type vmin: Optional[float]
    :param vmax: Max value for colormap.
    :type vmax: Optional[float]
    :param dpi: Figure DPI.
    :type dpi: int
    """
    grid = metric_values.reshape((num_grid_x, num_grid_y)).T
    grid_masked = np.ma.masked_invalid(grid)
    
    fig, ax = plt.subplots(figsize=(page_w, page_h))
    
    X, Y = np.meshgrid(x_grid_edges, y_grid_edges)
    pcm = ax.pcolormesh(X, Y, grid_masked, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax)
    
    cbar = fig.colorbar(pcm, ax=ax, orientation='horizontal', label=colorbar_label)
    cbar.ax.tick_params(labelsize=8)
    
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('X Coordinate', fontsize=8)
    ax.set_ylabel('Y Coordinate', fontsize=8)
    ax.tick_params(axis='both', which='both', labelsize=6)
    ax.set_aspect('equal')
    
    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close(fig)

