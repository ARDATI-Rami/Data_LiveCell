#!/usr/bin/env python3
"""
Lagrangian (cell-tracking) analysis script.

This script performs cell-tracking analysis where ROIs follow specific groups
of cells over time.
"""

import argparse
import sys
import os

# Add src to path if running from source
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import wing_disc_analysis as wda
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


def main():
    """Main entry point for Lagrangian analysis."""
    parser = argparse.ArgumentParser(
        description='Lagrangian ROI topology analysis of wing disc tissue'
    )
    parser.add_argument(
        '--data',
        type=str,
        default='Xls_Data/all_frames_wing_discs.xls',
        help='Path to XLS workbook (default: Xls_Data/all_frames_wing_discs.xls)'
    )
    parser.add_argument(
        '--target-cells',
        type=int,
        default=46,
        help='Target number of cells per ROI (default: 46)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='Results_wing_lagrangian',
        help='Output directory (default: Results_wing_lagrangian)'
    )
    parser.add_argument(
        '--alpha',
        type=float,
        default=0.5,
        help='Dirichlet prior parameter (default: 0.5)'
    )
    parser.add_argument(
        '--min-global-pct',
        type=float,
        default=0.2,
        help='Minimum global percentage to keep a polygon class (default: 0.2)'
    )
    parser.add_argument(
        '--random-seed',
        type=int,
        default=7,
        help='Random seed for reproducibility (default: 7)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Figure DPI (default: 300)'
    )
    parser.add_argument(
        '--map-style',
        choices=['scatter', 'polygon'],
        default='polygon',
        help='Lagrangian metric map style (scatter=centroids, polygon=Frame0 hulls)'
    )

    args = parser.parse_args()
    
    # Ensure output directory exists
    wda.utils.ensure_dir_exists(args.output)
    
    print(f"Loading data from: {args.data}")
    workbook = wda.io.load_workbook(args.data)
    
    print("Reading all frames...")
    frames = wda.io.read_all_frames(workbook, skip_boundary=False)
    print(f"Loaded {len(frames)} frames")
    
    # Initialize Lagrangian ROI manager
    print(f"Initializing Lagrangian ROI manager (target: {args.target_cells} cells/ROI)")
    lag_roi = wda.roi.LagrangianROI(
        target_cells_per_roi=args.target_cells,
        random_seed=args.random_seed
    )
    
    # Initialize from Frame 0
    print("Building initial ROIs from Frame 0...")
    ids0, x0, y0, p0 = frames[0]
    lag_roi.initialize(ids0, x0, y0)
    print(f"Created {lag_roi.num_rois} ROIs")
    
    # Collect Frame 0 data
    roi_indices_0 = lag_roi.track_frame(ids0, x0, y0, p0)
    lag_roi.collect(p0, roi_indices_0)
    # Build Frame 0 ROI membership index lists (non-boundary rows) for polygon maps
    roi_members_f0 = {r: [] for r in range(lag_roi.num_rois)}
    for i, cid in enumerate(ids0):
        if cid != -1:
            r = roi_indices_0[i]
            if r >= 0:
                roi_members_f0[int(r)].append(int(i))

    # Track through remaining frames
    print("Tracking cells through frames...")
    all_polygon_numbers = [p0]
    
    for t, (ids, x, y, polygons) in enumerate(frames[1:], start=1):
        roi_indices = lag_roi.track_frame(ids, x, y, polygons)
        lag_roi.collect(polygons, roi_indices)
        all_polygon_numbers.append(polygons)
        
        if (t + 1) % 50 == 0:
            wda.utils.print_progress(t + 1, len(frames), "Frames tracked")
    
    print(f"Tracked {len(frames)} frames")
    
    # Compute global distribution
    print("Computing global distribution...")
    filtered_polygons, mean_pct, std_pct = wda.topology.compute_global_distribution(
        all_polygon_numbers, threshold=0.5
    )
    
    print(f"Polygon classes found: {filtered_polygons}")
    
    # Plot global distribution
    output_path = os.path.join(args.output, 'global_distribution.png')
    print(f"Plotting global distribution: {output_path}")
    wda.visualization.plot_global_distribution(
        filtered_polygons, mean_pct, std_pct, output_path
    )
    
    # Compute topology metrics
    print("Computing topology deviation metrics...")
    
    # Convert roi_data to expected format
    roi_data_list = [lag_roi.roi_data[r] for r in range(lag_roi.num_rois)]
    
    metrics = wda.topology.compute_topology_deviation_metrics(
        roi_data_list=roi_data_list,
        all_polygon_numbers=all_polygon_numbers,
        polygons=None,
        alpha=args.alpha,
        min_global_pct=args.min_global_pct
    )
    print(f"Metrics computed for {len(metrics['support'])} polygon classes: {metrics['support']}")

    # --- Polygon-style map helpers (Frame 0 convex hull style) ---
    # Build ROI hulls from Frame 0
    roi_hulls = {}
    for r, members in roi_members_f0.items():
        pts = np.column_stack([x0[members], y0[members]])
        if pts.shape[0] < 3:
            continue
        try:
            h = ConvexHull(pts)
            roi_hulls[r] = pts[h.vertices]
        except Exception:
            roi_hulls[r] = pts

    # Global tissue hull (for outline)
    try:
        tissue_hull = ConvexHull(np.column_stack([x0, y0]))
        tissue_outline = np.column_stack([x0, y0])[tissue_hull.vertices]
    except Exception:
        tissue_outline = None

    def _plot_mean_sidedness():
        support = metrics['support'].astype(float)
        p = metrics['p_roi']
        mean_n = (p * support[None,:]).sum(axis=1)
        wda.visualization.plot_polygon_metric_map(
            mean_n,
            roi_hulls,
            'Mean polygon sidedness — Lagrangian ROIs',
            'E[n]',
            os.path.join(args.output, 'lag_map_mean_sidedness.png'),
            tissue_outline=tissue_outline,
            dpi=args.dpi
        )
    # --- Map generation branch ---
    if args.map_style == 'scatter':
        print("Using legacy scatter style maps.")
        num_rois = lag_roi.num_rois
        centroids_x = np.full(num_rois, np.nan); centroids_y = np.full(num_rois, np.nan)
        for r in range(num_rois):
            inds = np.where(roi_indices_0 == r)[0]
            if inds.size:
                centroids_x[r] = np.mean(x0[inds]); centroids_y[r] = np.mean(y0[inds])
        def _scatter(values, title, cmap, fname):
            v = np.asarray(values, dtype=float)
            valid = np.isfinite(v) & np.isfinite(centroids_x) & np.isfinite(centroids_y)
            if not np.any(valid):
                print(f"Skipping {fname}: no valid data"); return
            fig, ax = plt.subplots(figsize=(6,6))
            sc = ax.scatter(centroids_x[valid], centroids_y[valid], c=v[valid],
                            s=80, cmap=cmap, edgecolor='k', linewidth=0.3)
            ax.set_aspect('equal'); ax.set_title(title); ax.set_xlabel('X'); ax.set_ylabel('Y')
            cbar = fig.colorbar(sc, ax=ax); cbar.set_label(title, fontsize=8)
            plt.tight_layout()
            outpath = os.path.join(args.output, fname)
            plt.savefig(outpath, dpi=args.dpi); plt.close(fig); print(f"Saved {fname}")
        _scatter(metrics['jsd'],'JSD vs global','viridis','lag_map_topology_jsd.png')
        _scatter(metrics['w1'],'W1 distance','plasma','lag_map_topology_w1.png')
        _scatter(metrics['mean_abs_q'],'Mean |n-6|','magma','lag_map_mean_abs_charge.png')
        _plot_mean_sidedness()
    else:
        print("Generating polygon (Frame 0 hull) style Lagrangian maps.")
        wda.visualization.plot_polygon_metric_map(
            metrics['jsd'],
            roi_hulls,
            'Topology deviation (JSD vs global) — Lagrangian ROIs',
            'JSD (bits)',
            os.path.join(args.output, 'lag_map_topology_jsd.png'),
            tissue_outline=tissue_outline,
            dpi=args.dpi
        )
        wda.visualization.plot_polygon_metric_map(
            metrics['w1'],
            roi_hulls,
            'Topology deviation (W1) — Lagrangian ROIs',
            'W1 distance (sides)',
            os.path.join(args.output, 'lag_map_topology_w1.png'),
            tissue_outline=tissue_outline,
            cmap='plasma',
            dpi=args.dpi
        )
        wda.visualization.plot_polygon_metric_map(
            metrics['mean_abs_q'],
            roi_hulls,
            'Mean absolute topological charge — Lagrangian ROIs',
            'E[|n−6|]',
            os.path.join(args.output, 'lag_map_mean_abs_charge.png'),
            tissue_outline=tissue_outline,
            cmap='magma',
            dpi=args.dpi
        )
        # Hexagon z-residual if hexagons present
        if 6 in metrics['support']:
            k6 = int(np.where(metrics['support'] == 6)[0][0])
            z6 = metrics['z'][:, k6]
            lim = 4.0
            wda.visualization.plot_polygon_metric_map(
                z6,
                roi_hulls,
                'Hexagon residual (z) — Lagrangian ROIs',
                'z residual (hexagon)',
                os.path.join(args.output, 'lag_map_z_hexagon.png'),
                tissue_outline=tissue_outline,
                cmap='coolwarm',
                vmin=-lim,
                vmax=lim,
                dpi=args.dpi
            )
        _plot_mean_sidedness()

    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Global composition: {dict(zip(metrics['support'].tolist(), np.round(metrics['p_global'], 3)))}")
    print(f"Global E[|n-6|]: {np.sum(metrics['p_global'] * np.abs(metrics['support'] - 6)):.3f}")
    
    print(f"\nJSD (bits) - min/median/max: {np.nanmin(metrics['jsd']):.3f} / "
          f"{np.nanmedian(metrics['jsd']):.3f} / {np.nanmax(metrics['jsd']):.3f}")
    print(f"W1 (sides) - min/median/max: {np.nanmin(metrics['w1']):.3f} / "
          f"{np.nanmedian(metrics['w1']):.3f} / {np.nanmax(metrics['w1']):.3f}")
    print(f"E[|n-6|] - min/median/max: {np.nanmin(metrics['mean_abs_q']):.3f} / "
          f"{np.nanmedian(metrics['mean_abs_q']):.3f} / {np.nanmax(metrics['mean_abs_q']):.3f}")
    
    print(f"\nAnalysis complete! Results saved to: {args.output}")


if __name__ == '__main__':
    main()
