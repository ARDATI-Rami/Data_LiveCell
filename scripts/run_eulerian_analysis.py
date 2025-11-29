#!/usr/bin/env python3
"""
Eulerian (fixed-grid) analysis script.

This script performs spatial analysis of wing disc topology using a fixed
spatial grid to define regions of interest (ROIs).
"""

import argparse
import sys
import os

# Add src to path if running from source
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import wing_disc_analysis as wda


def main():
    """Main entry point for Eulerian analysis."""
    parser = argparse.ArgumentParser(
        description='Eulerian ROI topology analysis of wing disc tissue'
    )
    parser.add_argument(
        '--data',
        type=str,
        default='Xls_Data/all_frames_wing_discs.xls',
        help='Path to XLS workbook (default: Xls_Data/all_frames_wing_discs.xls)'
    )
    parser.add_argument(
        '--grid-x',
        type=int,
        default=20,
        help='Grid divisions along X axis (default: 20)'
    )
    parser.add_argument(
        '--grid-y',
        type=int,
        default=10,
        help='Grid divisions along Y axis (default: 10)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='Results_wing_eulerian',
        help='Output directory (default: Results_wing_eulerian)'
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
        '--dpi',
        type=int,
        default=300,
        help='Figure DPI (default: 300)'
    )
    parser.add_argument(
        '--skip-frames',
        action='store_true',
        help='Skip per-frame plots (faster)'
    )
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    wda.utils.ensure_dir_exists(args.output)
    
    print(f"Loading data from: {args.data}")
    workbook = wda.io.load_workbook(args.data)
    num_frames = wda.io.get_frame_count(workbook)
    print(f"Found {num_frames} frames")
    
    # Initialize Eulerian grid
    print(f"Creating Eulerian grid ({args.grid_x} x {args.grid_y} = {args.grid_x * args.grid_y} ROIs)")
    grid = wda.roi.EulerianGrid(num_grid_x=args.grid_x, num_grid_y=args.grid_y)
    
    # Collect data from all frames
    all_x_coords = []
    all_y_coords = []
    all_polygon_numbers = []
    
    print("Reading all frames...")
    for frame_idx in range(num_frames):
        x, y, polygons = wda.io.extract_data(workbook, frame_idx, skip_boundary=True)
        all_x_coords.append(x)
        all_y_coords.append(y)
        all_polygon_numbers.append(polygons)
        
        if (frame_idx + 1) % 50 == 0:
            wda.utils.print_progress(frame_idx + 1, num_frames, "Frames loaded")
    
    print(f"Loaded {num_frames} frames")
    
    # Fit grid to data
    print("Fitting grid to data bounds...")
    grid.fit(all_x_coords, all_y_coords)
    
    # Process each frame
    print("Processing frames and assigning to ROIs...")
    for frame_idx in range(num_frames):
        x = all_x_coords[frame_idx]
        y = all_y_coords[frame_idx]
        polygons = all_polygon_numbers[frame_idx]
        
        # Assign cells to grid
        roi_indices = grid.assign(x, y)
        
        # Collect data
        grid.collect(polygons, roi_indices)
        
        if (frame_idx + 1) % 50 == 0:
            wda.utils.print_progress(frame_idx + 1, num_frames, "Frames processed")
    
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
    metrics = wda.topology.compute_topology_deviation_metrics(
        roi_data_list=grid.roi_data_list,
        all_polygon_numbers=all_polygon_numbers,
        polygons=None,
        alpha=args.alpha,
        min_global_pct=args.min_global_pct
    )
    
    print(f"Metrics computed for {len(metrics['support'])} polygon classes: {metrics['support']}")
    
    # Plot metric maps
    page_w, page_h = wda.utils.get_page_dimensions()
    
    print("Generating metric maps...")
    
    # JSD map
    wda.visualization.plot_metric_map(
        metrics['jsd'],
        grid.x_grid_edges,
        grid.y_grid_edges,
        grid.num_grid_x,
        grid.num_grid_y,
        'Topology deviation map (JSD vs global)',
        'Jensen–Shannon divergence (bits)',
        os.path.join(args.output, 'map_topology_jsd.png'),
        page_w=page_w/2,
        page_h=page_h/6,
        dpi=args.dpi
    )
    
    # W1 map
    wda.visualization.plot_metric_map(
        metrics['w1'],
        grid.x_grid_edges,
        grid.y_grid_edges,
        grid.num_grid_x,
        grid.num_grid_y,
        'Topology deviation map (Wasserstein-1)',
        'Wasserstein-1 distance (sides)',
        os.path.join(args.output, 'map_topology_w1.png'),
        page_w=page_w/2,
        page_h=page_h/6,
        dpi=args.dpi
    )
    
    # Mean absolute charge map
    wda.visualization.plot_metric_map(
        metrics['mean_abs_q'],
        grid.x_grid_edges,
        grid.y_grid_edges,
        grid.num_grid_x,
        grid.num_grid_y,
        'Mean absolute topological charge per ROI',
        'E[|n−6|]',
        os.path.join(args.output, 'map_mean_abs_charge.png'),
        page_w=page_w/2,
        page_h=page_h/6,
        dpi=args.dpi
    )
    
    # Mean sidedness numbers map
    points0 = np.column_stack((all_x_coords[0], all_y_coords[0]))
    wda.visualization.plot_mean_sidedness_numbers(
        metrics=metrics,
        x_grid_edges=grid.x_grid_edges,
        y_grid_edges=grid.y_grid_edges,
        num_grid_x=grid.num_grid_x,
        num_grid_y=grid.num_grid_y,
        page_w=page_w,
        page_h=page_h,
        save_path=os.path.join(args.output, 'map_mean_sidedness_numbers.png'),
        hull_points=points0,
        posterior=True,
        empty_policy="nan",
        dpi=args.dpi
    )
    
    print(f"\nAnalysis complete! Results saved to: {args.output}")
    print(f"\nGenerated files:")
    print(f"  - global_distribution.png")
    print(f"  - map_topology_jsd.png")
    print(f"  - map_topology_w1.png")
    print(f"  - map_mean_abs_charge.png")
    print(f"  - map_mean_sidedness_numbers.png")


if __name__ == '__main__':
    main()

