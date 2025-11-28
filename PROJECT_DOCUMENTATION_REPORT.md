# Project Documentation Report: Live Cell Epithelial Topology Analysis

**Project:** Data_LiveCell - Wing Disc Epithelial Tissue Analysis  
**Date:** November 28, 2025  
**Analysis Framework:** Eulerian (Fixed-Grid) and Lagrangian (Cell-Tracking) Methods

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Project Overview](#project-overview)
3. [Scientific Workflow](#scientific-workflow)
4. [File Inventory and Relationships](#file-inventory-and-relationships)
5. [Core Analysis Scripts](#core-analysis-scripts)
6. [Utility Scripts](#utility-scripts)
7. [Data Flow Architecture](#data-flow-architecture)
8. [Dependencies and Environment](#dependencies-and-environment)
9. [External Tools and Libraries](#external-tools-and-libraries)
10. [Recommendations for Improvement](#recommendations-for-improvement)

---

## Executive Summary

This project performs quantitative analysis of epithelial cell topology in Drosophila wing disc tissue using time-lapse microscopy data. The analysis employs two complementary approaches:

- **Eulerian Analysis** (`Euler_wing_disc_cell_analysis.py`): Fixed spatial grid regions (ROIs) track how cell shapes evolve at specific tissue locations
- **Lagrangian Analysis** (`Lagrange_wing_disc_analysis.py`): Cell-centered ROIs that follow specific groups of cells over time

The project processes cell segmentation data (from EpiTools), computes topology metrics (polygon distributions, topological charge, Jensen-Shannon divergence), and generates spatial heatmaps showing tissue heterogeneity.

**Current Status:** Functional analysis scripts with outputs, but lacking project organization, documentation, and code reuse patterns.

---

## Project Overview

### Research Question
How does epithelial cell topology (3-gon, 4-gon, 5-gon, 6-gon, 7-gon distributions) vary spatially and temporally during tissue development, and how do cell division, death, and neighbor exchange events affect these patterns?

### Data Sources
- **Input:** XLS workbooks containing per-frame cell data (cell ID, X/Y coordinates, polygon count)
- **Source Tool:** EpiTools (MATLAB/Icy segmentation pipeline)
- **Format:** Each sheet = one time frame; rows = cells; columns = [ID, X, Y, Polygon_Number]

### Key Outputs
- Spatial topology deviation maps (PNG heatmaps)
- Cell shape distribution histograms with statistical metrics
- ROI-based analysis comparing local to global distributions
- Animations showing topology evolution over time
- Quantitative metrics: JSD, Wasserstein-1 distance, topological charge

---

## Scientific Workflow

```
┌─────────────────┐
│ Raw Microscopy  │
│  TIFF Images    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   EpiTools      │ (External: MATLAB + Icy)
│  Segmentation   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  XLS Workbooks  │ (Xls_Data/)
│ [ID, X, Y, N]   │ N = polygon sides
└────────┬────────┘
         │
         ├──────────────────┬──────────────────┐
         ▼                  ▼                  ▼
┌──────────────┐   ┌──────────────┐   ┌──────────────┐
│   Eulerian   │   │  Lagrangian  │   │   Utility    │
│   Analysis   │   │   Analysis   │   │   Scripts    │
└──────┬───────┘   └──────┬───────┘   └──────┬───────┘
       │                  │                  │
       ▼                  ▼                  ▼
┌─────────────────────────────────────────────────────┐
│              Visualization Outputs                  │
│  - Heatmaps (JSD, W1, charge, sidedness)           │
│  - Histograms (shape distributions)                 │
│  - Animations (topology evolution)                  │
│  - Statistical reports                              │
└─────────────────────────────────────────────────────┘
```

---

## File Inventory and Relationships

### Core Python Scripts (Root Level)

| File | Lines | Purpose | Outputs | Dependencies |
|------|-------|---------|---------|--------------|
| `Euler_wing_disc_cell_analysis.py` | ~900 | Fixed-grid ROI analysis | `Results_wing_eulerian/` | xlrd, numpy, matplotlib, scipy, shapely |
| `Lagrange_wing_disc_analysis.py` | ~2500 | Cell-tracking ROI analysis | `Results_wing_lagrangian/` | xlrd, numpy, matplotlib, scipy, shapely |
| `read_data.py` | ~55 | ODS file reader & neighbor distribution | Console output, histogram | ezodf, matplotlib |
| `plot_histo.py` | ~65 | Per-frame histogram plotter | PNG histograms | xlrd, numpy, matplotlib |
| `topology_variation_animation.py` | ~280 | Simulated topology animation | GIF animation | numpy, matplotlib, scipy.Voronoi |
| `topology_variation_by_division.py` | ~65 | Static hexagon division visualization | Plot display | matplotlib, numpy |
| `tiff_to_png_invert.py` | ~240 | Batch TIFF→PNG with inversion | PNG images | PIL, numpy, argparse |

### Utility Scripts (Subdirectories)

| File | Purpose |
|------|---------|
| `image_analysis/python_napari.py` | Video processing with napari viewer (Z-projection) |
| `Icy/ij/plugins/Scripts/Python_Animation.py` | ImageJ rotation animation demo |

### Data Files

| Type | Location | Description |
|------|----------|-------------|
| Input XLS | `Xls_Data/all_frames_wing_discs.xls` | Main analysis data (~201 frames) |
| Input ODS | `data.ods` | Alternative format data |
| Pickled | `all_frames_wing_discs.pkl` | Serialized data cache |
| NumPy | `Data.npy`, `wing_initial_contour.npy` | Cached arrays |
| WKT | `wing_initial_contour.wkt` | Tissue boundary geometry |

### Output Directories

| Directory | Producer | Contents |
|-----------|----------|----------|
| `Results_wing_eulerian/` | Euler script | Fixed-grid analysis PNGs, heatmaps |
| `Results_wing_lagrangian/` | Lagrange script | Cell-tracking analysis PNGs |
| `Real_wing_tissue_figs_png/` | tiff_to_png script | Converted images |
| `wing_grids/` | — | Grid overlay visualizations |

---

## Core Analysis Scripts

### 1. Euler_wing_disc_cell_analysis.py

**Purpose:** Spatial analysis using fixed Eulerian grid (20×10 = 200 ROIs)

**Key Components:**

```python
# Data extraction from XLS
def extract_data(workbook, sheet_index) -> (x_coords, y_coords, polygon_numbers)

# Statistical analysis
def compute_topology_deviation_metrics(roi_data_list, all_polygon_numbers, ...)
    → Dict[support, C, N, p_roi, p_global, jsd, w1, z, mean_q, mean_abs_q]

# Dirichlet posterior smoothing
def _dirichlet_posterior_mean(counts, alpha=0.5) → probabilities

# Divergence metrics
def _js_divergence(P, Q, base=2.0) → JSD per ROI
def _wasserstein1_discrete(P, Q, support) → W1 distance per ROI

# Topological charge
def _topological_charge(counts, support) → (mean_q, mean_abs_q)
    # mean_q = E[n-6], where n = polygon sides
```

**Analysis Flow:**

1. **Load Data:** Read all frames from XLS workbook (201 sheets)
2. **Define Grid:** Create 20×10 spatial grid over tissue extent
3. **Assign Cells:** Digitize cell centroids to grid ROIs
4. **Compute Metrics:** 
   - Global distribution (mean ± std per shape)
   - Per-ROI distributions with Bayesian smoothing (Jeffreys prior α=0.5)
   - Divergence: JSD (information-theoretic), W1 (geometric)
   - Topological charge: E[n-6] (deviation from hexagonal packing)
5. **Visualize:**
   - Per-frame ROI maps with convex hull boundaries
   - Deviation heatmaps (JSD, W1, charge, per-shape deviations)
   - Cell count maps
   - Mean sidedness numbers overlay

**Key Outputs:**
- `Results_R2_euler_roi_frame_{i}.png` - Per-frame spatial maps
- `Results_R2_Mean_topology_over_201_frames.png` - Global histogram
- `absolute_deviation_map_{n}-gon.png` - Per-shape deviation heatmaps
- `map_topology_jsd.png`, `map_topology_w1.png` - Divergence maps
- `map_mean_sidedness_numbers.png` - Numerical overlay

**Strengths:**
- Comprehensive statistical framework (Bayesian, divergence metrics)
- Well-documented functions with type hints
- Handles boundary cells and empty ROIs gracefully

**Relationships:**
- **Reads from:** `Xls_Data/all_frames_wing_discs.xls`
- **Writes to:** `Results_wing_eulerian/`
- **Shares logic with:** Lagrange script (many copied functions)

---

### 2. Lagrange_wing_disc_analysis.py

**Purpose:** Cell-centered analysis tracking ~46 cells per ROI over time

**Key Components:**

```python
# ROI construction at t=0
def build_lagrangian_rois(ids0, x0, y0, target=50, ...) 
    → (cell_id_to_roi, roi_seeds_xy)
    # Uses farthest-point sampling + adjacency growth
    # Circularity refinement with robust outlier removal

# Tracking over time
def track_and_collect(frames, id2roi_0, roi_seeds_xy)
    → (roi_data_list, roi_index_per_frame)
    # Assigns new cells to nearest member or centroid
    # Periodic "round-up" to maintain circularity

# Circularity metrics
def _compactness(idx) → 4πA/P²  # 1 = perfect circle
def _radial_cv(idx, cx, cy) → std(r)/mean(r)

# Multi-strategy assignment
assign_strategy options:
    - "nearest_member" (default)
    - "nearest_seed"
    - "nearest_centroid_prev"
    - "hybrid_member_then_centroid"
```

**Analysis Flow:**

1. **Frame 0 Initialization:**
   - Farthest-point seeding (reproducible with `random_seed=7`)
   - Adjacency-based growth to ~target cells
   - Circularity repair: eject radial outliers (MAD-based threshold)
   
2. **Tracking (Frames 1-200):**
   - Known cells retain ROI membership
   - New cells assigned via nearest-neighbor to existing members
   - Boundary cells (-1 ID) assigned to nearest ROI centroid
   - Periodic round-up (every 3 frames) to restore circularity
   
3. **Collection:** Aggregate polygon distributions per ROI across all frames

4. **Metrics:** Same as Eulerian (JSD, W1, charge, deviations)

5. **Visualization:** ROI-colored scatter plots, deviation maps

**Key Outputs:**
- Per-frame scatter plots with color-coded ROIs
- Lagrangian deviation heatmaps
- CSV debug logs for specific ROIs

**Strengths:**
- Sophisticated tracking with multiple fallback strategies
- Circularity maintenance prevents ROI fragmentation
- Cell lineage preservation (follows divisions/deaths)

**Challenges:**
- Complex 2500+ lines (needs modularization)
- Some code duplication with Euler script
- Limited documentation for tuning parameters

**Relationships:**
- **Reads from:** Same XLS as Euler
- **Writes to:** `Results_wing_lagrangian/`
- **Shares ~60% code** with Euler (topology metrics)

---

### 3. read_data.py

**Purpose:** Quick ODS file reader for neighbor distribution analysis

**Functionality:**
- Opens `data.ods` using `ezodf` library
- Extracts cell labels and neighbor counts
- Computes neighbor distribution percentages (1-13 neighbors)
- Plots histogram

**Key Code:**
```python
doc = ezodf.opendoc('data.ods')
sheet = doc.sheets[0]
# Extract 'label' and 'num_neighbours' columns
neighbors = [row[num_neighbours_column_index].value for row in rows[1:]]
# Plot distribution
plt.hist(neighbors, bins=14, edgecolor='black')
```

**Usage:** Standalone diagnostic tool (not integrated into main pipeline)

**Relationships:**
- **Independent script** - doesn't share code with others
- **Alternative data source** - ODS vs XLS

---

### 4. plot_histo.py

**Purpose:** Generate per-frame cell shape histograms from XLS

**Functionality:**
```python
def process_and_plot(workbook, sheet_index, base_name):
    # Extract polygon numbers (skip boundary cells with ID=-1)
    # Calculate percentages per polygon type
    # Plot bar chart with percentages labeled
    # Save as '{base_name}_goodframe_{sheet_index}.png'
```

**Usage:**
```python
workbook = xlrd.open_workbook('t_0.xls')
for sheet_index in range(1, workbook.nsheets):
    process_and_plot(workbook, sheet_index, 't_0')
```

**Output:** Individual frame histograms (not aggregated)

**Relationships:**
- **Subset of Euler functionality** - could be integrated
- **Similar logic** to global distribution plots

---

### 5. topology_variation_animation.py

**Purpose:** Animated visualization of cell topology changes via division/death/T1 exchange

**Key Components:**
```python
def generate_hexagonal_lattice(rows, cols, spacing) → centers
def divide_cell(cell_centers, index) → new_cell_centers
def remove_cell(cell_centers, index) → updated_centers
def neighbor_exchange_t1(cell_centers) → transformed_centers

# Creates 4-panel animation:
# - Panel 1: Initial state
# - Panel 2: Cell division (green)
# - Panel 3: Cell death (red)
# - Panel 4: T1 neighbor exchange (orange)
```

**Uses Voronoi Tessellation:**
- Generates Voronoi diagram for each state
- Colors cells by event type
- Maintains tissue centrality via recentering

**Output:** Animated GIF or interactive matplotlib window

**Relationships:**
- **Educational/illustrative** - not analysis
- **Demonstrates concepts** used in real data analysis
- **Independent** - no data dependencies

---

### 6. topology_variation_by_division.py

**Purpose:** Static hexagon diagram showing topology change after division

**Functionality:**
- Draws hexagonal cell lattice with labels (n-gons)
- Shows before/after states of cell division
- Color-codes affected neighbors

**Usage:** Figure generation for presentations/papers

**Relationships:**
- **Companion to animation script**
- **Static version** of conceptual visualization

---

### 7. tiff_to_png_invert.py

**Purpose:** Batch convert TIFF microscopy images to inverted PNG

**Features:**
- Recursive directory walking
- Intensity inversion (useful for fluorescence: bright→dark)
- Alpha channel preservation
- uint8/uint16 support
- Mirror directory structure option

**Usage Examples:**
```bash
# Convert directory in-place
python tiff_to_png_invert.py Real_wing_tissue_figs

# Recursive with output mirroring
python tiff_to_png_invert.py Real_wing_tissue_figs --recursive -o Real_wing_tissue_figs_png
```

**Key Functions:**
```python
def invert_ndarray(arr: np.ndarray) → inverted_array
    # Inverts color channels, preserves alpha
def _iter_tiffs(paths, recursive) → Iterator[Path]
    # Yields TIFF files from paths
```

**Relationships:**
- **Preprocessing utility** - prepares figures for publication
- **Independent** - no analysis integration

---

## Utility Scripts

### image_analysis/python_napari.py

**Purpose:** Video processing with napari viewer

**Functionality:**
- Loads video via imageio (`embj2013197-sup-0002.mov`)
- Applies Z-projection (max intensity)
- Saves processed frames as PNG screenshots
- Exports processed video

**Status:** Exploratory/experimental (hardcoded paths)

---

### Icy/ij/plugins/Scripts/Python_Animation.py

**Purpose:** ImageJ Jython demo for rotational animation

**Context:** Example script from ImageJ, not project-specific

---

## Data Flow Architecture

### Input Data Pipeline

```
Microscopy Images (TIFF)
    ↓
[EpiTools Segmentation]
    ↓
Cell Tracking Data (XLS)
    ├─ Sheet 0-200: Time frames
    │   ├─ Column 0: Cell ID (-1 = boundary)
    │   ├─ Column 1: X coordinate
    │   ├─ Column 2: Y coordinate
    │   └─ Column 3: Polygon number (3-7+)
    ↓
Analysis Scripts (Euler/Lagrange)
    ↓
Visualization Outputs (PNG)
```

### Shared Data Structures

Both Euler and Lagrange scripts use:

```python
# Per-frame extraction
x_coords: np.ndarray  # (N,) cell X positions
y_coords: np.ndarray  # (N,) cell Y positions
polygon_numbers: np.ndarray  # (N,) cell shape class

# ROI data structure
roi_data_list: List[Dict] = [
    {
        'polygon_numbers': [frame0_array, frame1_array, ...],
        'cell_counts': [count0, count1, ...]
    },
    ...  # one dict per ROI
]

# Topology metrics output
metrics: Dict = {
    'support': np.array([3, 4, 5, 6, 7]),  # polygon classes
    'C': np.ndarray,  # (n_roi, n_classes) counts
    'N': np.ndarray,  # (n_roi,) totals
    'p_roi': np.ndarray,  # (n_roi, n_classes) probabilities
    'p_global': np.ndarray,  # (n_classes,) global distribution
    'jsd': np.ndarray,  # (n_roi,) Jensen-Shannon divergence
    'w1': np.ndarray,  # (n_roi,) Wasserstein-1 distance
    'z': np.ndarray,  # (n_roi, n_classes) standardized residuals
    'mean_q': np.ndarray,  # (n_roi,) mean topological charge
    'mean_abs_q': np.ndarray  # (n_roi,) mean |charge|
}
```

---

## Dependencies and Environment

### Python Environment (environment.yml)

```yaml
name: data_live_cell
python: 3.11
dependencies:
  - numpy          # Array operations, statistics
  - matplotlib     # Plotting, visualization
  - scipy          # Spatial algorithms (Delaunay, ConvexHull, KDTree, Voronoi)
  - shapely        # Geometric operations (Polygon, alpha-shapes)
  - xlrd           # Excel file reading (.xls format)
```

### Additional Dependencies (not in environment.yml)

- `ezodf` - ODS file reading (read_data.py)
- `PIL` (Pillow) - Image processing (tiff_to_png_invert.py)
- `napari` - Image viewer (image_analysis/python_napari.py)
- `imageio` - Video I/O (image_analysis/python_napari.py)

### External Tools

| Tool | Version | Purpose | Location |
|------|---------|---------|----------|
| EpiTools (MATLAB) | 2.1.6 | Cell segmentation | `epitools-matlab-2.1.6/` |
| EpiTools (Icy) | 0.8.8 | Cell tracking plugins | `epitools_part2_icy_plugins_v0.8.8/` |
| Icy | Latest | Image analysis platform | `Icy/` |
| MorphoNet | Standalone | 3D tissue visualization | `MorphoNet/` |

---

## External Tools and Libraries

### EpiTools Suite

**Purpose:** Automated epithelial tissue segmentation and tracking

**Components:**
1. **MATLAB Part** (`epitools-matlab-2.1.6/`, `epitools_part1_matlab_v2.1.6/`)
   - Watershed segmentation
   - Cell boundary detection
   - Initial tracking

2. **Icy Plugins** (`epitools_part2_icy_plugins_v0.8.8/`, `epitools-icy-master/`)
   - Manual correction interface
   - Graph-based cell tracking
   - XLS export

**Integration:** Generates the XLS input files consumed by Python scripts

**Citation:** See LICENSE.txt in respective directories

---

### MorphoNet

**Purpose:** 3D tissue morphology visualization (not directly used in current scripts)

**Location:** `MorphoNet/libs/STANDALONE/`

**Status:** Standalone application, no Python integration

---

## Recommendations for Improvement

### 1. Code Organization & Modularity

**Current Issues:**
- Two 900+ and 2500+ line monolithic scripts
- ~60% code duplication between Euler and Lagrange
- No separation of concerns (I/O, analysis, visualization mixed)

**Recommended Structure:**
```
wing_disc_analysis/
├── README.md
├── environment.yml
├── requirements.txt
├── setup.py
├── data/
│   ├── raw/                    # Original XLS/ODS files
│   ├── processed/              # Cached NumPy arrays
│   └── contours/               # Boundary geometries
├── src/
│   ├── __init__.py
│   ├── io/
│   │   ├── __init__.py
│   │   ├── xls_reader.py       # extract_data(), read_all_frames()
│   │   ├── ods_reader.py       # read_data.py logic
│   │   └── cache.py            # pickle/npy save/load
│   ├── geometry/
│   │   ├── __init__.py
│   │   ├── tessellation.py     # Delaunay, Voronoi, ConvexHull
│   │   ├── boundaries.py       # Alpha-shapes, contour detection
│   │   └── spatial.py          # KDTree queries, centroids
│   ├── topology/
│   │   ├── __init__.py
│   │   ├── metrics.py          # JSD, W1, topological charge
│   │   ├── distributions.py    # Dirichlet smoothing, global stats
│   │   └── deviations.py       # Residuals, comparisons
│   ├── roi/
│   │   ├── __init__.py
│   │   ├── eulerian.py         # Fixed-grid ROI construction
│   │   └── lagrangian.py       # Cell-tracking ROI construction & tracking
│   ├── visualization/
│   │   ├── __init__.py
│   │   ├── heatmaps.py         # Deviation maps, colored grids
│   │   ├── histograms.py       # Distribution plots
│   │   └── animations.py       # topology_variation_animation logic
│   └── utils/
│       ├── __init__.py
│       └── image_processing.py # TIFF inversion, preprocessing
├── scripts/
│   ├── run_eulerian_analysis.py      # Main entry point
│   ├── run_lagrangian_analysis.py    # Main entry point
│   ├── plot_single_frame.py          # plot_histo.py refactored
│   └── convert_images.py             # tiff_to_png_invert.py
├── notebooks/
│   ├── 01_exploratory_analysis.ipynb
│   ├── 02_eulerian_workflow.ipynb
│   └── 03_lagrangian_workflow.ipynb
├── tests/
│   ├── test_io.py
│   ├── test_topology.py
│   └── test_roi.py
└── outputs/
    ├── eulerian/
    └── lagrangian/
```

**Benefits:**
- **Reusability:** Shared functions in modules, imported by scripts
- **Testability:** Unit tests for each module
- **Maintainability:** Changes isolated to relevant modules
- **Discoverability:** Clear hierarchy guides new users

---

### 2. Documentation

**Missing Elements:**
- Project README with scientific context
- Installation instructions
- Usage examples
- API documentation (docstrings are good but need aggregation)
- Parameter tuning guide

**Recommended Additions:**

**README.md:**
```markdown
# Wing Disc Epithelial Topology Analysis

## Overview
Quantitative analysis of cell shape distributions in Drosophila wing disc 
epithelium using Eulerian and Lagrangian frameworks.

## Installation
```bash
conda env create -f environment.yml
conda activate data_live_cell
pip install -r requirements.txt
```

## Quick Start
```bash
# Run Eulerian analysis
python scripts/run_eulerian_analysis.py --data data/raw/all_frames_wing_discs.xls

# Run Lagrangian analysis
python scripts/run_lagrangian_analysis.py --data data/raw/all_frames_wing_discs.xls
```

## Scientific Background
[Explain tissue topology, Euler characteristic, topological transitions...]

## Methods
- **Eulerian ROIs:** 20×10 fixed spatial grid
- **Lagrangian ROIs:** ~46 cells tracked over time
- **Metrics:** JSD, Wasserstein-1, topological charge E[n-6]

## Citation
[Publication info]
```

**requirements.txt:**
```
numpy>=1.24
matplotlib>=3.7
scipy>=1.10
shapely>=2.0
xlrd>=2.0
ezodf>=0.9
Pillow>=9.0
```

**API Docs:** Use Sphinx to generate from docstrings

---

### 3. Code Deduplication

**Shared Functions to Extract:**

From both Euler and Lagrange:
```python
# io/xls_reader.py
def extract_data(workbook, sheet_index)
def extract_frame(sheet)
def read_all_frames(wb)

# topology/metrics.py
def compute_topology_deviation_metrics(...)
def _infer_support(...)
def _tally_roi_counts(...)
def _dirichlet_posterior_mean(...)
def _js_divergence(...)
def _wasserstein1_discrete(...)
def _multinomial_residuals(...)
def _topological_charge(...)

# topology/distributions.py
def compute_global_distribution(all_polys, thresh=0.5)
def mean_sidedness_per_roi(metrics, ...)

# visualization/heatmaps.py
def plot_mean_sidedness_numbers(...)
def plot_roi_deviation_heatmaps(...)
```

**Estimated Reduction:** ~800 lines of duplicated code eliminated

---

### 4. Configuration Management

**Current Issues:**
- Hardcoded parameters scattered throughout scripts
- No easy way to run parameter sweeps
- Magic numbers (e.g., `num_grid_x=20`, `target_cells_per_roi=46`)

**Recommended Approach:**

**config/eulerian_config.yaml:**
```yaml
data:
  input_file: "data/raw/all_frames_wing_discs.xls"
  output_dir: "outputs/eulerian"

roi:
  num_grid_x: 20
  num_grid_y: 10

analysis:
  min_global_pct: 0.2
  dirichlet_alpha: 0.5

visualization:
  page_width_pt: 455.24411
  page_height_pt: 702.78308
  dpi: 300
  colormaps:
    deviation: "viridis"
    residual: "coolwarm"
```

**Usage in script:**
```python
import yaml

with open('config/eulerian_config.yaml') as f:
    config = yaml.safe_load(f)

num_grid_x = config['roi']['num_grid_x']
# ...
```

---

### 5. Testing Strategy

**Current State:** No automated tests

**Recommended Tests:**

```python
# tests/test_io.py
def test_extract_data_skips_boundary_cells():
    # Mock workbook with ID=-1 rows
    # Assert they're filtered out
    
def test_extract_frame_handles_empty_sheet():
    # Assert graceful handling

# tests/test_topology.py
def test_dirichlet_posterior_mean_symmetry():
    counts = np.array([10, 20, 30])
    alpha = 1.0
    result = _dirichlet_posterior_mean(counts, alpha)
    assert np.isclose(result.sum(), 1.0)
    
def test_topological_charge_hexagons():
    # All hexagons should give E[n-6]=0
    counts = np.array([0, 0, 0, 100, 0])  # all 6-gons
    support = np.array([3, 4, 5, 6, 7])
    mean_q, _ = _topological_charge(counts, support)
    assert np.isclose(mean_q, 0.0)

# tests/test_roi.py
def test_eulerian_grid_coverage():
    # All cells should be assigned to exactly one ROI
    
def test_lagrangian_tracking_persistence():
    # Cell IDs should maintain ROI membership if not reassigned
```

**Run with:** `pytest tests/`

---

### 6. Performance Optimization

**Current Bottlenecks:**
- No caching of intermediate results
- Redundant Delaunay triangulations per frame
- Inefficient grid assignment (digitize called per frame)

**Optimizations:**

```python
# Cache expensive computations
from functools import lru_cache

@lru_cache(maxsize=256)
def compute_delaunay(points_tuple):
    points = np.array(points_tuple)
    return Delaunay(points)

# Vectorize grid assignments
def assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges):
    # Use numpy broadcasting instead of loops
    x_idx = np.searchsorted(x_edges, x_coords, side='right') - 1
    y_idx = np.searchsorted(y_edges, y_coords, side='right') - 1
    return np.clip(x_idx, 0, len(x_edges)-2), np.clip(y_idx, 0, len(y_edges)-2)

# Parallelize frame processing
from multiprocessing import Pool

def process_frame(args):
    sheet_index, workbook_path = args
    # ...
    
with Pool(processes=4) as pool:
    results = pool.map(process_frame, frame_args)
```

---

### 7. Reproducibility

**Add:**
- Random seed control for all stochastic operations
- Version pinning in requirements.txt
- Docker container for complete environment isolation
- Data provenance tracking (log parameters, git commit, timestamps)

**Example logging:**
```python
import logging
import json
from datetime import datetime

logging.basicConfig(
    filename=f'outputs/run_{datetime.now():%Y%m%d_%H%M%S}.log',
    level=logging.INFO
)

config_log = {
    'timestamp': datetime.now().isoformat(),
    'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode().strip(),
    'config': config,
    'environment': {
        'numpy': np.__version__,
        'scipy': scipy.__version__,
        # ...
    }
}

logging.info(json.dumps(config_log, indent=2))
```

---

### 8. Interactive Analysis

**Add Jupyter Notebooks:**

```python
# notebooks/02_eulerian_workflow.ipynb

import wing_disc_analysis as wda

# Load data
data = wda.io.load_xls('data/raw/all_frames_wing_discs.xls')

# Create ROIs
rois = wda.roi.EulerianGrid(nx=20, ny=10)
rois.fit(data.frames[0])

# Compute metrics
metrics = wda.topology.compute_metrics(rois, data)

# Interactive visualization
import plotly.graph_objects as go

fig = wda.visualization.interactive_heatmap(
    metrics['jsd'], 
    rois.grid,
    title='Jensen-Shannon Divergence'
)
fig.show()
```

**Benefits:**
- Exploratory analysis without modifying scripts
- Parameter tuning with immediate feedback
- Educational tool for new collaborators

---

### 9. Version Control & Collaboration

**Current Issues:**
- Large binary files (TIFF, XLS, PKL) tracked in git
- No .gitignore (risks committing temp files)

**Recommendations:**

**.gitignore:**
```
# Data files
*.xls
*.xlsx
*.ods
*.pkl
*.npy
*.tif
*.tiff
*.png
*.gif
*.avi
*.mov
*.mp4

# Output directories
Results_*/
outputs/
Real_wing_tissue_figs_png/

# Python
__pycache__/
*.pyc
*.egg-info/
.pytest_cache/
.ipynb_checkpoints/

# IDEs
.idea/
.vscode/
.DS_Store
```

**Use Git LFS for essential data:**
```bash
git lfs track "*.xls"
git lfs track "*.pkl"
```

**Or use external data hosting:**
- Zenodo for published datasets
- Lab server with documented download instructions

---

### 10. Command-Line Interface

**Current:** Scripts require editing source code to change parameters

**Recommended:** Add argparse CLI

```python
# scripts/run_eulerian_analysis.py

import argparse
from wing_disc_analysis.roi import EulerianGrid
from wing_disc_analysis.io import load_xls

def main():
    parser = argparse.ArgumentParser(description='Eulerian ROI topology analysis')
    parser.add_argument('--data', required=True, help='Path to XLS workbook')
    parser.add_argument('--grid-x', type=int, default=20, help='Grid divisions (X)')
    parser.add_argument('--grid-y', type=int, default=10, help='Grid divisions (Y)')
    parser.add_argument('--output', default='outputs/eulerian', help='Output directory')
    parser.add_argument('--alpha', type=float, default=0.5, help='Dirichlet prior')
    parser.add_argument('--dpi', type=int, default=300, help='Figure DPI')
    
    args = parser.parse_args()
    
    # Run analysis with args
    data = load_xls(args.data)
    grid = EulerianGrid(nx=args.grid_x, ny=args.grid_y)
    # ...

if __name__ == '__main__':
    main()
```

**Usage:**
```bash
python scripts/run_eulerian_analysis.py \
    --data data/raw/wing_disc.xls \
    --grid-x 15 --grid-y 15 \
    --output outputs/test_15x15 \
    --alpha 1.0
```

---

## Summary of Current File Relationships

```
┌─────────────────────────────────────────────────────────┐
│                    Data Sources                         │
├─────────────────────────────────────────────────────────┤
│  Xls_Data/all_frames_wing_discs.xls  (201 sheets)      │
│  data.ods  (alternative format)                         │
│  TIFF images (Real_wing_tissue_figs/)                   │
└────────────┬────────────────────────────────────────────┘
             │
             ├──────────────────────┬─────────────────────┐
             │                      │                     │
             ▼                      ▼                     ▼
┌────────────────────┐  ┌────────────────────┐  ┌────────────────┐
│ Euler Analysis     │  │ Lagrange Analysis  │  │ Utilities      │
│ (900 lines)        │  │ (2500 lines)       │  │                │
├────────────────────┤  ├────────────────────┤  ├────────────────┤
│ • Fixed grid ROI   │  │ • Cell-tracking    │  │ • read_data    │
│ • 20×10 = 200 ROIs│  │ • ~46 cells/ROI    │  │ • plot_histo   │
│ • Spatial maps     │  │ • Lagrangian track │  │ • TIFF convert │
│ • JSD, W1, charge  │  │ • Circularity fix  │  │ • Animations   │
│                    │  │ • Same metrics     │  │                │
└────────┬───────────┘  └────────┬───────────┘  └────────┬───────┘
         │                       │                       │
         │    Shared Code (60% overlap):                 │
         │    • extract_data()                           │
         │    • compute_topology_deviation_metrics()     │
         │    • _js_divergence(), _wasserstein1_discrete()│
         │    • _topological_charge()                    │
         │    • Plotting functions                       │
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────────────────────────────────────────────┐
│                    Outputs                              │
├─────────────────────────────────────────────────────────┤
│  Results_wing_eulerian/  (spatial heatmaps, histograms) │
│  Results_wing_lagrangian/  (tracking plots)             │
│  PNG images, histograms, animations                     │
└─────────────────────────────────────────────────────────┘
```

---

## Conclusion

This project demonstrates sophisticated quantitative biology with robust statistical methods (Bayesian smoothing, information-theoretic divergence, geometric metrics). The analysis is scientifically sound but suffers from engineering debt:

**Strengths:**
✅ Comprehensive topology metrics  
✅ Dual analysis frameworks (Eulerian + Lagrangian)  
✅ Well-documented functions with type hints  
✅ Handles edge cases (boundary cells, empty ROIs)  

**Areas for Improvement:**
❌ Code duplication (~60% overlap between main scripts)  
❌ Monolithic structure (900-2500 line files)  
❌ No modularization or package structure  
❌ Missing tests, CI/CD, and reproducibility controls  
❌ Hardcoded parameters, no CLI  
❌ Incomplete dependency specification  

**Priority Recommendations:**
1. **Immediate:** Create requirements.txt, add .gitignore, write README
2. **Short-term:** Extract shared functions into modules, add CLI arguments
3. **Medium-term:** Restructure into package, add tests, create notebooks
4. **Long-term:** Performance optimization, Docker containerization, documentation site

With these improvements, the project would transition from "functional research code" to "maintainable scientific software" suitable for publication, collaboration, and long-term use.

---

**Report Generated:** November 28, 2025  
**Total Python Scripts Analyzed:** 8 core + 2 utility = 10 files  
**Total Lines of Code:** ~5,000+ lines  
**Recommended Refactor Target:** ~2,000 lines (60% reduction via deduplication)

