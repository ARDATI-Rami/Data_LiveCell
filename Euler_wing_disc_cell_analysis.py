import xlrd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull,Delaunay
from scipy.stats import entropy
import os
import random
from shapely.geometry import Polygon,MultiLineString
from shapely.ops import polygonize, unary_union
from matplotlib import colors  # Ensure this import is at the top of your script
from typing import Dict, List, Optional, Tuple



page_width_in = 455.24411 / 72.27   # ≈ 6.30 in for Pixels = (455.24411 / 72.27) × 300 = 1,889.55 pixels
page_height_in = 702.78308 / 72.27  # ≈ 9.73 in
page_w = page_width_in
page_h = page_height_in
print(f"page size in pixels at 300dpi: {page_w*300:.2f} x {page_h*300:.2f}")


# Define shape names
shape_names = {3: 'Triangle', 4: 'Quadrilateral', 5: 'Pentagon', 6: 'Hexagon', 7: 'Heptagon'}
stamp_prefix = 'Results_R2_'
save_dir = 'Results_wing_eulerian'
os.makedirs(save_dir, exist_ok=True)

# Function to extract data from a sheet in the workbook
def extract_data(workbook, sheet_index):
    # Load the sheet
    sheet = workbook.sheet_by_index(sheet_index)

    # Extract data
    x_coords = []
    y_coords = []
    polygon_numbers = []

    # Assuming the data starts from the second row (skip header)
    for row_idx in range(1, sheet.nrows):
        # Skip rows where the value in column 0 is -1
        if sheet.cell_value(row_idx, 0) == -1:
            continue

        x = sheet.cell_value(row_idx, 1)  # 2nd column (index 1) contains 'X coordinate'
        y = sheet.cell_value(row_idx, 2)  # 3rd column (index 2) contains 'Y coordinate'
        polygon_no = sheet.cell_value(row_idx, 3)  # 4th column (index 3) contains 'Polygon no'

        x_coords.append(float(x))
        y_coords.append(float(y))
        polygon_numbers.append(int(polygon_no))  # Convert to integer

    return np.array(x_coords), np.array(y_coords), np.array(polygon_numbers)

def _infer_support(all_data: List[np.ndarray],
                   min_global_pct: float = 0.0) -> np.ndarray:
    """
    Infer polygon classes present globally and optionally filter rare classes.

    :param all_data: List of 1D arrays (per frame) with polygon numbers (e.g., 3..10).
    :type  all_data: List[np.ndarray]
    :param min_global_pct: Minimum global percentage to keep a class (0..100).
    :type  min_global_pct: float
    :return: Sorted unique polygon classes retained.
    :rtype:  np.ndarray
    """
    g = np.concatenate(all_data) if len(all_data) else np.array([], dtype=int)
    if g.size == 0:
        return np.array([], dtype=int)
    classes, counts = np.unique(g, return_counts=True)
    pct = 100 * counts / counts.sum()
    keep = pct >= min_global_pct
    return classes[keep]

def _tally_roi_counts(roi_data_list: List[Dict],
                      support: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert ROI polygon lists to count matrices over a fixed support.

    :param roi_data_list: List of ROI dicts, each with 'polygon_numbers' as List[np.ndarray] over frames.
    :type  roi_data_list: List[Dict]
    :param support: Sorted polygon classes to count (e.g., np.array([3,4,5,6,7])).
    :type  support: np.ndarray
    :return: (C, N), where C[i, k] are counts of support[k] in ROI i across all frames,
             and N[i] is total counts in ROI i.
    :rtype:  Tuple[np.ndarray, np.ndarray]
    """
    num_rois = len(roi_data_list)
    K = len(support)
    C = np.zeros((num_rois, K), dtype=np.int64)
    for i, roi in enumerate(roi_data_list):
        # Concatenate all frames for this ROI
        arr = np.concatenate(roi.get('polygon_numbers', [])) if roi.get('polygon_numbers') else np.array([], dtype=int)
        if arr.size:
            # Vectorized counting per class
            C[i, :] = np.array([(arr == cls).sum() for cls in support], dtype=np.int64)
    N = C.sum(axis=1)
    return C, N


def _dirichlet_posterior_mean(counts: np.ndarray, alpha: float = 0.5) -> np.ndarray:
    """
    Compute Dirichlet posterior mean probabilities with symmetric prior.

    :param counts: Count vector or matrix [..., K] for K classes.
    :type  counts: np.ndarray
    :param alpha: Symmetric Dirichlet prior concentration per class (Jeffreys ~0.5).
    :type  alpha: float
    :return: Posterior mean probabilities with shape matching 'counts'.
    :rtype:  np.ndarray
    """
    counts = np.asarray(counts, dtype=np.float64)
    K = counts.shape[-1]
    num = counts + alpha
    den = counts.sum(axis=-1, keepdims=True) + alpha * K
    return num / np.clip(den, 1e-12, None)


def _js_divergence(P: np.ndarray, Q: np.ndarray, base: float = 2.0) -> np.ndarray:
    """
    Jensen–Shannon divergence between rows of P and a single Q.

    :param P: Array [M, K] of probability vectors (ROIs).
    :type  P: np.ndarray
    :param Q: Array [K] global probability vector.
    :type  Q: np.ndarray
    :param base: Log base (2 => bits).
    :type  base: float
    :return: JSD per row of P.
    :rtype:  np.ndarray
    """
    # Ensure proper shapes
    Q = np.asarray(Q, dtype=np.float64)
    P = np.asarray(P, dtype=np.float64)
    M = 0.5 * (P + Q[None, :])
    # KL(P||M) + KL(Q||M) averaged
    # scipy.stats.entropy broadcasts over rows if we iterate
    jsd = np.empty(P.shape[0], dtype=np.float64)
    for i in range(P.shape[0]):
        jsd[i] = 0.5 * entropy(P[i], M[i], base=base) + 0.5 * entropy(Q, M[i], base=base)
    return jsd


def _wasserstein1_discrete(P: np.ndarray, Q: np.ndarray, support: np.ndarray) -> np.ndarray:
    """
    1D Wasserstein-1 distance (Earth Mover's Distance) between rows of P and Q over ordered support.

    For discrete measures on ordered points x_1 < ... < x_K, W1 = Σ_{k=1}^{K-1} |CDF_P(k) - CDF_Q(k)| * (x_{k+1}-x_k).

    :param P: Array [M, K] of probability vectors (ROIs).
    :type  P: np.ndarray
    :param Q: Array [K] global probability vector.
    :type  Q: np.ndarray
    :param support: Sorted support values (polygon classes).
    :type  support: np.ndarray
    :return: Wasserstein-1 per row of P (same units as support, here “sides”).
    :rtype:  np.ndarray
    """
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    cdf_P = np.cumsum(P, axis=1)
    cdf_Q = np.cumsum(Q)
    gaps = np.diff(support.astype(np.float64))
    # |CDF diff| at internal cutpoints times gap sizes
    return np.sum(np.abs(cdf_P[:, :-1] - cdf_Q[:-1][None, :]) * gaps[None, :], axis=1)


def _multinomial_residuals(counts: np.ndarray, p0: np.ndarray) -> np.ndarray:
    """
    Per-class standardized residuals under multinomial(n, p0):
    z_k = (x_k - n p0_k) / sqrt(n p0_k (1 - p0_k))
    (Not independent across classes, but useful diagnostically.)

    :param counts: Array [M, K] of counts per ROI.
    :type  counts: np.ndarray
    :param p0: Array [K] of baseline probabilities.
    :type  p0: np.ndarray
    :return: Array [M, K] of z-scores; NaN for n=0 or p0_k≈0.
    :rtype:  np.ndarray
    """
    counts = counts.astype(np.float64)
    n = counts.sum(axis=1, keepdims=True)
    exp = n * p0[None, :]
    var = n * p0[None, :] * (1.0 - p0[None, :])
    z = (counts - exp) / np.sqrt(np.maximum(var, 1e-12))
    z[np.where(n == 0)[0], :] = np.nan
    # Avoid exploding where p0≈0
    z[:, p0 < 1e-9] = np.nan
    return z


def _topological_charge(counts: np.ndarray, support: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Mean charge E[n-6] and mean absolute charge E[|n-6|] per ROI.

    :param counts: Array [M, K] of counts per ROI.
    :type  counts: np.ndarray
    :param support: Sorted polygon classes (e.g., [3,4,5,6,7]).
    :type  support: np.ndarray
    :return: (mean_q, mean_abs_q) each [M].
    :rtype:  Tuple[np.ndarray, np.ndarray]
    """
    n = counts.sum(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        P = counts / np.clip(n, 1e-12, None)
    q = support.astype(np.float64) - 6.0
    mean_q = np.nansum(P * q[None, :], axis=1)
    mean_abs_q = np.nansum(P * np.abs(q)[None, :], axis=1)
    mean_q[n.ravel() == 0] = np.nan
    mean_abs_q[n.ravel() == 0] = np.nan
    return mean_q, mean_abs_q


def compute_topology_deviation_metrics(
    roi_data_list: List[Dict],
    all_polygon_numbers: List[np.ndarray],
    polygons: Optional[np.ndarray] = None,
    alpha: float = 0.5,
    min_global_pct: float = 0.0
) -> Dict[str, np.ndarray]:
    """
    Compute principled topology (cell-shape) deviation metrics per ROI.

    This function:
      1) pools counts across frames;
      2) applies Dirichlet(α) smoothing (Jeffreys α≈0.5 by default);
      3) compares ROI posterior means to the global posterior mean with:
         - Jensen–Shannon divergence (bits),
         - 1D Wasserstein distance over polygon classes (sides),
         - per-class standardized residuals (z),
         - mean charge E[n-6] and mean |n-6|.

    :param roi_data_list: Each ROI dict must have 'polygon_numbers': List[np.ndarray] over frames.
    :type  roi_data_list: List[Dict]
    :param all_polygon_numbers: Per-frame arrays of polygon numbers over the whole tissue.
    :type  all_polygon_numbers: List[np.ndarray]
    :param polygons: Optional sorted array of polygon classes to include; if None, inferred globally.
    :type  polygons: Optional[np.ndarray]
    :param alpha: Symmetric Dirichlet prior concentration per class for smoothing.
    :type  alpha: float
    :param min_global_pct: Drop classes with global percentage below this threshold.
    :type  min_global_pct: float
    :return: Dict with keys:
             'support'        -> [K] polygon classes,
             'C'              -> [R, K] counts per ROI,
             'N'              -> [R] totals per ROI,
             'p_roi'          -> [R, K] posterior mean per ROI,
             'p_global'       -> [K] global posterior mean,
             'jsd'            -> [R] Jensen–Shannon divergence (bits),
             'w1'             -> [R] Wasserstein-1 (sides),
             'z'              -> [R, K] per-class z residuals,
             'mean_q'         -> [R] E[n-6],
             'mean_abs_q'     -> [R] E[|n-6|].
    :rtype:  Dict[str, np.ndarray]
    """
    # 1) Support
    if polygons is None:
        support = _infer_support(all_polygon_numbers, min_global_pct=min_global_pct)
    else:
        support = np.asarray(polygons, dtype=int)
    if support.size == 0:
        raise ValueError("No polygon classes found/kept. Check inputs or min_global_pct.")

    # 2) Counts per ROI and global
    C, N = _tally_roi_counts(roi_data_list, support)
    G_counts = C.sum(axis=0)
    # If you prefer the frame-pooled global (not via ROIs), uncomment:
    # G_all = np.concatenate(all_polygon_numbers)
    # G_counts = np.array([(G_all == cls).sum() for cls in support], dtype=np.int64)

    # 3) Posterior means with Dirichlet smoothing
    p_roi = _dirichlet_posterior_mean(C, alpha=alpha)
    p_global = _dirichlet_posterior_mean(G_counts, alpha=alpha)

    # 4) Metrics
    jsd = _js_divergence(p_roi, p_global, base=2.0)
    w1 = _wasserstein1_discrete(p_roi, p_global, support)
    z = _multinomial_residuals(C, p_global)
    mean_q, mean_abs_q = _topological_charge(C, support)

    return dict(
        support=support, C=C, N=N,
        p_roi=p_roi, p_global=p_global,
        jsd=jsd, w1=w1, z=z,
        mean_q=mean_q, mean_abs_q=mean_abs_q
    )

def mean_sidedness_per_roi(metrics: dict,
                           posterior: bool = True,
                           empty_policy: str = "nan") -> np.ndarray:
    """
    Mean polygon sidedness E[n] per ROI with explicit handling of empty ROIs.

    For ROIs with zero total count (metrics['N'] == 0), the posterior mean would
    be the uniform average of the support (often exactly 6). Use `empty_policy`
    to override that.

    :param metrics: Output dict from `compute_topology_deviation_metrics`.
    :type  metrics: dict
    :param posterior: If True, use Dirichlet-smoothed p_roi; if False, use empirical
                      pooled counts via E[n] = 6 + mean_q.
    :type  posterior: bool
    :param empty_policy: How to handle ROIs with N==0. Options:
                         - "nan": return NaN (recommended; will be skipped in plotting).
                         - "global": use global mean sidedness sum(n * p_global(n)).
                         - "uniform": keep posterior uniform mean (usual behavior).
                         - a numeric string, e.g. "6.0": cast to float and use that.
    :type  empty_policy: str
    :return: Flat array [R] with mean sidedness per ROI.
    :rtype:  np.ndarray
    """
    support = metrics['support'].astype(float)
    if posterior:
        mean_n = metrics['p_roi'] @ support
    else:
        mean_n = 6.0 + metrics['mean_q']  # NaN when N==0

    # Identify empties using pooled counts
    empties = (metrics['N'] == 0)

    if empty_policy == "nan":
        mean_n = mean_n.astype(float)
        mean_n[empties] = np.nan
    elif empty_policy == "global":
        g_mean = float(metrics['p_global'] @ support)
        mean_n = mean_n.astype(float)
        mean_n[empties] = g_mean
    elif empty_policy == "uniform":
        # keep as-is (posterior uniform -> mean(support) ; empirical already NaN)
        pass
    else:
        # try numeric override
        try:
            fill = float(empty_policy)
            mean_n = mean_n.astype(float)
            mean_n[empties] = fill
        except Exception:
            # fallback to NaN if parsing failed
            mean_n = mean_n.astype(float)
            mean_n[empties] = np.nan
    return mean_n


def plot_mean_sidedness_numbers(
    metrics: dict,
    x_grid_edges: np.ndarray,
    y_grid_edges: np.ndarray,
    num_grid_x: int,
    num_grid_y: int,
    page_w: float,
    page_h: float,
    save_path: str,
    hull_points: np.ndarray | None = None,
    posterior: bool = True,
    empty_policy: str = "nan",
    text_fmt: str = "{:.2f}",
    na_label: str = "",  # e.g., "—" if you want to visibly mark empties
    cmap: str = "viridis",  # Colormap for background
    vmin: float = None,  # Min value for colormap (None = auto)
    vmax: float = None   # Max value for colormap (None = auto)
) -> None:
    """
    Write mean sidedness E[n] as numbers centered in each ROI with colored background.

    :param metrics: Output dict from `compute_topology_deviation_metrics`.
    :type  metrics: dict
    :param x_grid_edges: X edges of the grid (length = num_grid_x+1).
    :type  x_grid_edges: np.ndarray
    :param y_grid_edges: Y edges of the grid (length = num_grid_y+1).
    :type  y_grid_edges: np.ndarray
    :param num_grid_x: Number of grid divisions along x.
    :type  num_grid_x: int
    :param num_grid_y: Number of grid divisions along y.
    :type  num_grid_y: int
    :param page_w: Figure width in inches.
    :type  page_w: float
    :param page_h: Figure height in inches.
    :type  page_h: float
    :param save_path: Output PNG path.
    :type  save_path: str
    :param hull_points: Optional Nx2 array to draw the tissue hull contour.
    :type  hull_points: np.ndarray | None
    :param posterior: If True, use posterior-smoothed E[n]; else empirical.
    :type  posterior: bool
    :param empty_policy: Passed to `mean_sidedness_per_roi` for N==0 handling.
    :type  empty_policy: str
    :param text_fmt: Python format string for values, e.g., "{:.2f}".
    :type  text_fmt: str
    :param na_label: If not empty and a ROI is NaN, draw this label instead.
    :type  na_label: str
    :param cmap: Colormap name for the background colors.
    :type  cmap: str
    :param vmin: Minimum value for colormap scaling.
    :type  vmin: float
    :param vmax: Maximum value for colormap scaling.
    :type  vmax: float
    :return: None.
    :rtype:  None
    """
    mean_n = mean_sidedness_per_roi(metrics, posterior=posterior, empty_policy=empty_policy)
    mean_n_grid = mean_n.reshape((num_grid_x, num_grid_y))

    # Create a masked array to handle NaN values properly
    mean_n_grid_masked = np.ma.masked_invalid(mean_n_grid.T)  # Transpose for pcolormesh

    fig, ax = plt.subplots(figsize=(page_w, page_h/3))

    # Create the background colormap
    X, Y = np.meshgrid(x_grid_edges, y_grid_edges)
    pcm = ax.pcolormesh(X, Y, mean_n_grid_masked, cmap=cmap,
                       shading='auto', vmin=vmin, vmax=vmax)

    # Add a colorbar
    cbar = fig.colorbar(pcm, ax=ax,
                      label='Mean polygon sidedness E[n]', shrink=0.9)
    cbar.ax.tick_params(labelsize=8)

    # Add text labels with values
    dx = (x_grid_edges[1] - x_grid_edges[0]) / 2.0
    dy = (y_grid_edges[1] - y_grid_edges[0]) / 2.0
    for i in range(num_grid_x):
        for j in range(num_grid_y):
            val = mean_n_grid[i, j]
            cx = x_grid_edges[i] + dx
            cy = y_grid_edges[j] + dy
            if np.isfinite(val):
                # Use white text for dark backgrounds, black text for light backgrounds
                # This assumes the default colormap is sequential
                norm_val = (val - (vmin if vmin is not None else np.nanmin(mean_n_grid))) / \
                          ((vmax if vmax is not None else np.nanmax(mean_n_grid)) -
                           (vmin if vmin is not None else np.nanmin(mean_n_grid)))
                text_color = 'black' if norm_val > 0.5 else 'white'
                ax.text(cx, cy, text_fmt.format(val), ha='center', va='center',
                        fontsize=4.2, color=text_color, weight='bold')
            elif na_label:
                ax.text(cx, cy, na_label, ha='center', va='center',
                        fontsize=4.2, color='black')

    # Grid lines (dashed lines)
    for x in x_grid_edges:
        ax.plot([x, x], [y_grid_edges[0], y_grid_edges[-1]], 'k--', lw=0.5, alpha=0.7)
    for y in y_grid_edges:
        ax.plot([x_grid_edges[0], x_grid_edges[-1]], [y, y], 'k--', lw=0.5, alpha=0.7)

    # Draw hull if provided
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
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)


# Parameters to control which plots to generate
plot_contour = True
plot_mean_of_all_frames_cell_shapes = True
plot_roi_analysis = True
plot_roi_deviations = True  # New flag to control deviation plotting

# Parameters for Grid ROI definitions
num_grid_x = 20  # Number of grid divisions along x-axis
num_grid_y = 10  # Number of grid divisions along y-axis

# Generate ROI labels automatically
roi_labels = [f'ROI ({i}, {j})' for i in range(num_grid_x) for j in range(num_grid_y)]

# Generate ROI colors using a colormap with more colors
num_rois = num_grid_x * num_grid_y
cmap = plt.get_cmap('gist_ncar', num_rois)  # 'nipy_spectral' provides a wide range of colors
roi_colors = [cmap(i / num_rois) for i in range(num_rois)]  # Sample colors at regular intervals
random.shuffle(roi_colors)  # This shuffles the list in place


# Load the .xls file
file_path = "Xls_Data"
file_name = 'all_frames_wing_discs.xls'  # Replace with your actual file name
workbook = xlrd.open_workbook(os.path.join(file_path, file_name))
print(f'Number of sheets : {workbook.nsheets}')

# Initialize data storage
all_data = []
all_x_coords = []
all_y_coords = []
all_polygon_numbers = []

# Initialize ROI data storage
roi_data_list = [{} for _ in range(num_rois)]  # List of dictionaries for each ROI
for roi_data in roi_data_list:
    roi_data['polygon_numbers'] = []
    roi_data['cell_counts'] = []  # New list to store cell counts per frame





# Loop through each sheet in the workbook, starting from the first sheet (index 0)
for sheet_index in range(0, workbook.nsheets):
    x_coords, y_coords, polygon_numbers = extract_data(workbook, sheet_index)
    print(f'Sheet {sheet_index+1} : Number of cells = {len(x_coords)}')

    # Store data for mean percentage calculation
    all_data.append(polygon_numbers)
    all_x_coords.append(x_coords)
    all_y_coords.append(y_coords)
    all_polygon_numbers.append(polygon_numbers)

# Compute overall min and max x and y coordinates for consistent grid edges
min_x = min([x.min() for x in all_x_coords])
max_x = max([x.max() for x in all_x_coords])
min_y = min([y.min() for y in all_y_coords])
max_y = max([y.max() for y in all_y_coords])

# Define grid boundaries
x_grid_edges = np.linspace(min_x, max_x, num_grid_x + 1)
y_grid_edges = np.linspace(min_y, max_y, num_grid_y + 1)

# Loop through each sheet again for consistency
for sheet_index in range(1, workbook.nsheets+1):
    x_coords = all_x_coords[sheet_index - 1]
    y_coords = all_y_coords[sheet_index - 1]
    polygon_numbers = all_polygon_numbers[sheet_index - 1]

    # Assign cells to grid cells (ROIs)
    x_indices = np.digitize(x_coords, x_grid_edges) - 1  # Subtract 1 to get indices starting from 0
    y_indices = np.digitize(y_coords, y_grid_edges) - 1

    # Ensure indices are within bounds
    x_indices = np.clip(x_indices, 0, num_grid_x - 1)
    y_indices = np.clip(y_indices, 0, num_grid_y - 1)

    # Compute ROI indices
    roi_indices = x_indices * num_grid_y + y_indices  # Unique index for each grid cell

    # Collect cell counts for each ROI
    for i in range(num_rois):
        roi_mask = roi_indices == i
        cell_count = np.sum(roi_mask)
        roi_data_list[i]['cell_counts'].append(cell_count)  # Store cell count for this frame

    # Plot contour if enabled
    if plot_contour or plot_roi_analysis:
        # Combine x and y coordinates into a single array
        points = np.column_stack((x_coords, y_coords))

        try:
            # 1) Delaunay triangulation of the points
            tri = Delaunay(points)  # points: (N,2)
            T = points[tri.simplices]  # (M,3,2)

            # 2) Estimate a length scale from triangle edges
            #    (median of all edge lengths; robust to outliers)
            e1 = np.linalg.norm(T[:, 0] - T[:, 1], axis=1)
            e2 = np.linalg.norm(T[:, 1] - T[:, 2], axis=1)
            e3 = np.linalg.norm(T[:, 2] - T[:, 0], axis=1)
            edge_lengths = np.concatenate([e1, e2, e3])
            L = np.median(edge_lengths)

            # 3) Keep triangles whose circumradius R <= R_thresh.
            #    For an equilateral triangle, R = a / sqrt(3).
            #    scale > 1 => smoother; scale < 1 => tighter.
            scale = 15
            R_thresh = (scale * L) / np.sqrt(3)

            # Collect boundary edges from "kept" triangles
            edges = []
            for A, B, C in T:
                a = np.linalg.norm(B - C)
                b = np.linalg.norm(A - C)
                c = np.linalg.norm(A - B)
                s = 0.5 * (a + b + c)
                area_sq = s * (s - a) * (s - b) * (s - c)
                if area_sq <= 0:
                    continue
                area = np.sqrt(area_sq)
                R = (a * b * c) / (4.0 * area)
                if R <= R_thresh:
                    edges.append((tuple(A), tuple(B)))
                    edges.append((tuple(B), tuple(C)))
                    edges.append((tuple(C), tuple(A)))

            if edges:
                # 4) Polygonize edges to get one or more polygons; keep the largest
                mls = MultiLineString(edges)
                polys = list(polygonize(mls))
                if polys:
                    merged = unary_union(polys)
                    if merged.geom_type == "MultiPolygon":
                        poly = max(merged.geoms, key=lambda g: g.area).buffer(0)
                    else:
                        poly = merged.buffer(0)
                else:
                    print(f"falling back to convex hull for sheet {sheet_index + 1}")
                    # fallback to convex hull if polygonization failed
                    hull = ConvexHull(points)
                    poly = Polygon(points[hull.vertices]).buffer(0)
            else:
                print(f"falling back to convex hull for sheet {sheet_index + 1}")
                # fallback if no edges passed the filter (too small R_thresh)
                hull = ConvexHull(points)
                poly = Polygon(points[hull.vertices]).buffer(0)

        except Exception as e:
            print(f"Concave outline failed for sheet {sheet_index + 1}: {e}")
            hull = ConvexHull(points)
            poly = Polygon(points[hull.vertices]).buffer(0)
        # Plot the points and the convex hull
        fig, ax = plt.subplots(figsize=(page_w/3, (page_h/9)+0.25))

        # Overlay grid lines
        for x_edge in x_grid_edges:
            ax.plot([x_edge, x_edge], [min_y, max_y], 'k--', linewidth=0.25)
        for y_edge in y_grid_edges:
            ax.plot([min_x, max_x], [y_edge, y_edge], 'k--', linewidth=0.25)

        bx, by = poly.exterior.xy
        ax.plot(bx, by, 'k-', linewidth=0.5, label='Concave boundary (α-shape)')
        for hole in poly.interiors:
            hx, hy = hole.xy
            ax.plot(hx, hy, 'k--', linewidth=0.5)

        # Overlay cells colored by ROI
        for i in range(num_rois):
            roi_mask = roi_indices == i
            if np.any(roi_mask):
                ax.scatter(x_coords[roi_mask], y_coords[roi_mask], s=0.15, color=roi_colors[i])  # Remove label to prevent clutter

                # Collect cell shapes within this ROI
                roi_polygon_numbers = polygon_numbers[roi_mask]
                # Store the data for analysis
                roi_data_list[i]['polygon_numbers'].append(roi_polygon_numbers)

        # Add labels and title
        ax.set_title(f'Wing Disc Eulerian ROI - Frame {sheet_index}',fontsize=6)
        ax.set_xlabel('X Coordinate',fontsize=3)
        ax.set_ylabel('Y Coordinate',fontsize=3)
        ax.tick_params(axis='both', which='both', labelsize=3)
        for spine in ax.spines.values():
            spine.set_linewidth(0.25)  # optional: scale spine width
        ax.xaxis.set_tick_params(width=0.8, length=0.25)
        ax.yaxis.set_tick_params(width=0.8, length=0.25)
        ax.set_aspect('equal')  # Equal scaling for x and y axes

        # Save the plot as a .png file
        plt.savefig(f'{save_dir}/{stamp_prefix}euler_roi_frame_{sheet_index}.png',dpi=300)
        plt.close()

# Calculate and plot mean percentages if enabled
if plot_mean_of_all_frames_cell_shapes:
    # Determine unique polygon types across all data
    all_data_flat = np.concatenate(all_data)
    unique_polygons = np.unique(all_data_flat)

    # Calculate mean and standard deviation for each polygon type
    mean_percentages = []
    std_percentages = []

    for polygon in unique_polygons:
        percentages = [(data == polygon).sum() / len(data) * 100 for data in all_data]
        mean_percentages.append(np.mean(percentages))
        std_percentages.append(np.std(percentages))

    # Convert to numpy arrays for indexing
    unique_polygons = np.array(unique_polygons)
    mean_percentages = np.array(mean_percentages)
    std_percentages = np.array(std_percentages)

    # Filter out polygons with mean percentage less than or equal to 0.5%
    threshold = 0.5
    mask = mean_percentages > threshold
    filtered_polygons = unique_polygons[mask]
    filtered_mean_percentages = mean_percentages[mask]
    filtered_std_percentages = std_percentages[mask]

    # After calculating 'filtered_polygons', initialize deviations array
    num_rois = len(roi_data_list)
    num_shapes = len(filtered_polygons)
    deviations = np.zeros((num_rois, num_shapes))

    # Create a dictionary for easy lookup
    global_percentages_dict = dict(zip(filtered_polygons, filtered_mean_percentages))

    # Plot the mean with standard deviation
    fig, ax = plt.subplots(figsize=(page_w, page_h/3))

    # Create a bar plot with error bars
    bars = ax.bar(filtered_polygons, filtered_mean_percentages, yerr=filtered_std_percentages,
                 color=plt.cm.tab20.colors, capsize=2)

    # Adding titles and labels
    ax.set_title('Cell shape distribution across all frames', weight="bold")
    ax.set_xlabel('Shape')
    ax.set_ylabel('Percentage of Cells')
    ax.set_ylim(0, 65)
    ax.set_xticks(filtered_polygons)
    ax.set_xticklabels([shape_names.get(i, f'{i}-gon') for i in filtered_polygons])
    ax.tick_params(axis='y')
    ax.grid(True)

    # Annotate the bars with mean and std
    for bar, mean, std in zip(bars, filtered_mean_percentages, filtered_std_percentages):
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, yval + std + 0.5, f'{mean:.2f}±{std:.2f}%', ha='center', va='bottom', fontsize=8)

    # Save the plot as a .png file
    plt.tight_layout()
    plt.savefig(f'{save_dir}/{stamp_prefix}Mean_topology_over_201_frames.png',dpi=300)
    plt.close(fig)


# Keep all observed polygon classes, but drop globally ultra-rare ones (<0.2%)
metrics = compute_topology_deviation_metrics(
    roi_data_list=roi_data_list,
    all_polygon_numbers=all_polygon_numbers,
    polygons=None,          # or np.array([3,4,5,6,7])
    alpha=0.5,              # Jeffreys prior; set to 1.0 for Laplace
    min_global_pct=0.2
)

print(f"Computed metrics for {len(metrics['support'])} polygon classes: {metrics['support']}")
# Analyze cell shapes within ROIs if enabled
if plot_roi_analysis or plot_roi_deviations:
    roi_mean_deviations = []  # List to store mean deviations per ROI

    for i in range(num_rois):
        roi_polygon_numbers_all_frames = np.concatenate(roi_data_list[i]['polygon_numbers']) if roi_data_list[i][
            'polygon_numbers'] else np.array([])
        if roi_polygon_numbers_all_frames.size == 0:
            deviations[i, :] = 0  # Set deviations to zero for empty ROIs
            continue

        # Calculate mean percentages for each polygon type within the ROI
        mean_percentages_roi = []
        for polygon in filtered_polygons:
            percentages = []
            for data in roi_data_list[i]['polygon_numbers']:
                if len(data) > 0:
                    percentages.append((data == polygon).sum() / len(data) * 100)
            mean_percentage = np.mean(percentages) if percentages else 0.0
            mean_percentages_roi.append(mean_percentage)

        roi_percentages_vector = np.array(mean_percentages_roi)
        global_percentages_vector = np.array(
            [global_percentages_dict.get(polygon, 0.0) for polygon in filtered_polygons])

        # Compute absolute deviation vector
        deviation_vector = np.abs(roi_percentages_vector - global_percentages_vector)

        # Store deviations
        deviations[i, :] = deviation_vector

# Plot deviation maps for each shape with custom colormap normalization using levels
for polygon_index, polygon in enumerate(filtered_polygons):
    shape_deviations = deviations[:, polygon_index]
    shape_deviation_grid = shape_deviations.reshape((num_grid_x, num_grid_y)).T
    X, Y = np.meshgrid(x_grid_edges, y_grid_edges)

    # Define custom levels
    levels = [0, 5, 10, 15, 20]

    # Create custom colormap normalization
    norm = colors.BoundaryNorm(boundaries=levels, ncolors=256, clip=False)

    fig, ax = plt.subplots(figsize=(page_w/2, page_h/6))
    pcm = ax.pcolormesh(X, Y, shape_deviation_grid, cmap='viridis', norm=norm, shading='auto')

    # Add color bar with extended ends (arrow-like shapes)
    cbar = fig.colorbar(pcm, ax=ax, extend='max', orientation='horizontal', ticks=levels, label=f'Absolute Deviation in percentage of {polygon}-gons')

    # Set tick labels for the intervals
    tick_labels = ['0-5', '5-10', '10-15', '15-20', '>20']
    cbar.set_ticklabels(tick_labels)

    # Plot the convex hull
    points = np.column_stack((all_x_coords[0], all_y_coords[0]))
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]

    # Save the initial contour as a NumPy array
    np.save('wing_initial_contour.npy', hull_points)

    # Or create a Shapely polygon
    initial_contour = Polygon(hull_points)

    # Optionally, save the Shapely polygon to a WKT file
    with open('wing_initial_contour.wkt', 'w') as f:
        f.write(initial_contour.wkt)

    for simplex in hull.simplices:
        ax.plot(points[simplex, 0], points[simplex, 1], 'red', linewidth=3)

    ax.set_title(
        f'Absolute Deviation of {shape_names.get(polygon, f"{polygon}-gon")} Percentage per ROI from Global Distribution')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal')

    # Save the plot
    plt.savefig(f'{save_dir}/absolute_deviation_map_{polygon}-gon.png', dpi=300)
    plt.close()


# After processing all frames, compute average cell counts per ROI
average_cell_counts = []
for roi_data in roi_data_list:
    if roi_data['cell_counts']:
        avg_count = np.mean(roi_data['cell_counts'])
    else:
        avg_count = 0
    average_cell_counts.append(avg_count)

# Reshape average cell counts into grid
average_cell_count_grid = np.array(average_cell_counts).reshape((num_grid_x, num_grid_y))
average_cell_count_grid_T = average_cell_count_grid.T  # Transpose for plotting

# Create grid edges
# Create grid edges
X, Y = np.meshgrid(x_grid_edges, y_grid_edges)

# Calculate and print the average of average cell counts across all ROIs
mean_cell_count_all_rois = np.mean(average_cell_counts)
print(f"Average number of cells per ROI across the entire wing disc: {mean_cell_count_all_rois:.2f}")
# Create figure with proper dimensions
fig, ax = plt.subplots(figsize=(page_w, page_h/3))

# Plot the average cell count grid
pcm = ax.pcolormesh(X, Y, average_cell_count_grid_T, cmap='viridis', shading='auto')

# Optionally, overlay cell count numbers on each grid cell
for i in range(num_grid_x):
    for j in range(num_grid_y):
        cell_count = average_cell_count_grid[i, j]
        if cell_count > 0:
            ax.text(
                x_grid_edges[i] + (x_grid_edges[1] - x_grid_edges[0]) / 2,
                y_grid_edges[j] + (y_grid_edges[1] - y_grid_edges[0]) / 2,
                f'{cell_count:.1f}',
                ha='center',
                va='center',
                fontsize=4,
                color='white'
            )

# Plot the convex hull
points = np.column_stack((all_x_coords[0], all_y_coords[0]))
hull = ConvexHull(points)
for simplex in hull.simplices:
    ax.plot(points[simplex, 0], points[simplex, 1], 'k-', linewidth=2)

# Adjust font sizes for labels and title
cbar = fig.colorbar(pcm, ax=ax, label='Average Number of Cells')
ax.set_title('Average Number of Cells per ROI')
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.tick_params(axis='both')
ax.set_aspect('equal')

# Save the plot with higher dpi
plt.savefig(f'{save_dir}/{stamp_prefix}average_cell_count_map.png', dpi=300)
plt.close(fig)


def _to_grid(vals, nx, ny):
    g = vals.reshape((nx, ny)).T  # note transpose to match your pcolormesh
    return np.ma.masked_invalid(g)  # mask NaNs (empty ROIs)

# Example 1: JSD map (bits, 0 = identical, ~0.1–0.2 indicates noticeable shift)
jsd_grid = _to_grid(metrics['jsd'], num_grid_x, num_grid_y)
fig, ax = plt.subplots(figsize=(page_w/2, page_h/6))
# Overlay grid lines
for x_edge in x_grid_edges:
    ax.plot([x_edge, x_edge], [min_y, max_y], 'k--', linewidth=0.25)
for y_edge in y_grid_edges:
    ax.plot([min_x, max_x], [y_edge, y_edge], 'k--', linewidth=0.25)
pcm = ax.pcolormesh(X, Y, jsd_grid, cmap='viridis', shading='auto')
cbar = fig.colorbar(pcm, ax=ax, orientation='horizontal', label='Jensen–Shannon divergence (bits)')
cbar.ax.tick_params(labelsize=3)
cbar.set_label('Jensen–Shannon divergence (bits)', fontsize=3)
ax.set_title('Topology deviation map (JSD vs global)', fontsize=6)

ax.tick_params(axis='both', which='both', labelsize=3)
ax.set_aspect('equal')
plt.savefig(f'{save_dir}/{stamp_prefix}map_topology_jsd.png', dpi=300)
plt.close(fig)

# Example 2: Wasserstein-1 map (units = "polygon sides" moved)
w1_grid = _to_grid(metrics['w1'], num_grid_x, num_grid_y)
fig, ax = plt.subplots(figsize=(page_w/2, page_h/6))
# Overlay grid lines
for x_edge in x_grid_edges:
    ax.plot([x_edge, x_edge], [min_y, max_y], 'k--', linewidth=0.25)
for y_edge in y_grid_edges:
    ax.plot([min_x, max_x], [y_edge, y_edge], 'k--', linewidth=0.25)
pcm = ax.pcolormesh(X, Y, w1_grid, cmap='viridis', shading='auto')
cbar = fig.colorbar(pcm, ax=ax, orientation='horizontal', label='Wasserstein-1 distance (sides)')
cbar.ax.tick_params(labelsize=3)
cbar.set_label('Wasserstein-1 distance (sides)', fontsize=3)
ax.set_title('Topology deviation map (adjacency-aware W1)', fontsize=6)

ax.tick_params(axis='both', which='both', labelsize=3)
ax.set_aspect('equal')
plt.savefig(f'{save_dir}/{stamp_prefix}map_topology_w1.png', dpi=300)
plt.close(fig)

# Example 3: Per-class z residuals for hexagons (diagnostic: positive = enrichment)
support = metrics['support']
if 6 in support:
    k6 = int(np.where(support == 6)[0][0])
    z6_grid = _to_grid(metrics['z'][:, k6], num_grid_x, num_grid_y)
    fig, ax = plt.subplots(figsize=(page_w/2, page_h/6))
    # Overlay grid lines
    for x_edge in x_grid_edges:
        ax.plot([x_edge, x_edge], [min_y, max_y], 'k--', linewidth=0.25)
    for y_edge in y_grid_edges:
        ax.plot([min_x, max_x], [y_edge, y_edge], 'k--', linewidth=0.25)
    # Diverging colormap around zero is appropriate; you can pick one you already use
    pcm = ax.pcolormesh(X, Y, z6_grid, cmap='coolwarm', shading='auto', vmin=-4, vmax=4)
    cbar = fig.colorbar(pcm, ax=ax, orientation='horizontal', label='Hexagon residual z')
    cbar.ax.tick_params(labelsize=3)
    cbar.set_label('Hexagon residual z', fontsize=3)
    ax.set_title('Hexagon enrichment/depletion (standardized residual)', fontsize=6)
    ax.tick_params(axis='both', which='both', labelsize=3)
    ax.set_aspect('equal')
    plt.savefig(f'{save_dir}/{stamp_prefix}map_z_hexagon.png', dpi=300)
    plt.close(fig)

# Example 4: Mean absolute topological charge (packing irregularity)
absq_grid = _to_grid(metrics['mean_abs_q'], num_grid_x, num_grid_y)
fig, ax = plt.subplots(figsize=(page_w/2, page_h/6))
# Overlay grid lines
for x_edge in x_grid_edges:
    ax.plot([x_edge, x_edge], [min_y, max_y], 'k--', linewidth=0.25)
for y_edge in y_grid_edges:
    ax.plot([min_x, max_x], [y_edge, y_edge], 'k--', linewidth=0.25)
pcm = ax.pcolormesh(X, Y, absq_grid, cmap='viridis', shading='auto')
cbar = fig.colorbar(pcm, ax=ax, orientation='horizontal', label='E[|n−6|]')
cbar.ax.tick_params(labelsize=3)
cbar.set_label('E[|n−6|]', fontsize=3)
ax.set_title('Mean absolute topological charge per ROI', fontsize=6)

ax.tick_params(axis='both', which='both', labelsize=3)
ax.set_aspect('equal')
plt.savefig(f'{save_dir}/{stamp_prefix}map_mean_abs_charge.png', dpi=300)
plt.close(fig)
# Global composition and charge
support = metrics['support']
p_global = metrics['p_global']      # length K
Eabs_global = np.sum(p_global * np.abs(support - 6))
print("Global p(n):", dict(zip(support.tolist(), np.round(p_global, 3))))
print("Global E[|n-6|]:", round(Eabs_global, 3))

# Extremes across ROIs
for name, arr in [("JSD (bits)", metrics['jsd']),
                  ("W1 (sides)", metrics['w1']),
                  ("E[|n-6|]", metrics['mean_abs_q'])]:
    print(name, "min/median/max =",
          np.nanmin(arr), np.nanmedian(arr), np.nanmax(arr))


# Optional hull from your first frame (same as you already do)
points0 = np.column_stack((all_x_coords[0], all_y_coords[0]))

plot_mean_sidedness_numbers(
    metrics=metrics,
    x_grid_edges=x_grid_edges,
    y_grid_edges=y_grid_edges,
    num_grid_x=num_grid_x,
    num_grid_y=num_grid_y,
    page_w=page_w,
    page_h=page_h,
    save_path=f'{save_dir}/{stamp_prefix}map_mean_sidedness_numbers.png',
    hull_points=points0,
    posterior=True,          # Jeffreys-smoothed
    empty_policy="nan",      # <- NEW: blanks out N==0 tiles
    text_fmt="{:.2f}",
    na_label=""              # or "—" if you want to show a mark
)
