import os
import csv
import random
import xlrd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,MultiLineString
from shapely.ops import polygonize, unary_union
from matplotlib import colors
from matplotlib.colors import Normalize
from typing import Dict, List, Tuple, Set, Optional
from scipy.spatial import Delaunay, KDTree, ConvexHull

# ====== page and IO params (kept from your script) ======
page_width_in = 455.24411 / 72.27   # ≈ 6.30 in
page_height_in = 702.78308 / 72.27  # ≈ 9.73 in
page_w = page_width_in
page_h = page_height_in
stamp_prefix = 'Results_R2_'

shape_names = {3: 'Triangle', 4: 'Quadrilateral', 5: 'Pentagon', 6: 'Hexagon', 7: 'Heptagon'}
save_dir = '../Results_wing_lagrangian'
os.makedirs(save_dir, exist_ok=True)

# ====== switches (kept/extended) ======
plot_contour = True
plot_mean_of_all_frames_cell_shapes = True
plot_roi_analysis = True
plot_roi_deviations = True

# ====== ROI parameters ======
target_cells_per_roi = 46  # Lagrangian ROI size (non-boundary cells)
random_seed = 7            # for reproducibility of seed selection if distances tie


# Initialize data storage
all_data = []
all_x_coords = []
all_y_coords = []
all_polygon_numbers = []




# -------------------------------------------------------------------------
#                             DATA I/O
# -------------------------------------------------------------------------

def load_workbook(xls_dir: str, filename: str) -> xlrd.book.Book:
    """
    Open the workbook.

    :param xls_dir: Directory containing the .xls file.
    :type  xls_dir: str
    :param filename: Workbook filename.
    :type  filename: str
    :return: Opened workbook.
    :rtype:  xlrd.book.Book
    """
    return xlrd.open_workbook(os.path.join(xls_dir, filename))

# ---------------------------------------------------------------------
# I/O helper (kept for completeness; identical name/signature)
# ---------------------------------------------------------------------
def extract_data(workbook, sheet_index):
    """
    Extract per-row (id, x, y, polygon_no) from an xlrd sheet.

    :param workbook: Opened xlrd workbook.
    :type  workbook: xlrd.book.Book
    :param sheet_index: Sheet index to read.
    :type  sheet_index: int
    :return: Arrays (x_coords, y_coords, polygon_numbers).
    :rtype:  Tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    sheet = workbook.sheet_by_index(sheet_index)
    x_coords, y_coords, polygon_numbers = [], [], []
    for row_idx in range(1, sheet.nrows):
        # keep boundary rows if you wish; caller can drop them earlier
        x = float(sheet.cell_value(row_idx, 1))
        y = float(sheet.cell_value(row_idx, 2))
        pg = int(sheet.cell_value(row_idx, 3))
        x_coords.append(x); y_coords.append(y); polygon_numbers.append(pg)
    return np.array(x_coords), np.array(y_coords), np.array(polygon_numbers)

def extract_frame(sheet) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract a frame from an xlrd sheet.

    :param sheet: xlrd sheet object.
    :type  sheet: xlrd.sheet.Sheet
    :return: Arrays (ids, x, y, polygons). ids may contain -1 for boundary.
    :rtype:  Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    ids, xs, ys, polys = [], [], [], []
    for row_idx in range(1, sheet.nrows):
        if sheet.cell_value(row_idx, 0) == -1:
            continue
        cid = int(sheet.cell_value(row_idx, 0))
        x = float(sheet.cell_value(row_idx, 1))
        y = float(sheet.cell_value(row_idx, 2))
        pg = int(sheet.cell_value(row_idx, 3))
        ids.append(cid); xs.append(x); ys.append(y); polys.append(pg)
    return np.array(ids), np.array(xs), np.array(ys), np.array(polys)


def read_all_frames(wb: xlrd.book.Book) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Read all frames (sheets 1..nsheets-1). Sheet 0 is often metadata; adapt if needed.

    :param wb: Opened workbook.
    :type  wb: xlrd.book.Book
    :return: List of (ids, x, y, polygons) for each frame (index aligned with sheet_index-1).
    :rtype:  List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    """
    frames = []
    for si in range(0, wb.nsheets):
        sheet = wb.sheet_by_index(si)
        frames.append(extract_frame(sheet))
    return frames

# -------------------------------------------------------------------------
#                  LAGRANGIAN ROI (CONSTRUCTION & TRACKING)
# -------------------------------------------------------------------------

def delaunay_adjacency(x: np.ndarray, y: np.ndarray) -> List[Set[int]]:
    """
    Build adjacency from Delaunay triangulation.

    :param x: X coordinates of points.
    :type  x: np.ndarray
    :param y: Y coordinates of points.
    :type  y: np.ndarray
    :return: List of neighbor index sets for each point (same length as x).
    :rtype:  List[Set[int]]
    """
    pts = np.column_stack([x, y])
    tri = Delaunay(pts)
    n = len(pts)
    nbrs = [set() for _ in range(n)]
    for simplex in tri.simplices:
        a, b, c = simplex
        nbrs[a].update([b, c])
        nbrs[b].update([a, c])
        nbrs[c].update([a, b])
    return nbrs


def farthest_point_seeds(pts: np.ndarray, k: int, rng: np.random.Generator) -> List[int]:
    """
    Farthest-point sampling for k seeds.

    :param pts: (N,2) array of coordinates.
    :type  pts: np.ndarray
    :param k: Number of seeds to generate.
    :type  k: int
    :param rng: NumPy Generator for reproducibility.
    :type  rng: np.random.Generator
    :return: Indices of seed points.
    :rtype:  List[int]
    """
    n = pts.shape[0]
    if k >= n:
        return list(range(n))
    seeds = [rng.integers(0, n)]
    d2 = np.full(n, np.inf)
    for _ in range(1, k):
        last = seeds[-1]
        diff = pts - pts[last]
        d2 = np.minimum(d2, np.einsum('ij,ij->i', diff, diff))
        d2[seeds] = -1.0
        nxt = int(np.argmax(d2))
        seeds.append(nxt)
    return seeds

def compute_first_seen(frames):
    """
    Compute the first frame in which each cell id (> -1) appears.

    :param frames: List of (ids, x, y, polygons) per frame.
    :type  frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    :return: Mapping cell_id -> first_seen_frame_index (1-based).
    :rtype:  Dict[int, int]
    """
    first_seen = {}
    for t, (ids, _, _, _) in enumerate(frames, start=1):
        for cid in ids:
            cid = int(cid)
            if cid >= 0 and cid not in first_seen:
                first_seen[cid] = t
    return first_seen

def build_lagrangian_rois(
    ids0: np.ndarray,
    x0: np.ndarray,
    y0: np.ndarray,
    target: int = 50,
    *,
    # circularity controls (defaults match what worked for you later in time)
    circ_k_outlier: float = 1.0,          # MAD threshold (stricter = smaller)
    circ_cv_max: float = 0.35,            # max allowed std(r)/mean(r)
    circ_compactness_min: float = 0.65,   # min allowed 4πA/P^2 of ROI hull
    circ_min_members_eval: int = 20,      # skip very small ROIs in repair
    circ_min_members_keep: Optional[int] = None,  # don't shrink below this
    circ_max_iter: int = 4,               # max repair iterations
) -> Tuple[Dict[int, int], Dict[int, np.ndarray]]:
    """
    Construct Lagrangian ROIs at Frame 0 by adjacency growth to ~target cells;
    then improve circularity by reassigning peripheral outliers to neighbouring ROIs.

    Boundary rows (id=-1) are ignored for growth/repair but can be attached later.

    :param ids0: Cell IDs at Frame 0 (may include -1).
    :type  ids0: np.ndarray
    :param x0: X coordinates at Frame 0.
    :type  x0: np.ndarray
    :param y0: Y coordinates at Frame 0.
    :type  y0: np.ndarray
    :param target: Target number of non-boundary cells per ROI.
    :type  target: int
    :param circ_k_outlier: Robust radial outlier threshold multiplier (MAD).
    :type  circ_k_outlier: float
    :param circ_cv_max: Max allowed coefficient of variation of radii.
    :type  circ_cv_max: float
    :param circ_compactness_min: Min allowed hull compactness 4πA/P² (∈(0,1]).
    :type  circ_compactness_min: float
    :param circ_min_members_eval: Minimum members to attempt circularity repair.
    :type  circ_min_members_eval: int
    :param circ_min_members_keep: Do not shrink an ROI below this many members
                                  during repair; default = max(10, floor(0.5*target)).
    :type  circ_min_members_keep: Optional[int]
    :param circ_max_iter: Maximum circularity repair iterations.
    :type  circ_max_iter: int
    :return: (cell_id_to_roi, roi_seeds_xy) where roi_seeds_xy holds post-repair
             ROI centers (centroids) for robust seeding/fallback in later frames.
    :rtype:  Tuple[Dict[int, int], Dict[int, np.ndarray]]
    """
    rng = np.random.default_rng(random_seed)
    if circ_min_members_keep is None:
        circ_min_members_keep = max(10, int(np.floor(0.5 * target)))

    # ------------------------- seed & grow (your code) -------------------------
    mask_nb = ids0 != -1
    idx_nb = np.where(mask_nb)[0]
    pts_nb = np.column_stack([x0[mask_nb], y0[mask_nb]])
    n_nb = len(idx_nb)
    if n_nb == 0:
        raise ValueError("No non-boundary cells in Frame 0; cannot seed ROIs.")

    n_roi = max(1, n_nb // target + (1 if (n_nb % target) >= target // 2 else 0))

    nbrs_full = delaunay_adjacency(x0[mask_nb], y0[mask_nb])

    seed_local_indices = farthest_point_seeds(pts_nb, n_roi, rng)
    seed_global_indices = idx_nb[seed_local_indices]

    assigned = np.full(len(ids0), fill_value=-1, dtype=int)
    quota = np.full(n_roi, target, dtype=int)
    queues: List[List[int]] = [[s] for s in seed_local_indices]

    for r, s_loc in enumerate(seed_local_indices):
        g = idx_nb[s_loc]
        assigned[g] = r
        quota[r] = max(0, quota[r] - 1)

    frontier_levels = 0
    while (assigned[idx_nb] == -1).any() and (quota > 0).any():
        frontier_levels += 1
        next_queues: List[List[int]] = [[] for _ in range(n_roi)]
        for r in range(n_roi):
            if quota[r] <= 0 or not queues[r]:
                continue
            for u_loc in queues[r]:
                for v_loc in nbrs_full[u_loc]:
                    g = idx_nb[v_loc]
                    if assigned[g] == -1:
                        assigned[g] = r
                        quota[r] -= 1
                        if quota[r] > 0:
                            next_queues[r].append(v_loc)
                if quota[r] <= 0:
                    break
        queues = next_queues
        if frontier_levels > 5 * n_roi:
            break

    remaining_nb = idx_nb[assigned[idx_nb] == -1]
    if len(remaining_nb) > 0:
        seed_pts = np.column_stack([x0[seed_global_indices], y0[seed_global_indices]])
        kd = KDTree(seed_pts)
        rem_pts = np.column_stack([x0[remaining_nb], y0[remaining_nb]])
        _, nn = kd.query(rem_pts, k=1)
        for g, r in zip(remaining_nb, nn):
            assigned[g] = int(r)

    idx_bnd = np.where(ids0 == -1)[0]
    if len(idx_bnd) > 0:
        seed_pts = np.column_stack([x0[seed_global_indices], y0[seed_global_indices]])
        kd = KDTree(seed_pts)
        bnd_pts = np.column_stack([x0[idx_bnd], y0[idx_bnd]])
        _, nn = kd.query(bnd_pts, k=1)
        for g, r in zip(idx_bnd, nn):
            assigned[g] = int(r)

    # -------------------- circularity repair at t0 (new) ----------------------
    def _roi_members():
        """Return dict roi -> np.array of NON-boundary frame-0 indices."""
        out = {r: [] for r in range(n_roi)}
        nb_idx = np.where(ids0 != -1)[0]
        for g in nb_idx:
            r = int(assigned[g])
            if r >= 0:
                out[r].append(int(g))
        return {r: np.array(v, dtype=int) for r, v in out.items()}

    def _centroids(members: Dict[int, np.ndarray]) -> np.ndarray:
        C = np.full((n_roi, 2), np.nan, dtype=float)
        for r, idx in members.items():
            if idx.size > 0:
                C[r, 0] = float(np.mean(x0[idx]))
                C[r, 1] = float(np.mean(y0[idx]))
        return C

    def _compactness(idx: np.ndarray) -> float:
        """4πA/P² of convex hull; 1 is a perfect circle (upper bound)."""
        if idx.size < 3:
            return 0.0
        pts = np.column_stack([x0[idx], y0[idx]])
        try:
            hull = ConvexHull(pts)
            hv = pts[hull.vertices]
        except Exception:
            hv = pts
        # polygon area & perimeter
        x, y = hv[:, 0], hv[:, 1]
        area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
        per = np.sum(np.hypot(np.diff(x, append=x[0]), np.diff(y, append=y[0])))
        if per <= 1e-12:
            return 0.0
        return float(4.0 * np.pi * area / (per * per))

    def _radial_cv(idx: np.ndarray, cx: float, cy: float) -> float:
        if idx.size < 3:
            return 1.0
        r = np.hypot(x0[idx] - cx, y0[idx] - cy)
        mu = np.mean(r)
        sd = np.std(r)
        return float(sd / max(mu, 1e-12))

    members = _roi_members()
    C = _centroids(members)

    for it in range(circ_max_iter):
        moved = 0
        # KD over current centroids (valid only)
        valid = ~np.isnan(C).any(axis=1)
        valid_idx = np.where(valid)[0]
        if not np.any(valid):
            break
        kdC = KDTree(C[valid])

        for r in range(n_roi):
            idx = members[r]
            if idx.size < circ_min_members_eval or not valid[r]:
                continue

            cx, cy = C[r]
            d = np.hypot(x0[idx] - cx, y0[idx] - cy)

            # robust radial threshold: median + k * MAD (k=1.0 default)
            med = np.median(d)
            mad = 1.4826 * np.median(np.abs(d - med)) if d.size > 0 else 0.0
            scale = mad if mad > 1e-12 else (np.std(d) if d.size > 1 else 0.0)
            thr = med + circ_k_outlier * scale

            # circularity scores BEFORE
            cv_before = _radial_cv(idx, cx, cy)
            comp_before = _compactness(idx)
            needs_fix = (cv_before > circ_cv_max) or (comp_before < circ_compactness_min)

            if not needs_fix:
                continue

            # candidate outliers (sorted farthest first to make big gains early)
            out_mask = d > thr
            if not np.any(out_mask):
                continue
            out_rows = idx[out_mask]
            order = np.argsort(d[out_mask])[::-1]
            out_rows = out_rows[order]

            for g in out_rows:
                # do not shrink below keep-threshold
                if members[r].size <= circ_min_members_keep:
                    break

                # choose nearest OTHER centroid
                dist, nn = kdC.query([x0[g], y0[g]], k=min(2, len(valid_idx)))
                # if nearest is itself and we have another option, pick 2nd
                pick = int(nn[1]) if (isinstance(nn, np.ndarray) and valid_idx[int(nn[0])] == r and len(valid_idx) > 1) \
                       else (int(nn) if not isinstance(nn, np.ndarray) else int(nn[0]))
                to_roi = int(valid_idx[pick])
                if to_roi == r:
                    # nothing to do
                    continue

                # move g: update assigned + members
                assigned[g] = to_roi
                members[r] = members[r][members[r] != g]
                members[to_roi] = np.append(members[to_roi], g)
                moved += 1

            if moved > 0:
                # refresh centroid of r (and target to_roi will refresh in outer recompute)
                pass

        if moved == 0:
            break
        # recompute centroids after this sweep
        C = _centroids(members)

    # After repair, set seeds to (post-repair) centroids for more stable later assignment
    roi_seeds_xy = {}
    for r in range(n_roi):
        if np.isnan(C[r]).any():
            # fallback to original seed
            sg = seed_global_indices[r]
            roi_seeds_xy[r] = np.array([x0[sg], y0[sg]])
        else:
            roi_seeds_xy[r] = C[r].copy()

    # ------------------------- finalize mapping -------------------------
    cell_id_to_roi: Dict[int, int] = {}
    for g, roi in enumerate(assigned):
        cid = int(ids0[g])
        cell_id_to_roi[cid] = int(roi)

    return cell_id_to_roi, roi_seeds_xy


def debug_print_lagrangian_roi(
    frames,
    roi_index_per_frame,
    roi_colors,
    roi_id,
    first_seen,
    csv_path: str | None = None,
):
    """
    Print per-frame diagnostics for one Lagrangian ROI:
    center (mean x,y), ROI color, number of cells, and age stats.

    :param frames: List of (ids, x, y, polygons) per frame.
    :type  frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    :param roi_index_per_frame: Per-frame vector of ROI indices (aligned to 'frames').
    :type  roi_index_per_frame: List[np.ndarray]
    :param roi_colors: RGBA color list (index = ROI id).
    :type  roi_colors: List[Tuple[float, float, float, float]]
    :param roi_id: ROI index to debug.
    :type  roi_id: int
    :param first_seen: Map cell_id -> first-seen frame (1-based).
    :type  first_seen: Dict[int, int]
    :param csv_path: Optional output CSV to log the same info.
    :type  csv_path: Optional[str]
    :return: None.
    :rtype:  None
    """
    color = roi_colors[roi_id]
    color_str = "(" + ", ".join(f"{c:.3f}" for c in color) + ")"

    # baseline center at frame 1 for drift
    cx0 = cy0 = None

    print(f"\n=== DEBUG ROI {roi_id} ===")
    print(f"Color RGBA: {color_str}")

    # Convert color_str (e\.g\., "(0\.500, 0\.200, 0\.700, 1\.000)") to RGBA tuple
    rgba = tuple(float(c) for c in color_str.strip("()").split(","))
    r, g, b = [int(255 * c) for c in rgba[:3]]

    # ANSI escape code for truecolor foreground
    print(f"\033[38;2;{r};{g};{b}mROI color here\033[0m")

    writer = None
    fh = None
    if csv_path is not None:
        fh = open(csv_path, "w", newline="")
        writer = csv.writer(fh)
        writer.writerow(["frame", "roi_id", "cx", "cy", "drift_from_f1",
                         "n_cells", "n_aged", "age_min", "age_med", "age_max"])

    for t, (ids, xs, ys, _) in enumerate(frames, start=1):
        roi_idx = roi_index_per_frame[t-1]
        mask = (roi_idx == roi_id)
        n = int(np.sum(mask))
        if n == 0:
            print(f"Frame {t:3d}: EMPTY")
            if writer:
                writer.writerow([t, roi_id, "", "", "", 0, 0, "", "", ""])
            continue

        cx = float(np.mean(xs[mask]))
        cy = float(np.mean(ys[mask]))

        if cx0 is None:
            cx0, cy0 = cx, cy
        drift = float(np.hypot(cx - cx0, cy - cy0))

        # ages in frames since first sighting
        ages = []
        for cid in ids[mask]:
            cid = int(cid)
            if cid >= 0:
                ages.append(t - first_seen.get(cid, t))
            else:
                # boundary entries: skip from age stats
                pass

        if len(ages) > 0:
            arr = np.array(ages, dtype=int)
            age_min = int(np.min(arr))
            age_med = float(np.median(arr))
            age_max = int(np.max(arr))
            age_txt = f"ages n={len(arr)} min/med/max={age_min}/{age_med:.1f}/{age_max}"
            n_aged = len(arr)
        else:
            age_min = age_med = age_max = ""
            age_txt = "ages n=0"
            n_aged = 0

        print(f"Frame {t:3d}: center=({cx:.2f}, {cy:.2f})  drift={drift:.2f}  "
              f"n_cells={n:4d}  {age_txt}")

        if writer:
            writer.writerow([t, roi_id, f"{cx:.6f}", f"{cy:.6f}", f"{drift:.6f}",
                             n, n_aged, age_min, f"{age_med}", age_max])

    if fh is not None:
        fh.close()
        print(f"[wrote] {csv_path}")

def track_and_collect(
    frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    id2roi_0: Dict[int, int],
    roi_seeds_xy: Dict[int, np.ndarray]
) -> Tuple[Dict[int, Dict[str, List[np.ndarray]]], np.ndarray, np.ndarray, np.ndarray, List[np.ndarray]]:
    """
    Track Lagrangian ROIs; append new cells; collect polygons; and
    return per-frame ROI index arrays aligned with each frame's rows.

    :param frames: List of (ids, x, y, polygons) per frame.
    :type  frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    :param id2roi_0: Mapping from cell id to ROI at Frame 0 (includes -1).
    :type  id2roi_0: Dict[int, int]
    :param roi_seeds_xy: ROI seed coordinates (Frame 0) as fallback for orphan cells.
    :type  roi_seeds_xy: Dict[int, np.ndarray]
    :return: (roi_data_list, x0, y0, p0, roi_index_per_frame)
    :rtype:  Tuple[Dict[int, Dict[str, List[np.ndarray]]], np.ndarray, np.ndarray, np.ndarray, List[np.ndarray]]
    """

    def _roi_centroids_current(xs: np.ndarray, ys: np.ndarray, roi_idx_vec: np.ndarray, n_roi: int) -> Tuple[
        np.ndarray, np.ndarray]:
        """
        Compute per-ROI centroids for the current frame.

        :param xs: X coordinates (N,).
        :type  xs: np.ndarray
        :param ys: Y coordinates (N,).
        :type  ys: np.ndarray
        :param roi_idx_vec: ROI index per row (N,).
        :type  roi_idx_vec: np.ndarray
        :param n_roi: Number of ROIs.
        :type  n_roi: int
        :return: (centroids [n_roi,2], counts [n_roi,])
        :rtype:  Tuple[np.ndarray, np.ndarray]
        """
        C = np.full((n_roi, 2), np.nan, dtype=float)
        counts = np.zeros(n_roi, dtype=int)
        for r in range(n_roi):
            m = roi_idx_vec == r
            counts[r] = int(np.sum(m))
            if counts[r] > 0:
                C[r, 0] = float(np.mean(xs[m]))
                C[r, 1] = float(np.mean(ys[m]))
        return C, counts

    def round_up_clusters_once(
            ids: np.ndarray,
            xs: np.ndarray,
            ys: np.ndarray,
            roi_idx_vec: np.ndarray,
            id2roi: Dict[int, int],
            *,
            k_outlier: float = 2.5,
            min_members_for_roundup: int = 15
    ) -> Tuple[np.ndarray, Dict[int, int], Dict[str, object]]:
        """
        One "round-up" pass: make clusters more circular by ejecting radial outliers
        (beyond a robust threshold) to the nearest other ROI centroid.

        :param ids: Cell ids for this frame (N,).
        :type  ids: np.ndarray
        :param xs: X coordinates (N,).
        :type  xs: np.ndarray
        :param ys: Y coordinates (N,).
        :type  ys: np.ndarray
        :param roi_idx_vec: Current ROI index per row; will be updated.
        :type  roi_idx_vec: np.ndarray
        :param id2roi: Evolving id->roi map; will be updated for non-boundary ids.
        :type  id2roi: Dict[int, int]
        :param k_outlier: Threshold multiplier on MAD (robust outlier detection).
        :type  k_outlier: float
        :param min_members_for_roundup: Skip tiny clusters below this size.
        :type  min_members_for_roundup: int
        :return: (updated_roi_idx_vec, updated_id2roi, debug dict)
        :rtype:  Tuple[np.ndarray, Dict[int, int], Dict[str, object]]
        """
        n_roi = int(np.max(roi_idx_vec)) + 1
        C, counts = _roi_centroids_current(xs, ys, roi_idx_vec, n_roi)

        # KDTree over centroids (only where defined)
        valid = ~np.isnan(C).any(axis=1)
        valid_idx = np.where(valid)[0]
        kdC = KDTree(C[valid]) if np.any(valid) else None

        moved = 0
        moved_records: List[Tuple[int, int, int, float, float]] = []
        # (row_i, from_roi, to_roi, d_from, d_to)

        # For each ROI, compute robust radius and eject outliers
        for r in range(n_roi):
            if not valid[r]:
                continue
            mask = roi_idx_vec == r
            if np.sum(mask) < min_members_for_roundup:
                continue
            cx, cy = C[r]
            dx = xs[mask] - cx
            dy = ys[mask] - cy
            d = np.hypot(dx, dy)

            med = np.median(d)
            mad = 1.4826 * np.median(np.abs(d - med)) if d.size > 0 else 0.0
            thr = med + k_outlier * (mad if mad > 1e-12 else (np.std(d) if d.size > 1 else 0.0))

            # indices (global rows) that are outliers for ROI r
            roi_rows = np.where(mask)[0]
            out_rows = roi_rows[d > thr]

            if kdC is None or len(out_rows) == 0:
                continue

            # For each outlier: send to nearest OTHER centroid
            for i in out_rows:
                if ids[i] == -1:
                    # boundary: only change the per-frame view
                    _, idx = kdC.query([xs[i], ys[i]], k=1)
                    to_roi = int(valid_idx[int(idx)])
                    if to_roi == r and len(valid_idx) > 1:
                        # if nearest is itself and we have options, take the second nearest
                        _, idx2 = kdC.query([xs[i], ys[i]], k=2)
                        to_roi = int(valid_idx[int(idx2[1])])
                    if to_roi != r:
                        d_from = float(np.hypot(xs[i] - C[r, 0], ys[i] - C[r, 1]))
                        d_to = float(np.hypot(xs[i] - C[to_roi, 0], ys[i] - C[to_roi, 1]))
                        roi_idx_vec[i] = to_roi
                        moved += 1
                        moved_records.append((int(i), int(r), int(to_roi), d_from, d_to))
                else:
                    # non-boundary: also update global id2roi
                    _, idx = kdC.query([xs[i], ys[i]], k=1)
                    to_roi = int(valid_idx[int(idx)])
                    if to_roi == r and len(valid_idx) > 1:
                        _, idx2 = kdC.query([xs[i], ys[i]], k=2)
                        to_roi = int(valid_idx[int(idx2[1])])
                    if to_roi != r:
                        d_from = float(np.hypot(xs[i] - C[r, 0], ys[i] - C[r, 1]))
                        d_to = float(np.hypot(xs[i] - C[to_roi, 0], ys[i] - C[to_roi, 1]))
                        roi_idx_vec[i] = to_roi
                        id2roi[int(ids[i])] = to_roi
                        moved += 1
                        moved_records.append((int(i), int(r), int(to_roi), d_from, d_to))

        dbg = {
            "moved": moved,
            "moved_records": moved_records,
            "centroids": C,
            "counts_before": counts,
        }
        return roi_idx_vec, id2roi, dbg

    def _centroids_from_prev(
            prev_ids: np.ndarray,
            prev_x: np.ndarray,
            prev_y: np.ndarray,
            prev_roi_idx: np.ndarray,
            n_roi: int
    ) -> Optional[np.ndarray]:
        """
        Compute per-ROI centroids from the previous frame.

        :param prev_ids: IDs of previous frame (not used except for alignment).
        :type  prev_ids: np.ndarray
        :param prev_x: X of previous frame.
        :type  prev_x: np.ndarray
        :param prev_y: Y of previous frame.
        :type  prev_y: np.ndarray
        :param prev_roi_idx: ROI index per row (previous frame).
        :type  prev_roi_idx: np.ndarray
        :param n_roi: Number of ROIs.
        :type  n_roi: int
        :return: Array (n_roi, 2) of centroids or None if any centroid cannot be computed.
        :rtype:  Optional[np.ndarray]
        """
        C = np.full((n_roi, 2), np.nan, dtype=float)
        for r in range(n_roi):
            mask = prev_roi_idx == r
            if np.any(mask):
                C[r, 0] = float(np.mean(prev_x[mask]))
                C[r, 1] = float(np.mean(prev_y[mask]))
        if np.isnan(C).any():
            return None
        return C

    def _assign_roi_for_frame(
            ids: np.ndarray,
            xs: np.ndarray,
            ys: np.ndarray,
            id2roi: Dict[int, int],
            kd_seeds: 'KDTree',
            max_id0: int,
            *,
            frame_index: int,
            debug: bool = True,
            debug_sample_n: int = 12,
            assign_strategy: str = "nearest_member",
            prev_frame: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = None,
            prev_roi_idx: Optional[np.ndarray] = None,
            n_roi: Optional[int] = None,
    ) -> Tuple[np.ndarray, Dict[int, int], Dict[str, object]]:
        """
        Assign ROI to every row in this frame with a selectable strategy.

        :param assign_strategy: One of {"nearest_member", "nearest_seed",
                                 "nearest_centroid_prev", "hybrid_member_then_centroid"}.
        :type  assign_strategy: str
        :param prev_frame: Previous frame (ids, x, y, polys) if using centroid-based strategies.
        :type  prev_frame: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
        :param prev_roi_idx: ROI index vector for previous frame.
        :type  prev_roi_idx: Optional[np.ndarray]
        :param n_roi: Number of ROIs (required for centroid-based strategies).
        :type  n_roi: Optional[int]
        """
        N = len(ids)
        roi_idx_vec = np.full(N, -1, dtype=int)

        # Known (already mapped) non-boundary ids present in this frame
        known_mask = np.array([(int(cid) in id2roi) and (int(cid) != -1) for cid in ids], dtype=bool)
        known_idx = np.where(known_mask)[0]
        kd_known = KDTree(np.column_stack([xs[known_idx], ys[known_idx]])) if len(known_idx) > 0 else None

        # Optional centroids from previous frame
        kd_centroids = None
        if assign_strategy in ("nearest_centroid_prev", "hybrid_member_then_centroid"):
            if prev_frame is not None and prev_roi_idx is not None and n_roi is not None:
                _, px, py, _ = prev_frame
                C = _centroids_from_prev(prev_frame[0], px, py, prev_roi_idx, n_roi)
                if C is not None:
                    kd_centroids = KDTree(C)

        # Pass A: copy known ids (keep their ROI)
        for i in known_idx:
            roi_idx_vec[i] = int(id2roi[int(ids[i])])

        # Helper: decide ROI for a single unseen, non-boundary row i
        def assign_unseen(i: int) -> int:
            """
            Decide ROI for row i using the selected strategy.

            :param i: Row index in current frame.
            :type  i: int
            :return: ROI index.
            :rtype: int
            """
            if assign_strategy == "nearest_member":
                if kd_known is not None and len(known_idx) > 0:
                    _, nn = kd_known.query([xs[i], ys[i]], k=1)
                    j = known_idx[int(nn)]
                    return int(id2roi[int(ids[j])])
                _, seed_i = kd_seeds.query([xs[i], ys[i]], k=1)
                return int(seed_i)

            elif assign_strategy == "nearest_seed":
                _, seed_i = kd_seeds.query([xs[i], ys[i]], k=1)
                return int(seed_i)

            elif assign_strategy == "nearest_centroid_prev":
                if kd_centroids is not None:
                    _, r = kd_centroids.query([xs[i], ys[i]], k=1)
                    return int(r)
                _, seed_i = kd_seeds.query([xs[i], ys[i]], k=1)
                return int(seed_i)

            elif assign_strategy == "hybrid_member_then_centroid":
                if kd_known is not None and len(known_idx) > 0:
                    _, nn = kd_known.query([xs[i], ys[i]], k=1)
                    j = known_idx[int(nn)]
                    return int(id2roi[int(ids[j])])
                if kd_centroids is not None:
                    _, r = kd_centroids.query([xs[i], ys[i]], k=1)
                    return int(r)
                _, seed_i = kd_seeds.query([xs[i], ys[i]], k=1)
                return int(seed_i)

            else:
                raise ValueError(f"Unknown assign_strategy='{assign_strategy}'")

        # Pass B: unseen/non-boundary ids
        for i, cid_raw in enumerate(ids):
            cid = int(cid_raw)
            if cid == -1 or roi_idx_vec[i] != -1:
                continue
            if (cid not in id2roi) or (cid > max_id0):
                roi = assign_unseen(i)
                id2roi[cid] = roi
                roi_idx_vec[i] = roi

        # Pass C: boundary rows → nearest seed (per row; no id2roi mutation)
        for i, cid_raw in enumerate(ids):
            if int(cid_raw) == -1:
                _, seed_i = kd_seeds.query([xs[i], ys[i]], k=1)
                roi_idx_vec[i] = int(seed_i)

        # Pass D: late orphans (seed fallback)
        for i in range(N):
            if roi_idx_vec[i] < 0:
                _, seed_i = kd_seeds.query([xs[i], ys[i]], k=1)
                roi_idx_vec[i] = int(seed_i)

        # (Optional: your existing debug prints/metrics can stay here)
        debug_info = {"frame": frame_index, "assign_strategy": assign_strategy}
        return roi_idx_vec, id2roi, debug_info

    roi_data: Dict[int, Dict[str, List[np.ndarray]]] = {}
    n_roi = len(roi_seeds_xy)
    for r in range(n_roi):
        roi_data[r] = {"polygon_numbers": [], "cell_counts": []}

    ids0, x0, y0, p0 = frames[0]
    max_id0 = np.max(ids0[ids0 >= 0]) if np.any(ids0 >= 0) else -1

    # evolving id->roi map
    id2roi = dict(id2roi_0)

    # fallback KD on seeds
    seed_mat = np.vstack([roi_seeds_xy[r] for r in range(n_roi)])
    kd_fallback = KDTree(seed_mat)

    roi_index_per_frame: List[np.ndarray] = []

    for t, (ids, xs, ys, polys) in enumerate(frames):
        # KD on already-assigned cells present in this frame
        known_mask = np.array([cid in id2roi and cid != -1 for cid in ids], dtype=bool)
        known_idx = np.where(known_mask)[0]
        kd_current = KDTree(np.column_stack([xs[known_idx], ys[known_idx]])) if len(known_idx) > 0 else None
        prev_fr = frames[t - 2] if t > 1 else None
        prev_roi = roi_index_per_frame[t - 2] if t > 1 else None
        roi_idx_vec, id2roi, dbg = _assign_roi_for_frame(
            ids=ids, xs=xs, ys=ys, id2roi=id2roi, kd_seeds=kd_fallback, max_id0=max_id0,
            frame_index=t, debug=True, debug_sample_n=12,
            assign_strategy="nearest_member",
            # or "nearest_seed" / "nearest_centroid_prev" / "hybrid_member_then_centroid"
            prev_frame=prev_fr, prev_roi_idx=prev_roi, n_roi=len(roi_seeds_xy)
        )
        roi_index_per_frame.append(roi_idx_vec)
        # ---- end assignment ----
        n_orphans = int(np.sum(roi_idx_vec < 0))
        if n_orphans != 0:
            raise RuntimeError(f"Found {n_orphans} unassigned rows after robust pass — should never happen.")
        # Periodic round-up: every 10th frame (t is 1-based)
        if (t % 3) == 0:
            roi_idx_vec, id2roi, dbg_round = round_up_clusters_once(
                ids=ids, xs=xs, ys=ys,
                roi_idx_vec=roi_idx_vec,
                id2roi=id2roi,
                k_outlier=1,  # tune if you want tighter/looser circles
                min_members_for_roundup=15  # skip tiny clusters to avoid instability
            )
            if dbg_round["moved"] > 0:
                print(f"[round-up] Frame {t}: moved {dbg_round['moved']} members out of ROIs")

        # collect polygons per ROI for this frame
        roi_to_polys = {r: [] for r in range(n_roi)}
        for i, cid in enumerate(ids):
            roi = roi_idx_vec[i]
            if roi is not None and roi >= 0:
                roi_to_polys[int(roi)].append(int(polys[i]))
        for r in range(n_roi):
            arr = np.array(roi_to_polys[r], dtype=int)
            roi_data[r]["polygon_numbers"].append(arr)
            roi_data[r]["cell_counts"].append(len(arr))

    return roi_data, x0, y0, p0, roi_index_per_frame


# -------------------------------------------------------------------------
#                         ANALYSIS & PLOTTING
# -------------------------------------------------------------------------

def compute_global_distribution(all_polys: List[np.ndarray], thresh: float = 0.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute global mean/std polygon distribution across frames and filter rare shapes.

    :param all_polys: List of per-frame polygon arrays.
    :type  all_polys: List[np.ndarray]
    :param thresh: Minimum mean percentage to keep a polygon type.
    :type  thresh: float
    :return: (polygons, mean_pct, std_pct) filtered.
    :rtype:  Tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    unique = np.unique(np.concatenate(all_polys))
    mean_pct, std_pct = [], []
    for pg in unique:
        percentages = [(arr == pg).sum() / len(arr) * 100 if len(arr) else 0.0 for arr in all_polys]
        mean_pct.append(np.mean(percentages))
        std_pct.append(np.std(percentages))
    unique = np.array(unique, dtype=int)
    mean_pct = np.array(mean_pct)
    std_pct = np.array(std_pct)
    mask = mean_pct > thresh
    return unique[mask], mean_pct[mask], std_pct[mask]


def plot_global_distribution(filtered_polys, mean_pct, std_pct, outpath: str) -> None:
    """
    Plot global polygon distribution with error bars.

    :param filtered_polys: Polygon categories kept.
    :type  filtered_polys: np.ndarray
    :param mean_pct: Mean percentages.
    :type  mean_pct: np.ndarray
    :param std_pct: Std deviations.
    :type  std_pct: np.ndarray
    :param outpath: Output filepath (PNG).
    :type  outpath: str
    :return: None.
    :rtype:  None
    """
    fig, ax = plt.subplots(figsize=(page_w, page_h/3))
    bars = ax.bar(filtered_polys, mean_pct, yerr=std_pct, color=plt.cm.tab20.colors, capsize=2)
    ax.set_title('Cell shape distribution across all frames', weight="bold")
    ax.set_xlabel('Shape'); ax.set_ylabel('Percentage of Cells')
    ax.set_ylim(0, 65); ax.set_xticks(filtered_polys)
    ax.set_xticklabels([shape_names.get(i, f'{i}-gon') for i in filtered_polys])
    ax.grid(True)
    for bar, m, s in zip(bars, mean_pct, std_pct):
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, yval + s + 0.5, f'{m:.2f}±{s:.2f}%', ha='center', va='bottom', fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300); plt.close(fig)


def compute_roi_deviations(
    roi_data: Dict[int, Dict[str, List[np.ndarray]]],
    filtered_polys: np.ndarray,
    global_mean: np.ndarray
) -> np.ndarray:
    """
    Absolute deviations of each ROI's mean polygon distribution from global.

    :param roi_data: Per-ROI collected polygons over time.
    :type  roi_data: Dict[int, Dict[str, List[np.ndarray]]]
    :param filtered_polys: Polygon categories used.
    :type  filtered_polys: np.ndarray
    :param global_mean: Global mean percentages aligned with filtered_polys.
    :type  global_mean: np.ndarray
    :return: deviations[roi, shape_index] absolute percentage differences.
    :rtype:  np.ndarray
    """
    n_roi = len(roi_data)
    n_shapes = len(filtered_polys)
    deviations = np.zeros((n_roi, n_shapes), dtype=float)
    global_dict = {int(pg): float(mu) for pg, mu in zip(filtered_polys, global_mean)}

    for r in range(n_roi):
        # mean percent within ROI over frames
        means = []
        for pg in filtered_polys:
            percs = []
            for arr in roi_data[r]["polygon_numbers"]:
                if arr.size > 0:
                    percs.append((arr == pg).sum() / len(arr) * 100)
            means.append(np.mean(percs) if percs else 0.0)
        roi_vec = np.array(means, dtype=float)
        glob_vec = np.array([global_dict[int(pg)] for pg in filtered_polys], dtype=float)
        deviations[r, :] = np.abs(roi_vec - glob_vec)
    return deviations


def plot_roi_deviation_heatmaps(
    deviations: np.ndarray,
    filtered_polys: np.ndarray,
    x0: np.ndarray,
    y0: np.ndarray,
    outdir: str
) -> None:
    """
    Make one heatmap per polygon type across ROIs (Lagrangian set). For visualization,
    we embed ROIs back on Frame 0 space by coloring each ROI’s convex hull.

    :param deviations: Array [n_roi, n_shapes] of absolute percentage differences.
    :type  deviations: np.ndarray
    :param filtered_polys: Polygon categories (aligned with deviations columns).
    :type  filtered_polys: np.ndarray
    :param x0: Frame 0 x-coordinates (for overall hull overlay).
    :type  x0: np.ndarray
    :param y0: Frame 0 y-coordinates (for overall hull overlay).
    :type  y0: np.ndarray
    :param outdir: Output directory.
    :type  outdir: str
    :return: None.
    :rtype:  None
    """
    # color levels as in your script
    levels = [0, 5, 10, 15, 20]
    norm = colors.BoundaryNorm(boundaries=levels, ncolors=256, clip=False)

    # Overall wing hull for reference
    pts = np.column_stack((x0, y0))
    hull = ConvexHull(pts)

    # For spatial depiction we need a per-ROI representative polygon at Frame 0.
    # We approximate each ROI’s spatial extent by the convex hull of its Frame 0 members (computed earlier).
    # Since we didn’t carry explicit per-ROI membership here, this plotting helper expects a precomputed
    # file ‘roi_members_frame0.npy’ produced below (roi -> array of point indices). This keeps interfaces clean.

    if not os.path.exists(os.path.join(outdir, "roi_members_frame0.npy")):
        # Skip spatial hull-based plots if membership file missing.
        return

    roi_members = np.load(os.path.join(outdir, "roi_members_frame0.npy"), allow_pickle=True).item()

    for j, pg in enumerate(filtered_polys):
        fig, ax = plt.subplots(figsize=(page_w, page_h/3))
        # paint each ROI hull with a single color (deviation magnitude)
        for r, members in roi_members.items():
            if len(members) < 3:
                continue
            poly_pts = np.column_stack([x0[members], y0[members]])
            try:
                roi_hull = ConvexHull(poly_pts)
                hull_xy = poly_pts[roi_hull.vertices]
                poly = Polygon(hull_xy)
                x, y = poly.exterior.xy
                # filled polygon with value-coded color
                val = deviations[int(r), j]
                pc = plt.Polygon(np.column_stack([x, y]), closed=True, fill=True, linewidth=0,
                                 facecolor=plt.cm.viridis(norm(val)))
                ax.add_patch(pc)
            except Exception:
                continue

        # overlay global hull outline
        for simplex in hull.simplices:
            ax.plot(pts[simplex, 0], pts[simplex, 1], 'red', linewidth=2)

        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, extend='max', ticks=levels,
                            label=f'|Δ| % of {pg}-gons')
        cbar.set_ticklabels(['0-5', '5-10', '10-15', '15-20', '>20'])

        ax.set_title(f'Absolute deviation for {shape_names.get(int(pg), f"{pg}-gon")}')
        ax.set_aspect('equal')
        plt.savefig(os.path.join(outdir, f'{stamp_prefix}lagrangian_absolute_deviation_{pg}-gon.png'),
                    bbox_inches='tight', dpi=300)
        plt.close(fig)

def plot_lagrangian_rois_like_before(
    frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    roi_index_per_frame: List[np.ndarray],
    roi_colors: List[Tuple[float, float, float, float]],
    save_dir: str,
    page_w: float,
    page_h: float,
    plot_contour: bool = True,
    delaunay_polygonize: bool =False,
) -> None:
    """
    Make per-frame plots identical in style to the Eulerian version, but color by Lagrangian ROI.

    :param frames: List of frame tuples (ids, x, y, polygons).
    :type  frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    :param roi_index_per_frame: Per-frame array of ROI index per row (aligned to frames).
    :type  roi_index_per_frame: List[np.ndarray]
    :param roi_colors: RGBA colors, one per ROI.
    :type  roi_colors: List[Tuple[float, float, float, float]]
    :param save_dir: Output directory.
    :type  save_dir: str
    :param page_w: Figure width (inches).
    :type  page_w: float
    :param page_h: Figure height (inches).
    :type  page_h: float
    :param plot_contour: Whether to overlay convex hull of the points.
    :type  plot_contour: bool
    :return: None.
    :rtype:  None
    """
    for sheet_index, (ids, x, y, polys) in enumerate(frames, start=1):
        print(f"Plotting Lagrangian ROIs for Frame {sheet_index}/{len(frames)} with {len(ids)} cells...")
        x_coords, y_coords = x, y
        points = np.column_stack((x_coords, y_coords))

        # Plot contour if enabled
        if plot_contour or plot_roi_analysis:
            # Combine x and y coordinates into a single array
            points = np.column_stack((x_coords, y_coords))

            if delaunay_polygonize:
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
        fig, ax = plt.subplots(figsize=(page_w/3, (page_h/9)+0.25))

        bx, by = poly.exterior.xy
        ax.plot(bx, by, 'k-', linewidth=0.5, label='Concave boundary (α-shape)')
        for hole in poly.interiors:
            hx, hy = hole.xy
            ax.plot(hx, hy, 'k--', linewidth=0.5)

        # overlay: color by Lagrangian ROI
        roi_idx_vec = roi_index_per_frame[sheet_index-1]
        n_roi = len(roi_colors)
        for r in range(n_roi):
            mask = roi_idx_vec == r
            if np.any(mask):
                ax.scatter(x_coords[mask], y_coords[mask], s=0.25, color=roi_colors[r])

        ax.set_title(f'Wing Disc Lagrangian ROI - Frame {sheet_index}',fontsize=6)
        ax.set_xlabel('X Coordinate', fontsize=3)
        ax.set_ylabel('Y Coordinate', fontsize=3)
        ax.tick_params(axis='both', which='both', labelsize=3)
        for spine in ax.spines.values():
            spine.set_linewidth(0.25)  # optional: scale spine width
        ax.xaxis.set_tick_params(width=0.8, length=0.25)
        ax.yaxis.set_tick_params(width=0.8, length=0.25)
        ax.set_aspect('equal')  # Equal scaling for x and y axes


        plt.savefig(os.path.join(save_dir, f'{stamp_prefix}lagrange_roi_frame_{sheet_index}.png'), dpi=300)
        plt.close(fig)

# ---------------------------------------------------------------------
# Support inference (naming kept: _infer_support)
# ---------------------------------------------------------------------
def _infer_support(
    all_data: List[np.ndarray],
    min_global_pct: float = 0.0
) -> np.ndarray:
    """
    Infer polygon classes present globally and optionally drop rare ones.

    :param all_data: List of per-frame polygon arrays over the whole tissue.
    :type  all_data: List[np.ndarray]
    :param min_global_pct: Minimum global percentage to keep a class.
    :type  min_global_pct: float
    :return: Sorted unique polygon classes.
    :rtype: np.ndarray
    """
    if len(all_data) == 0:
        return np.array([], dtype=int)
    G = np.concatenate([a.ravel().astype(int) for a in all_data if a.size > 0])
    if G.size == 0:
        return np.array([], dtype=int)
    classes, counts = np.unique(G, return_counts=True)
    pct = counts / counts.sum() * 100.0
    keep = pct >= float(min_global_pct)
    return classes[keep].astype(int)

# ---------------------------------------------------------------------
# Counts tally (naming kept: _tally_roi_counts)
# ---------------------------------------------------------------------
def _tally_roi_counts(
    roi_data_list: List[Dict],
    support: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pool counts across frames for each ROI on the given support.

    :param roi_data_list: Each ROI dict has 'polygon_numbers': List[np.ndarray].
    :type  roi_data_list: List[Dict]
    :param support: Polygon classes (K,).
    :type  support: np.ndarray
    :return: (C, N) with C [R,K] counts and N [R] totals.
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    R, K = len(roi_data_list), len(support)
    C = np.zeros((R, K), dtype=np.int64)
    for r, roi in enumerate(roi_data_list):
        pooled = (np.concatenate([a.ravel().astype(int) for a in roi.get("polygon_numbers", []) if a.size > 0])
                  if roi.get("polygon_numbers") else np.array([], dtype=int))
        if pooled.size == 0:
            continue
        for j, cls in enumerate(support):
            C[r, j] = int(np.sum(pooled == cls))
    N = C.sum(axis=1)
    return C, N

# ---------------------------------------------------------------------
# Dirichlet posterior mean (naming kept)
# ---------------------------------------------------------------------
def _dirichlet_posterior_mean(
    counts: np.ndarray,
    alpha: float = 0.5
) -> np.ndarray:
    """
    Symmetric-Dirichlet posterior mean(s) for counts.

    :param counts: Either [K] or [R,K] integer counts.
    :type  counts: np.ndarray
    :param alpha: Symmetric prior concentration per class.
    :type  alpha: float
    :return: Posterior mean(s) (sum to 1 along last axis).
    :rtype: np.ndarray
    """
    arr = np.asarray(counts, dtype=float)
    single = (arr.ndim == 1)
    if single:
        arr = arr[None, :]
    post = arr + float(alpha)
    post /= post.sum(axis=1, keepdims=True).clip(min=np.finfo(float).eps)
    return post[0] if single else post

# ---------------------------------------------------------------------
# JSD (naming kept: _js_divergence)
# ---------------------------------------------------------------------
def _js_divergence(
    P: np.ndarray,
    Q: np.ndarray,
    base: float = 2.0
) -> np.ndarray:
    """
    Jensen–Shannon divergence between each row of P and Q (bits by default).

    :param P: [R,K] distributions.
    :type  P: np.ndarray
    :param Q: [K] reference distribution.
    :type  Q: np.ndarray
    :param base: Log base (2 → bits).
    :type  base: float
    :return: [R] JSD values.
    :rtype: np.ndarray
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float).ravel()
    Q = Q / Q.sum()
    M = 0.5 * (P + Q[None, :])

    def _kl(A, B):
        A = np.clip(A, np.finfo(float).eps, 1.0)
        B = np.clip(B, np.finfo(float).eps, 1.0)
        return np.sum(A * (np.log(A) - np.log(B)), axis=1)

    js = 0.5 * _kl(P, M) + 0.5 * _kl(Q[None, :], M)
    if base != np.e:
        js /= np.log(base)
    return js

# ---------------------------------------------------------------------
# W1 on discrete support (naming kept: _wasserstein1_discrete)
# ---------------------------------------------------------------------
def _wasserstein1_discrete(
    P: np.ndarray,
    Q: np.ndarray,
    support: np.ndarray
) -> np.ndarray:
    """
    1D Wasserstein-1 distance on the ordered discrete support.

    :param P: [R,K] distributions.
    :type  P: np.ndarray
    :param Q: [K] reference distribution.
    :type  Q: np.ndarray
    :param support: Sorted classes (length K).
    :type  support: np.ndarray
    :return: [R] distances (units: “sides”).
    :rtype: np.ndarray
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float).ravel()
    cP = np.cumsum(P, axis=1)
    cQ = np.cumsum(Q, axis=0)
    s = support.astype(float)
    ds = np.diff(s, prepend=s[0])
    return np.sum(np.abs(cP - cQ[None, :]) * ds[None, :], axis=1)

# ---------------------------------------------------------------------
# Multinomial residuals (naming kept: _multinomial_residuals)
# ---------------------------------------------------------------------
def _multinomial_residuals(
    counts: np.ndarray,
    p0: np.ndarray
) -> np.ndarray:
    """
    Standardized residuals per class under a multinomial with reference p0.

    :param counts: [R, K] integer counts per ROI.
    :type  counts: np.ndarray
    :param p0: [K] reference probabilities (sum to 1).
    :type  p0: np.ndarray
    :return: [R, K] z-residuals; rows with N==0 are set to 0.
    :rtype: np.ndarray
    """
    counts = np.asarray(counts, dtype=float)
    p0 = np.asarray(p0, dtype=float).ravel()
    R, K = counts.shape

    N = counts.sum(axis=1, keepdims=True)            # [R,1]
    mu = N * p0[None, :]                             # [R,K]
    var = N * p0[None, :] * (1.0 - p0[None, :])      # [R,K]

    z = (counts - mu) / np.sqrt(np.clip(var, np.finfo(float).eps, None))

    # rows with zero total -> set entire row to 0
    row_zero = (N.squeeze() == 0)                    # [R]
    if np.any(row_zero):
        z[row_zero, :] = 0.0

    return z


# ---------------------------------------------------------------------
# Topological charge (naming kept: _topological_charge)
# ---------------------------------------------------------------------
def _topological_charge(
    counts: np.ndarray,
    support: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Mean charge E[n-6] and absolute |n-6| using empirical probabilities.

    :param counts: [R,K] counts per ROI.
    :type  counts: np.ndarray
    :param support: [K] polygon classes.
    :type  support: np.ndarray
    :return: (mean_q, mean_abs_q) arrays [R].
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    N = counts.sum(axis=1, keepdims=True).astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        p = counts / np.clip(N, np.finfo(float).eps, None)
    q = (support.astype(float) - 6.0)[None, :]
    mean_q = np.sum(p * q, axis=1)
    mean_abs_q = np.sum(p * np.abs(q), axis=1)
    mean_q[np.squeeze(N) == 0] = np.nan
    mean_abs_q[np.squeeze(N) == 0] = np.nan
    return mean_q, mean_abs_q
def build_roi_data_list_lagrangian(
    frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    roi_index_per_frame: List[np.ndarray],
    n_roi: int,
    include_boundary: bool = True
) -> List[Dict]:
    """
    Construct `roi_data_list` for Lagrangian ROIs from tracked memberships.

    :param frames: List of (ids, x, y, polygons) per frame (your `frames` list).
    :type  frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    :param roi_index_per_frame: Per-frame arrays of ROI indices aligned with `frames` rows.
    :type  roi_index_per_frame: List[np.ndarray]
    :param n_roi: Number of ROIs.
    :type  n_roi: int
    :param include_boundary: If False, drop rows with id == -1 from the ROI tallies.
    :type  include_boundary: bool
    :return: roi_data_list where each item has:
             - 'roi_id'           : int
             - 'polygon_numbers'  : List[np.ndarray] (one per frame; may be empty arrays)
             - 'cell_counts'      : List[int]        (same length as polygon_numbers)
    :rtype:  List[Dict]
    """
    roi_data_list: List[Dict] = []
    for r in range(n_roi):
        roi_data_list.append({
            "roi_id": r,
            "polygon_numbers": [],
            "cell_counts": []
        })

    for t, (ids, xs, ys, polys) in enumerate(frames):
        roi_idx = roi_index_per_frame[t]
        if roi_idx.shape[0] != ids.shape[0]:
            raise ValueError(f"Frame {t+1}: roi_index_per_frame length {roi_idx.shape[0]} "
                             f"!= ids length {ids.shape[0]}")

        # optional boundary filtering
        mask_valid = np.ones_like(ids, dtype=bool)
        if not include_boundary:
            mask_valid &= (ids != -1)

        for r in range(n_roi):
            mask = (roi_idx == r) & mask_valid
            pg = np.array(polys[mask], dtype=int)  # may be empty
            roi_data_list[r]["polygon_numbers"].append(pg)
            roi_data_list[r]["cell_counts"].append(int(mask.sum()))

    return roi_data_list
# ---------------------------------------------------------------------
# MAIN: compute_topology_deviation_metrics (naming kept exactly)
# ---------------------------------------------------------------------
def compute_topology_deviation_metrics(
    roi_data_list: List[Dict],
    all_polygon_numbers: List[np.ndarray],
    polygons: Optional[np.ndarray] = None,
    alpha: float = 0.5,
    min_global_pct: float = 0.0
) -> Dict[str, np.ndarray]:
    """
    Compute topology deviation metrics per ROI (Lagrangian/Eulerian agnostic).

    :param roi_data_list: Each ROI dict has 'polygon_numbers': List[np.ndarray] over frames.
    :type  roi_data_list: List[Dict]
    :param all_polygon_numbers: Per-frame arrays over the whole tissue.
    :type  all_polygon_numbers: List[np.ndarray]
    :param polygons: Optional sorted array of classes to include; if None, inferred.
    :type  polygons: Optional[np.ndarray]
    :param alpha: Symmetric Dirichlet prior concentration per class.
    :type  alpha: float
    :param min_global_pct: Drop classes with global percentage below this threshold.
    :type  min_global_pct: float
    :return: Dict of arrays (support, C, N, p_roi, p_global, jsd, w1, z, mean_q, mean_abs_q).
    :rtype: Dict[str, np.ndarray]
    """
    # 1) support
    support = (np.asarray(polygons, dtype=int)
               if polygons is not None else
               _infer_support(all_polygon_numbers, min_global_pct=min_global_pct))
    if support.size == 0:
        raise ValueError("No polygon classes found/kept. Check inputs or min_global_pct.")

    # 2) counts
    C, N = _tally_roi_counts(roi_data_list, support)
    G_counts = C.sum(axis=0)  # ROI-pooled global (OK if ROIs cover all cells)

    # 3) posteriors
    p_roi = _dirichlet_posterior_mean(C, alpha=alpha)           # [R,K]
    p_global = _dirichlet_posterior_mean(G_counts, alpha=alpha) # [K]

    # 4) metrics
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

# ---------------------------------------------------------------------
# Convenience: mean_sidedness_per_roi (naming kept)
# ---------------------------------------------------------------------
def mean_sidedness_per_roi(
    metrics: Dict[str, np.ndarray],
    posterior: bool = True,
    empty_policy: str = "nan"
) -> np.ndarray:
    """
    Compute E[n] per ROI from `metrics`.

    :param metrics: Dict returned by `compute_topology_deviation_metrics`.
    :type  metrics: Dict[str, np.ndarray]
    :param posterior: If True use posterior p_roi; else empirical C/N.
    :type  posterior: bool
    :param empty_policy: For N==0, one of {"nan","zero"}.
    :type  empty_policy: str
    :return: Array [R] of mean sidedness per ROI.
    :rtype: np.ndarray
    """
    support = metrics["support"].astype(float)  # [K]
    if posterior:
        p = metrics["p_roi"]                    # [R,K]
        out = (p * support[None, :]).sum(axis=1)
    else:
        C = metrics["C"]                        # [R,K]
        N = metrics["N"].astype(float)          # [R]
        with np.errstate(invalid='ignore', divide='ignore'):
            p_emp = (C.T / np.clip(N, np.finfo(float).eps, None)).T
        out = (p_emp * support[None, :]).sum(axis=1)
        out[N == 0] = np.nan

    if empty_policy == "zero":
        out = np.nan_to_num(out, nan=0.0)
    else:
        # "nan" (default): leave NaNs as-is
        pass
    return out




def plot_mean_sidedness_numbers(
    metrics: Dict[str, np.ndarray],
    # --- Eulerian grid inputs (set all 4 or none) ---
    x_grid_edges: Optional[np.ndarray] = None,
    y_grid_edges: Optional[np.ndarray] = None,
    num_grid_x: Optional[int] = None,
    num_grid_y: Optional[int] = None,
    # --- Lagrangian outlines inputs (set all 3 or none) ---
    roi_members_f0: Optional[Dict[int, List[int]]] = None,
    x0: Optional[np.ndarray] = None,
    y0: Optional[np.ndarray] = None,
    # --- common / styling ---
    page_w: float = 6.3,
    page_h: float = 9.73,
    save_path: str = "mean_sidedness.png",
    hull_points: Optional[np.ndarray] = None,
    posterior: bool = True,
    empty_policy: str = "nan",
    text_fmt: str = "{:.2f}",
    na_label: str = "",
    cmap_name: str = "viridis",
    face_alpha: float = 0.25,
    edge_lw: float = 1.0,
) -> None:
    """
    Plot mean polygon sidedness E[n] per ROI either on an Eulerian grid (old behavior)
    or on Lagrangian Frame-0 ROI outlines (new behavior), using the same function name.

    The mode is auto-detected:
      • If grid arguments are provided (x_grid_edges, y_grid_edges, num_grid_x, num_grid_y),
        draw the Eulerian grid heatmap with numbers centered in each cell.
      • Else, if (roi_members_f0, x0, y0) are provided, draw filled polygons for each
        Frame-0 Lagrangian ROI (convex hull of its members) and place the numbers at
        polygon centroids.

    :param metrics: Dict returned by `compute_topology_deviation_metrics`.
    :type  metrics: Dict[str, np.ndarray]
    :param x_grid_edges: X edges of Eulerian grid (len = num_grid_x+1). If None, use Lagrangian mode.
    :type  x_grid_edges: Optional[np.ndarray]
    :param y_grid_edges: Y edges of Eulerian grid (len = num_grid_y+1). If None, use Lagrangian mode.
    :type  y_grid_edges: Optional[np.ndarray]
    :param num_grid_x: Number of grid divisions along x (Eulerian mode).
    :type  num_grid_x: Optional[int]
    :param num_grid_y: Number of grid divisions along y (Eulerian mode).
    :type  num_grid_y: Optional[int]
    :param roi_members_f0: Mapping roi_id -> list of Frame-0 row indices (Lagrangian mode).
    :type  roi_members_f0: Optional[Dict[int, List[int]]]
    :param x0: Frame-0 x coordinates aligned with roi_members_f0 indices (Lagrangian mode).
    :type  x0: Optional[np.ndarray]
    :param y0: Frame-0 y coordinates aligned with roi_members_f0 indices (Lagrangian mode).
    :type  y0: Optional[np.ndarray]
    :param page_w: Figure width (inches).
    :type  page_w: float
    :param page_h: Figure height (inches).
    :type  page_h: float
    :param save_path: Output PNG path.
    :type  save_path: str
    :param hull_points: Optional Nx2 array to draw tissue hull.
    :type  hull_points: Optional[np.ndarray]
    :param posterior: If True, use posterior p_roi; else empirical C/N for E[n].
    :type  posterior: bool
    :param empty_policy: For N==0 in empirical mode: {"nan","zero"}.
    :type  empty_policy: str
    :param text_fmt: Format string for numbers drawn inside cells/polygons.
    :type  text_fmt: str
    :param na_label: Label to draw when a value is NaN.
    :type  na_label: str
    :param cmap_name: Matplotlib colormap name.
    :type  cmap_name: str
    :param face_alpha: Fill transparency for Lagrangian polygons.
    :type  face_alpha: float
    :param edge_lw: Edge linewidth for Lagrangian polygons.
    :type  edge_lw: float
    :return: None.
    :rtype: None
    """
    # --- compute E[n] per ROI from metrics (uses your naming) ---
    support = metrics["support"].astype(float)
    if posterior:
        p = metrics["p_roi"]                    # [R,K]
        mean_n = (p * support[None, :]).sum(axis=1)  # [R]
    else:
        C = metrics["C"]; N = metrics["N"].astype(float)
        with np.errstate(invalid='ignore', divide='ignore'):
            p_emp = (C.T / np.clip(N, np.finfo(float).eps, None)).T
        mean_n = (p_emp * support[None, :]).sum(axis=1)
        mean_n[N == 0] = np.nan
    # -------------------------------------------------------------

    # ---- MODE SELECTION ----
    grid_mode = (x_grid_edges is not None and y_grid_edges is not None
                 and num_grid_x is not None and num_grid_y is not None)
    lagr_mode = (roi_members_f0 is not None and x0 is not None and y0 is not None)

    if grid_mode and lagr_mode:
        raise ValueError("Provide either grid arguments OR Lagrangian outlines, not both.")
    if not grid_mode and not lagr_mode:
        raise ValueError("Missing inputs: provide grid args for Eulerian mode or (roi_members_f0, x0, y0) for Lagrangian mode.")

    # common color normalization
    finite_vals = mean_n[np.isfinite(mean_n)]
    if finite_vals.size == 0:
        vmin, vmax = 0.0, 1.0
    else:
        vmin, vmax = float(np.min(finite_vals)), float(np.max(finite_vals))
        if np.isclose(vmin, vmax):
            eps = 1e-6
            vmin, vmax = vmin - eps, vmax + eps
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    fig, ax = plt.subplots(figsize=(page_w, page_h/3))

    if grid_mode:
        # ----------------- EULERIAN: your original behavior -----------------
        # faint background to improve legibility
        Z = np.zeros((num_grid_y, num_grid_x))
        ax.pcolormesh(x_grid_edges, y_grid_edges, Z, cmap='Greys', shading='auto', alpha=0.12)

        # write numbers at grid centers
        dx = (x_grid_edges[1] - x_grid_edges[0]) / 2.0
        dy = (y_grid_edges[1] - y_grid_edges[0]) / 2.0
        mean_n_grid = mean_n.reshape((num_grid_x, num_grid_y))
        for i in range(num_grid_x):
            for j in range(num_grid_y):
                val = mean_n_grid[i, j]
                cx = x_grid_edges[i] + dx
                cy = y_grid_edges[j] + dy
                if np.isfinite(val):
                    ax.text(cx, cy, text_fmt.format(val), ha='center', va='center', fontsize=4.2, color='black')
                elif na_label:
                    ax.text(cx, cy, na_label, ha='center', va='center', fontsize=4.2, color='black')

        # grid lines
        for x in x_grid_edges:
            ax.plot([x, x], [y_grid_edges[0], y_grid_edges[-1]], 'k--', lw=1)
        for y in y_grid_edges:
            ax.plot([x_grid_edges[0], x_grid_edges[-1]], [y, y], 'k--', lw=1)

        title = 'Mean polygon sidedness per ROI (Eulerian grid)'

    else:
        # ----------------- LAGRANGIAN: Frame-0 outlines -----------------
        roi_ids = sorted(roi_members_f0.keys())
        # draw filled polygons with colormap + numbers at centroids
        for r in roi_ids:
            idxs = np.asarray(roi_members_f0[r], dtype=int)
            if idxs.size < 3:
                continue
            pts = np.column_stack([x0[idxs], y0[idxs]])
            try:
                hull = ConvexHull(pts)
                hull_xy = pts[hull.vertices]
            except Exception:
                hull_xy = pts
            poly = Polygon(hull_xy)
            if (not poly.is_valid) or poly.area == 0:
                continue
            # Use negative buffer to shrink polygon
            buffer_distance = -3.0  # Adjust based on your data scale
            shrunk_poly = poly.buffer(buffer_distance)
            # Only use if buffer operation created a valid polygon
            if shrunk_poly.is_valid and not shrunk_poly.is_empty:
                poly = shrunk_poly
            val = mean_n[r] if r < len(mean_n) else np.nan
            face_color = (0.85, 0.85, 0.85, face_alpha) if (not np.isfinite(val)) else (*cmap(norm(val))[:3], face_alpha)
            patch = plt.Polygon(np.asarray(poly.exterior.xy).T, closed=True,
                                facecolor=face_color, edgecolor='k', lw=edge_lw)
            ax.add_patch(patch)

            # label
            cx, cy = poly.centroid.x, poly.centroid.y
            if np.isfinite(val):
                ax.text(cx, cy, text_fmt.format(val), ha='center', va='center', fontsize=4.0, color='black')
            elif na_label:
                ax.text(cx, cy, na_label, ha='center', va='center', fontsize=4.0, color='black')

        title = 'Mean polygon sidedness per Lagrangian ROI'

    # tissue hull overlay (optional)
    if hull_points is not None and len(hull_points) >= 3:
        try:
            h = ConvexHull(hull_points)
            for simplex in h.simplices:
                ax.plot(hull_points[simplex, 0], hull_points[simplex, 1], 'k-', lw=2)
        except Exception:
            pass

    # colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.9)
    cbar.set_label('Mean polygon sidedness E[n]')

    ax.set_title(title)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)

def average_cell_count_per_roi(roi_data_list: List[Dict]) -> np.ndarray:
    """
    Compute the average number of cells per ROI from roi_data_list['cell_counts'].

    :param roi_data_list: List of ROI dicts; each has 'cell_counts': List[int] per frame.
    :type  roi_data_list: List[Dict]
    :return: Array [R] with mean cell count per ROI (NaN if no frames).
    :rtype: np.ndarray
    """
    avgs = []
    for roi in roi_data_list:
        counts = roi.get("cell_counts", [])
        if counts and len(counts) > 0:
            avgs.append(float(np.mean(counts)))
        else:
            avgs.append(np.nan)
    return np.asarray(avgs, dtype=float)


def plot_average_cell_count(
    roi_data_list: List[Dict],
    # --- Eulerian grid inputs (set all 4 or none) ---
    x_grid_edges: Optional[np.ndarray] = None,
    y_grid_edges: Optional[np.ndarray] = None,
    num_grid_x: Optional[int] = None,
    num_grid_y: Optional[int] = None,
    # --- Lagrangian outlines inputs (set all 3 or none) ---
    roi_members_f0: Optional[Dict[int, List[int]]] = None,
    x0: Optional[np.ndarray] = None,
    y0: Optional[np.ndarray] = None,
    # --- shared/plotting ---
    page_w: float = 6.3,
    page_h: float = 9.73,
    save_path: str = "average_cell_count_map.png",
    hull_points: Optional[np.ndarray] = None,
    cmap_name: str = "viridis",
    face_alpha: float = 0.30,
    edge_lw: float = 1.0,
    text_fmt: Optional[str] = "{:d}",
    na_label: str = "",
    title: Optional[str] = None,
) -> None:
    """
    Plot the *average number of cells per ROI* either on an Eulerian grid or on
    Lagrangian Frame-0 outlines, using a unified function signature.

    Mode selection:
      • If grid args are provided (x_grid_edges, y_grid_edges, num_grid_x, num_grid_y),
        render an Eulerian heatmap with numbers at grid centers.
      • Else, if (roi_members_f0, x0, y0) are provided, render Frame-0 polygons
        (convex hulls of ROI members) colored by the metric with optional labels.

    :param roi_data_list: Per-ROI data; each dict has 'cell_counts': List[int] per frame.
    :type  roi_data_list: List[Dict]
    :param x_grid_edges: X edges of grid (len=num_grid_x+1); set with other grid args for Eulerian mode.
    :type  x_grid_edges: Optional[np.ndarray]
    :param y_grid_edges: Y edges of grid (len=num_grid_y+1); set with other grid args for Eulerian mode.
    :type  y_grid_edges: Optional[np.ndarray]
    :param num_grid_x: Number of grid cells in x (Eulerian mode).
    :type  num_grid_x: Optional[int]
    :param num_grid_y: Number of grid cells in y (Eulerian mode).
    :type  num_grid_y: Optional[int]
    :param roi_members_f0: roi_id -> list of Frame-0 row indices (Lagrangian mode).
    :type  roi_members_f0: Optional[Dict[int, List[int]]]
    :param x0: Frame-0 x coordinates aligned with roi_members_f0 indices.
    :type  x0: Optional[np.ndarray]
    :param y0: Frame-0 y coordinates aligned with roi_members_f0 indices.
    :type  y0: Optional[np.ndarray]
    :param page_w: Figure width (inches).
    :type  page_w: float
    :param page_h: Figure height (inches).
    :type  page_h: float
    :param save_path: Output PNG path.
    :type  save_path: str
    :param hull_points: Optional Nx2 array; if provided, draw tissue hull.
    :type  hull_points: Optional[np.ndarray]
    :param cmap_name: Matplotlib colormap name.
    :type  cmap_name: str
    :param face_alpha: Fill transparency for Lagrangian polygons.
    :type  face_alpha: float
    :param edge_lw: Edge linewidth for Lagrangian polygons.
    :type  edge_lw: float
    :param text_fmt: If not None, draw numeric values; e.g., "{:.1f}".
    :type  text_fmt: Optional[str]
    :param na_label: Label drawn when a value is NaN (used if text_fmt is not None).
    :type  na_label: str
    :param title: Optional title override.
    :type  title: Optional[str]
    :return: None.
    :rtype: None
    """
    values = average_cell_count_per_roi(roi_data_list)  # [R]

    # mode detection
    grid_mode = (x_grid_edges is not None and y_grid_edges is not None
                 and num_grid_x is not None and num_grid_y is not None)
    lagr_mode = (roi_members_f0 is not None and x0 is not None and y0 is not None)

    if grid_mode and lagr_mode:
        raise ValueError("Provide either grid arguments OR Lagrangian outlines, not both.")
    if not grid_mode and not lagr_mode:
        raise ValueError("Missing inputs: grid args for Eulerian mode or (roi_members_f0, x0, y0) for Lagrangian mode.")

    finite_vals = values[np.isfinite(values)]
    if finite_vals.size == 0:
        vmin, vmax = 0.0, 1.0
    else:
        vmin, vmax = float(np.min(finite_vals)), float(np.max(finite_vals))
        if np.isclose(vmin, vmax):
            eps = 1e-6
            vmin, vmax = vmin - eps, vmax + eps

    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    fig, ax = plt.subplots(figsize=(page_w, page_h/3))

    if grid_mode:
        # Eulerian: reshape and draw colormap + numbers
        grid_vals = values.reshape((num_grid_x, num_grid_y))  # assumes ROI order matches grid order
        X, Y = np.meshgrid(x_grid_edges, y_grid_edges)
        pcm = ax.pcolormesh(X, Y, grid_vals.T, cmap=cmap, shading='auto', norm=norm)

        dx = (x_grid_edges[1] - x_grid_edges[0]) / 2.0
        dy = (y_grid_edges[1] - y_grid_edges[0]) / 2.0
        for i in range(num_grid_x):
            for j in range(num_grid_y):
                val = grid_vals[i, j]
                if text_fmt is None:
                    continue
                cx = x_grid_edges[i] + dx
                cy = y_grid_edges[j] + dy
                if np.isfinite(val):
                    ax.text(cx, cy, text_fmt.format(val), ha='center', va='center',
                            fontsize=4.2, color='white')
                elif na_label:
                    ax.text(cx, cy, na_label, ha='center', va='center',
                            fontsize=4.2, color='white')

        for x in x_grid_edges:
            ax.plot([x, x], [y_grid_edges[0], y_grid_edges[-1]], 'k--', lw=1)
        for y in y_grid_edges:
            ax.plot([x_grid_edges[0], x_grid_edges[-1]], [y, y], 'k--', lw=1)

        cbar_label = "Average number of cells per ROI"
        default_title = "Average number of cells per ROI (Eulerian grid)"

    else:
        # Lagrangian: draw Frame-0 convex hull for each ROI and fill by value
        roi_ids = sorted(roi_members_f0.keys())
        for r in roi_ids:
            idxs = np.asarray(roi_members_f0[r], dtype=int)
            if idxs.size < 3:
                continue
            pts = np.column_stack([x0[idxs], y0[idxs]])
            try:
                hull = ConvexHull(pts)
                hull_xy = pts[hull.vertices]
            except Exception:
                hull_xy = pts
            poly = Polygon(hull_xy)
            if (not poly.is_valid) or poly.area == 0.0:
                continue

            val = values[r] if r < len(values) else np.nan
            face_color = (0.85, 0.85, 0.85, face_alpha) if (not np.isfinite(val)) \
                         else (*cmap(norm(val))[:3], face_alpha)

            patch = plt.Polygon(np.asarray(poly.exterior.xy).T, closed=True,
                                facecolor=face_color, edgecolor='k', lw=edge_lw)
            ax.add_patch(patch)

            if text_fmt is not None:
                cx, cy = poly.centroid.x, poly.centroid.y
                if np.isfinite(val):
                    ax.text(cx, cy, text_fmt.format(val), ha='center', va='center',
                            fontsize=5.0, color='black')
                elif na_label:
                    ax.text(cx, cy, na_label, ha='center', va='center',
                            fontsize=5.0, color='black')

        cbar_label = "Average number of cells per ROI"
        default_title = "Average number of cells per Lagrangian ROI"

    # tissue hull overlay (optional)
    if hull_points is not None and len(hull_points) >= 3:
        try:
            h = ConvexHull(hull_points)
            for simplex in h.simplices:
                ax.plot(hull_points[simplex, 0], hull_points[simplex, 1], 'k-', lw=2)
        except Exception:
            pass

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.9)
    cbar.set_label(cbar_label)

    ax.set_title(default_title if title is None else title)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.spatial import ConvexHull
from typing import Dict, List, Optional


def plot_lagrangian_metric_on_frame0_outlines(
    values: np.ndarray,
    roi_members_f0: Dict[int, List[int]],
    x0: np.ndarray,
    y0: np.ndarray,
    *,
    page_w: float,
    page_h: float,
    save_path: str,
    hull_points: Optional[np.ndarray] = None,
    cmap_name: str = "viridis",
    face_alpha: float = 0.30,
    edge_lw: float = 1.0,
    text_fmt: Optional[str] = None,
    na_label: str = "",
    cbar_label: str = "",
    title: str = "",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None
) -> None:
    """
    Paint a scalar metric per Lagrangian ROI onto Frame-0 ROI outlines.

    Each ROI polygon is the convex hull of its Frame-0 member cells. The polygon
    is filled using a colormap scaled across all finite ROI values. Optionally,
    draw the numeric value at the polygon centroid.

    :param values: Per-ROI scalar values, shape [R] (NaNs allowed).
    :type  values: np.ndarray
    :param roi_members_f0: Mapping roi_id -> list of Frame-0 row indices (non-boundary).
    :type  roi_members_f0: Dict[int, List[int]]
    :param x0: Frame-0 x coordinates (aligned to those row indices).
    :type  x0: np.ndarray
    :param y0: Frame-0 y coordinates (aligned to those row indices).
    :type  y0: np.ndarray
    :param page_w: Figure width in inches.
    :type  page_w: float
    :param page_h: Figure height in inches.
    :type  page_h: float
    :param save_path: Output PNG path.
    :type  save_path: str
    :param hull_points: Optional Nx2 array to draw the global tissue hull.
    :type  hull_points: Optional[np.ndarray]
    :param cmap_name: Matplotlib colormap name.
    :type  cmap_name: str
    :param face_alpha: Fill transparency for ROI polygons.
    :type  face_alpha: float
    :param edge_lw: Polygon edge line width.
    :type  edge_lw: float
    :param text_fmt: If provided (e.g., "{:.2f}"), draw numbers at centroids.
    :type  text_fmt: Optional[str]
    :param na_label: Label to draw when value is NaN (used if text_fmt is not None).
    :type  na_label: str
    :param cbar_label: Colorbar label.
    :type  cbar_label: str
    :param title: Figure title.
    :type  title: str
    :param vmin: Optional lower color limit (else auto from finite values).
    :type  vmin: Optional[float]
    :param vmax: Optional upper color limit (else auto from finite values).
    :type  vmax: Optional[float]
    :return: None.
    :rtype: None
    """
    # Ensure values length covers all roi ids (pad with NaN if needed)
    roi_ids = sorted(roi_members_f0.keys())
    R_needed = (max(roi_ids) + 1) if roi_ids else 0
    if values.shape[0] < R_needed:
        vv = np.full(R_needed, np.nan, dtype=float)
        vv[:values.shape[0]] = values
        values = vv

    # Color normalization
    finite_vals = values[np.isfinite(values)]
    if vmin is None or vmax is None:
        if finite_vals.size == 0:
            vmin_auto, vmax_auto = 0.0, 1.0
        else:
            vmin_auto, vmax_auto = float(np.min(finite_vals)), float(np.max(finite_vals))
            if np.isclose(vmin_auto, vmax_auto):
                eps = 1e-6
                vmin_auto, vmax_auto = vmin_auto - eps, vmax_auto + eps
        if vmin is None: vmin = vmin_auto
        if vmax is None: vmax = vmax_auto

    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    # Simple polygon centroid (shoelace) to avoid heavy deps
    def _poly_centroid(xy: np.ndarray) -> np.ndarray:
        """
        :param xy: Array [M,2] of polygon vertices (closed or open).
        :type  xy: np.ndarray
        :return: [2] centroid (cx, cy).
        :rtype: np.ndarray
        """
        x = xy[:, 0]; y = xy[:, 1]
        # close if not closed
        if not (np.isclose(x[0], x[-1]) and np.isclose(y[0], y[-1])):
            x = np.append(x, x[0]); y = np.append(y, y[0])
        a = x[:-1]*y[1:] - x[1:]*y[:-1]
        A = np.sum(a) / 2.0
        if np.isclose(A, 0.0):
            return np.array([np.mean(x[:-1]), np.mean(y[:-1])], dtype=float)
        cx = np.sum((x[:-1] + x[1:]) * a) / (6.0 * A)
        cy = np.sum((y[:-1] + y[1:]) * a) / (6.0 * A)
        return np.array([cx, cy], dtype=float)

    fig, ax = plt.subplots(figsize=(page_w, page_h/3))

    # Draw each ROI polygon
    for r in roi_ids:
        idxs = np.asarray(roi_members_f0[r], dtype=int)
        if idxs.size < 3:
            continue
        pts = np.column_stack([x0[idxs], y0[idxs]])
        try:
            hull = ConvexHull(pts)
            hull_xy = pts[hull.vertices]
        except Exception:
            hull_xy = pts  # fallback: use all points

        val = values[r] if r < len(values) else np.nan
        rgba = (0.85, 0.85, 0.85, face_alpha) if (not np.isfinite(val)) else (*cmap(norm(val))[:3], face_alpha)

        # filled polygon + edge
        patch = plt.Polygon(hull_xy, closed=True, facecolor=rgba, edgecolor='k', lw=edge_lw)
        ax.add_patch(patch)

        # label at centroid
        if text_fmt is not None:
            cx, cy = _poly_centroid(hull_xy)
            if np.isfinite(val):
                ax.text(cx, cy, text_fmt.format(val), ha='center', va='center', fontsize=5.0, color='black')
            elif na_label:
                ax.text(cx, cy, na_label, ha='center', va='center', fontsize=5.0, color='black')

    # Optional global hull overlay
    if hull_points is not None and len(hull_points) >= 3:
        try:
            h = ConvexHull(hull_points)
            for simplex in h.simplices:
                ax.plot(hull_points[simplex, 0], hull_points[simplex, 1], 'k-', lw=2)
        except Exception:
            pass

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.9)
    if cbar_label:
        cbar.set_label(cbar_label)

    if title:
        ax.set_title(title)
    ax.set_xlabel('X Coordinate'); ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)

def plot_lagrangian_topology_maps(
    metrics: Dict[str, np.ndarray],
    roi_members_f0: Dict[int, List[int]],
    x0: np.ndarray,
    y0: np.ndarray,
    *,
    page_w: float,
    page_h: float,
    save_dir: str,
    stamp_prefix: str = "",
    hull_points: Optional[np.ndarray] = None,
    z_sym_limit: float = 4.0
) -> None:
    """
    Render Lagrangian topology maps (Frame-0 outlines) for:
    JSD, W1, hexagon z-residual, and mean absolute charge E[|n-6|].

    :param metrics: Output of `compute_topology_deviation_metrics`.
    :type  metrics: Dict[str, np.ndarray]
    :param roi_members_f0: roi_id -> list of Frame-0 row indices (non-boundary).
    :type  roi_members_f0: Dict[int, List[int]]
    :param x0: Frame-0 x coordinates aligned with roi_members_f0 indices.
    :type  x0: np.ndarray
    :param y0: Frame-0 y coordinates aligned with roi_members_f0 indices.
    :type  y0: np.ndarray
    :param page_w: Figure width (inches).
    :type  page_w: float
    :param page_h: Figure height (inches).
    :type  page_h: float
    :param save_dir: Directory to write PNGs.
    :type  save_dir: str
    :param stamp_prefix: Optional filename prefix.
    :type  stamp_prefix: str
    :param hull_points: Optional Nx2 array for tissue hull overlay.
    :type  hull_points: Optional[np.ndarray]
    :param z_sym_limit: Symmetric color limits for z residual plots (±value).
    :type  z_sym_limit: float
    :return: None.
    :rtype: None
    """
    os.makedirs(save_dir, exist_ok=True)

    # ---- 1) JSD (bits) ----
    jsd = metrics["jsd"]  # [R]
    plot_lagrangian_metric_on_frame0_outlines(
        values=jsd,
        roi_members_f0=roi_members_f0,
        x0=x0, y0=y0,
        page_w=page_w, page_h=page_h,
        save_path=os.path.join(save_dir, f"{stamp_prefix}lagrangian_map_topology_jsd.png"),
        hull_points=hull_points,
        cmap_name="viridis",
        face_alpha=0.30, edge_lw=1.0,
        text_fmt=None,  # cmap only
        cbar_label="Jensen–Shannon divergence (bits)",
        title="Topology deviation (JSD vs global) — Lagrangian ROIs"
    )

    # ---- 2) W1 (sides) ----
    w1 = metrics["w1"]
    plot_lagrangian_metric_on_frame0_outlines(
        values=w1,
        roi_members_f0=roi_members_f0,
        x0=x0, y0=y0,
        page_w=page_w, page_h=page_h,
        save_path=os.path.join(save_dir, f"{stamp_prefix}lagrangian_map_topology_w1.png"),
        hull_points=hull_points,
        cmap_name="viridis",
        face_alpha=0.30, edge_lw=1.0,
        text_fmt=None,
        cbar_label='Wasserstein-1 distance (sides)',
        title="Topology deviation (W1) — Lagrangian ROIs"
    )

    # ---- 3) z residual for hexagons ----
    support = metrics["support"]
    if 6 in set(int(k) for k in support):
        k6 = int(np.where(support == 6)[0][0])
        z6 = metrics["z"][:, k6]
        # We'll use the generic plotter but with symmetric limits via a wrapper
        # Since the generic function doesn't expose vmin/vmax, we provide a small shim:
        _plot_lagrangian_metric_with_limits(
            values=z6,
            roi_members_f0=roi_members_f0,
            x0=x0, y0=y0,
            page_w=page_w, page_h=page_h,
            save_path=os.path.join(save_dir, f"{stamp_prefix}lagrangian_map_z_hexagon.png"),
            hull_points=hull_points,
            cmap_name="coolwarm",
            face_alpha=0.30, edge_lw=1.0,
            text_fmt=None,
            cbar_label="Hexagon residual z",
            title="Hexagon enrichment/depletion (z residual) — Lagrangian ROIs",
            vmin=-abs(z_sym_limit), vmax=abs(z_sym_limit)
        )

    # ---- 4) Mean absolute charge E[|n-6|] ----
    absq = metrics["mean_abs_q"]
    plot_lagrangian_metric_on_frame0_outlines(
        values=absq,
        roi_members_f0=roi_members_f0,
        x0=x0, y0=y0,
        page_w=page_w, page_h=page_h,
        save_path=os.path.join(save_dir, f"{stamp_prefix}lagrangian_map_mean_abs_charge.png"),
        hull_points=hull_points,
        cmap_name="viridis",
        face_alpha=0.30, edge_lw=1.0,
        text_fmt=None,
        cbar_label="E[|n−6|]",
        title="Mean absolute topological charge — Lagrangian ROIs"
    )

    # ---- Console summary (global composition & charge) ----
    p_global = metrics["p_global"]
    Eabs_global = float(np.sum(p_global * np.abs(metrics["support"] - 6)))
    print("Global p(n):", dict(zip(metrics["support"].tolist(), np.round(p_global, 3))))
    print("Global E[|n-6|]:", round(Eabs_global, 3))


def _plot_lagrangian_metric_with_limits(
    values: np.ndarray,
    roi_members_f0: Dict[int, List[int]],
    x0: np.ndarray,
    y0: np.ndarray,
    *,
    page_w: float,
    page_h: float,
    save_path: str,
    hull_points: Optional[np.ndarray],
    cmap_name: str,
    face_alpha: float,
    edge_lw: float,
    text_fmt: Optional[str],
    cbar_label: str,
    title: str,
    vmin: float,
    vmax: float
) -> None:
    """
    Same as plot_lagrangian_metric_on_frame0_outlines but with explicit vmin/vmax
    (used for the symmetric z-residual map).
    """
    from shapely.geometry import Polygon  # local import to avoid issues if Shapely is optional

    roi_ids = sorted(roi_members_f0.keys())
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    fig, ax = plt.subplots(figsize=(page_w, page_h/3))
    for r in roi_ids:
        idxs = np.asarray(roi_members_f0[r], dtype=int)
        if idxs.size < 3:
            continue
        pts = np.column_stack([x0[idxs], y0[idxs]])
        try:
            hull = ConvexHull(pts)
            hull_xy = pts[hull.vertices]
        except Exception:
            hull_xy = pts
        poly = Polygon(hull_xy)
        if (not poly.is_valid) or poly.area == 0.0:
            continue

        val = values[r] if r < len(values) else np.nan
        if np.isfinite(val):
            rgba = (*cmap(norm(val))[:3], face_alpha)
        else:
            rgba = (0.85, 0.85, 0.85, face_alpha)

        patch = plt.Polygon(np.asarray(poly.exterior.xy).T, closed=True,
                            facecolor=rgba, edgecolor='k', lw=edge_lw)
        ax.add_patch(patch)

        if text_fmt is not None:
            cx, cy = poly.centroid.x, poly.centroid.y
            if np.isfinite(val):
                ax.text(cx, cy, text_fmt.format(val), ha='center', va='center', fontsize=5.0, color='black')

    # Hull overlay (optional)
    if hull_points is not None and len(hull_points) >= 3:
        try:
            h = ConvexHull(hull_points)
            for simplex in h.simplices:
                ax.plot(hull_points[simplex, 0], hull_points[simplex, 1], 'k-', lw=2)
        except Exception:
            pass

    # Colorbar with fixed limits
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.9)
    cbar.set_label(cbar_label)

    ax.set_title(title)
    ax.set_xlabel('X Coordinate'); ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close(fig)

# -------------------------------------------------------------------------
#                          MAIN PIPELINE
# -------------------------------------------------------------------------

def main():
    """
    Execute the Lagrangian ROI pipeline on the workbook and produce figures.

    :return: None.
    :rtype:  None
    """
    # 1) Load workbook and frames
    file_path = "../Xls_Data"
    file_name = "all_frames_wing_discs.xls"
    wb = load_workbook(file_path, file_name)
    frames = read_all_frames(wb)

    # # Loop through each sheet in the workbook, starting from the first sheet (index 0)
    # for sheet_index in range(0, wb.nsheets):
    #     x_coords, y_coords, polygon_numbers = extract_data(wb, sheet_index)
    #     print(f'Sheet {sheet_index + 1} : Number of cells = {len(x_coords)}')
    #
    #     # Store data for mean percentage calculation
    #     all_data.append(polygon_numbers)
    #     all_x_coords.append(x_coords)
    #     all_y_coords.append(y_coords)
    #     all_polygon_numbers.append(polygon_numbers)
    print(f"how many framse? : {len(frames)}")
    print(f"Number of sheets (including header): {wb.nsheets}; processed frames: {len(frames)}")

    # 2) Build Lagrangian ROIs at Frame 0
    ids0, x0, y0, polys0 = frames[0]
    id2roi_0, roi_seeds_xy = build_lagrangian_rois(ids0, x0, y0, target=target_cells_per_roi)
    # Colors for Lagrangian ROIs (stable across frames)
    num_lagr_rois = len(roi_seeds_xy)
    cmap = plt.get_cmap('gist_ncar', num_lagr_rois)
    roi_colors = [cmap(i / max(1, num_lagr_rois)) for i in range(num_lagr_rois)]
    random.Random(7).shuffle(roi_colors)  # optional shuffle for visual separation

    # Save membership indices for visualization (Frame 0)
    # We store indices of entries in Frame 0 (not IDs) per ROI (non-boundary members only)
    mask_nb0 = ids0 != -1
    idx_nb0 = np.where(mask_nb0)[0]
    # reconstruct assignment over frame 0 indices:
    # map each non-boundary *id* to roi, collect their frame indices
    roi_members_f0: Dict[int, List[int]] = {r: [] for r in range(len(roi_seeds_xy))}
    for i in idx_nb0:
        r = id2roi_0.get(int(ids0[i]), None)
        if r is not None and r >= 0:
            roi_members_f0[int(r)].append(int(i))
    np.save(os.path.join(save_dir, "roi_members_frame0.npy"), roi_members_f0, allow_pickle=True)

    # 3) Track and collect per-ROI data across frames
    roi_data, x0_for_plot, y0_for_plot, p0_dummy, roi_index_per_frame = track_and_collect(
        frames, id2roi_0, roi_seeds_xy
    )
    # First-seen frame per cell id (for ages)
    first_seen = compute_first_seen(frames)

    # Pick the ROI you want to inspect (change this index as needed)
    roi_to_debug = 24  # e.g., start with 0, then try others that seem to intersect

    # Print to console and also write a CSV (optional)
    debug_csv = os.path.join(save_dir, f"debug_roi_{roi_to_debug}.csv")
    debug_print_lagrangian_roi(
        frames=frames,
        roi_index_per_frame=roi_index_per_frame,
        roi_colors=roi_colors,
        roi_id=roi_to_debug,
        first_seen=first_seen,
        csv_path=debug_csv,  # set to None to skip file
    )

    # plot_lagrangian_rois_like_before(
    #     frames=frames,
    #     roi_index_per_frame=roi_index_per_frame,
    #     roi_colors=roi_colors,
    #     save_dir=save_dir,
    #     page_w=page_w,
    #     page_h=page_h,
    #     plot_contour=True
    # )

    # 4) Global topology across frames (for deviation baseline)
    all_frame_polys = [fr[3] for fr in frames]  # take polygons array per frame
    filtered_polys, mean_pct, std_pct = compute_global_distribution(all_frame_polys, thresh=0.5)

    if plot_mean_of_all_frames_cell_shapes:
        plot_global_distribution(filtered_polys, mean_pct, std_pct,
                                 outpath=os.path.join(save_dir, f"{stamp_prefix}Mean_topology_over_all_frames.png"))

    # 5) Per-ROI deviations from global
    deviations = compute_roi_deviations(roi_data, filtered_polys, mean_pct)

    # 6) Heatmaps overlaid in space via ROI hulls at Frame 0
    if plot_roi_deviations:
        plot_roi_deviation_heatmaps(
            deviations, filtered_polys, x0_for_plot, y0_for_plot, outdir=save_dir
        )

    # 7) Diagnostics akin to your average count map (Lagrangian)
    # Here: average #entries per ROI across frames (not spatial grid)
    avg_counts = [np.mean(roi_data[r]["cell_counts"]) if roi_data[r]["cell_counts"] else 0.0
                  for r in range(len(roi_seeds_xy))]
    # Calculate the mean of means
    mean_of_means = np.mean(avg_counts)

    fig, ax = plt.subplots(figsize=(page_w, page_h/3))
    ax.bar(np.arange(len(avg_counts)), avg_counts)
    # Add horizontal line for the mean of means
    ax.axhline(y=mean_of_means, color='r', linestyle='--', linewidth=1.5)
    ax.text(len(avg_counts)*0.02, mean_of_means*1.05, f'Mean: {mean_of_means:.2f} cells',
            color='r', fontweight='bold')

    ax.set_xlabel('Lagrangian ROI index'); ax.set_ylabel('Average #cells (incl. boundary where present)')
    ax.set_title('Average number of cells per Lagrangian ROI')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f"{stamp_prefix}avg_cell_count_per_lagrangian_roi.png"), dpi=300)
    plt.close(fig)

    # 8) Save initial tissue outline (as you did)
    pts0 = np.column_stack((x0_for_plot, y0_for_plot))
    hull0 = ConvexHull(pts0)
    hull_points = pts0[hull0.vertices]
    np.save(os.path.join(save_dir, 'wing_initial_contour.npy'), hull_points)
    initial_contour = Polygon(hull_points)
    with open(os.path.join(save_dir, 'wing_initial_contour.wkt'), 'w') as f:
        f.write(initial_contour.wkt)

    roi_data_list = build_roi_data_list_lagrangian(
        frames=frames,
        roi_index_per_frame=roi_index_per_frame,
        n_roi=len(roi_seeds_xy),
        include_boundary=False  # set False if you want to exclude id == -1 from metrics
    )

    # compute metrics using your shared API
    all_polygon_numbers = [fr[3] for fr in frames]
    # Keep all observed polygon classes, but drop globally ultra-rare ones (<0.2%)
    metrics = compute_topology_deviation_metrics(
        roi_data_list=roi_data_list,
        all_polygon_numbers=all_polygon_numbers,
        polygons=None,  # or np.array([3,4,5,6,7])
        alpha=0.5,  # Jeffreys prior; set to 1.0 for Laplace
        min_global_pct=0.2
    )
    roi_members_f0 = np.load(os.path.join(save_dir, "roi_members_frame0.npy"), allow_pickle=True).item()
    plot_mean_sidedness_numbers(
        metrics,
        roi_members_f0=roi_members_f0,
        x0=x0_for_plot,
        y0=y0_for_plot,
        page_w=page_w, page_h=page_h,
        save_path=f"{save_dir}/{stamp_prefix}mean_sidedness_lagrangian.png",
        hull_points=np.column_stack((x0_for_plot, y0_for_plot)),
        posterior=True, empty_policy="nan"
    )
    roi_members_f0 = np.load(os.path.join(save_dir, "roi_members_frame0.npy"), allow_pickle=True).item()
    plot_average_cell_count(
        roi_data_list=roi_data_list,
        roi_members_f0=roi_members_f0,
        x0=x0_for_plot, y0=y0_for_plot,
        page_w=page_w, page_h=page_h,
        save_path=f"{save_dir}/{stamp_prefix}average_cell_count_map_lagrangian.png",
        hull_points=np.column_stack((x0_for_plot, y0_for_plot)),
        text_fmt="{:.0f}"
    )
    roi_members_f0 = np.load(os.path.join(save_dir, "roi_members_frame0.npy"), allow_pickle=True).item()
    hull_points = np.column_stack((x0_for_plot, y0_for_plot))  # or your precomputed hull vertices

    plot_lagrangian_topology_maps(
        metrics=metrics,
        roi_members_f0=roi_members_f0,
        x0=x0_for_plot,
        y0=y0_for_plot,
        page_w=page_w,
        page_h=page_h,
        save_dir=save_dir,
        stamp_prefix=stamp_prefix,
        hull_points=hull_points,
        z_sym_limit=4.0  # matches your Eulerian choice vmin/vmax = ±4
    )

    print("Done.")


if __name__ == "__main__":
    main()
