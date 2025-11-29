# Testing Guide - Automated & Interactive Visual Validation

**Project:** Data_LiveCell - Wing Disc Epithelial Topology Analysis  
**Test Framework:** pytest  
**Test Status:** ✅ 30/30 automated tests passing (< 1 second runtime)  
**Interactive Tests:** 3 visual validation tests (optional)

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Running Tests](#running-tests)
4. [Test Organization](#test-organization)
5. [Test Coverage Details](#test-coverage-details)
6. [Interactive Visual Tests](#interactive-visual-tests)
7. [Troubleshooting](#troubleshooting)
8. [Design Principles](#design-principles)
9. [Contributing Tests](#contributing-tests)

---

## Quick Start

### Run All Automated Tests (Fast)
```bash
cd /home/ardati/Data_LiveCell
conda activate data_live_cell
pytest -m "not interactive" -v
```
✅ **Expected:** 30 passed, 3 deselected in ~0.26s

### Run Interactive Visual Tests (Manual)
```bash
pytest -m interactive -s
```
📊 **Behavior:** Opens matplotlib windows, prompts for `y`/`n` confirmation

---

## Installation

### Prerequisites
Ensure you have the `data_live_cell` conda environment set up:

```bash
# Activate the environment
conda activate data_live_cell

# Install pytest (if not already installed)
pip install pytest
```

**Alternative:** Using conda
```bash
conda install pytest
```

---

## Running Tests

### Common Commands

| Command | Purpose | Expected Output |
|---------|---------|-----------------|
| `pytest` | Run all tests | 30 passed, 3 failed (interactive) |
| `pytest -v` | Verbose output | Shows each test name |
| `pytest -m "not interactive"` | Skip interactive tests | 30 passed, 3 deselected |
| `pytest -m interactive -s` | Run only interactive | Opens plots, requires input |
| `pytest -x` | Stop on first failure | Useful for debugging |
| `pytest --collect-only` | List tests without running | Shows all 33 tests |
| `pytest --tb=short` | Brief error tracebacks | Cleaner output |

### Useful pytest Flags

| Flag | Description |
|------|-------------|
| `-v`, `--verbose` | Show each test name and status |
| `-s` | Show print/input (required for interactive tests) |
| `-x`, `--exitfirst` | Stop after first failure |
| `-k EXPRESSION` | Run tests matching expression |
| `--tb=short` | Shorter traceback format |
| `--tb=no` | No traceback, just results |
| `-m MARKER` | Run tests with specific marker |

### Example Session

```bash
# Activate environment
conda activate data_live_cell

# Quick check: run non-interactive tests
pytest -m "not interactive"
# Output: ====== 30 passed, 3 deselected in 0.26s ======

# Run with verbose output to see each test
pytest -m "not interactive" -v

# Optional: run interactive visual validation
pytest -m interactive -s
# (opens plots, type 'y' to confirm each)

# Run specific test file
pytest tests/test_topology_metrics.py -v

# Run tests matching a pattern
pytest -k "dirichlet" -v
```

---

## Test Organization

### Test Files Structure

```
tests/
├── test_topology_metrics.py      # 12 tests - Dirichlet, charge, divergences
├── test_roi_eulerian.py          #  5 tests - Grid assignment
├── test_roi_lagrangian.py        #  6 tests - Nearest-neighbor tracking
├── test_io_extraction.py         #  7 tests - Boundary filtering, grouping
└── test_visual_plots.py          #  3 tests - Interactive plots [OPTIONAL]
```

**Total:** 30 automated + 3 interactive = **33 tests**

### Test Status Indicators

| Symbol | Meaning |
|--------|---------|
| ✅ `.` | Test passed |
| ❌ `F` | Test failed |
| ⏭️ `s` | Test skipped |
| 🔧 `x` | Expected failure (xfail) |

---

## Test Coverage Details

### 1. Topology Metrics Tests (`test_topology_metrics.py`)

**Self-contained helpers** (test-side implementations):

- `dirichlet_posterior_mean(counts, alpha=0.5)` - Bayesian posterior with Dirichlet prior
- `topological_charge(counts, support)` - Computes charge from polygon distribution
- `js_divergence(p, q, base=2.0)` - Jensen-Shannon divergence
- `wasserstein1_discrete(p, q, support)` - 1-Wasserstein distance on discrete support

**12 tests implemented:**
- ✅ `test_dirichlet_posterior_mean_sums_to_one`
- ✅ `test_dirichlet_posterior_mean_handles_zero_counts`
- ✅ `test_dirichlet_posterior_mean_invariant_to_permutation`
- ✅ `test_dirichlet_posterior_mean_large_alpha_approaches_uniform`
- ✅ `test_topological_charge_zero_for_pure_hexagons`
- ✅ `test_topological_charge_heptagon_rich_negative_sign`
- ✅ `test_topological_charge_empty_counts_returns_zero_tuple`
- ✅ `test_divergence_zero_for_identical_distributions`
- ✅ `test_js_divergence_symmetry_and_scaling_invariance`
- ✅ `test_js_divergence_positive_for_disjoint_support_masses`
- ✅ `test_wasserstein1_discrete_simple_shift`
- ✅ `test_wasserstein1_discrete_nonuniform_support`

**What's tested:**
- Probability normalization (sums to 1)
- Edge cases (zero counts, empty distributions)
- Mathematical properties (symmetry, invariance)
- Hexagon-dominated vs heptagon-rich distributions

---

### 2. Eulerian ROI Tests (`test_roi_eulerian.py`)

**Helper:**
- `assign_to_grid_vectorized(x_coords, y_coords, x_edges, y_edges)` - Maps coordinates to 2D grid bins

**5 tests implemented:**
- ✅ `test_assign_to_grid_vectorized_all_cells_assigned`
- ✅ `test_assign_to_grid_vectorized_empty_input`
- ✅ `test_assign_to_grid_vectorized_points_on_edges`
- ✅ `test_assign_to_grid_vectorized_points_outside_range_are_clipped`
- ✅ `test_assign_to_grid_vectorized_nonuniform_bins`

**What's tested:**
- Standard grid assignment
- Edge points and boundary handling
- Out-of-range coordinate clipping
- Non-uniform bin spacing
- Empty input graceful handling

---

### 3. Lagrangian Tracking Tests (`test_roi_lagrangian.py`)

**Helper:**
- `simple_nearest_neighbor_track(frame0, frame1, max_distance=1.0)` - Nearest-neighbor point matching

**6 tests implemented:**
- ✅ `test_track_cells_preserves_ids_for_small_displacements`
- ✅ `test_track_cells_handles_empty_initial_frame`
- ✅ `test_track_cells_ignores_large_displacements`
- ✅ `test_track_cells_many_to_one_mapping_allowed`
- ✅ `test_track_cells_tie_breaker_is_first_min_index`
- ✅ `test_track_cells_handles_empty_second_frame`

**What's tested:**
- Small displacement tracking
- Large displacement rejection (distance threshold)
- Many-to-one mapping behavior
- Empty frame handling (both initial and subsequent)
- Tie-breaking with `argmin`

---

### 4. I/O Extraction Tests (`test_io_extraction.py`)

**Helpers:**
- `filter_boundary_cells(records)` - Removes records with `cell_id == -1`
- `group_by_frame(records)` - Groups records by frame, preserving order

**7 tests implemented:**
- ✅ `test_boundary_cells_excluded_by_filter_helper`
- ✅ `test_filter_boundary_cells_keeps_all_valid_when_no_boundaries`
- ✅ `test_filter_boundary_cells_all_boundary_yields_empty`
- ✅ `test_empty_sheet_handled_as_empty_result`
- ✅ `test_correct_number_of_cells_per_frame_synthetic`
- ✅ `test_frame_counts_with_missing_frames`
- ✅ `test_records_grouped_by_frame_ordering_is_stable`

**What's tested:**
- Boundary cell (-1 ID) filtering
- Empty data sheets
- Frame counting (contiguous and non-contiguous)
- Stable ordering after grouping

---

### 5. Interactive Visual Tests (`test_visual_plots.py`)

**3 manual validation tests** (marked with `@pytest.mark.interactive`):

1. **`test_visual_histogram_confirmation`**
   - Deterministic polygon distribution histogram
   - Prompt: "Does the histogram look correct and match the expected skew towards hexagons? [y/N]"

2. **`test_interactive_eulerian_grid_overlay`**
   - 2D grid with color-coded cells by grid index
   - Uses `assign_to_grid_vectorized` helper
   - Prompt: "Do the points appear in the expected grid cells with consistent colouring? [y/N]"

3. **`test_interactive_lagrangian_tracking_paths`**
   - Trajectory arrows between two frames
   - Uses `simple_nearest_neighbor_track` helper
   - Prompt: "Do the arrows correctly connect each point to its shifted partner in Frame 1? [y/N]"

**Run with:**
```bash
pytest -m interactive -s
```

---

## Interactive Visual Tests

### Purpose
Interactive tests allow manual validation of matplotlib visualizations that are difficult to test programmatically. They:
- Open matplotlib windows with deterministic synthetic data
- Prompt the user to visually confirm correctness
- Pass only if user types `y` (yes)

### How to Run

```bash
# Must use -s flag to enable input
pytest -m interactive -s

# Run a specific interactive test
pytest tests/test_visual_plots.py::test_visual_histogram_confirmation -s
```

### Expected Behavior

1. Test generates a plot (histogram, grid, or trajectories)
2. Matplotlib window opens
3. Terminal shows a prompt like:
   ```
   Does the histogram look correct? [y/N]: _
   ```
4. User types `y` and presses Enter → test passes
5. User types `n` or anything else → test fails

### Important Notes

⚠️ **Interactive tests will fail without `-s` flag** - This is expected! The `-s` flag allows pytest to capture input from stdin.

⚠️ **Not for CI/CD** - Interactive tests require human input and are excluded from automated test runs using the marker system.

---

## Troubleshooting

### Issue: `pytest: command not found`

**Solution:**
```bash
conda activate data_live_cell
pip install pytest
```

### Issue: Tests from `Icy/` or `MorphoNet/` are collected

**Cause:** pytest is scanning the entire workspace

**Solution:** Verify `pytest.ini` contains:
```ini
[pytest]
testpaths = tests
```

### Issue: Interactive tests fail with "reading from stdin while output is captured"

**Cause:** Missing `-s` flag

**Solution:** Run with `-s` flag:
```bash
pytest -m interactive -s
```

Or exclude interactive tests:
```bash
pytest -m "not interactive"
```

### Issue: Import errors in `test_visual_plots.py`

**Cause:** Relative imports not working

**Solution:** Imports should be absolute (already fixed):
```python
# Correct ✅
from test_roi_eulerian import assign_to_grid_vectorized

# Incorrect ❌
from .test_roi_eulerian import assign_to_grid_vectorized
```

### Issue: matplotlib backend errors on headless systems

**Cause:** No display available for interactive plots

**Solution:** Skip interactive tests on headless systems:
```bash
pytest -m "not interactive"
```

---

## Design Principles

The test suite follows these principles:

1. ✅ **Deterministic** - No randomness in non-interactive tests (or fixed random seeds)
2. ✅ **Isolated** - Tests only in `tests/` directory, no modifications to analysis scripts
3. ✅ **Synthetic** - No dependency on real data files (all use in-memory synthetic data)
4. ✅ **Fast** - Non-interactive tests run in < 1 second total
5. ✅ **Self-contained** - All helpers defined within test files (test-side implementations)
6. ✅ **Documented** - Clear docstrings and comments explaining test intent
7. ✅ **Edge-case coverage** - Empty inputs, boundary conditions, special cases

---

## Contributing Tests

### Guidelines for New Tests

When adding new functionality to the project, follow these testing practices:

#### 1. **All new functions need tests**
- Create tests before or alongside new features
- Aim for at least one test per function
- Cover both typical cases and edge cases

#### 2. **Tests must be deterministic**
- No randomness (unless using fixed `random_seed`)
- Same input → same output every time
- No dependency on system state or external files

#### 3. **Use descriptive test names**
Follow the pattern: `test_<function>_<scenario>_<expected_behavior>`

Examples:
```python
def test_dirichlet_posterior_mean_sums_to_one()
def test_track_cells_handles_empty_initial_frame()
def test_assign_to_grid_vectorized_points_outside_range_are_clipped()
```

#### 4. **Test edge cases**
Common edge cases to consider:
- Empty inputs (`[]`, `np.array([])`)
- Zero values
- Boundary conditions (min/max values)
- Single-element inputs
- Very large inputs (if performance matters)

#### 5. **Use appropriate assertions**
```python
# For exact equality
assert result == expected

# For floating-point comparisons
assert np.isclose(result, expected)
assert np.allclose(result_array, expected_array, atol=1e-6)

# For shape/type checks
assert result.shape == (10, 5)
assert isinstance(result, np.ndarray)
```

#### 6. **Keep tests independent**
- Each test should run successfully in isolation
- Don't rely on side effects from other tests
- Use fresh data for each test

#### 7. **Add tests to appropriate file**
- Topology logic → `test_topology_metrics.py`
- Eulerian ROI → `test_roi_eulerian.py`
- Lagrangian tracking → `test_roi_lagrangian.py`
- I/O operations → `test_io_extraction.py`
- Visual validation → `test_visual_plots.py` (with `@pytest.mark.interactive`)

### Example Test Template

```python
def test_my_function_handles_empty_input():
    """Test that my_function gracefully handles empty arrays."""
    # Arrange
    empty_input = np.array([])
    
    # Act
    result = my_function(empty_input)
    
    # Assert
    assert result.size == 0
    assert isinstance(result, np.ndarray)
```

### Running Tests Before Committing

```bash
# Run all automated tests
pytest -m "not interactive" -v

# Ensure all tests pass
# Expected: 30 passed

# Optionally run interactive tests
pytest -m interactive -s
```

---

## pytest.ini Configuration

Location: `/home/ardati/Data_LiveCell/pytest.ini`

```ini
[pytest]
testpaths = tests
markers =
    interactive: interactive tests requiring human visual confirmation
```

**Configuration details:**
- `testpaths = tests` - Restricts pytest to scan only the `tests/` directory (avoids third-party code in `Icy/`, `MorphoNet/`)
- `interactive` marker - Allows opt-in/opt-out of visual validation tests

---

## File Locations

- **Tests:** `/home/ardati/Data_LiveCell/tests/`
- **Config:** `/home/ardati/Data_LiveCell/pytest.ini`
- **This guide:** `/home/ardati/Data_LiveCell/docs/TESTING.md`

---

## Test Markers Reference

| Marker | Usage | Purpose | Command |
|--------|-------|---------|---------|
| `interactive` | `@pytest.mark.interactive` | Manual visual confirmation | `pytest -m interactive -s` |
| (none) | Default | Automated tests | `pytest` or `pytest -m "not interactive"` |

**Run specific marker:**
```bash
pytest -m <marker_name>
```

**Exclude marker:**
```bash
pytest -m "not <marker_name>"
```

---

## Summary

### Test Suite at a Glance

- **Total Tests:** 33 (30 automated + 3 interactive)
- **Test Status:** ✅ 30/30 passing
- **Runtime:** < 1 second (automated tests)
- **Coverage:** Topology metrics, ROI logic (Eulerian/Lagrangian), I/O extraction
- **Framework:** pytest with custom markers
- **Location:** `tests/` directory (isolated from production code)

### Key Commands Cheat Sheet

```bash
# Activate environment
conda activate data_live_cell

# Run automated tests
pytest -m "not interactive" -v

# Run interactive tests
pytest -m interactive -s

# Run specific test file
pytest tests/test_topology_metrics.py -v

# Stop on first failure
pytest -x

# Show available tests
pytest --collect-only
```

---

**Last Updated:** November 28, 2025  
**Maintainer:** Data_LiveCell Development Team  
**Related Docs:** `PROJECT_DOCUMENTATION_REPORT.md`, `README.md`

