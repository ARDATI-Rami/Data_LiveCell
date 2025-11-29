"""
XLS/Excel file reading utilities for cell segmentation data.

This module provides functions to read cell tracking data from EpiTools-generated
XLS workbooks. Each sheet represents one time frame with columns:
- Column 0: Cell ID (-1 for boundary cells)
- Column 1: X coordinate
- Column 2: Y coordinate
- Column 3: Polygon number (cell shape class)
"""

import xlrd
import numpy as np
from typing import Tuple, List
import os


def load_workbook(xls_path: str) -> xlrd.book.Book:
    """
    Open an XLS workbook.

    :param xls_path: Path to the .xls file (can be absolute or relative).
    :type xls_path: str
    :return: Opened workbook.
    :rtype: xlrd.book.Book
    :raises FileNotFoundError: If the file does not exist.
    :raises xlrd.XLRDError: If the file cannot be opened.
    """
    if not os.path.exists(xls_path):
        raise FileNotFoundError(f"XLS file not found: {xls_path}")
    return xlrd.open_workbook(xls_path)


def extract_data(workbook: xlrd.book.Book, sheet_index: int,
                 skip_boundary: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract cell data from a single sheet in the workbook.

    :param workbook: Opened xlrd workbook.
    :type workbook: xlrd.book.Book
    :param sheet_index: Sheet index to read (0-based).
    :type sheet_index: int
    :param skip_boundary: If True, skip rows where cell ID == -1.
    :type skip_boundary: bool
    :return: Tuple of (x_coords, y_coords, polygon_numbers) as numpy arrays.
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray]
    """
    sheet = workbook.sheet_by_index(sheet_index)

    x_coords = []
    y_coords = []
    polygon_numbers = []

    # Assuming data starts from row 1 (row 0 is header)
    for row_idx in range(1, sheet.nrows):
        cell_id = sheet.cell_value(row_idx, 0)
        
        # Skip boundary cells if requested
        if skip_boundary and cell_id == -1:
            continue

        x = sheet.cell_value(row_idx, 1)
        y = sheet.cell_value(row_idx, 2)
        polygon_no = sheet.cell_value(row_idx, 3)

        x_coords.append(float(x))
        y_coords.append(float(y))
        polygon_numbers.append(int(polygon_no))

    return np.array(x_coords), np.array(y_coords), np.array(polygon_numbers)


def extract_frame(sheet: xlrd.sheet.Sheet,
                  skip_boundary: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract a complete frame from an xlrd sheet (including cell IDs).

    :param sheet: xlrd sheet object.
    :type sheet: xlrd.sheet.Sheet
    :param skip_boundary: If True, skip rows where cell ID == -1.
    :type skip_boundary: bool
    :return: Tuple of (ids, x_coords, y_coords, polygon_numbers).
             ids may contain -1 for boundary cells if skip_boundary=False.
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
    """
    ids, xs, ys, polys = [], [], [], []
    
    for row_idx in range(1, sheet.nrows):
        cell_id = sheet.cell_value(row_idx, 0)
        
        if skip_boundary and cell_id == -1:
            continue
        
        cid = int(cell_id)
        x = float(sheet.cell_value(row_idx, 1))
        y = float(sheet.cell_value(row_idx, 2))
        pg = int(sheet.cell_value(row_idx, 3))
        
        ids.append(cid)
        xs.append(x)
        ys.append(y)
        polys.append(pg)
    
    return np.array(ids), np.array(xs), np.array(ys), np.array(polys)


def read_all_frames(workbook: xlrd.book.Book,
                    skip_boundary: bool = False) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Read all frames from a workbook.

    :param workbook: Opened xlrd workbook.
    :type workbook: xlrd.book.Book
    :param skip_boundary: If True, skip boundary cells (ID == -1) in all frames.
    :type skip_boundary: bool
    :return: List of (ids, x, y, polygons) tuples, one per sheet/frame.
    :rtype: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]
    """
    frames = []
    for sheet_index in range(workbook.nsheets):
        sheet = workbook.sheet_by_index(sheet_index)
        frames.append(extract_frame(sheet, skip_boundary=skip_boundary))
    return frames


def get_frame_count(workbook: xlrd.book.Book) -> int:
    """
    Get the number of frames (sheets) in a workbook.

    :param workbook: Opened xlrd workbook.
    :type workbook: xlrd.book.Book
    :return: Number of sheets.
    :rtype: int
    """
    return workbook.nsheets

