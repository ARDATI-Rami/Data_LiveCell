"""
ODS (OpenDocument Spreadsheet) file reading utilities.

This module provides functions to read data from ODS files,
primarily for diagnostic analysis of neighbor distributions.
"""

import ezodf
import numpy as np
from typing import List, Tuple, Optional


def load_ods(ods_path: str) -> ezodf.document.PackagedDocument:
    """
    Open an ODS file.

    :param ods_path: Path to the .ods file.
    :type ods_path: str
    :return: Opened ODS document.
    :rtype: ezodf.document.PackagedDocument
    :raises FileNotFoundError: If the file does not exist.
    """
    return ezodf.opendoc(ods_path)


def extract_column_data(doc: ezodf.document.PackagedDocument,
                       sheet_index: int,
                       column_name: str,
                       skip_header: bool = True) -> List:
    """
    Extract data from a specific column in an ODS sheet.

    :param doc: Opened ODS document.
    :type doc: ezodf.document.PackagedDocument
    :param sheet_index: Index of the sheet to read (0-based).
    :type sheet_index: int
    :param column_name: Name of the column to extract (must match header).
    :type column_name: str
    :param skip_header: If True, skip the first row (header).
    :type skip_header: bool
    :return: List of values from the specified column.
    :rtype: List
    :raises ValueError: If column name not found in header.
    """
    sheet = doc.sheets[sheet_index]
    rows = list(sheet.rows())
    
    # Get column headers
    if skip_header and len(rows) > 0:
        headers = [cell.value for cell in rows[0]]
        
        if column_name not in headers:
            raise ValueError(f"Column '{column_name}' not found. Available: {headers}")
        
        column_index = headers.index(column_name)
        
        # Extract data from subsequent rows
        data = []
        for row in rows[1:]:
            if column_index < len(row):
                data.append(row[column_index].value)
        
        return data
    else:
        raise ValueError("Sheet is empty or skip_header=False not supported without column index")


def extract_neighbor_distribution(doc: ezodf.document.PackagedDocument,
                                 sheet_index: int = 0,
                                 label_column: str = "label",
                                 neighbors_column: str = "num_neighbours") -> Tuple[List, List]:
    """
    Extract cell IDs and neighbor counts from an ODS file.

    :param doc: Opened ODS document.
    :type doc: ezodf.document.PackagedDocument
    :param sheet_index: Index of the sheet to read (default 0).
    :type sheet_index: int
    :param label_column: Name of the cell ID/label column.
    :type label_column: str
    :param neighbors_column: Name of the neighbor count column.
    :type neighbors_column: str
    :return: Tuple of (cell_ids, neighbor_counts) as lists.
    :rtype: Tuple[List, List]
    """
    cell_ids = extract_column_data(doc, sheet_index, label_column)
    neighbors = extract_column_data(doc, sheet_index, neighbors_column)
    
    return cell_ids, neighbors


def compute_neighbor_statistics(neighbors: List[int],
                               max_neighbors: int = 13) -> dict:
    """
    Compute statistics about neighbor distribution.

    :param neighbors: List of neighbor counts per cell.
    :type neighbors: List[int]
    :param max_neighbors: Maximum number of neighbors to report.
    :type max_neighbors: int
    :return: Dictionary with statistics (counts, percentages per neighbor count).
    :rtype: dict
    """
    total_cells = len(neighbors)
    
    stats = {
        'total_cells': total_cells,
        'distribution': {}
    }
    
    for n in range(1, max_neighbors + 1):
        count = neighbors.count(n)
        percentage = (count / total_cells * 100) if total_cells > 0 else 0.0
        stats['distribution'][n] = {
            'count': count,
            'percentage': percentage
        }
    
    return stats

