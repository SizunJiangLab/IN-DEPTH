"""I/O utility functions for the TMA pipeline."""

import os
import re
from typing import List, Optional, Tuple


def list_files(directory: str, extension: Optional[str] = None, 
               contains: Optional[str] = None) -> List[str]:
    """
    Recursively list all files in a directory.
    
    Parameters
    ----------
    directory : str
        Root directory to search
    extension : str, optional
        Filter by file extension (e.g., '.tiff')
    contains : str, optional
        Filter by substring in filename
        
    Returns
    -------
    List[str]
        List of file paths
    """
    paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if extension and not file.endswith(extension):
                continue
            if contains and contains not in file:
                continue
            paths.append(os.path.join(root, file))
    return sorted(paths)


def load_marker_list(filepath: str) -> List[str]:
    """
    Load marker names from a text file.
    
    Parameters
    ----------
    filepath : str
        Path to marker list file (one marker per line)
        
    Returns
    -------
    List[str]
        List of marker names
    """
    with open(filepath, 'r') as f:
        markers = [line.strip() for line in f if line.strip()]
    return markers


def get_cycle_number(filename: str) -> int:
    """
    Extract cycle number from filename.
    
    Parameters
    ----------
    filename : str
        Filename containing cycle number (e.g., 'cycle_03.tiff')
        
    Returns
    -------
    int
        Cycle number
    """
    match = re.search(r'cycle[_\s]*(\d+)', filename, re.IGNORECASE)
    if match:
        return int(match.group(1))
    
    # Try alternative patterns
    match = re.search(r'_(\d+)\.', filename)
    if match:
        return int(match.group(1))
    
    return 0


def get_core_name(filepath: str) -> str:
    """
    Extract core name from filepath.
    
    Parameters
    ----------
    filepath : str
        Path to core file
        
    Returns
    -------
    str
        Core name (e.g., 'A1', 'B2')
    """
    basename = os.path.basename(filepath)
    # Match patterns like 'A1', 'B02', 'core_A1'
    match = re.search(r'([A-Z]\d+)', basename, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    return os.path.splitext(basename)[0]


def ensure_dir(path: str) -> str:
    """
    Create directory if it doesn't exist.
    
    Parameters
    ----------
    path : str
        Directory path
        
    Returns
    -------
    str
        The same path (for chaining)
    """
    os.makedirs(path, exist_ok=True)
    return path


def parse_channel_name(channel_str: str) -> Tuple[int, str, int]:
    """
    Parse channel string like '3_Cy5_200' into components.
    
    Parameters
    ----------
    channel_str : str
        Channel string in format 'cycle_filter_exposure'
        
    Returns
    -------
    Tuple[int, str, int]
        (cycle_number, filter_name, exposure_time)
    """
    parts = channel_str.split('_')
    if len(parts) >= 3:
        return int(parts[0]), parts[1], int(parts[2])
    return 0, '', 0
