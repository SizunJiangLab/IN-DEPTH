"""
Step 01: Format Raw QPTIFF to OME-TIFF

Converts raw QPTIFF files from the microscope to OME-TIFF format,
separating cycle 0 (reference/DST) from subsequent cycles (SRC).

Edit the CONFIGURATION section below before running.
"""

import os
import numpy as np
import tifffile
from tqdm import tqdm

# Try to import pyqupath, fall back to basic implementation
try:
    from pyqupath.tiff import TiffZarrReader, PyramidWriter
    HAS_PYQUPATH = True
except ImportError:
    HAS_PYQUPATH = False
    print("Warning: pyqupath not available, using basic tifffile implementation")

# =============================================================================
# CONFIGURATION - EDIT THESE PATHS
# =============================================================================

TMA_NAMES = ['Sample1', 'Sample2']

INPUT_QPTIFF = {
    'Sample1': '/path/to/raw/Sample1.qptiff',
    'Sample2': '/path/to/raw/Sample2.qptiff',
}

MARKER_LIST = {
    'Sample1': '/path/to/markers/Sample1_markers.txt',
    'Sample2': '/path/to/markers/Sample2_markers.txt',
}

OUTPUT_DIR = '/path/to/output/01_ometiff'

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_marker_list(filepath):
    """Load marker names from text file."""
    with open(filepath, 'r') as f:
        markers = [line.strip() for line in f if line.strip()]
    return markers


def format_qptiff_basic(input_path, output_dir, tma_name, marker_list):
    """
    Basic QPTIFF to OME-TIFF conversion using tifffile.
    
    Separates cycle 0 (DST) from other cycles (SRC).
    """
    os.makedirs(output_dir, exist_ok=True)
    
    with tifffile.TiffFile(input_path) as tif:
        # Get image shape and metadata
        series = tif.series[0]
        shape = series.shape
        dtype = series.dtype
        
        print(f"  Image shape: {shape}")
        print(f"  Data type: {dtype}")
        
        # Assume shape is (C, Y, X) or (Z, C, Y, X)
        if len(shape) == 3:
            n_channels, height, width = shape
        elif len(shape) == 4:
            n_channels = shape[0] * shape[1]
            height, width = shape[2], shape[3]
        else:
            raise ValueError(f"Unexpected shape: {shape}")
        
        print(f"  Total channels: {n_channels}")
        
        # Read all data
        data = series.asarray()
        if len(data.shape) == 4:
            data = data.reshape(-1, height, width)
        
        # Channels per cycle (DAPI + Cy3 + Cy5 = 3 typically)
        channels_per_cycle = 3
        n_cycles = n_channels // channels_per_cycle
        
        # Cycle 0 is DST (reference)
        dst_data = data[:channels_per_cycle]
        dst_path = os.path.join(output_dir, f'{tma_name}_dst.ome.tiff')
        tifffile.imwrite(dst_path, dst_data, ome=True, photometric='minisblack')
        print(f"  Saved DST: {dst_path}")
        
        # Remaining cycles are SRC
        src_data = data[channels_per_cycle:]
        src_path = os.path.join(output_dir, f'{tma_name}_src.ome.tiff')
        
        # Add channel names if marker list provided
        metadata = {}
        if marker_list:
            metadata['Channel'] = {'Name': marker_list[:len(src_data)]}
        
        tifffile.imwrite(src_path, src_data, ome=True, photometric='minisblack',
                        metadata=metadata)
        print(f"  Saved SRC: {src_path}")


def format_qptiff_pyqupath(input_path, output_dir, tma_name, marker_list):
    """
    QPTIFF to OME-TIFF conversion using pyqupath (if available).
    """
    os.makedirs(output_dir, exist_ok=True)
    
    reader = TiffZarrReader(input_path)
    
    # Get dimensions
    shape = reader.shape
    print(f"  Image shape: {shape}")
    
    # Similar logic but using pyqupath's optimized readers
    # ... (implement based on pyqupath API)
    
    # Fallback to basic for now
    format_qptiff_basic(input_path, output_dir, tma_name, marker_list)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    for tma in TMA_NAMES:
        print(f"\nProcessing {tma}...")
        
        input_path = INPUT_QPTIFF[tma]
        markers = load_marker_list(MARKER_LIST[tma]) if tma in MARKER_LIST else None
        
        if HAS_PYQUPATH:
            format_qptiff_pyqupath(input_path, OUTPUT_DIR, tma, markers)
        else:
            format_qptiff_basic(input_path, OUTPUT_DIR, tma, markers)
    
    print("\nDone!")
