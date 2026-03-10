"""
Step 02: Generate Core GeoJSON

Automated TMA dearraying - detects circular TMA cores from DAPI channel
and outputs their boundaries as GeoJSON for downstream cropping.

Edit the CONFIGURATION section below before running.

Note: After running, manually inspect the GeoJSON and correct if needed.
"""

import os
import numpy as np
import tifffile
import cv2 as cv
from tqdm import tqdm

# Try to import pyqupath for optimized dearraying
try:
    from pyqupath.tiff import TiffZarrReader, PyramidWriter
    from pyqupath.tma import plot_contours, tma_dearrayer
    HAS_PYQUPATH = True
except ImportError:
    HAS_PYQUPATH = False
    print("Warning: pyqupath not available, using OpenCV-based fallback")

try:
    import geopandas as gpd
    from shapely.geometry import Polygon
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False
    print("Warning: geopandas not available, will output JSON instead")

# =============================================================================
# CONFIGURATION - EDIT THESE PATHS
# =============================================================================

TMA_NAMES = ['Sample1', 'Sample2']

# Path to DST OME-TIFF (output from step 01) containing DAPI
INPUT_DST = {
    'Sample1': '/path/to/output/01_ometiff/Sample1_dst.ome.tiff',
    'Sample2': '/path/to/output/01_ometiff/Sample2_dst.ome.tiff',
}

OUTPUT_DIR = '/path/to/output/02_geojson'

# Dearraying parameters
DOWNSAMPLE_STEP = 8      # Downsampling for faster processing
KERNEL_RADIUS = 5        # Morphological kernel radius
AREA_THRESHOLD = 1000    # Minimum core area (in downsampled pixels)
RADIUS_EXPAND = 0        # Expand detected radius by this amount

# =============================================================================
# FALLBACK IMPLEMENTATION (without pyqupath)
# =============================================================================

def detect_cores_opencv(dapi_image, downsample=8, min_area=1000, circularity_thresh=0.7):
    """
    Detect TMA cores using OpenCV circle detection.
    
    Parameters
    ----------
    dapi_image : np.ndarray
        DAPI channel image
    downsample : int
        Downsampling factor for faster processing
    min_area : int
        Minimum contour area
    circularity_thresh : float
        Minimum circularity (4*pi*area/perimeter^2)
        
    Returns
    -------
    list
        List of (center_x, center_y, radius) tuples
    """
    # Downsample
    h, w = dapi_image.shape[:2]
    small = cv.resize(dapi_image, (w // downsample, h // downsample))
    
    # Normalize to 8-bit
    if small.dtype == np.uint16:
        small = (small / 256).astype(np.uint8)
    elif small.dtype != np.uint8:
        small = ((small - small.min()) / (small.max() - small.min()) * 255).astype(np.uint8)
    
    # Threshold
    _, binary = cv.threshold(small, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
    
    # Morphological operations
    kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (5, 5))
    binary = cv.morphologyEx(binary, cv.MORPH_CLOSE, kernel)
    binary = cv.morphologyEx(binary, cv.MORPH_OPEN, kernel)
    
    # Find contours
    contours, _ = cv.findContours(binary, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
    
    cores = []
    for cnt in contours:
        area = cv.contourArea(cnt)
        if area < min_area:
            continue
        
        perimeter = cv.arcLength(cnt, True)
        if perimeter == 0:
            continue
            
        circularity = 4 * np.pi * area / (perimeter ** 2)
        if circularity < circularity_thresh:
            continue
        
        # Fit minimum enclosing circle
        (cx, cy), radius = cv.minEnclosingCircle(cnt)
        
        # Scale back to original coordinates
        cx_orig = cx * downsample
        cy_orig = cy * downsample
        radius_orig = radius * downsample
        
        cores.append((cx_orig, cy_orig, radius_orig))
    
    return cores


def cores_to_geojson(cores, output_path, n_points=64):
    """
    Convert detected cores to GeoJSON format.
    
    Parameters
    ----------
    cores : list
        List of (center_x, center_y, radius) tuples
    output_path : str
        Path to save GeoJSON
    n_points : int
        Number of points for circle approximation
    """
    import json
    
    features = []
    for i, (cx, cy, r) in enumerate(cores):
        # Generate circle points
        angles = np.linspace(0, 2*np.pi, n_points, endpoint=False)
        coords = [(cx + r * np.cos(a), cy + r * np.sin(a)) for a in angles]
        coords.append(coords[0])  # Close the polygon
        
        feature = {
            "type": "Feature",
            "properties": {
                "name": f"core_{i:02d}",
                "objectType": "annotation"
            },
            "geometry": {
                "type": "Polygon",
                "coordinates": [coords]
            }
        }
        features.append(feature)
    
    geojson = {
        "type": "FeatureCollection",
        "features": features
    }
    
    with open(output_path, 'w') as f:
        json.dump(geojson, f, indent=2)


# =============================================================================
# MAIN
# =============================================================================

def process_tma(tma_name, input_path, output_dir):
    """Process a single TMA to generate core GeoJSON."""
    
    print(f"\nProcessing {tma_name}...")
    
    # Load DAPI channel
    if HAS_PYQUPATH:
        reader = TiffZarrReader.from_ometiff(input_path)
        # Try to get DAPI channel
        if 'DAPI_0' in reader.zimg_dict:
            dapi = reader.zimg_dict['DAPI_0'][:]
        elif 'DAPI' in reader.zimg_dict:
            dapi = reader.zimg_dict['DAPI'][:]
        else:
            # Assume first channel is DAPI
            dapi = reader.zimg_dict[list(reader.zimg_dict.keys())[0]][:]
    else:
        with tifffile.TiffFile(input_path) as tif:
            # Assume first channel is DAPI
            dapi = tif.asarray()[0] if len(tif.asarray().shape) > 2 else tif.asarray()
    
    print(f"  DAPI shape: {dapi.shape}")
    
    # Detect cores
    if HAS_PYQUPATH:
        result = tma_dearrayer(
            dapi,
            downsample_step=DOWNSAMPLE_STEP,
            kernel_radius=KERNEL_RADIUS,
            area_threshold=AREA_THRESHOLD,
            radius_expand=RADIUS_EXPAND,
            merge_list=None
        )
        gdf = result[-1]  # Last element is the GeoDataFrame
        
        # Save
        output_path = os.path.join(output_dir, f'{tma_name}_cores.geojson')
        gdf.to_file(output_path, driver='GeoJSON')
        print(f"  Saved: {output_path}")
        
        # Optional: save visualization
        try:
            fig, ax = plot_contours(gdf, merge_list=result[1])
            fig.savefig(os.path.join(output_dir, f'{tma_name}_cores_preview.png'), dpi=150)
            print(f"  Saved preview image")
        except Exception as e:
            print(f"  Could not save preview: {e}")
    else:
        # Fallback to OpenCV
        cores = detect_cores_opencv(
            dapi, 
            downsample=DOWNSAMPLE_STEP,
            min_area=AREA_THRESHOLD
        )
        print(f"  Detected {len(cores)} cores")
        
        output_path = os.path.join(output_dir, f'{tma_name}_cores.geojson')
        cores_to_geojson(cores, output_path)
        print(f"  Saved: {output_path}")


if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    for tma in TMA_NAMES:
        if tma in INPUT_DST:
            process_tma(tma, INPUT_DST[tma], OUTPUT_DIR)
        else:
            print(f"Warning: No input path configured for {tma}")
    
    print("\nDone!")
    print("\nIMPORTANT: Manually inspect the generated GeoJSON files and correct if needed.")
