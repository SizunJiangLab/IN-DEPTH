import os
import json
import logging
import numpy as np
import tifffile
from tqdm import tqdm
from shapely.geometry import shape

# Third-party / Custom libraries found in notebook
from pyqupath.geojson import GeojsonProcessor
from pyqupath.tiff import TiffZarrReader, PyramidWriter

# Setup Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

def stitch_grid_to_image(geojson_path, markerList, tiles_dir, output_path, img_width, img_height):
    """
    Stitches tile images back into a single OME-TIFF.
    Assumes tile filenames correspond to the index of features in the GeoJSON (0.ome.tiff, 1.ome.tiff...).
    """
    
    # ---------------------------------------------------------
    # 1. LOAD FEATURES & PREPARE METADATA
    # ---------------------------------------------------------
    logger.info(f"Loading GeoJSON from {geojson_path}...")
    with open(geojson_path, 'r') as f:
        data = json.load(f)

    features = []
    
    # Iterate using index 'i' to match the file naming convention
    for i, f in enumerate(data['features']):
        poly = shape(f['geometry'])
        bounds = poly.bounds # (minx, miny, maxx, maxy)

        # Round float coordinates to nearest integer for pixel indices
        int_bounds = [int(round(b)) for b in bounds]

        features.append({
            'index': i,
            'filename': f"{i}.ome.tiff",
            'bounds': int_bounds
        })

    logger.info(f"Loaded {len(features)} features. Canvas size set to: {img_width}x{img_height}")

    with open(markerList, 'r') as f:
        MARKERLIST = [line.strip() for line in f if line.strip() != 'blank']

    cycle_dir = [i.path for i in os.scandir(tiles_dir)]
    
    cycle_num = len(cycle_dir)

    # ---------------------------------------------------------
    # 2. INFER IMAGE PROPERTIES FROM FIRST TILE
    # ---------------------------------------------------------
    # We check '0.ome.tiff' specifically
    first_tile_path = os.path.join(cycle_dir[0], "0.ome.tiff")
    
    if not os.path.exists(first_tile_path):
        raise FileNotFoundError(f"First tile not found at {first_tile_path}. Ensure tiles are named by index (e.g., 0.ome.tiff).")

    with tifffile.TiffFile(first_tile_path) as tf:
        sample_series = tf.series[0]
        sample_shape = sample_series.shape
        sample_dtype = sample_series.dtype

    # Determine shape for the full canvas
    # Check if image is (C, Y, X) or (Y, X)
    if len(sample_shape) == 3:
        n_channels = int(sample_shape[0] / 2 * cycle_num + 3) # this is a one-off. Since the dst is cycle 0 and contains dapi, cy3-first-blank, cy5-first-blank. 
        logger.info(f'Channels to be created: {n_channels}...')
        if n_channels != len(MARKERLIST):
            raise ValueError(f"Number of channels to be created, {n_channels}, is not equal to the number of entries, {len(MARKERLIST)}, in MarkerList.txt without blank.")
        mmap_shape = (n_channels, img_height, img_width)
        logger.info(f"Output shape: {mmap_shape}")
    else:
        n_channels = 1
        mmap_shape = (cycle_num, img_height, img_width)
        logger.info(f"Output shape: {mmap_shape}")

    # ---------------------------------------------------------
    # 3. CREATE MEMMAP & STITCH
    # ---------------------------------------------------------
    
    # Create a file-backed array to avoid RAM overflow
    stitched_arr = np.zeros(mmap_shape, dtype=sample_dtype)

    logger.info("Starting stitch process...")
    for cycle in cycle_dir:
        cycle_idx = int(os.path.basename(cycle)) - 1
        logger.info(f"Processing cycle {cycle_idx}...")
        
        for feat in tqdm(features, desc="Pasting tiles"):
            tile_filename = feat['filename']
            logger.info(f'Pasting {tile_filename}...')
            x_min, y_min, x_max, y_max = feat['bounds']
            
            tile_path = os.path.join(cycle, tile_filename)
            
            if not os.path.exists(tile_path):
                logger.warning(f"Tile {tile_filename} missing. Skipping.")
                continue
                
            # Read tile
            tile_data = tifffile.imread(tile_path)
            
            # Calculate expected slot dimensions based on GeoJSON
            slot_h = y_max - y_min
            slot_w = x_max - x_min
            
            # --- Safety Checks for Bounds ---

            if cycle_idx == 0:

            
                # 1. Check if the GeoJSON box exceeds the manual canvas size
                if x_max > img_width or y_max > img_height:
                    target_x_max = min(x_max, img_width)
                    target_y_max = min(y_max, img_height)
                    
                    # Recalculate how much of the tile we can actually paste
                    paste_w = target_x_max - x_min
                    paste_h = target_y_max - y_min
                    
                    if paste_w <= 0 or paste_h <= 0:
                        continue

                    # Slice the source tile to match the smaller paste area
                    if n_channels > 1:
                        idx_start = 0 
                        idx_end = 6
                        stitched_arr[idx_start:idx_end, y_min:target_y_max, x_min:target_x_max] = tile_data[[0,2,3,1,4,5], :paste_h, :paste_w]
                    else:
                        stitched_arr[cycle_idx, y_min:target_y_max, x_min:target_x_max] = tile_data[:paste_h, :paste_w]
                        
                # 2. Standard Paste
                else:
                    # We slice tile_data[:slot_h, :slot_w] to handle cases where 
                    # the tile image might be slightly larger than the GeoJSON box due to compression blocks
                    if n_channels > 1:
                        idx_start = 0
                        idx_end = 6
                        stitched_arr[idx_start:idx_end, y_min:y_max, x_min:x_max] = tile_data[[0,2,3,1,4,5], :slot_h, :slot_w]
                    else:
                        stitched_arr[cycle_idx, y_min:y_max, x_min:x_max] = tile_data[:slot_h, :slot_w]
            else: 
                # 1. Check if the GeoJSON box exceeds the manual canvas size
                if x_max > img_width or y_max > img_height:
                    target_x_max = min(x_max, img_width)
                    target_y_max = min(y_max, img_height)
                    
                    # Recalculate how much of the tile we can actually paste
                    paste_w = target_x_max - x_min
                    paste_h = target_y_max - y_min
                    
                    if paste_w <= 0 or paste_h <= 0:
                        continue

                    # Slice the source tile to match the smaller paste area
                    if n_channels > 1:
                        idx_start = (cycle_idx + 1) * 3 
                        idx_end = idx_start + 3
                        stitched_arr[idx_start:idx_end, y_min:target_y_max, x_min:target_x_max] = tile_data[[1,4,5], :paste_h, :paste_w]
                    else:
                        stitched_arr[cycle_idx, y_min:target_y_max, x_min:target_x_max] = tile_data[:paste_h, :paste_w]
                        
                # 2. Standard Paste
                else:
                    # We slice tile_data[:slot_h, :slot_w] to handle cases where 
                    # the tile image might be slightly larger than the GeoJSON box due to compression blocks
                    if n_channels > 1:
                        idx_start = (cycle_idx + 1) * 3 
                        idx_end = idx_start + 3
                        stitched_arr[idx_start:idx_end, y_min:y_max, x_min:x_max] = tile_data[[1,4,5], :slot_h, :slot_w]
                    else:
                        stitched_arr[cycle_idx, y_min:y_max, x_min:x_max] = tile_data[:slot_h, :slot_w]


    # ---------------------------------------------------------
    # 4. SAVE FINAL OUTPUT
    # ---------------------------------------------------------
        
        
    logger.info("Saving final OME-TIFF...")
    writer = PyramidWriter.from_array(stitched_arr, MARKERLIST)
    writer.export_ometiff_pyramid(output_path)
    # Clean up
    del stitched_arr
        
    logger.info(f"Success! Image saved to {output_path}")


# =========================================================
# CONFIGURATION
# =========================================================

if __name__ == "__main__":
    
    # 1. Dimensions (From your previous notebook)
    ORIGINAL_WIDTH  = 29760
    ORIGINAL_HEIGHT = 31680
    # 2. Paths
    GEOJSON_INPUT = "/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/02.5_crop_cores/Tub97.2_grid_boxes.geojson"
    TILES_FOLDER  = "/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/03_warping/Tub97.2"
    OUTPUT_FILE   = "/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/04_stitched_cycles/Tub972.ome.tiff"
    MARKER = "/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/raw/rawMarkerList.txt"
    # 3. Execution
    try:
        stitch_grid_to_image(
            geojson_path=GEOJSON_INPUT, 
            markerList= MARKER,
            tiles_dir=TILES_FOLDER, 
            output_path=OUTPUT_FILE,
            img_width=ORIGINAL_WIDTH,
            img_height=ORIGINAL_HEIGHT
        )
    except KeyboardInterrupt:
        logger.warning("Process interrupted by user.")
    except Exception as e:
        logger.exception("An error occurred during stitching.")