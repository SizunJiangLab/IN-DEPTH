import os
import numpy as np
from tqdm import tqdm
from pyqupath.geojson import GeojsonProcessor
from pyqupath.tiff import TiffZarrReader, PyramidWriter
import gc

# --- Configuration ---
# Input Directories
TMA = ['Tub971', 'Tub972']
# MARKER_DIR = '../data/02_AF_corrected_markers_DFCI/marker_array'
# OMEGA_DIR = '../data/02_AF_corrected_markers_DFCI/omega'
# DAPI_IMAGE_PATH = '../data/image_registration/03_warping_DFCI/DFCI.ome.tiff'
# GEOJSON_PATH = '../data/02.5_crop_cores/dearrayer_DFCI_correct.geojson'

# # Output Directories
# OUT_CROP_IMG_DIR = '../data/02.5_crop_cores/cropped_arrays_with_dapi_DFCI'
# OUT_CROP_OMEGA_DIR = '../data/02.5_crop_cores/cropped_omega_DFCI'
# OUT_OMETIFF_DIR = '../data/02.5_crop_cores/AFcorr_ometiff_DFCI'

# # Output Metadata Lists
# MARKER_LIST_FILE = '../data/02.5_crop_cores/markerList_with_dapi.txt'
# OMEGA_LIST_FILE = '../data/02.5_crop_cores/omegaList_DFCI.txt'

# # Create output dirs
# os.makedirs(OUT_CROP_IMG_DIR, exist_ok=True)
# os.makedirs(OUT_CROP_OMEGA_DIR, exist_ok=True)
# os.makedirs(OUT_OMETIFF_DIR, exist_ok=True)
# os.makedirs(os.path.dirname(MARKER_LIST_FILE), exist_ok=True)

def list_files(directory: str):
    paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            paths.append(os.path.join(root, file))
    return paths

def main():
    for tma in TMA:

        print(f'Processing {tma}...')

        # ---------------------------------------------------------
        # STEP 0: Make directories
        # ---------------------------------------------------------

        MARKER_DIR = f'../data/05_AF_corrected_markers/{tma}/marker_array'
        OMEGA_DIR = f'../data/05_AF_corrected_markers/{tma}/omega'
        DAPI_IMAGE_PATH = f'../data/04_stitched_cycles/{tma}.ome.tiff'
        GEOJSON_PATH = f'../data/02.5_crop_cores/dearrayer_{tma}_correct.geojson'

        # Output Directories
        OUT_CROP_IMG_DIR = f'../data/06_crop_cores/{tma}/cropped_arrays_with_dapi'
        OUT_CROP_OMEGA_DIR = f'../data/06_crop_cores/{tma}/cropped_omega'
        OUT_OMETIFF_DIR = f'../data/06_crop_cores/{tma}/AFcorr_ometiff'

        # Output Metadata Lists
        MARKER_LIST_FILE = f'../data/06_crop_cores/{tma}/markerList_with_dapi.txt'
        OMEGA_LIST_FILE = f'../data/06_crop_cores/{tma}/omegaList.txt'

        # Create output dirs
        os.makedirs(OUT_CROP_IMG_DIR, exist_ok=True)
        os.makedirs(OUT_CROP_OMEGA_DIR, exist_ok=True)
        os.makedirs(OUT_OMETIFF_DIR, exist_ok=True)

        # ---------------------------------------------------------
        # STEP 1: Prepare GeoJSON Processor
        # ---------------------------------------------------------
        print(f"Loading GeoJSON from {GEOJSON_PATH}...")
        geojson_processor = GeojsonProcessor.from_path(GEOJSON_PATH)

        # ---------------------------------------------------------
        # STEP 2: Process MARKERS + DAPI (Sequential Memmap)
        # ---------------------------------------------------------
        print("\n--- Phase 1: Processing Image Channels ---")
        
        # 1. Build List
        marker_files = [f for f in list_files(MARKER_DIR) if f.endswith('.npy')]
        marker_names = [os.path.basename(f).split('.npy')[0] for f in marker_files]
        
        try:
            tiff_reader = TiffZarrReader.from_ometiff(DAPI_IMAGE_PATH)
            dapi_names = [k for k in tiff_reader.zimg_dict.keys() if 'DAPI' in k]
        except Exception as e:
            print(f"Error reading DAPI metadata: {e}")
            dapi_names = []

        full_channel_list = marker_names + dapi_names
        total_channels = len(full_channel_list)
        
        # Save List
        with open(MARKER_LIST_FILE, 'w') as f:
            for name in full_channel_list:
                f.write(f'{name}\n')

        # 2. Loop Channels
        for channel_idx, channel_name in enumerate(tqdm(full_channel_list, desc="Image Channels")):
            
            # Load single channel into memory
            current_image = None
            if channel_name in marker_names:
                path = next(f for f in marker_files if channel_name in f)
                current_image = np.load(path).astype('float32')
            elif channel_name in dapi_names:
                tiff = TiffZarrReader.from_ometiff(DAPI_IMAGE_PATH)
                current_image = tiff.zimg_dict[channel_name][:].astype('float32')

            if current_image is None: continue

            # Crop for all cores
            single_dict = {channel_name: current_image}
            
            for core_name, cropped_dict in geojson_processor.crop_dict_by_polygons(single_dict):
                crop = cropped_dict[channel_name]
                h, w = crop.shape
                out_path = f'{OUT_CROP_IMG_DIR}/{core_name}.npy'

                # Write to disk (Memmap)
                if channel_idx == 0:
                    # First pass: Create file
                    mm = np.lib.format.open_memmap(
                        out_path, mode='w+', dtype='float32', shape=(total_channels, h, w)
                    )
                    mm[channel_idx] = crop
                    mm.flush()
                else:
                    # Subsequent passes: Open and insert
                    if os.path.exists(out_path):
                        mm = np.lib.format.open_memmap(out_path, mode='r+')
                        mm[channel_idx] = crop
                        mm.flush()
                del mm

            # Clean Memory
            del current_image
            del single_dict
            gc.collect()

        # ---------------------------------------------------------
        # STEP 3: Process OMEGA MASKS (Sequential Memmap)
        # ---------------------------------------------------------
        print("\n--- Phase 2: Processing Omega Masks ---")

        # 1. Build List
        omega_files = [f for f in list_files(OMEGA_DIR) if f.endswith('.npy')]
        omega_names = [os.path.basename(f).split('.npy')[0] for f in omega_files]
        total_omegas = len(omega_names)

        # Save List
        with open(OMEGA_LIST_FILE, 'w') as f:
            for name in omega_names:
                f.write(f'{name}\n')

        # 2. Loop Masks
        for omega_idx, omega_name in enumerate(tqdm(omega_names, desc="Omega Masks")):
            
            # Load single mask
            path = next(f for f in omega_files if omega_name in f)
            current_mask = np.load(path) # Keeps original dtype (likely uint16 or bool)
            
            single_dict = {omega_name: current_mask}

            for core_name, cropped_dict in geojson_processor.crop_dict_by_polygons(single_dict):
                crop = cropped_dict[omega_name]
                h, w = crop.shape
                out_path = f'{OUT_CROP_OMEGA_DIR}/{core_name}.npy'

                # Write to disk (Memmap)
                if omega_idx == 0:
                    # First pass: Create file
                    # NOTE: We assume all omega masks share the dtype of the first one
                    mm = np.lib.format.open_memmap(
                        out_path, mode='w+', dtype=crop.dtype, shape=(total_omegas, h, w)
                    )
                    mm[omega_idx] = crop
                    mm.flush()
                else:
                    # Subsequent passes
                    if os.path.exists(out_path):
                        mm = np.lib.format.open_memmap(out_path, mode='r+')
                        mm[omega_idx] = crop
                        mm.flush()
                del mm

            # Clean Memory
            del current_mask
            del single_dict
            gc.collect()

        print("\nAll Processing Complete.")

if __name__ == "__main__":
    main()