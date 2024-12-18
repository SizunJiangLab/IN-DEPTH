import cv2 
import tifffile 
import numpy as np
import math
import pandas as pd
import json

from typing import List, Dict

import logging 

import skimage 
import skimage.io
import skimage.measure
import skimage.morphology
from scipy.io import loadmat

import os 


import skimage.io as io
from datetime import datetime
from matplotlib import pyplot as plt
date = datetime.today().strftime('%Y-%m-%d')

from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')



# Constants

PATH_TO_MARKER = '../data/intermediate_channels_sub_cycle15_processed'
PATH_TO_SEGMENTATION = '../output/seg_results'
OUTPUT_DIR = '../output'
MARKERS_TO_EXTRACT = ['PD-1', 'GZMB', 'Ki67', 'CD45RA', 'CD45RO', 'LAG3', 'HLA1', 'HLA-DR',
                      'DAPI_cycle4', 'DAPI_cycle6', 'DAPI_cycle7', 'DAPI_cycle9', 'DAPI_cycle10', 'DAPI_cycle11']
CHANNELS = ['PD-1', 'GZMB', 'Ki67', 'CD45RA', 'CD45RO', 'LAG3', 'HLA1', 'HLA-DR']
LOWERBOUND_DF = pd.read_csv('../data/threshold_csv/lowerBound.csv')
UPPERBOUND_DF = pd.read_csv('../data/threshold_csv/upperBound.csv')

SELECTED_CORES = ["Rochester_4", "Rochester_6",
                    "Rochester_7", "Rochester_9", "Rochester_11", "Rochester_12",
                    "Rochester_13", "Rochester_14",
                    "Rochester_15", "Rochester_16", "Rochester_17", "Rochester_18",
                    "Rochester_19", "Rochester_21", "Rochester_23",
                    "Rochester_25", "DFCI_2.2", "DFCI_3.2",
                    "DFCI_4.1", "DFCI_7.1", "DFCI_8.1",
                    "DFCI_12.1", "DFCI_13.2", "DFCI_14.1", "DFCI_15.2", "DFCI_17.1",
                    "DFCI_18.2", "DFCI_19.2", "DFCI_22.2", "DFCI_23.2"]



#MARKERS_TO_EXTRACT = [i + '_filtered' for i in MARKERS_TO_EXTRACT]

def list_files(directory: str):
    paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            paths.append(os.path.join(root, file))
    return paths

def organize_core_data(tma:str, path: str) -> Dict:

    core_dict = {}

    core_path = [i.path for i in os.scandir(f'{PATH_TO_MARKER}/{tma}')]

    for path in tqdm(core_path):
        core_name = path.split('/')[-1]
        if not f'{tma}_{core_name}' in SELECTED_CORES:
            logging.info(f'Skipping {tma}_{core_name}. Not selected')
            continue
        else: 
            core_dict[core_name] = {}
            markers_path = list_files(path)
            for marker in MARKERS_TO_EXTRACT:
                marker_path = [j for j in markers_path if f'{marker}.tiff' in j][0]
                core_dict[core_name][marker] = tifffile.imread(marker_path)

    return core_dict

def extract_single_cell_info(core_img: Dict[str, np.ndarray], segmentation_mask: np.ndarray, 
                             interested_markers: List[str], output_path: str):
    """Extract single cell information from a core."""
    array_list = [core_img[channel] for channel in interested_markers]
    counts_no_noise = np.stack(array_list, axis=2)

    stats = skimage.measure.regionprops(segmentation_mask)
    label_num = len(stats)

    channel_num = len(array_list)
    data = np.zeros((label_num, channel_num))
    data_scale_size = np.zeros((label_num, channel_num))
    cell_sizes = np.zeros((label_num, 1))
    cell_props = np.zeros((label_num, 3))

    for i, region in enumerate(stats):
        cell_label = region.label
        label_counts = [counts_no_noise[coord[0], coord[1], :] for coord in region.coords]
        data[i] = np.sum(label_counts, axis=0)
        data_scale_size[i] = data[i] / region.area
        cell_sizes[i] = region.area
        cell_props[i] = [cell_label, region.centroid[0], region.centroid[1]]

    col_names = [marker for marker in interested_markers if marker != 'Empty']

    data_df = pd.DataFrame(data, columns=col_names)
    data_full = pd.concat([
        pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
        pd.DataFrame(cell_sizes, columns=["cellSize"]),
        data_df
    ], axis=1)

    data_scale_size_df = pd.DataFrame(data_scale_size, columns=col_names)
    data_scale_size_full = pd.concat([
        pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
        pd.DataFrame(cell_sizes, columns=["cellSize"]),
        data_scale_size_df
    ], axis=1)

    os.makedirs(output_path, exist_ok=True)
    data_full.to_csv(os.path.join(output_path, "data.csv"), index=False)
    data_scale_size_full.to_csv(os.path.join(output_path, "dataScaleSize.csv"), index=False)



def gate_markers(core_img: Dict[str, np.ndarray], coreName, channels, lowerBound_df, upperBound_df):

    for marker in channels:
        channel_img = core_img[marker]
        lowerBound = lowerBound_df[lowerBound_df['Core'] == coreName][marker].values.item()
        upperBound = upperBound_df[upperBound_df['Core'] == coreName][marker].values.item()
        if not upperBound == 0:
            img_filtered = np.copy(channel_img)
            img_filtered[img_filtered <= lowerBound] = 0
            img_filtered[img_filtered >= upperBound] = 0
        else:
            img_filtered = np.copy(channel_img)
            img_filtered[img_filtered <= lowerBound] = 0
        core_img[marker] = img_filtered
    return core_img



def process_tma(tma: str):

    """Process a single TMA."""

    logging.info(f'Processing {tma}')

    core_dict = organize_core_data(tma, f'{PATH_TO_MARKER}/{tma}')
    logging.info('Core data organized')

    for core_num, core_img in tqdm(core_dict.items(), desc="Processing cores"):
        # Gate markers
        core_img_filtered = gate_markers(core_img, f'{tma}_{core_num}', CHANNELS, LOWERBOUND_DF, UPPERBOUND_DF)

        # Load segmentation mask
        seg_mask_path = f'{PATH_TO_SEGMENTATION}/{tma}/core_{core_num}/MESMER_mask.tiff'
        try:
            segmentation_mask = tifffile.imread(seg_mask_path)
        except FileNotFoundError:
            logging.warning(f"Segmentation mask not found for core {core_num}. Skipping.")
            continue

        # Extract single cell information
        info_output_path = f'{OUTPUT_DIR}/extracted_info_16bit_sub_cycle15_processed_gated/{tma}/{core_num}'
        extract_single_cell_info(core_img_filtered, segmentation_mask, MARKERS_TO_EXTRACT, info_output_path)
        logging.info(f'Single cell information extracted for core {core_num}')

            
def main():
    tma_to_process = ['DFCI', 'Rochester']

    for tma in tma_to_process:
        try:
            process_tma(tma)
        except Exception as e:
            logging.error(f"Error processing {tma}: {str(e)}")

if __name__ == "__main__":
    main()
            




