import tifffile
import os 
import numpy as np
import zarr 
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cv2
import math
from xml.etree import ElementTree
import re
from tqdm import tqdm

import skimage.io as io
from datetime import datetime
from matplotlib import pyplot as plt
date = datetime.today().strftime('%Y-%m-%d')

from tqdm.autonotebook import tqdm 
from pyqupath.geojson import GeojsonProcessor
from pyqupath.tiff import TiffZarrReader, PyramidWriter
# from pyqupath.tma import plot_contours, tma_dearrayer

import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

import pandas as pd 

from typing import List, Dict

import skimage 

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

TMA = ['DFCI', 'Rochester', 'Tub971', 'Tub972']

IMG_PATH = "/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/07_Unspecific_staining_corrected_moderate"

with open('../data/06_crop_cores/DFCI/markerList_with_dapi.txt') as f:
    MARKER_LIST = [line.strip() for line in f]


# CHANNEL_NAME = ['CD3',
#  'CD4',
#  'CD8',
#  'Pax5',
#  'CD31',
#  'CD4',
#  'CD45RA',
#  'CD45RO',
#  'CD68',
#  'CD8',
#  'DAPI_0',
#  'DAPI_12',
#  'DAPI_16',
#  'DAPI_20',
#  'DAPI_24',
#  'DAPI_28',
#  'DAPI_32',
#  'DAPI_36',
#  'DAPI_4',
#  'DAPI_40',
#  'DAPI_44',
#  'DAPI_48',
#  'DAPI_52',
#  'DAPI_56',
#  'DAPI_8',
#  'DC-SIGN',
#  'FoxP3',
#  'GLUT1',
#  'GzmB',
#  'HLA-1',
#  'HLA-DR',
#  'HLA-E',
#  'IDO1',
#  'Ki-67',
#  'LAG3',
#  'LMP1',
#  'PD-1',
#  'Pax5',
#  'RelA',
#  'TCF1',
#  'Tox']

def extract_single_cell_info(core_img: Dict[str, np.ndarray], segmentation_mask: np.ndarray, 
                             interested_markers: List[str], output_path: str):
    """Extract single cell information from a core."""
    array_list = [core_img[channel] for channel in interested_markers]
    counts_no_noise = np.stack(array_list, axis=2)

    stats = skimage.measure.regionprops(segmentation_mask)
    label_num = len(stats)

    channel_num = len(array_list)
    data = np.zeros((label_num, channel_num))
    segment_mean = np.zeros((label_num, channel_num)) # average over all pixels in a segment
    marker_mean = np.zeros((label_num, channel_num)) # average over all non-zero pixels for a marker in a segment
    cell_sizes = np.zeros((label_num, 1)) # number of pixels in the segment
    cell_props = np.zeros((label_num, 3))
    # effective_sizes = np.zeros((label_num, channel_num)) # number of non-zero pixels for a marker in a segment 
    # variance_array = np.zeros((label_num, channel_num)) # variance over all pixels in a segment
    # IQR_array = np.zeros((label_num, channel_num)) # IQR over all pixels in a segment
    # range_array = np.zeros((label_num, channel_num)) # range over all pixels in a segment 

    for i, region in enumerate(stats):
        cell_label = region.label
        label_counts = np.array([counts_no_noise[coord[0], coord[1], :] for coord in region.coords])
        effective_area = np.sum(np.where(label_counts > 0, 1, 0), axis = 0)
        # variance = np.var(label_counts, axis = 0)
        # q95 = np.percentile(label_counts, 95, axis = 0)
        # q5 = np.percentile(label_counts, 5, axis = 0)
        # IQR = q95 - q5
        # range = np.ptp(label_counts, axis = 0)
        data[i] = np.sum(label_counts, axis=0)
        segment_mean[i] = data[i] / region.area
        marker_mean[i] = data[i] / effective_area
        marker_mean[i][np.isnan(marker_mean[i])] = 0
        cell_sizes[i] = region.area
        cell_props[i] = [cell_label, region.centroid[0], region.centroid[1]]
        # effective_sizes[i] = effective_area
        # variance_array[i] = variance
        # IQR_array[i] = IQR
        # range_array[i] = range

    col_names = []

    col_names.extend([marker + '_segmentMean' for marker in interested_markers if marker != 'Empty'])

    col_names.extend([marker + '_markerMean' for marker in interested_markers if marker != 'Empty'])

    # morph_names = []

    # morph_names.extend([marker + '_effectiveArea' for marker in interested_markers if marker != 'Empty'])

    # morph_names.extend([marker + '_var' for marker in interested_markers if marker != 'Empty'])

    # morph_names.extend([marker + '_IQR' for marker in interested_markers if marker != 'Empty'])

    # morph_names.extend([marker + '_range' for marker in interested_markers if marker != 'Empty'])

    # morph_array = np.column_stack([effective_sizes, variance_array, IQR_array, range_array])

    # morph_df = pd.DataFrame(morph_array, columns=morph_names)

    # data_df = pd.DataFrame(data, columns=col_names)
    
    # data_full = pd.concat([
    #     pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
    #     pd.DataFrame(cell_sizes, columns=["cellSize"]),
    #     morph_df,
    #     data_df,
    # ], axis=1)

    marker_array = np.column_stack([segment_mean, marker_mean])

    data_scale_size_df = pd.DataFrame(marker_array, columns=col_names)
    data_scale_size_full = pd.concat([
        pd.DataFrame(cell_props, columns=["cellLabel", "Y_cent", "X_cent"]),
        pd.DataFrame(cell_sizes, columns=["cellSize"]),
        data_scale_size_df
    ], axis=1)

    os.makedirs(output_path, exist_ok=True)
    # data_full.to_csv(os.path.join(output_path, "data.csv"), index=False)
    data_scale_size_full.to_csv(os.path.join(output_path, "dataScaleSize.csv"), index=False)


def main():

    for tma in TMA:
        tma_path = f'{IMG_PATH}/{tma}/by_core_auto_U_hat_20'
        for core in [i.path for i in os.scandir(tma_path)]:
            core_name = os.path.basename(core).split('.npy')[0]
            print(f"Extract info for {tma}-{core_name}")
            seg_path = f"../data/09_segmentation/{tma}/seg_mask_0.075_0.2/{core_name}.tiff"
            out_path = f"../data/10_extracted_info/{tma}/{core_name}"
            os.makedirs(out_path, exist_ok=True)
            seg_mask = tifffile.imread(seg_path)
            input_dict = {}
            img = np.load(core)
            for marker in MARKER_LIST:
                idx = MARKER_LIST.index(marker)
                input_dict[marker] = np.where(img[idx] < 0, 0, img[idx])
            extract_single_cell_info(input_dict, seg_mask, MARKER_LIST, out_path)


    

if __name__ == "__main__":
    main()




