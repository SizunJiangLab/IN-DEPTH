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

os.environ["CUDA_VISIBLE_DEVICES"]="0" # specify visible gpu(s) 0-3
import tensorflow as tf
try:
    tf_gpus = tf.config.list_physical_devices('GPU')
    for gpu in tf_gpus:
        tf.config.experimental.set_memory_growth(gpu, True) # allow memory growth on visible gpu(s)
except:
    pass 

import skimage.io as io
from datetime import datetime
from matplotlib import pyplot as plt
date = datetime.today().strftime('%Y-%m-%d')

from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay


from tqdm.autonotebook import tqdm 
from pyqupath.geojson import GeojsonProcessor
from pyqupath.tiff import TiffZarrReader, PyramidWriter
# from pyqupath.tma import plot_contours, tma_dearrayer

import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

TMA = ['DFCI', 'Rochester', 'Tub971', 'Tub972']

IMG_PATH = "/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/07_Unspecific_staining_corrected_harsh"

with open('../data/06_crop_cores/DFCI/markerList_with_dapi.txt') as f:
    MARKER_LIST = [line.strip() for line in f]

# CORENAME_PATTERN = r"/([^/]+)\.ome\.tiff"

NUCLEAR = ['DAPI_0']
MEMBRANE = ['HLADR', 'CD3', 'CD68']




MPP = 0.5068
MAXIMA_THRESHOLD = 0.075 # MAXIMA_THRESHOLD = 0.075
INTERIOR_THRESHOLD = 0.2 # INTERIOR_THRESHOLD = 0.2

def list_files(directory: str):
    paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            paths.append(os.path.join(root, file))
    return paths



def main():


    for tma in TMA:
        logging.info(f'Processing {tma}...')
        core_path = f'{IMG_PATH}/{tma}/by_core_auto_U_hat_20_ometiff'
        for core in list_files(core_path):
            core_name = os.path.basename(core).split('.ome.tiff')[0]
            tiff = TiffZarrReader.from_ometiff(core)
            nuclear_channel = tiff.zimg_dict['DAPI_0'][:]
            membrane_list = [tiff.zimg_dict[i][:] for i in MEMBRANE]
            membrane_channel = np.clip(np.sum(membrane_list, axis = 0), 0, 65535)
            seg_stack = np.stack((nuclear_channel, membrane_channel), axis = -1)
            seg_stack = np.expand_dims(seg_stack, 0)

            mesmer = Mesmer()

            predictions = mesmer.predict(seg_stack, image_mpp = MPP,
                                postprocess_kwargs_whole_cell={
                                "maxima_threshold": MAXIMA_THRESHOLD,
                                "interior_threshold": INTERIOR_THRESHOLD})

            rgb_image = create_rgb_image(seg_stack, channel_colors = ["green", "blue"])

            overlay = make_outline_overlay(rgb_data = rgb_image, predictions = predictions)

            segmask_output = f"../data/09_segmentation/{tma}/seg_mask_{MAXIMA_THRESHOLD}_{INTERIOR_THRESHOLD}"

            os.makedirs(segmask_output, exist_ok=True)

            tifffile.imwrite(f"{segmask_output}/{core_name}.tiff", predictions[0,...,0])

            tifffile.imwrite(f"{segmask_output}/{core_name}_overlay.tiff", overlay[0,...])

            logging.info(f'Finished segmenting {tma}: {core}.')


    
    

if __name__ == "__main__":
    main()



