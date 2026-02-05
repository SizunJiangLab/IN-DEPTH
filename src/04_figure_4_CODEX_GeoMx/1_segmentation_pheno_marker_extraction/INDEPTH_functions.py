import cv2 
import tifffile 
import numpy as np
import math
import pandas as pd
import skimage.io as io
import json
from skimage import img_as_ubyte
from datetime import datetime
from matplotlib import pyplot as plt
date = datetime.today().strftime('%Y-%m-%d')

import os 
os.environ["CUDA_VISIBLE_DEVICES"]="0" # specify visible gpu(s) 0-3
import tensorflow as tf
try:
    tf_gpus = tf.config.list_physical_devices('GPU')
    for gpu in tf_gpus:
        tf.config.experimental.set_memory_growth(gpu, True) # allow memory growth on visible gpu(s)
except:
    pass 

from deepcell.applications import Mesmer
from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay


def process_fusion_image(fusion_img, channelNameFile, configFile, coreName = None):
    """
    Store the fusion image into a dictionary and process the channels specified by the user.

    Args:
    - fusion_img (np.array): numpy array of rotated fusion image.
    - channelNameFile (str): Path to your channelName.txt
    - configFile (str): Path to your configuration json.
    - coreName (str): Name of the core. For example, "1"

    Returns:
    -  img_dict (dict): A dictionary of numpy arrays that represents each channel of your input image.
    """
    # Create an empty dictionary to store the image 
    img_dict = {}

    # read in the marker list 

    txt_file = open(channelNameFile, 'r')

    MarkerList = txt_file.read().split('\n')

    txt_file.close()

    # read in the config file

    with open(configFile) as json_file:
        config = json.load(json_file)

    # Check if global threshold 
    
    is_global = config['is_global']

    if is_global:
        channels = config['global']['channels']
        thresholds = config['global']['lowerBound']
        print(f'Processing {channels} with lower bound {thresholds}')
    
        # iterate through each channels and process the channel if specified by the user

        for i, marker in enumerate(MarkerList):
            #print(f'Channel: {marker}')
            if marker in channels:
                # get the corresponding cut off
                cutoff = thresholds[channels.index(marker)]
                #print(f'Processing {marker} with lower bound {cutoff}')
                # get the corresponding image for the marker 
                img = fusion_img[i]
                # process the image 
                img_filtered = np.copy(img)
                #print(f'{marker}: {np.max(img) - cutoff}')
                img_filtered[img_filtered <= cutoff] = 0
                # write the array to the dictionary 
                #print(f'Write {marker} to dictionary')
                img_dict[marker] = img_filtered
            else:
                img = fusion_img[i]
                # rotate the image
                #img_rotate = cv2.rotate(img, cv2.ROTATE_180)
                #print(f'Write {marker} to dictionary')
                img_dict[marker] = img
    else:

        if coreName is None:
            raise ValueError("Need to input a list of core names.")
        
        core_dict = config['cores'][coreName]
        channels = core_dict['channels']
        thresholds = core_dict['lowerBound']
        for i, marker in enumerate(MarkerList):
            #print(f'Channel: {marker}')
            if marker in channels:
                # get the corresponding cut off
                cutoff = thresholds[channels.index(marker)]
                #print(f'Processing {marker} with lower bound {cutoff}')
                # get the corresponding image for the marker 
                img = fusion_img[i]
                # process the image 
                img_filtered = np.copy(img)
                img_filtered[img_filtered <= cutoff] = 0
                # write the array to the dictionary 
                #print(f'Write {marker} to dictionary')
                img_dict[marker] = img_filtered
            else:
                img = fusion_img[i]
                # rotate the image
                #img_rotate = cv2.rotate(img, cv2.ROTATE_180)
                #print(f'Write {marker} to dictionary')
                img_dict[marker] = img

    return img_dict

def segment_image(img_dict, nuclear_channel, membrane_channel, 
                  mpp, output_path, maxima_threshold = 0.075, interior_threshold = 0.2,
                  search_mode = True, run_mode = False):
    """
    Segment image.

    Args:
    - img_dict (dictionary): A dictionary storing the channels of the image.
    - nuclear_channel (list): A list of strings specifying the nuclear channels.
    - membrane_channel (list): A list of strings specifying the membrane channels.
    - mpp (numeric): Micron per pixel of the instrument.
    - maxima_threshold (numeric): maxima_threshold controls what is considered a unique cell (lower values = more separate cells, higher values = fewer cells).
    - interior_threshold (numeric): interior_threshold determines what is considered background/not part of a cell (lower value = larger cells).
    - search_mode (boolean): If true, a search will be performed for the user to tune maxima_threshold and interior_threshold.
    - run_mode (boolean): If true, the user needs to input the maxima_threshold and interior_threshold.

    Returns:
    -  mesmer_mask (array): segmentation mask, a 2d numpy array.
    """

    nuclear_img = np.clip(np.sum([img_dict[channel] for channel in nuclear_channel], axis = 0),0,255).astype('uint8')

    membrane_img = np.clip(np.sum([img_dict[channel] for channel in membrane_channel], axis = 0),0,255).astype('uint8')

    seg_stack = np.stack((nuclear_img, membrane_img), axis = -1)

    seg_stack = np.expand_dims(seg_stack, 0)

    mesmer = Mesmer()

    if search_mode:
        img_centroid_y = math.ceil(seg_stack.shape[1]/2)
        img_centroid_x = math.ceil(seg_stack.shape[2]/2)
        img_start_y = img_centroid_y - 250
        img_end_y = img_centroid_y + 250
        img_start_x = img_centroid_x - 250
        img_end_x = img_centroid_x + 250
        stack_crop = seg_stack[:,img_start_y:img_end_y, img_start_x:img_end_x, :]
        maxima_lower = np.arange(0.015, maxima_threshold, 0.01).round(3)
        maxima_upper = np.arange(maxima_threshold + 0.01, 0.2, 0.01).round(3)
        interior_lower = np.arange(0.05, interior_threshold, 0.05).round(3)
        interior_upper = np.arange(interior_threshold + 0.05, 0.4, 0.05).round(3)
        maxima_range = np.concatenate((maxima_lower, maxima_threshold, maxima_upper))
        interior_range = np.concatenate((interior_lower, interior_threshold, interior_upper))
        for i in maxima_range:
            for j in interior_range:
                print(f"Performing segmentation with maxima_threshold:{i} and interior_threshold:{j}")
                predictions = mesmer.predict(stack_crop, image_mpp = mpp,
                                             postprocess_kwargs_whole_cell={"maxima_threshold" : i,
                                                                            "interior_threshold": j})
                rgb_image = create_rgb_image(stack_crop, channel_colors = ["green", "blue"])
                overlay = make_outline_overlay(rgb_data = rgb_image, predictions = predictions)
                if os.path.exists(output_path) == False:
                    os.makedirs(output_path)
                tifffile.imwrite(os.path.join(output_path, f"maxima_{i}_interior_{j}_overlay.tiff"), overlay[0,...])
    elif run_mode:
        predictions = mesmer.predict(seg_stack, image_mpp = mpp,
                                    postprocess_kwargs_whole_cell={
                                        "maxima_threshold": maxima_threshold,
                                        "interior_threshold": interior_threshold
                                    })
        rgb_image = create_rgb_image(seg_stack, channel_colors = ["green", "blue"])
        overlay = make_outline_overlay(rgb_data = rgb_image, predictions = predictions)
        if os.path.exists(output_path) == False:
            os.makedirs(output_path)
        tifffile.imwrite(os.path.join(output_path, 'MESMER_mask.tiff'), predictions[0,...,0])
        tifffile.imwrite(os.path.join(output_path, 'overlay.tiff'), overlay[0,...])
        tifffile.imwrite(os.path.join(output_path, 'overlay_downsample.tiff'), overlay[0,::4,::4])
        
        
        


        


