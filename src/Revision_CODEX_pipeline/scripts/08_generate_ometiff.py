import numpy as np
import os 
import tifffile 
import re
from xml.etree import ElementTree
import pandas as pd 
import json

from tqdm.autonotebook import tqdm 
from pyqupath.geojson import GeojsonProcessor
from pyqupath.tiff import TiffZarrReader, PyramidWriter
from pyqupath.tma import plot_contours, tma_dearrayer
from pathlib import Path


def list_files(directory: str):
    paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            paths.append(os.path.join(root, file))
    return paths


if __name__ == '__main__':

    for tma in [i.path for i in os.scandir('../data/07_Unspecific_staining_corrected_moderate_new') if '_ometiff' not in i.path]:
        print(f'Processing {tma}...')
        tma_name = os.path.basename(tma)
        with open(f'/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/06_crop_cores/{tma_name}/markerList_with_dapi.txt', 'r') as f:
            markerList = [line.strip() for line in f.readlines()]
        for path in [j.path for j in os.scandir(tma) if '_ometiff' not in j.path]:
            for core in list_files(path):
                core_name = os.path.basename(core).split('.npy')[0]
                print(f'Working on {path} : {core_name}...')
                img = np.load(core)
                out = {}
                out_dir = f'{path}_ometiff'
                os.makedirs(out_dir, exist_ok=True)
                for marker in markerList:
                    idx = markerList.index(marker)
                    out[marker] = np.where(img[idx] < 0, 0, img[idx]).astype(np.uint16)
                writer = PyramidWriter.from_dict(out)
                writer.export_ometiff_pyramid(f'{out_dir}/{core_name}.ome.tiff')

    #  for tma in [i.path for i in os.scandir('../data/07_Unspecific_staining_corrected_harsh') if '_ometiff' not in i.path]:
    #     print(f'Processing {tma}...')
    #     tma_name = os.path.basename(tma)
    #     if tma_name == 'DFCI':

    #         with open(f'/mnt/nfs/home/huayingqiu/for_steph/INDEPTH_NOV_26/data/06_crop_cores/{tma_name}/markerList_with_dapi.txt', 'r') as f:
    #             markerList = [line.strip() for line in f.readlines()]
    #         for path in [j.path for j in os.scandir(tma)]:
    #             for core in list_files(path):
    #                 core_name = os.path.basename(core).split('.npy')[0]
    #                 print(f'Working on {path} : {core_name}...')
    #                 img = np.load(core)
    #                 out = {}
    #                 out_dir = f'{path}_ometiff'
    #                 os.makedirs(out_dir, exist_ok=True)
    #                 for marker in markerList:
    #                     idx = markerList.index(marker)
    #                     out[marker] = np.where(img[idx] < 0, 0, img[idx]).astype(np.uint16)
    #                 writer = PyramidWriter.from_dict(out)
    #                 writer.export_ometiff_pyramid(f'{out_dir}/{core_name}.ome.tiff')
    #     else:
    #         continue

