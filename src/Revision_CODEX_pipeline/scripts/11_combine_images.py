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

rule_dict = {
    'DAPI': 'clean',
    'Myc-Cy5-200': 'clean',
    'Pax5': 'qptiff',
    'FoxP3': 'clean',
    'CD8': 'clean',
    'LMP1': 'clean',
    'PD1': 'clean',
    'RelB': 'qptiff',
    'IDO1': 'clean',
    'C1q': 'clean',
    'CD45RO': 'clean',
    'HLA1': 'clean',
    'GZMB': 'clean',
    'Podoplanin': 'clean',
    'CD3': 'clean',
    'CD4': 'clean',
    'CD163': 'clean',
    'CD206': 'clean',
    'CD68': 'clean',
    'BCL2': 'qptiff',
    'BCL6': 'qptiff',
    'CD79B': 'clean',
    'PDL1': 'clean',
    'LAG3': 'clean',
    'HLADR': 'clean',
    'COL1A1': 'clean',
    'Ki67': 'clean'
}

core_to_toss = {
    'DFCI': ['c00', 'c01', 'c04', 'C07', 'c07', 'c08', 'c15', 'c16', 'C18', 'C19', 'c19', 'c22', 'c25', 'c30', 
    'c35', 'c36', 'c37', 'c38', 'c39', 'c40', 'c42', 'c44', 'c43', 'c45', 'C46', 'C47', 'C48', 'C49', 'c46', 
    'c48', 'c49', 'C50', 'c50', 'C51', 'c52', 'c51', 'c53'],
    
    'Rochester': ['c24', 'C19', 'c05', 'c27', 'c10', 'c26', 'c12', 'c09', 'c00', 'c01', 'c02', 'C03'],
    
    'Tub971': ['c20', 'c04', 'c09', 'C14', 'c23', 'c01', 'c22', 'c16', 'c21'],
    
    'Tub972': ['None_4', 'None_5', 'c19', 'c20', 'None_2', 'None_3', 'c08']
}

if __name__ == '__main__':
    for tma in ['DFCI', 'Rochester', 'Tub971', 'Tub972']:
        for core in os.scandir(f'../data/07_Unspecific_staining_corrected_moderate/{tma}/by_core_auto_U_hat_20_ometiff'):
            core_name = os.path.basename(core).split('.ome.tiff')[0]
            if core_name in core_to_toss[tma]:
                print(f'Skip {tma}: {core_name}...')
                continue
            qptiff_name = os.path.basename(core).split('.ome.tiff')[0]
            if core_name == 'None_1':
                qptiff_name = 'NA_0'
            if core_name == 'None_2':
                qptiff_name = 'NA_1'
            if core_name == 'None_3':
                qptiff_name = 'NA_2'
            if core_name == 'None_4':
                qptiff_name = 'NA_3'
            if core_name == 'None_5':
                qptiff_name = 'NA_4'
            clean_img = TiffZarrReader.from_ometiff(core)
            qptiff_img = TiffZarrReader.from_ometiff(f'../data/crop_final_qptiff/{tma}/{qptiff_name}.ome.tiff')
            out_dict = {}
            for key, item in rule_dict.items():
                if item == 'clean':
                    if key == 'DAPI':
                        clean_key = 'DAPI_0'
                        out_dict[key] = clean_img.zimg_dict[clean_key]
                    else:
                        out_dict[key] = clean_img.zimg_dict[key]
                elif item == 'qptiff':
                    out_dict[key] = qptiff_img.zimg_dict[key]
            writer = PyramidWriter.from_dict(out_dict)
            out_dir = f'../data/combined_img/{tma}'
            os.makedirs(out_dir, exist_ok=True)
            writer.export_ometiff_pyramid(f'../data/combined_img/{tma}/{core_name}.ome.tiff')
            