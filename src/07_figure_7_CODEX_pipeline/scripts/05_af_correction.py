import os
import re
import math
import pickle
import numpy as np
import cv2 as cv
import tifffile
from skimage import morphology
from sklearn.linear_model import HuberRegressor
from tqdm import tqdm
from pyqupath.tiff import TiffZarrReader
import gc

# =============================================================================
# CONFIGURATION & CONSTANTS
# =============================================================================

# TMA = ['DFCI', 'Rochester', 'Tub971', 'Tub972']
TMA = ['Tub971', 'Tub972']

INPUT_OME_TIFF = {
    'DFCI': '../data/04_stitched_cycles/DFCI.ome.tiff',
    'Rochester': '../data/04_stitched_cycles/Rochester.ome.tiff',
    'Tub971': '../data/04_stitched_cycles/Tub971.ome.tiff',
    'Tub972': '../data/04_stitched_cycles/Tub972.ome.tiff'
}

OUTPUT_DIR = {
    'DFCI': '../data/05_AF_corrected_markers/DFCI',
    'Rochester': '../data/05_AF_corrected_markers/Rochester',
    'Tub971': '../data/05_AF_corrected_markers/Tub971',
    'Tub972': '../data/05_AF_corrected_markers/Tub972'
}


# Ensure output directory exists
for key, item in OUTPUT_DIR.items():
    if not os.path.exists(item):
        os.makedirs(item)

CYCLE_CHANNEL_EXPOSURE_DICT = {
    'Myc-Cy3-100': '3_Cy3_100',
    'CD3': '3_Cy5_200',
    'Pax5': '4_Cy3_200',
    'Myc-Cy5-200': '4_Cy5_200',
    'FoxP3': '5_Cy3_200',
    'CD4': '5_Cy5_200',
    'CD8': '6_Cy3_200',
    'CD206': '6_Cy5_200',
    'CD163': '7_Cy3_200',
    'CD68': '7_Cy5_200',
    'LMP1': '8_Cy3_200',
    'BCL2': '8_Cy5_200',
    'PD1': '9_Cy3_200',
    'BCL6': '9_Cy5_200',
    'RelB': '10_Cy3_200',
    'CD79B': '10_Cy5_200',
    'TCF1': '11_Cy3_200',
    'PDL1': '11_Cy5_300',
    'IDO1': '12_Cy3_200',
    'Tox12': '12_Cy5_300',
    'C1q': '13_Cy3_200',
    'TIM3': '13_Cy5_300',
    'SPP1': '14_Cy3_200',
    'LAG3': '14_Cy5_300',
    'CD45RO': '15_Cy3_200',
    'CD45RA': '15_Cy5_200',
    'Geminin': '16_Cy3_200',
    'CD31': '16_Cy5_200',
    'HLA1': '17_Cy3_200',
    'HLADR': '17_Cy5_200',
    'GZMB': '18_Cy3_200',
    'COL1A1': '18_Cy5_200',
    'Podoplanin': '19_Cy3_200',
    'Ki67': '19_Cy5_200',
}

COMPARTMENT_DICT = {
    "Cy3-first-blank-100": "extranuclear",
    "Cy5-first-blank-300": "extranuclear",

    "Myc-Cy3-100": "nuclear",
    "CD3": "extranuclear",

    "Pax5": "nuclear",
    "Myc-Cy5-200": "nuclear",

    "FoxP3": "nuclear",
    "CD4": "extranuclear",

    "CD8": "extranuclear",
    "CD206": "extranuclear",

    "CD163": "extranuclear",
    "CD68": "extranuclear",

    "LMP1": "extranuclear",
    "BCL2": "extranuclear",

    "PD1": "extranuclear",
    "BCL6": "nuclear",

    "RelB": "nuclear",
    "CD79B": "extranuclear",

    "TCF1": "nuclear",
    "PDL1": "extranuclear",

    "IDO1": "extranuclear",
    "Tox12": "nuclear",

    "C1q": "extranuclear",
    "TIM3": "extranuclear",

    "SPP1": "extranuclear",
    "LAG3": "extranuclear",

    "CD45RO": "extranuclear",
    "CD45RA": "extranuclear",

    "Geminin": "nuclear",
    "CD31": "extranuclear",

    "HLA1": "extranuclear",
    "HLADR": "extranuclear",

    "GZMB": "extranuclear",
    "COL1A1": "extranuclear",

    "Podoplanin": "extranuclear",
    "Ki67": "nuclear"
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def moment_threshold(image, num_bins=65536):
    """
    Generalized Tsai moment-preserving thresholding method.
    """
    img = image.astype(np.float64)
    img_min, img_max = img.min(), img.max()

    if img_max == img_min:
        return img_min 

    hist, bin_edges = np.histogram(img, bins=num_bins, range=(img_min, img_max))
    hist = hist.astype(np.float64)
    total = hist.sum()

    if total == 0:
        return img_min

    p = hist / total
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    m0 = 1.0
    m1 = np.sum(centers * p)
    m2 = np.sum((centers**2) * p)
    m3 = np.sum((centers**3) * p)

    cd = m0 * m2 - m1 * m1
    if cd == 0:
        return m1

    c0 = (-m2 * m2 + m1 * m3) / cd
    c1 = (m0 * (-m3) + m2 * m1) / cd

    disc = c1 * c1 - 4 * c0
    if disc < 0:
        return m1

    sqrt_disc = math.sqrt(disc)
    z0 = 0.5 * (-c1 - sqrt_disc)
    z1 = 0.5 * (-c1 + sqrt_disc)

    p0 = 0.5 if z1 == z0 else (z1 - m1) / (z1 - z0)
    p0 = max(0.0, min(1.0, p0))

    cumsum = np.cumsum(p)
    idx = np.searchsorted(cumsum, p0)

    if idx >= len(centers):
        return centers[-1]
    return centers[idx]

def yen_threshold(image, show_process=False):
    """
    Implements Yen's thresholding method from ImageJ.
    
    Based on:
    1) Yen J.C., Chang F.J., and Chang S. (1995) "A New Criterion 
       for Automatic Multilevel Thresholding" IEEE Trans. on Image 
       Processing, 4(3): 370-378
    2) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
       Techniques and Quantitative Performance Evaluation" Journal of 
       Electronic Imaging, 13(1): 146-165
    
    The algorithm finds the threshold that maximizes a criterion based on
    the discrepancy between the original and thresholded images using
    information theory measures.
    
    Parameters:
    -----------
    image : numpy.ndarray
        Input grayscale image (8-bit or 16-bit)
    show_process : bool
        Whether to visualize the criterion function
        
    Returns:
    --------
    threshold : int
        Optimal threshold value (0-255 for 8-bit, 0-65535 for 16-bit)
    """
    
    # Determine bit depth and set parameters accordingly
    if image.dtype == np.uint8:
        # 8-bit image
        bins = 256
        max_val = 255
        hist_range = (0, 255)
    elif image.dtype == np.uint16:
        # 16-bit image
        bins = 65536
        max_val = 65535
        hist_range = (0, 65535)
    else:
        # Convert other dtypes to uint8 (existing behavior for float, etc.)
        image = (image * 255).astype(np.uint8)
        bins = 256
        max_val = 255
        hist_range = (0, 255)
    
    # Calculate histogram
    hist, _ = np.histogram(image.flatten(), bins=bins, range=hist_range)
    total = np.sum(hist)
    
    if total == 0:
        return 0
    
    # Normalize histogram to get probabilities
    norm_histo = hist.astype(np.float64) / total
    
    # Calculate cumulative probabilities P1[t] = sum of probabilities from 0 to t
    P1 = np.zeros(bins)
    P1[0] = norm_histo[0]
    for i in range(1, bins):
        P1[i] = P1[i-1] + norm_histo[i]
    
    # Calculate cumulative squared probabilities
    P1_sq = np.zeros(bins)
    P1_sq[0] = norm_histo[0] * norm_histo[0]
    for i in range(1, bins):
        P1_sq[i] = P1_sq[i-1] + norm_histo[i] * norm_histo[i]
    
    # Calculate reverse cumulative squared probabilities
    P2_sq = np.zeros(bins)
    P2_sq[bins-1] = 0.0
    for i in range(bins-2, -1, -1):
        P2_sq[i] = P2_sq[i + 1] + norm_histo[i + 1] * norm_histo[i + 1]
    
    # Find threshold that maximizes the criterion
    threshold = -1
    max_crit = -np.inf
    criterion_values = []
    
    for t in range(bins):
        # Yen's criterion formula
        # This measures the discrepancy between original and binarized image
        term1 = P1_sq[t] * P2_sq[t]
        term2 = P1[t] * (1.0 - P1[t])
        
        # Calculate criterion (with safe logarithm)
        log_term1 = math.log(term1) if term1 > 0 else 0.0
        log_term2 = math.log(term2) if term2 > 0 else 0.0
        
        crit = -1.0 * log_term1 + 2.0 * log_term2
        criterion_values.append(crit)
        
        if crit > max_crit:
            max_crit = crit
            threshold = t
    
    return threshold

def blank_estimator(start, end, total_cycle, current_cycle):
    """Linearly interpolate blank image between start and end cycles."""
    w_t = (total_cycle - 1 - current_cycle) / max(total_cycle - 1, 1)
    # estimator blank for current cycle
    # Note: Loading images inside here if they are passed as paths, 
    # but assuming passed as arrays for speed if memory allows for just 2 reference images.
    B_t = w_t * start.astype(np.float32) + (1 - w_t) * end.astype(np.float32)
    return B_t

def cell_area_mask(dapi: np.ndarray, rbc_mask: np.ndarray):
    """Generate nucleus and extranucleus masks."""
    good_areas = 1 - rbc_mask
    dapi_filter = (dapi * good_areas).astype(np.uint16)
    _, otsu = cv.threshold(dapi_filter, 0, 2**16 - 1, cv.THRESH_BINARY + cv.THRESH_OTSU)
    otsu[otsu != 0] = 1
    nuclear_bkg_mask = 1 - otsu
    
    neucleus_dilate = morphology.binary_dilation(otsu > 0, morphology.disk(2)).astype(np.uint16) * 1
    extraneucleus_dilate = morphology.binary_dilation(nuclear_bkg_mask > 0, morphology.disk(2)).astype(np.uint16) * 1
    return neucleus_dilate, extraneucleus_dilate

def AF_regressor(I_c: np.ndarray, B_t: np.ndarray, Omega: np.ndarray, use_huber: bool = True):
    """Perform background subtraction regression."""
    I_c = I_c.astype(np.float32)
    B_t = B_t.astype(np.float32)

    mask = Omega.astype(bool)
    
    # Sample a subset if image is huge to save memory during fit? 
    # For now, using full data as per notebook logic but be mindful.
    mask_flat = mask.ravel()
    I_flat = I_c.ravel()
    B_flat = B_t.ravel()

    y = I_flat[mask_flat]
    x = B_flat[mask_flat]

    if y.size < 100:
        raise ValueError('Not enough pixels to fit AF regression')
    
    X = np.stack([np.ones_like(x, dtype=np.float32), x], axis=1)
    eps = 1e-6
    w = 1 / np.maximum(y, eps)

    if use_huber:
        hub = HuberRegressor(fit_intercept=False, epsilon=1.35, alpha=0.0)
        hub.fit(X, y, sample_weight=w)
        beta0, beta1 = map(float, hub.coef_)
    else:
        WX = X * np.sqrt(w[:, None])
        Wy = y * np.sqrt(w)
        beta_hat, *_ = np.linalg.lstsq(WX, Wy, rcond=None)
        beta0, beta1 = map(float, beta_hat)
    
    I_afcorr = I_c - (beta0 + beta1 * B_t)
    return I_afcorr.astype(np.float32), beta0, beta1

def get_blank_lists(tiff):
    """Identify blank channels."""
    pattern_cy3 = re.compile(r'^Cy3-')
    cy3_list = [key for key in tiff.zimg_dict.keys() if pattern_cy3.search(key)]
    
    pattern_cy5 = re.compile(r'^Cy5-')
    cy5_list = [key for key in tiff.zimg_dict.keys() if pattern_cy5.search(key)]
    
    return cy3_list, cy5_list

# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    for tma in TMA:
        print(f'Processing {tma}...')
        BETA_LOG_FILE = os.path.join(OUTPUT_DIR[tma], 'beta','beta_estimate.txt')
        beta_dir = os.path.dirname(BETA_LOG_FILE)
        os.makedirs(beta_dir, exist_ok=True)
        print("Loading TIFF reader...")
        tiff = TiffZarrReader.from_ometiff(INPUT_OME_TIFF[tma])
        
        # 1. Identify Blanks
        cy3_blank_list, cy5_blank_list = get_blank_lists(tiff)
        
        # Load reference blanks for masking and scaling calculation
        # We keep these in memory as they are needed frequently, but they are few.
        print("Loading reference blanks...")
        cy3_blank_dict = {k: tiff.zimg_dict[k][:] for k in cy3_blank_list}
        cy5_blank_dict = {k: tiff.zimg_dict[k][:] for k in cy5_blank_list}

        # 2. Generate RBC/Artifact Masks
        print("Generating Cy3 Masks...")
        mask_dict = {}
        for blank_key, img in cy3_blank_dict.items():
            thres = np.percentile(img, 99)
            print(f'{blank_key}: {thres}')
            mask = (img > thres).astype(img.dtype) * 1
            mask_dict[f'{blank_key}_mask'] = mask

        print("Generating Cy5 Masks...")
        cy5_mask_dict = {}
        for blank_key, img in cy5_blank_dict.items():
            # Use 65535 as maxval for 16-bit images
            thres = np.percentile(img, 99)
            print(f'{blank_key}: {thres}')
            mask = (img > thres).astype(img.dtype) * 1
            cy5_mask_dict[f'{blank_key}_mask'] = mask

        # # 3. Calculate Relative Scaling Factors
        # print("Calculating exposure scaling factors...")
        # # Cy3 Scaling
        # dilate_cy3_150 = morphology.binary_dilation(mask_dict['Cy3-Blank-150_2_moment_mask'] > 0, morphology.disk(2)).astype(np.uint16)
        # dilate_cy3_250 = morphology.binary_dilation(mask_dict['Cy3-Blank-250_moment_mask'] > 0, morphology.disk(2)).astype(np.uint16)
        
        # cy3_150_rbc = cy3_blank_dict['Cy3-Blank-150_2'] * dilate_cy3_150
        # cy3_250_rbc = cy3_blank_dict['Cy3-Blank-250'] * dilate_cy3_250
        
        # cy3_m_150 = np.mean(cy3_150_rbc[cy3_150_rbc > 0])
        # cy3_m_250 = np.mean(cy3_250_rbc[cy3_250_rbc > 0])
        # cy3_s_250 = cy3_m_250 / cy3_m_150

        # # Cy5 Scaling
        # dilate_cy5_150 = morphology.binary_dilation(cy5_mask_dict['Cy5-Blank-150_2_otsu_mask'] > 0, morphology.disk(2)).astype(np.uint16)
        # dilate_cy5_250 = morphology.binary_dilation(cy5_mask_dict['Cy5-Blank-250_otsu_mask'] > 0, morphology.disk(2)).astype(np.uint16)
        
        # cy5_150_rbc = cy5_blank_dict['Cy5-Blank-150_2'] * dilate_cy5_150
        # cy5_250_rbc = cy5_blank_dict['Cy5-Blank-250'] * dilate_cy5_250
        
        # cy5_m_150 = np.mean(cy5_150_rbc[cy5_150_rbc > 0])
        # cy5_m_250 = np.mean(cy5_250_rbc[cy5_250_rbc > 0])
        # cy5_s_250 = cy5_m_250 / cy5_m_150
        
        # s_200 = 200 / 150

        # print(f"Cy3 Scale 250ms: {cy3_s_250}")
        # print(f"Cy5 Scale 250ms: {cy5_s_250}")

        # 4. Create Union Masks for Global Exclusion
        # Reduce and Dilate
        print("Creating Union Masks...")
        cy3_mask_list = [v for k, v in mask_dict.items() if 'Cy3' in k]
        cy5_mask_list = [v for k, v in cy5_mask_dict.items() if 'Cy5' in k ]

        cy3_union_mask = morphology.binary_dilation(
            np.logical_or.reduce(cy3_mask_list).astype(np.uint16) > 0, 
            morphology.disk(2)
        ).astype(np.uint16)

        cy5_union_mask = morphology.binary_dilation(
            np.logical_or.reduce(cy5_mask_list).astype(np.uint16), 
            morphology.disk(2)
        ).astype(np.uint16)


        # # Define Reference Images for Interpolation to avoid re-accessing dictionary
        # # NOTE: We keep these in memory. If they are too large, load them inside the loop too.
        # ref_cy3_start = cy3_blank_dict['Cy3-Blank-150_1']
        # ref_cy3_end = cy3_blank_dict['Cy3-Blank-150_2']
        # ref_cy5_start = cy5_blank_dict['Cy5-Blank-150_1']
        # ref_cy5_end = cy5_blank_dict['Cy5-Blank-150_2']
        # # Handle rescue blanks if they exist
        # if 'Cy5-Blank-250_1_rescue' in cy5_blank_dict:
        #     ref_cy5_res_start = cy5_blank_dict['Cy5-Blank-250_1_rescue']
        #     ref_cy5_res_end = cy5_blank_dict['Cy5-Blank-250_2_rescue']
        
        # # Clear temporary dictionaries to free memory before heavy processing
        # del mask_dict, cy5_mask_dict, cy3_blank_dict, cy5_blank_dict
        # gc.collect()

        # 5. Processing Markers
        print("Starting Marker Processing...")
        
        # Clear log file
        with open(BETA_LOG_FILE, 'w') as f:
            f.write("Marker: Beta0, Beta1\n")

        for key, item in tqdm(CYCLE_CHANNEL_EXPOSURE_DICT.items()):
            print(f'Processing {key}...')
            
            parser = item.split('_')
            cycle = int(parser[0]) - 1
            channel = parser[1]
            exposure = parser[2]
            dapi_key = f'DAPI_{cycle}'  

            # Load Images
            try:
                dapi = tiff.zimg_dict[dapi_key][:]
                marker = tiff.zimg_dict[key][:]
            except KeyError as e:
                print(f"Skipping {key}: Missing image data {e}")
                continue

            compartment = COMPARTMENT_DICT.get(key, 'nuclear') # Default to nucleus if missing
            
            # Logic per channel
            if channel == 'Cy3':
                rbc_mask = cy3_union_mask
                total_cycles = 21 # From notebook logic
                
                # On-the-fly blank estimation
                start_blank = cy3_blank_dict[f'Cy3-first-blank-{exposure}']
                end_blank = cy3_blank_dict[f'Cy3-last-blank-{exposure}']
                blank = blank_estimator(start_blank, end_blank, total_cycles, cycle)
                    
            elif channel == 'Cy5':
                rbc_mask = cy5_union_mask
                total_cycles = 21
                start_blank = cy5_blank_dict[f'Cy5-first-blank-{exposure}']
                end_blank = cy5_blank_dict[f'Cy5-last-blank-{exposure}']
                blank = blank_estimator(start_blank, end_blank, total_cycles, cycle)


            else:
                print(f"Unknown channel {channel} for {key}, skipping.")
                continue

            # Filter Images
            good_area = 1 - rbc_mask
            marker_filtered = marker * good_area
            blank_filtered = blank * good_area

            # Determine Omega (background set)
            nuclear, extranuclear = cell_area_mask(dapi, rbc_mask)
            
            if compartment == 'nuclear':
                Omega = extranuclear
            else:
                Omega = nuclear
            
            if not os.path.exists(f'{OUTPUT_DIR[tma]}/omega'):
                os.makedirs(f'{OUTPUT_DIR[tma]}/omega')

            np.save(f'{OUTPUT_DIR[tma]}/omega/{key}_{compartment}_omega.npy', Omega)

            # Perform Regression
            try:
                marker_corrected, beta0, beta1 = AF_regressor(
                    marker_filtered, 
                    blank_filtered, 
                    Omega, 
                    use_huber=False
                )
                
                # Save Result
                if not os.path.exists(f'{OUTPUT_DIR[tma]}/marker_array'):
                    os.makedirs(f'{OUTPUT_DIR[tma]}/marker_array')
                save_path = os.path.join(OUTPUT_DIR[tma], 'marker_array', f'{key}.npy')
                np.save(save_path, marker_corrected)
                
                with open(BETA_LOG_FILE, 'a') as f:
                    f.write(f'{key}: {beta0}, {beta1}\n')
                    
            except Exception as e:
                print(f"Error processing {key}: {e}")

            # Memory Cleanup
            del marker, blank, dapi, marker_filtered, blank_filtered, Omega, nuclear, extranuclear, marker_corrected
            gc.collect()

        print("Processing Complete.")

if __name__ == "__main__":
    main()