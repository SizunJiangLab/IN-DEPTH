import os
import numpy as np
import re
import logging

# --- Configuration ---
# Inputs from Script 1
TMA = ['DFCI', 'Rochester', 'Tub971', 'Tub972']


def list_files(directory: str):
    paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            paths.append(os.path.join(root, file))
    return paths

def decide_U_hat(vals, 
                 trim=(1, 99),
                 pct_lowband=(5, 25),
                 pos_pct_thresholds=(10, 70),
                 mid_thres = 20): 
    """
    Adaptive estimation of U_c based on AF-corrected values in Omega.
    """
    vals = np.asarray(vals)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.0, None

    # 1. Trim extremes
    lo, hi = np.percentile(vals, trim)
    trimmed = vals[(vals >= lo) & (vals <= hi)]
    if trimmed.size == 0:
        return 0.0, None

    # 2. Find first positive percentile within trimmed vals
    percentiles = np.linspace(0, 100, 1001)  # 0.1% resolution
    q_vals = np.percentile(trimmed, percentiles)
    pos_mask = q_vals > 0

    if not pos_mask.any():
        return 0.0, None

    first_pos_pct = percentiles[pos_mask.argmax()]

    # Decision thresholds
    medium_th, bad_th = pos_pct_thresholds 
    
    # 3. Case 1: Bad channels
    if first_pos_pct > bad_th:
        print(f'Bad channel. First positive percentile: {first_pos_pct}')
        return 0.0, float(first_pos_pct)

    # 4. Case 2: Medium-quality channels
    if first_pos_pct > medium_th:
        print(f'Moderate channel. First positive percentile: {first_pos_pct}')
        U_raw = np.percentile(q_vals[pos_mask], mid_thres) 
        return max(U_raw, 0.0), float(first_pos_pct)

    # 5. Case 3: Good channels
    print(f'Good channel. First positive percentile: {first_pos_pct}')
    p1, p2 = pct_lowband
    band_lo, band_hi = np.percentile(trimmed, [p1, p2])
    band = trimmed[(trimmed >= band_lo) & (trimmed <= band_hi)]
    if band.size == 0:
        return 0.0, float(first_pos_pct)

    U_raw = np.median(band)
    return max(U_raw, 0.0), float(first_pos_pct)

def main():
    for tma in TMA:

        CROPPED_CORE_DIR = f'../data/06_crop_cores/{tma}/cropped_arrays_with_dapi'
        CROPPED_OMEGA_DIR = f'../data/06_crop_cores/{tma}/cropped_omega'
        MARKER_LIST_FILE = f'../data/06_crop_cores/{tma}/markerList_with_dapi.txt'
        OMEGA_LIST_FILE = f'../data/06_crop_cores/{tma}/omegaList.txt'

        # Output
        OUT_BASE_DIR = f'../data/07_Unspecific_staining_corrected_moderate_new/{tma}'

        os.makedirs(OUT_BASE_DIR, exist_ok=True)


        # 1. Load Metadata Lists
        if not os.path.exists(MARKER_LIST_FILE) or not os.path.exists(OMEGA_LIST_FILE):
            raise FileNotFoundError("Metadata lists not found. Run Script 1 first.")

        with open(MARKER_LIST_FILE, 'r') as f:
            markerList = [line.strip() for line in f.readlines()]

        with open(OMEGA_LIST_FILE, 'r') as f:
            omegaList = [line.strip() for line in f.readlines()]

        # 2. Get list of cores
        core_paths = list_files(CROPPED_CORE_DIR)
        if not core_paths:
            print("No cropped cores found.")
            return

        # 3. Iterate through parameters (as per notebook)
        for mid_thres in [20]:
            print(f"--- Processing with mid_thres = {mid_thres} ---")
            
            # Create output directory for this parameter
            current_out_dir = os.path.join(OUT_BASE_DIR, f'by_core_auto_U_hat_{mid_thres}')
            if not os.path.exists(current_out_dir):
                os.makedirs(current_out_dir)

            for core_path in core_paths:
                core_name = os.path.basename(core_path).split('.npy')[0]
                
                # Load Data
                img_stack = np.load(core_path) # Shape: (C, H, W)
                omega_path = os.path.join(CROPPED_OMEGA_DIR, f'{core_name}.npy')
                
                if not os.path.exists(omega_path):
                    print(f"Warning: Omega mask not found for {core_name}, skipping.")
                    continue
                    
                omega_stack = np.load(omega_path) # Shape: (M, H, W)

                output_list = []

                for idx, marker in enumerate(markerList):
                    if idx >= img_stack.shape[0]: break # Safety check

                    # Pass DAPI through unchanged
                    if 'DAPI' in marker:
                        output_list.append(img_stack[idx])
                        continue

                    I_afcorr = img_stack[idx]

                    # Find corresponding Omega Mask index
                    # Matches marker name substring in omega filename
                    omega_idx = next((i for i, s in enumerate(omegaList) if f"{marker}" in s), None)

                    omega_name = next((s for i, s in enumerate(omegaList) if f"{marker}" in s), None)

                    compartment = omega_name.split('_')[1]

                    if compartment == 'nuclear':
                        pct_lowband = (50, 75)
                    elif compartment == 'extranuclear':
                        pct_lowband = (75,90)
                    
                    if omega_idx is None:
                        print(f"Warning: No omega mask found for marker {marker} in core {core_name}. Passing unchanged.")
                        output_list.append(I_afcorr)
                        continue

                    omega_mask = omega_stack[omega_idx].astype(bool)

                    # Calculate U_hat
                    print(f'Calculating U_hat for {core_name} - {marker} - {compartment} - {pct_lowband}')
                    U_hat, pct = decide_U_hat(I_afcorr[omega_mask], pct_lowband = pct_lowband, mid_thres=mid_thres)
                    print(f'Result: {U_hat}')

                    # Apply Correction
                    I_af_U_corr = (I_afcorr - U_hat).astype(I_afcorr.dtype)
                    output_list.append(I_af_U_corr)

                # Stack and Save
                output_array = np.stack(output_list, axis=0)
                np.save(f'{current_out_dir}/{core_name}.npy', output_array)

if __name__ == "__main__":
    main()