"""Image processing utility functions for the TMA pipeline."""

import numpy as np
import cv2 as cv
from scipy import ndimage
from skimage import morphology
from sklearn.linear_model import HuberRegressor
from typing import Tuple, Optional


def blank_estimator(blanks: np.ndarray, blank_cycles: list, target_cycle: int) -> np.ndarray:
    """
    Estimate blank channel at target cycle using linear interpolation.
    
    Parameters
    ----------
    blanks : np.ndarray
        Array of blank images, shape (n_blanks, H, W)
    blank_cycles : list
        Cycle numbers corresponding to each blank
    target_cycle : int
        Target cycle to estimate
        
    Returns
    -------
    np.ndarray
        Interpolated blank image
    """
    blank_cycles = np.array(blank_cycles)
    
    # Find bracketing blanks
    before_mask = blank_cycles <= target_cycle
    after_mask = blank_cycles >= target_cycle
    
    if not before_mask.any():
        return blanks[0]
    if not after_mask.any():
        return blanks[-1]
    
    before_idx = np.where(before_mask)[0][-1]
    after_idx = np.where(after_mask)[0][0]
    
    if before_idx == after_idx:
        return blanks[before_idx]
    
    # Linear interpolation
    t = (target_cycle - blank_cycles[before_idx]) / (blank_cycles[after_idx] - blank_cycles[before_idx])
    interpolated = (1 - t) * blanks[before_idx] + t * blanks[after_idx]
    
    return interpolated.astype(blanks.dtype)


def cell_area_mask(dapi_image: np.ndarray, 
                   method: str = 'otsu') -> Tuple[np.ndarray, np.ndarray]:
    """
    Create nuclear and extranuclear masks from DAPI image.
    
    Parameters
    ----------
    dapi_image : np.ndarray
        DAPI channel image
    method : str
        Thresholding method ('otsu', 'yen', 'moment')
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (nuclear_mask, extranuclear_mask) as boolean arrays
    """
    from .thresholding import otsu_threshold, yen_threshold, moment_threshold
    
    if method == 'otsu':
        thresh = otsu_threshold(dapi_image)
    elif method == 'yen':
        thresh = yen_threshold(dapi_image)
    elif method == 'moment':
        thresh = moment_threshold(dapi_image)
    else:
        thresh = otsu_threshold(dapi_image)
    
    nuclear_mask = dapi_image > thresh
    
    # Dilate nuclear mask to get cell area, then subtract for extranuclear
    dilated = morphology.binary_dilation(nuclear_mask, morphology.disk(5))
    extranuclear_mask = dilated & ~nuclear_mask
    
    return nuclear_mask, extranuclear_mask


def autofluorescence_regressor(marker_image: np.ndarray,
                                blank_image: np.ndarray,
                                omega_mask: np.ndarray,
                                use_huber: bool = False,
                                huber_epsilon: float = 1.35) -> Tuple[np.ndarray, float]:
    """
    Perform autofluorescence regression correction.
    
    Parameters
    ----------
    marker_image : np.ndarray
        Marker channel to correct
    blank_image : np.ndarray
        Estimated blank at same cycle
    omega_mask : np.ndarray
        Boolean mask of pixels to use for regression
    use_huber : bool
        Use robust Huber regression instead of WLS
    huber_epsilon : float
        Huber regression epsilon parameter
        
    Returns
    -------
    Tuple[np.ndarray, float]
        (corrected_image, beta_coefficient)
    """
    # Get pixels under mask
    y = marker_image[omega_mask].astype(np.float64)
    x = blank_image[omega_mask].astype(np.float64)
    
    if len(x) == 0 or x.max() == 0:
        return marker_image.astype(np.float32), 0.0
    
    if use_huber:
        # Robust regression
        reg = HuberRegressor(epsilon=huber_epsilon, fit_intercept=False)
        reg.fit(x.reshape(-1, 1), y)
        beta = reg.coef_[0]
    else:
        # Weighted least squares (weight by blank intensity)
        weights = x / x.max()
        weights = np.clip(weights, 0.01, 1.0)
        
        # Solve weighted least squares: beta = sum(w*x*y) / sum(w*x*x)
        numerator = np.sum(weights * x * y)
        denominator = np.sum(weights * x * x)
        
        if denominator == 0:
            beta = 0.0
        else:
            beta = numerator / denominator
    
    # Correct the image
    corrected = marker_image.astype(np.float64) - beta * blank_image.astype(np.float64)
    
    return corrected.astype(np.float32), float(beta)


def create_artifact_mask(image: np.ndarray,
                         percentile: float = 99.0,
                         dilation_radius: int = 2) -> np.ndarray:
    """
    Create mask for high-intensity artifacts.
    
    Parameters
    ----------
    image : np.ndarray
        Input image
    percentile : float
        Percentile threshold for artifact detection
    dilation_radius : int
        Radius for morphological dilation
        
    Returns
    -------
    np.ndarray
        Boolean mask (True = artifact)
    """
    threshold = np.percentile(image, percentile)
    mask = image > threshold
    
    if dilation_radius > 0:
        mask = morphology.binary_dilation(mask, morphology.disk(dilation_radius))
    
    return mask


def create_union_mask(*masks: np.ndarray) -> np.ndarray:
    """
    Create union of multiple boolean masks.
    
    Parameters
    ----------
    *masks : np.ndarray
        Boolean masks to combine
        
    Returns
    -------
    np.ndarray
        Union mask (True where any input is True)
    """
    if not masks:
        raise ValueError("At least one mask required")
    
    result = masks[0].copy()
    for mask in masks[1:]:
        result = result | mask
    
    return result


def decide_unspecific_threshold(image: np.ndarray,
                                omega_mask: np.ndarray,
                                compartment: str = 'extranuclear',
                                mid_threshold: float = 20,
                                nuclear_lowband: Tuple[int, int] = (50, 75),
                                extranuclear_lowband: Tuple[int, int] = (75, 90),
                                positive_thresholds: Tuple[int, int] = (10, 70)
                                ) -> Tuple[float, str]:
    """
    Decide unspecific staining threshold adaptively.
    
    Parameters
    ----------
    image : np.ndarray
        Marker image
    omega_mask : np.ndarray
        Compartment mask
    compartment : str
        'nuclear' or 'extranuclear'
    mid_threshold : float
        Threshold for quality classification
    nuclear_lowband : Tuple[int, int]
        Percentile range for nuclear markers
    extranuclear_lowband : Tuple[int, int]
        Percentile range for extranuclear markers
    positive_thresholds : Tuple[int, int]
        Percentile range for positive signal estimation
        
    Returns
    -------
    Tuple[float, str]
        (threshold_value, quality_label)
    """
    values = image[omega_mask]
    
    if len(values) == 0:
        return 0.0, 'bad'
    
    # Get percentile bands
    if compartment == 'nuclear':
        lowband = nuclear_lowband
    else:
        lowband = extranuclear_lowband
    
    low_pct = np.percentile(values, lowband[0])
    high_pct = np.percentile(values, lowband[1])
    
    # Estimate positive signal
    pos_low = np.percentile(values, positive_thresholds[0])
    pos_high = np.percentile(values, positive_thresholds[1])
    
    # Quality classification
    mid_val = np.percentile(values, mid_threshold)
    
    if mid_val < low_pct:
        quality = 'good'
        threshold = (low_pct + high_pct) / 2
    elif mid_val < high_pct:
        quality = 'moderate'
        threshold = mid_val
    else:
        quality = 'bad'
        threshold = high_pct
    
    return float(threshold), quality


def normalize_image(image: np.ndarray,
                    percentile_low: float = 1,
                    percentile_high: float = 99) -> np.ndarray:
    """
    Normalize image to 0-1 range using percentile clipping.
    
    Parameters
    ----------
    image : np.ndarray
        Input image
    percentile_low : float
        Lower percentile for clipping
    percentile_high : float
        Upper percentile for clipping
        
    Returns
    -------
    np.ndarray
        Normalized image (float32, 0-1 range)
    """
    low = np.percentile(image, percentile_low)
    high = np.percentile(image, percentile_high)
    
    if high == low:
        return np.zeros_like(image, dtype=np.float32)
    
    normalized = (image.astype(np.float32) - low) / (high - low)
    normalized = np.clip(normalized, 0, 1)
    
    return normalized
