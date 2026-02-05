"""Thresholding algorithms for the TMA pipeline."""

import math
import numpy as np
import cv2 as cv


def moment_threshold(image: np.ndarray, num_bins: int = 65536) -> float:
    """
    Tsai moment-preserving thresholding method.
    
    Reference:
    Tsai, W. H. (1985). Moment-preserving thresholding: A new approach.
    Computer Vision, Graphics, and Image Processing, 29(3), 377-393.
    
    Parameters
    ----------
    image : np.ndarray
        Input grayscale image
    num_bins : int
        Number of histogram bins
        
    Returns
    -------
    float
        Threshold value
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

    # Calculate moments
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


def yen_threshold(image: np.ndarray) -> float:
    """
    Yen's thresholding method (ImageJ implementation).
    
    References:
    1) Yen J.C., Chang F.J., and Chang S. (1995) "A New Criterion 
       for Automatic Multilevel Thresholding" IEEE Trans. on Image 
       Processing, 4(3): 370-378
    2) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
       Techniques and Quantitative Performance Evaluation" Journal of 
       Electronic Imaging, 13(1): 146-165
    
    Parameters
    ----------
    image : np.ndarray
        Input grayscale image (8-bit or 16-bit)
        
    Returns
    -------
    float
        Optimal threshold value
    """
    # Determine bit depth
    if image.dtype == np.uint8:
        bins = 256
        max_val = 255
        hist_range = (0, 255)
    elif image.dtype == np.uint16:
        bins = 65536
        max_val = 65535
        hist_range = (0, 65535)
    else:
        image = (image * 255).astype(np.uint8)
        bins = 256
        max_val = 255
        hist_range = (0, 255)
    
    # Calculate histogram
    hist, _ = np.histogram(image.flatten(), bins=bins, range=hist_range)
    total = np.sum(hist)
    
    if total == 0:
        return 0
    
    # Normalize histogram
    norm_histo = hist.astype(np.float64) / total
    
    # Calculate cumulative probabilities
    P1 = np.zeros(bins)
    P1[0] = norm_histo[0]
    for i in range(1, bins):
        P1[i] = P1[i-1] + norm_histo[i]
    
    # Calculate cumulative entropy sums
    P1_sq = np.zeros(bins)
    P2_sq = np.zeros(bins)
    
    # Forward cumulative sum of squared probabilities
    P1_sq[0] = norm_histo[0] ** 2 if norm_histo[0] > 0 else 0
    for i in range(1, bins):
        P1_sq[i] = P1_sq[i-1] + (norm_histo[i] ** 2 if norm_histo[i] > 0 else 0)
    
    # Backward cumulative sum
    P2_sq[bins-1] = norm_histo[bins-1] ** 2 if norm_histo[bins-1] > 0 else 0
    for i in range(bins-2, -1, -1):
        P2_sq[i] = P2_sq[i+1] + (norm_histo[i] ** 2 if norm_histo[i] > 0 else 0)
    
    # Find threshold maximizing criterion
    max_crit = -np.inf
    threshold = 0
    
    for t in range(bins):
        if P1[t] > 0 and P1[t] < 1:
            # Yen's criterion
            crit = -np.log(P1_sq[t] * P2_sq[t]) + 2 * np.log(P1[t] * (1 - P1[t]))
            
            if crit > max_crit:
                max_crit = crit
                threshold = t
    
    return threshold


def otsu_threshold(image: np.ndarray) -> float:
    """
    Otsu's thresholding method using OpenCV.
    
    Parameters
    ----------
    image : np.ndarray
        Input grayscale image
        
    Returns
    -------
    float
        Optimal threshold value
    """
    if image.dtype == np.uint16:
        # OpenCV Otsu works with 8-bit, so scale
        img_8bit = (image / 256).astype(np.uint8)
        thresh, _ = cv.threshold(img_8bit, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
        return thresh * 256
    else:
        if image.dtype != np.uint8:
            image = (image * 255).astype(np.uint8)
        thresh, _ = cv.threshold(image, 0, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
        return thresh


def percentile_threshold(image: np.ndarray, percentile: float = 95) -> float:
    """
    Simple percentile-based threshold.
    
    Parameters
    ----------
    image : np.ndarray
        Input image
    percentile : float
        Percentile value (0-100)
        
    Returns
    -------
    float
        Threshold value
    """
    return np.percentile(image, percentile)
