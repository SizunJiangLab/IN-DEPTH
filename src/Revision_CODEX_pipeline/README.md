# TMA Multiplexed Imaging Pipeline

Code for processing multiplexed immunofluorescence tissue microarray (TMA) images, as described in [Your Paper Title].

## Overview

This pipeline processes cyclic immunofluorescence (CycIF) TMA data through the following steps:

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_format_qptiff.py` | Convert raw QPTIFF to OME-TIFF format |
| 2 | `02_generate_geojson.py` | Automated TMA core detection (dearraying) |
| 3 | `03_crop_grid.py` | Tile images for parallel processing |
| 4 | `04_stitch_grid.py` | Reassemble processed tiles |
| 5 | `05_af_correction.py` | Autofluorescence correction using blank channels |
| 6 | `06_crop_cores.py` | Extract individual TMA cores |
| 7 | `07_staining_correction.py` | Unspecific staining background subtraction |
| 8 | `08_generate_ometiff.py` | Convert corrected arrays to OME-TIFF |
| 9 | `09_segmentation.py` | Cell segmentation using Mesmer |
| 10 | `10_extract_features.py` | Single-cell feature extraction |
| 11 | `11_combine_images.py` | Combine channels from multiple sources |

## Requirements

```bash
pip install -r requirements.txt
```

### pyqupath

This pipeline relies on [pyqupath](https://github.com/wuwenrui555/pyqupath) for OME-TIFF handling and TMA dearraying. Install from GitHub:

```bash
pip install git+https://github.com/wuwenrui555/pyqupath.git
```

**Note:** Step 9 (segmentation) requires TensorFlow and the Mesmer model. GPU is recommended.

## Usage

Each script is standalone. Edit the **CONFIGURATION** section at the top of each script to set your input/output paths, then run:

```bash
python scripts/01_format_qptiff.py
python scripts/02_generate_geojson.py
# ... continue through pipeline
```

### Configuration

Each script has a configuration block at the top:

```python
# =============================================================================
# CONFIGURATION - EDIT THESE PATHS
# =============================================================================
TMA_NAMES = ['Sample1', 'Sample2']
INPUT_DIR = '/path/to/raw/data'
OUTPUT_DIR = '/path/to/output'
```

### Shared Utilities

Common functions are in `utils/`:
- `utils/io_utils.py` - File listing, marker parsing
- `utils/image_utils.py` - Image processing, regression, masking
- `utils/thresholding.py` - Otsu, Yen, moment-preserving thresholds

## Data Structure

Expected input structure:
```
data/
├── raw/
│   ├── Sample1.qptiff
│   └── Sample2.qptiff
└── markers/
    ├── Sample1_markers.txt
    └── Sample2_markers.txt
```

Output structure:
```
output/
├── 01_ometiff/
├── 02_geojson/
├── 03_grid/
├── 04_stitched/
├── 05_af_corrected/
├── 06_cropped_cores/
├── 07_staining_corrected/
├── 08_core_ometiff/
├── 09_segmentation/
└── 10_features/
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `IMAGE_MPP` | 0.5068 | Microns per pixel |
| `ARTIFACT_PERCENTILE` | 99.0 | Percentile for artifact masking |
| `DILATION_RADIUS` | 2 | Morphological dilation for masks |
| `MAXIMA_THRESHOLD` | 0.075 | Mesmer cell detection sensitivity |
| `INTERIOR_THRESHOLD` | 0.2 | Mesmer cell boundary threshold |



## License

MIT License - see [LICENSE](LICENSE)
