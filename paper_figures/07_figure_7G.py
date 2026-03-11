"""
Spatial Transcriptomics Visualization for LMP1+/- Tumor-Macrophage Interactions

This script provides 5 main visualization functions:
1. Cell type annotation overview
2. Macrophage-centered L-R expression
3. AUCell scores per cell
4. Tumor-macrophage L-R pairs (circles)
5. Tumor-macrophage L-R pairs (lines)

env: wrplot
"""

# %% Import Libraries ==========
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import tifffile
from matplotlib.patches import Circle
from skimage.segmentation import find_boundaries
from tqdm import tqdm

# %% Constants ==========
TQDM_FORMAT = "{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
RADIUS = 40

BACKGROUND_COLOR_RGB = tuple(int(192) for _ in range(3))
BOUNDARY_COLOR_RGB = tuple(int(128) for _ in range(3))

# Color mapping for cell annotations
COLOR_MAP = {
    "lmp1p_tumor_hopped": "#cd2829",
    "lmp1p_macrophage": "#dc6eb7",
    "lmp1n_tumor_hopped": "#226ba7",
    "lmp1n_macrophage": "#885db0",
}

# Colormap for heatmap visualizations (e.g., expression levels, scores)
# Options: 'viridis', 'Reds', 'Blues', 'RdYlBu_r', 'coolwarm', 'plasma', 'inferno', etc.
HEATMAP_COLORMAP = "viridis"

UP_QUANTILE = 0.95
DOWN_QUANTILE = 0


# %% Data Paths ==========

data_dir = Path(__file__).parent.parent / "src" / "07_figure_7_neighborhood_analysis"

# the directory should be the "06_figure_6_Segmentation_mask" folder in 10.5281/zenodo.18379156
seg_dir = ""

output_root = Path(__file__).parent / "figures" / "03_squidpy"
output_root.mkdir(parents=True, exist_ok=True)


# %% Load All Data ==========

# Load spatial transcriptomics data
print("Loading spatial transcriptomics data...")
adata = sc.read_h5ad(data_dir / "01_adata_with_spatial.h5ad")

# Add combined sample ID column
adata.obs["smp_id_fov"] = (
    adata.obs["smp_id"].astype(str) + "_" + adata.obs["fov_id"].astype(str)
)

# Add LMP1 status column
adata.obs["lmp1"] = adata.obs["annotation_3"].apply(
    lambda x: "lmp1p" if "lmp1p" in x else ("lmp1n" if "lmp1n" in x else None)
)

# Load macrophage-centered L-R pair mean expression
print("Loading macrophage-centered L-R pair data...")
lr_pairs = [
    "CCL22::CCR4",
    "EBI3::IL6ST",
    "EBI3::IL27RA",
    "IL6::IL6ST",
    "EBI3::STAT3",
    "ICOSLG::ICOS",
    "ICAM3::CD209",
    "DCN::TLR4",
    "CCN1::TLR4",
]

df_mac_lr_mean = pd.read_feather(data_dir / "06_macrophage_center_lr_mean.feather")[
    ["cell_id_mac", "x_mac", "y_mac", "lmp1"] + lr_pairs
]
df_mac_lr_mean.rename(
    columns={"cell_id_mac": "center", "x_mac": "x", "y_mac": "y"},
    inplace=True,
)

# Load macrophage-centered AUCell scores
print("Loading macrophage-centered AUCell scores...")
df_aucell_mac_centered = pd.read_feather(
    data_dir / "07_macrophage_centered_aucell_scores.feather"
)

# Merge macrophage-centered data
df_mac_center = df_mac_lr_mean.merge(df_aucell_mac_centered, on=["center"])
df_mac_center = df_mac_center.merge(
    adata.obs[["smp_id_fov"]], left_on="center", right_index=True
)

# Load AUCell scores for all cells
print("Loading AUCell scores for all cells...")
df_aucell_scores = pd.read_feather(data_dir / "07_aucell_scores.feather")

# Add combined sample ID column
df_aucell_scores["smp_id_fov"] = (
    df_aucell_scores["smp_id"].astype(str)
    + "_"
    + df_aucell_scores["fov_id"].astype(str)
)

# Add LMP1 status column
df_aucell_scores["lmp1"] = df_aucell_scores["annotation_3"].apply(
    lambda x: "lmp1p" if "lmp1p" in x else ("lmp1n" if "lmp1n" in x else None)
)

# Subset for macrophages and tumor cells
df_aucell_mac = df_aucell_scores[
    df_aucell_scores["annotation_3"].isin(["lmp1p_macrophage", "lmp1n_macrophage"])
]
df_aucell_tumor = df_aucell_scores[
    df_aucell_scores["annotation_3"].isin(["lmp1p_tumor_hopped", "lmp1n_tumor_hopped"])
]

# Load tumor-macrophage L-R pair mean expression
print("Loading tumor-macrophage L-R pair data...")
meta_cols = [
    "smp_id",
    "smp_id_fov",
    "cell_id_mac",
    "cell_id_tumor",
    "annotation_mac",
    "annotation_tumor",
    "x_mac",
    "y_mac",
    "x_tumor",
    "y_tumor",
    "x_lr",
    "y_lr",
    "lmp1",
]
df_lr_mean_tumor_mac = pd.read_feather(data_dir / "06_tumor_macrophage_lr_mean.feather")

# Add sample ID with FOV column
df_lr_mean_tumor_mac = df_lr_mean_tumor_mac.merge(
    adata.obs["smp_id_fov"], left_on="cell_id_mac", right_index=True
)
df_lr_mean_tumor_mac = df_lr_mean_tumor_mac[meta_cols + lr_pairs]

# Load tumor-macrophage ligand and receptor expression
print("Loading tumor-macrophage ligand and receptor expression data...")
df_lr_expr = pd.read_feather(data_dir / "06_tumor_macrophage_lr_expr.feather")

# Add sample ID with FOV column
df_lr_expr = df_lr_expr.merge(
    adata.obs["smp_id_fov"], left_on="cell_id_mac", right_index=True
)

lr_cols = []
for lr_pair in lr_pairs:
    ligand, receptor = lr_pair.split("::")
    lr_cols.append(ligand)
    lr_cols.append(receptor)
lr_cols = list(set(lr_cols))

df_lr_expr = df_lr_expr[meta_cols + lr_cols]

print("All data loaded successfully!")


# %% Get All Samples for Visualization ==========

n_annotation = (
    adata.obs.groupby(["smp_id_fov", "annotation_3"], observed=True)
    .size()
    .reset_index(name="counts")
    .pivot_table(
        index="smp_id_fov",
        columns="annotation_3",
        values="counts",
        fill_value=0,
        observed=True,
    )
)

mask = (
    (n_annotation["lmp1n_macrophage"] > 0)
    & (n_annotation["lmp1p_macrophage"] > 0)
    & (n_annotation["lmp1n_tumor_hopped"] > 0)
    & (n_annotation["lmp1p_tumor_hopped"] > 0)
)

smp_id_fov_list = n_annotation.index[mask].tolist()


# %% Functions ==========


def ax_plot_legend(
    color_dict,
    marker="s",
    marker_size=10,
    marker_edgecolor="black",
    loc="center left",
    bbox_to_anchor=(1.05, 0.5),
    text_size=10,
    ax=None,
):
    """
    Plot a legend for a given color dictionary on a specified axes object.

    Parameters
    ----------
    color_dict : dict
        A dictionary where keys are labels and values are colors (hex codes).
    marker : str, optional
        The marker style to use in the legend. Default is "s" (square).
    marker_size : int, optional
        The size of the markers in the legend. Default is 10.
    marker_edgecolor : str, optional
        The edge color of the markers in the legend. Default is "black".
    loc : str, optional
        The location of the legend. Default is "center left".
    bbox_to_anchor : tuple, optional
        The bounding box anchor for the legend. Default is (1.05, 0.5).
    text_size : int, optional
        The font size of the legend text. Default is 10.
    ax : matplotlib.axes.Axes, optional
        The axes object to plot the legend on. If None, uses the current axes.
        Default is None.
    """
    if ax is None:
        ax = plt.gca()

    # Create legend handles with custom markers for each label
    handles = [
        mlines.Line2D(
            [0],
            [0],
            marker=marker,
            markerfacecolor=color_dict[label],
            markeredgecolor=marker_edgecolor,
            markersize=marker_size,
            label=label,
            linestyle="None",  # No line, only marker
        )
        for label in color_dict.keys()
    ]

    # Plot the legend with specified parameters
    ax.legend(
        handles=handles,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
        fontsize=text_size,
    )


def create_rgb_annotation(
    segmentation_mask: np.ndarray,
    annotation_dict: dict[int, str],
    color_dict: dict[str, str],
    boundary: bool = False,
    boundary_color_rgb: tuple[int, int, int] = (0, 0, 0),
    background_color_rgb: tuple[int, int, int] = (0, 0, 0),
) -> np.ndarray:
    """
    Create an RGB image for segmentation annotations.

    Parameters
    ----------
    segmentation_mask : np.ndarray
        2D array representing the segmentation mask.
    annotation_dict : dict[int, str]
        A dictionary mapping label integers to annotation names.
    color_dict : dict[str, str]
        A dictionary mapping annotation names to colors.
    outline : bool, optional
        If True, outlines the segmentation mask in white, by default False.
    color_bg : tuple[int, int, int], optional
        Background color for the RGB image, by default (255, 255, 255).

    Returns
    -------
    np.ndarray
        RGB image of annotations and outline (if enabled).
    """
    # Map label integers to RGB colors based on annotations
    unique_labels = np.unique(segmentation_mask[segmentation_mask != 0])
    color_mapping = np.zeros((np.max(segmentation_mask) + 1, 3), dtype=np.uint8)

    # Set background color for label 0 and unmapped labels
    color_mapping[:] = background_color_rgb

    for label in tqdm(
        unique_labels, desc="Mapping labels to colors", bar_format=TQDM_FORMAT
    ):
        annotation = annotation_dict.get(label, None)
        color_hex = color_dict.get(annotation, None)
        if color_hex is not None:
            color_rgb = tuple(int(color_hex[i : i + 2], 16) for i in (1, 3, 5))
            color_mapping[label] = color_rgb

    # Create RGB image
    rgb_image = color_mapping[segmentation_mask]

    # Add cell boundaries if outline option is enabled
    if boundary:
        segmentation_boundaries = find_segmentation_boundaries(segmentation_mask)
        rgb_image[segmentation_boundaries.astype(bool)] = boundary_color_rgb

    return rgb_image


def find_segmentation_boundaries(
    segmentation_mask: np.ndarray,
    mode: str = "inner",
    connectivity: int = 1,
    background: int = 0,
) -> np.ndarray:
    """
    Convert a segmentation mask to a mask for boundaries.

    Parameters
    ----------
    segmentation_mask : np.ndarray
        An array in which different regions are labeled with different integers.
    mode : str, optional
        The mode of boundary detection. Options are 'inner', 'outer', 'thick',
        and 'subpixel'. Default is 'inner'.
    connectivity : int, optional
        The connectivity defining the neighborhood of a pixel. Default is 1.
    background : int, optional
        The value representing the background in the segmentation mask. Default
        is 0.

    Returns
    -------
    np.ndarray
        An array with same shape as `segmentation_mask` and same dtype, where
        values of 0 represent background pixels and other values represent
        the boundaries of labels.
    """
    boundaries = find_boundaries(
        segmentation_mask,
        mode=mode,
        connectivity=connectivity,
        background=background,
    )
    segmentation_boundary = segmentation_mask * boundaries

    return segmentation_boundary


def load_and_crop_segmentation(
    smp_id_fov,
    seg_dir,
    adata,
    expand_ratio=0.05,
    background_color_rgb=None,
    boundary_color_rgb=None,
):
    """
    Load segmentation mask and crop to bounding box of cells in the sample.

    Parameters
    ----------
    smp_id_fov : str
        Sample ID with FOV (format: tma_id_core_id_fov_id)
    seg_dir : Path
        Directory containing segmentation masks
    adata : AnnData
        Annotated data object with spatial information
    expand_ratio : float, default=0.05
        Fraction to expand bounding box on each side
    background_color_rgb : tuple, optional
        RGB color for background (defaults to BACKGROUND_COLOR_RGB)
    boundary_color_rgb : tuple, optional
        RGB color for cell boundaries (defaults to BOUNDARY_COLOR_RGB)

    Returns
    -------
    seg_sm : np.ndarray
        Cropped segmentation mask
    rgb_img_bg : np.ndarray
        Background RGB image with boundaries
    x_min : int
        Minimum x coordinate of the cropped region
    y_min : int
        Minimum y coordinate of the cropped region
    """
    if background_color_rgb is None:
        background_color_rgb = BACKGROUND_COLOR_RGB
    if boundary_color_rgb is None:
        boundary_color_rgb = BOUNDARY_COLOR_RGB

    tma_id, core_id, fov_id = smp_id_fov.split("_")

    # Load segmentation mask
    print("Loading segmentation mask...")
    seg_f = seg_dir / tma_id / "seg_mask_0.075_0.2" / f"{core_id}.tiff"
    seg = tifffile.imread(seg_f)

    # Find bounding box for all cells in the sample
    cell_ids = adata.obs[adata.obs.smp_id_fov == smp_id_fov].cell_id.tolist()
    seg_cell = np.isin(seg, cell_ids)

    ys, xs = np.where(seg_cell)
    x_min, x_max = xs.min(), xs.max()
    y_min, y_max = ys.min(), ys.max()

    # Expand bounding box by expand_ratio
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_min = max(0, int(x_min - expand_ratio * x_range))
    x_max = min(seg.shape[1] - 1, int(x_max + expand_ratio * x_range))
    y_min = max(0, int(y_min - expand_ratio * y_range))
    y_max = min(seg.shape[0] - 1, int(y_max + expand_ratio * y_range))

    # Crop segmentation to bounding box
    seg_sm = seg[y_min : y_max + 1, x_min : x_max + 1]

    print(f"Cropped segmentation shape: {seg_sm.shape}")

    # Create background image for later visualizations
    rgb_img_bg = create_rgb_annotation(
        seg_sm,
        {},
        {},
        boundary=True,
        boundary_color_rgb=boundary_color_rgb,
        background_color_rgb=background_color_rgb,
    )

    return seg_sm, rgb_img_bg, x_min, y_min


def plot_cell_annotation_overview(
    smp_id_fov,
    adata,
    seg_sm,
    x_min,
    y_min,
    line_width=1,
    line_alpha=0.5,
    color_map=None,
    radius=None,
    boundary_color_rgb=None,
    background_color_rgb=None,
):
    """
    Plot cell type annotation overview with LMP1- and LMP1+ cells side by side.

    Parameters
    ----------
    smp_id_fov : str
        Sample ID with FOV to visualize
    adata : AnnData
        Annotated data object with spatial information
    seg_sm : np.ndarray
        Cropped segmentation mask
    x_min : int
        Minimum x coordinate of the cropped region
    y_min : int
        Minimum y coordinate of the cropped region
    line_width : float, default=1
        Width of circle edge lines
    line_alpha : float, default=0.5
        Transparency of circle edge lines
    color_map : dict, optional
        Color mapping for cell annotations (defaults to COLOR_MAP)
    radius : int, optional
        Circle radius for macrophage visualization (defaults to RADIUS)
    boundary_color_rgb : tuple, optional
        RGB color for cell boundaries (defaults to BOUNDARY_COLOR_RGB)
    background_color_rgb : tuple, optional
        RGB color for background (defaults to BACKGROUND_COLOR_RGB)

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    if color_map is None:
        color_map = COLOR_MAP
    if radius is None:
        radius = RADIUS
    if boundary_color_rgb is None:
        boundary_color_rgb = BOUNDARY_COLOR_RGB
    if background_color_rgb is None:
        background_color_rgb = BACKGROUND_COLOR_RGB

    df_plot = adata.obs[adata.obs.smp_id_fov == smp_id_fov]
    annotation_dict = df_plot.set_index("cell_id")["annotation_3"].to_dict()

    # Create RGB annotation image for different cell groups
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    axes = axes.flatten()

    # Plot 1: LMP1- cells (tumor and macrophages) - no legend
    ax = axes[0]
    df_plot_sm = df_plot[
        df_plot.annotation_3.isin(["lmp1n_tumor_hopped", "lmp1n_macrophage"])
    ]
    cell_ids_sm = df_plot_sm.cell_id.values
    annotation_dict_sm = {
        cid: annotation_dict[cid] for cid in cell_ids_sm if cid in annotation_dict
    }
    color_dict_sm = {
        k: v
        for k, v in color_map.items()
        if k in ["lmp1n_tumor_hopped", "lmp1n_macrophage"]
    }

    rgb_img = create_rgb_annotation(
        seg_sm,
        annotation_dict_sm,
        color_dict_sm,
        boundary=True,
        boundary_color_rgb=boundary_color_rgb,
        background_color_rgb=background_color_rgb,
    )

    ax.imshow(rgb_img)
    # Draw circles around macrophages
    df_mac = df_plot_sm[df_plot_sm.annotation_3 == "lmp1n_macrophage"]
    for _, row in df_mac.iterrows():
        circle = Circle(
            (row["x_cent"] - x_min, row["y_cent"] - y_min),
            radius,
            fill=False,
            edgecolor=color_map["lmp1n_macrophage"],
            linewidth=line_width,
            alpha=line_alpha,
        )
        ax.add_patch(circle)
    ax.set_title(f"LMP1- cells in sample {smp_id_fov}")
    ax.axis("off")

    # Plot 2: LMP1+ cells (tumor and macrophages) - no legend
    ax = axes[1]
    df_plot_sm = df_plot[
        df_plot.annotation_3.isin(["lmp1p_tumor_hopped", "lmp1p_macrophage"])
    ]
    cell_ids_sm = df_plot_sm.cell_id.values
    annotation_dict_sm = {
        cid: annotation_dict[cid] for cid in cell_ids_sm if cid in annotation_dict
    }
    color_dict_sm = {
        k: v
        for k, v in color_map.items()
        if k in ["lmp1p_tumor_hopped", "lmp1p_macrophage"]
    }

    rgb_img = create_rgb_annotation(
        seg_sm,
        annotation_dict_sm,
        color_dict_sm,
        boundary=True,
        boundary_color_rgb=boundary_color_rgb,
        background_color_rgb=background_color_rgb,
    )

    ax.imshow(rgb_img)
    # Draw circles around macrophages
    df_mac = df_plot_sm[df_plot_sm.annotation_3 == "lmp1p_macrophage"]
    for _, row in df_mac.iterrows():
        circle = Circle(
            (row["x_cent"] - x_min, row["y_cent"] - y_min),
            radius,
            fill=False,
            edgecolor=color_map["lmp1p_macrophage"],
            linewidth=line_width,
            alpha=line_alpha,
        )
        ax.add_patch(circle)
    ax.set_title(f"LMP1+ cells in sample {smp_id_fov}")
    ax.axis("off")

    # Plot 3: Empty panel for legend only
    ax = axes[2]
    ax.axis("off")

    cbar_ax = ax.inset_axes([0.3, 0.1, 0.05, 0.8])  # [x, y, width, height]
    cbar_ax.axis("off")

    ax_plot_legend(color_map, ax=cbar_ax)

    plt.tight_layout()

    return fig

def plot_tumor_macrophage_lr_expression_2(
    smp_id_fov,
    score,
    df_mac_center,
    rgb_img_bg,
    x_min,
    y_min,
    down_quantile=None,
    up_quantile=None,
    sigma=20,
    alpha=0.6,
    cmap=None,
    radius=None,
):
    """
    Plot tumor-macrophage L-R pair expression using PSF method

    Creates a heatmap by overlaying circular regions (point spread functions)
    centered at each macrophage position. Each center contributes a weighted
    disk to the heatmap, where the weight is determined by its L-R pair
    expression score. Multiple overlapping disks accumulate to create a
    continuous intensity map showing both spatial distribution and L-R activity.

    The visualization workflow:
    1. For each center, create a soft-edged circular disk (PSF kernel)
    2. Weight the disk by the raw L-R expression score of the center
    3. Accumulate all weighted disks to form a continuous heatmap
    4. Normalize heatmap intensity using quantile-based display range
    5. Map normalized intensities to colors and overlay on background image

    Parameters
    ----------
    smp_id_fov : str
        Sample ID with FOV to visualize (format: tma_id_core_id_fov_id)
    score : str
        L-R pair name to visualize (e.g., "CCL22::CCR4")
    df_mac_center : pd.DataFrame
        DataFrame containing macrophage-centered L-R pair mean expression data.
        Must include columns: 'smp_id_fov', 'x', 'y', 'lmp1', and the
        specified score column.
    rgb_img_bg : np.ndarray
        Background RGB image (typically showing cell boundaries). Shape: (H, W, 3)
    x_min, y_min : int
        Minimum x and y coordinates of the cropped region, used to adjust
        macrophage positions to the cropped coordinate system
    down_quantile : float, optional
        Lower quantile for colorbar display range. Values below this quantile
        will be mapped to the darkest color. Default: DOWN_QUANTILE (0.0)
    up_quantile : float, optional
        Upper quantile for colorbar display range. Values above this quantile
        will be mapped to the brightest color. Default: UP_QUANTILE (0.95)
    sigma : float, default=20
        Controls the width of the soft edge falloff in pixels. Determines how
        gradually the disk intensity decreases from center to edge:
        - Small sigma (e.g., 5): Sharp edge, disk intensity drops quickly
        - Medium sigma (e.g., 20): Moderate edge softness
        - Large sigma (e.g., 100): Very soft edge, disk appears diffuse
        - Very large sigma (e.g., 1000): Nearly uniform disk, minimal falloff
    alpha : float, default=0.6
        Transparency of the heatmap overlay (0=fully transparent, 1=fully opaque)
    cmap : str, optional
        Matplotlib colormap name for heatmap visualization.
        Default: HEATMAP_COLORMAP (typically 'viridis')
    radius : int, optional
        Radius of the circular disk kernel in pixels. Defines the spatial extent
        of each macrophage's influence. Default: RADIUS (40 pixels)

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure object containing three subplots:
        - Left: LMP1- macrophage heatmap
        - Middle: LMP1+ macrophage heatmap
        - Right: Colorbar showing accumulated intensity scale
    """
    # Set defaults
    if down_quantile is None:
        down_quantile = DOWN_QUANTILE
    if up_quantile is None:
        up_quantile = UP_QUANTILE
    if cmap is None:
        cmap = HEATMAP_COLORMAP
    if radius is None:
        radius = RADIUS

    # Prepare and filter data
    df_plot = df_mac_center[df_mac_center["smp_id_fov"] == smp_id_fov].copy()
    df_plot = df_plot[["x", "y", "lmp1", score]].dropna(subset=["lmp1"])
    df_plot.rename(columns={score: "score"}, inplace=True)

    # Adjust coordinates to cropped region
    df_plot["x"] -= x_min
    df_plot["y"] -= y_min

    # Calculate score range
    vmax = df_plot["score"].max()
    vmin = df_plot["score"].min()

    # Create soft-edged circular disk kernel
    half_size = radius
    y_grid, x_grid = np.ogrid[-half_size : half_size + 1, -half_size : half_size + 1]
    radius_map = np.sqrt(x_grid**2 + y_grid**2)

    # Disk with linear edge falloff
    disk_kernel = np.clip(
        (half_size - radius_map) / sigma if sigma > 0 else (radius_map <= half_size),
        0,
        1,
    )

    # Generate heatmaps for both LMP1 groups
    def create_heatmap(df_data):
        heatmap = np.zeros(rgb_img_bg.shape[:2], dtype=float)

        for _, row in df_data.iterrows():
            x_center = int(round(row["x"]))
            y_center = int(round(row["y"]))

            # Normalize score with minimum visibility
            if vmax > vmin:
                norm_score = np.clip((row["score"] - vmin) / (vmax - vmin), 10e-5, 1)
            else:
                norm_score = 10e-5

            # Calculate region bounds
            y_start = max(0, y_center - half_size)
            y_end = min(heatmap.shape[0], y_center + half_size + 1)
            x_start = max(0, x_center - half_size)
            x_end = min(heatmap.shape[1], x_center + half_size + 1)

            # Corresponding kernel bounds
            ky_start = half_size - (y_center - y_start)
            ky_end = half_size + (y_end - y_center)
            kx_start = half_size - (x_center - x_start)
            kx_end = half_size + (x_end - x_center)

            # Accumulate weighted disk
            heatmap[y_start:y_end, x_start:x_end] += (
                disk_kernel[ky_start:ky_end, kx_start:kx_end] * norm_score
            )

        return heatmap

    # Split by LMP1 status and create heatmaps
    heatmaps = {}
    for lmp1_status in ["lmp1n", "lmp1p"]:
        df_group = df_plot[df_plot["lmp1"] == lmp1_status]
        if len(df_group) > 0:
            heatmaps[lmp1_status] = create_heatmap(df_group)

    # Calculate display range based on combined heatmaps
    if heatmaps:
        all_values = np.concatenate([hm[hm > 0].flatten() for hm in heatmaps.values()])
        vmin_display = np.quantile(all_values, down_quantile)
        vmax_display = np.quantile(all_values, up_quantile)
    else:
        vmin_display = 0
        vmax_display = 1

    # Create visualization
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    cmap_obj = plt.get_cmap(cmap)

    # Plot both LMP1 groups
    for idx, (lmp1_status, title) in enumerate(
        [("lmp1n", "LMP1-"), ("lmp1p", "LMP1+")]
    ):
        ax = axes[idx]
        ax.imshow(rgb_img_bg)

        if lmp1_status in heatmaps:
            heatmap = heatmaps[lmp1_status]
            # Use display range to normalize heatmap
            normalized_heatmap = np.clip(
                (heatmap - vmin_display) / (vmax_display - vmin_display), 0, 1
            )
            rgba = cmap_obj(normalized_heatmap)
            rgba[..., 3] = np.where(heatmap > 0, alpha, 0)
            ax.imshow(rgba)

        ax.set_title(f"{title} Macrophage {score} PSF")
        ax.axis("off")
        ax.set_aspect("equal")

    # Add colorbar to third panel
    ax = axes[2]
    ax.axis("off")
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin_display, vmax=vmax_display)
    )
    sm.set_array([])
    cbar_ax = ax.inset_axes([0.3, 0.1, 0.05, 0.8])
    fig.colorbar(sm, cax=cbar_ax, label=f"{score} accumulated intensity")

    plt.tight_layout()
    return fig

def plot_ligand_receptor_expression(
    smp_id_fov,
    score,
    df_lr_expr,
    df_lr_mean_tumor_mac,
    seg_sm,
    x_min,
    y_min,
    line_width=2,
    line_alpha=0.75,
    line_exclude_quantile=0.05,
    down_quantile=None,
    up_quantile=None,
    cmap=None,
    boundary_color_rgb=None,
    background_color_rgb=None,
    plot_connections=False,
):
    """
    Plot ligand and receptor expression for tumor and macrophage cells.

    This function visualizes the expression of ligand-receptor pairs by coloring
    cells according to their expression levels. Optionally, it can also plot
    connections between tumor and macrophage cells.

    Parameters
    ----------
    smp_id_fov : str
        Sample ID with FOV to visualize (format: tma_id_core_id_fov_id)
    score : str
        L-R pair name to visualize (e.g., "EBI3::IL27RA")
    df_lr_expr : pd.DataFrame
        DataFrame with ligand and receptor expression data per cell
    df_lr_mean_tumor_mac : pd.DataFrame
        DataFrame with tumor-macrophage L-R pair mean expression
    seg_sm : np.ndarray
        Cropped segmentation mask
    x_min : int
        Minimum x coordinate of the cropped region
    y_min : int
        Minimum y coordinate of the cropped region
    line_width : float, default=2
        Width of connection lines
    line_alpha : float, default=0.75
        Transparency of connection lines
    line_exclude_quantile : float, default=0.05
        Quantile threshold to exclude low-expression connections
    down_quantile : float, optional
        Lower quantile for colorbar range (defaults to DOWN_QUANTILE)
    up_quantile : float, optional
        Upper quantile for colorbar range (defaults to UP_QUANTILE)
    cmap : str, optional
        Colormap name (not used in current implementation, uses 'Reds' and 'Blues')
    boundary_color_rgb : tuple, optional
        RGB color for cell boundaries (defaults to BOUNDARY_COLOR_RGB)
    background_color_rgb : tuple, optional
        RGB color for background (defaults to BACKGROUND_COLOR_RGB)
    plot_connections : bool, default=False
        Whether to plot tumor-macrophage connection lines

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    if down_quantile is None:
        down_quantile = DOWN_QUANTILE
    if up_quantile is None:
        up_quantile = UP_QUANTILE
    if boundary_color_rgb is None:
        boundary_color_rgb = BOUNDARY_COLOR_RGB
    if background_color_rgb is None:
        background_color_rgb = BACKGROUND_COLOR_RGB

    ligand, receptor = score.split("::")

    color_dict_lmp1n = {}
    color_dict_lmp1p = {}

    annotation_dict_lmp1n = {}
    annotation_dict_lmp1p = {}

    sm_list = []

    for score_lr, cell_id_col, cmap in [
        (ligand, "cell_id_tumor", "Reds"),
        (receptor, "cell_id_mac", "Blues"),
    ]:
        # Prepare data for the selected sample
        df_plot = df_lr_expr.loc[
            df_lr_expr["smp_id_fov"] == smp_id_fov, [cell_id_col, "lmp1", score_lr]
        ].rename(columns={score_lr: "score", cell_id_col: "cell_id"})
        df_plot["cell_id"] = (
            df_plot["cell_id"].str.extract(r"[^_]+_[^_]+_(\d+)$")[0]
        ).astype(int)

        # Calculate score range using quantiles
        vmin = df_plot["score"].quantile(down_quantile)
        vmax = df_plot["score"].quantile(up_quantile)
        sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=plt.Normalize(vmin=vmin, vmax=vmax),
        )
        sm.set_array([])
        sm_list.append(sm)

        # Create cell_id to score mapping
        cell_to_score = df_plot.set_index("cell_id")["score"].to_dict()

        # Split by LMP1 status
        df_lmp1n = df_plot[df_plot["lmp1"] == "lmp1n"]
        df_lmp1p = df_plot[df_plot["lmp1"] == "lmp1p"]

        # Create annotation_dict (cell_id -> cell_id) for each LMP1 group
        annotation_dict_lmp1n_lr = {cell_id: cell_id for cell_id in df_lmp1n["cell_id"]}
        annotation_dict_lmp1p_lr = {cell_id: cell_id for cell_id in df_lmp1p["cell_id"]}
        annotation_dict_lmp1n.update(annotation_dict_lmp1n_lr)
        annotation_dict_lmp1p.update(annotation_dict_lmp1p_lr)

        # Create color_dict (cell_id -> hex color) based on scores
        for cell_id, score_val in cell_to_score.items():
            # Normalize score to [0, 1] range
            normalized_score = (score_val - vmin) / (vmax - vmin) if vmax > vmin else 0
            normalized_score = np.clip(normalized_score, 0, 1)  # Clip to valid range

            # Get RGB color from colormap
            color_rgba = plt.get_cmap(cmap)(normalized_score)
            color_rgb = tuple(int(255 * c) for c in color_rgba[:3])

            # Convert RGB to hex
            color_hex = f"#{color_rgb[0]:02x}{color_rgb[1]:02x}{color_rgb[2]:02x}"

            # Assign to appropriate dict based on LMP1 status
            if cell_id in annotation_dict_lmp1n_lr:
                color_dict_lmp1n[cell_id] = color_hex
            elif cell_id in annotation_dict_lmp1p_lr:
                color_dict_lmp1p[cell_id] = color_hex

    # Create RGB images using create_rgb_annotation
    print("  Creating RGB image for LMP1- cells...")
    rgb_img_lmp1n = create_rgb_annotation(
        seg_sm,
        annotation_dict_lmp1n,
        color_dict_lmp1n,
        boundary=True,
        boundary_color_rgb=boundary_color_rgb,
        background_color_rgb=background_color_rgb,
    )

    print("  Creating RGB image for LMP1+ cells...")
    rgb_img_lmp1p = create_rgb_annotation(
        seg_sm,
        annotation_dict_lmp1p,
        color_dict_lmp1p,
        boundary=True,
        boundary_color_rgb=boundary_color_rgb,
        background_color_rgb=background_color_rgb,
    )

    fig, axes = plt.subplots(1, 3, figsize=(24, 8))

    # Plot LMP1- cells
    ax = axes[0]
    ax.imshow(rgb_img_lmp1n)
    ax.set_title(f"LMP1- {score}, Ligand-Receptor-Connections")
    ax.axis("off")
    ax.set_aspect("equal")

    # Plot LMP1+ cells
    ax = axes[1]
    ax.imshow(rgb_img_lmp1p)
    ax.set_title(f"LMP1+ {score}, Ligand-Receptor-Connections")
    ax.axis("off")
    ax.set_aspect("equal")

    # Plot legends
    ax = axes[2]
    ax.axis("off")
    cbar_ax_l = ax.inset_axes([0.3, 0.7, 0.05, 0.2])
    fig.colorbar(sm_list[0], cax=cbar_ax_l, label="Ligand expression")
    cbar_ax_r = ax.inset_axes([0.3, 0.4, 0.05, 0.2])
    fig.colorbar(sm_list[1], cax=cbar_ax_r, label="Receptor expression")

    if not plot_connections:
        plt.tight_layout()
        return fig

    # L-R connections if needed
    cmap = "viridis"

    # Prepare data for the selected sample
    df_plot = (
        df_lr_mean_tumor_mac.loc[
            df_lr_mean_tumor_mac["smp_id_fov"] == smp_id_fov,
            ["x_lr", "y_lr", "x_mac", "y_mac", "x_tumor", "y_tumor", "lmp1", score],
        ]
        .copy()
        .rename(columns={score: "score"})
    )
    df_plot.dropna(subset=["lmp1"], inplace=True)
    df_plot.sort_values(by="score", ascending=True, inplace=True)

    # Calculate colorbar limits
    vmin = df_plot["score"].quantile(down_quantile)
    vmax = df_plot["score"].quantile(up_quantile)
    line_exclude_threshold = df_plot["score"].quantile(line_exclude_quantile)

    df_plot_lmp1n = df_plot[df_plot.lmp1 == "lmp1n"]
    df_plot_lmp1p = df_plot[df_plot.lmp1 == "lmp1p"]

    # Plot LMP1- tumor-macrophage connections (left) - no legend
    ax = axes[0]
    for _, row in df_plot_lmp1n.iterrows():
        if row["score"] < line_exclude_threshold:
            continue

        # Calculate color based on score
        normalized_score = (row["score"] - vmin) / (vmax - vmin)
        color = plt.get_cmap(cmap)(normalized_score)

        # Draw line connecting macrophage and tumor
        ax.plot(
            [row["x_mac"] - x_min, row["x_tumor"] - x_min],
            [row["y_mac"] - y_min, row["y_tumor"] - y_min],
            color=color,
            linewidth=line_width,
            alpha=line_alpha,
        )

    # Plot LMP1+ tumor-macrophage connections (middle) - no legend
    ax = axes[1]
    for _, row in df_plot_lmp1p.iterrows():
        if row["score"] < line_exclude_threshold:
            continue

        # Calculate color based on score
        normalized_score = (row["score"] - vmin) / (vmax - vmin)
        color = plt.get_cmap(cmap)(normalized_score)

        # Draw line connecting macrophage and tumor
        ax.plot(
            [row["x_mac"] - x_min, row["x_tumor"] - x_min],
            [row["y_mac"] - y_min, row["y_tumor"] - y_min],
            color=color,
            linewidth=line_width,
            alpha=line_alpha,
        )

    # Plot 3: Empty panel for legend only
    ax = axes[2]
    sm = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=plt.Normalize(vmin=vmin, vmax=vmax),
    )
    sm.set_array([])
    cbar_ax = ax.inset_axes([0.3, 0.1, 0.05, 0.2])  # [x, y, width, height]
    fig.colorbar(sm, cax=cbar_ax, label="Mean expression")

    plt.tight_layout()

    return fig


# %% Batch Processing: Generate All Visualizations ==========

for idx, smp_id_fov in enumerate(smp_id_fov_list, 1):
    print(f"\n{'=' * 60}")
    print(f"Processing sample {idx}/{len(smp_id_fov_list)}: {smp_id_fov}")
    print(f"{'=' * 60}")

    output_dir = output_root / smp_id_fov
    output_dir.mkdir(parents=True, exist_ok=True)

    # Visualization 1: Cell Type Annotation Overview
    print("  Generating cell type annotation overview...")

    seg_sm, rgb_img_bg, x_min, y_min = load_and_crop_segmentation(
        smp_id_fov,
        seg_dir,
        adata,
    )

    rgb_img_lmp1 = create_rgb_annotation(
        seg_sm,
        annotation_dict={
            row["cell_id"]: row["annotation_3"]
            for _, row in adata.obs.loc[
                adata.obs["smp_id_fov"] == smp_id_fov, ["cell_id", "annotation_3"]
            ].iterrows()
        },
        color_dict=COLOR_MAP,
        boundary=True,
        boundary_color_rgb=BOUNDARY_COLOR_RGB,
        background_color_rgb=BACKGROUND_COLOR_RGB,
    )
    # tifffile.imshow(rgb_img_lmp1)

    fig = plot_cell_annotation_overview(smp_id_fov, adata, seg_sm, x_min, y_min)
    output_f = output_dir / "01_cell_annotation.pdf"
    fig.savefig(output_f, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved: {output_f.name}")



    # Visualization 2: Ligand-Receptor Expression with Connections
    print("  Generating ligand-receptor expression with connections plots...")

    for score in lr_pairs:
        print(f"    Processing L-R pair: {score}...")
        fig = plot_ligand_receptor_expression(
            smp_id_fov,
            score,
            df_lr_expr,
            df_lr_mean_tumor_mac,
            seg_sm,
            x_min,
            y_min,
            plot_connections=True,
            line_exclude_quantile=0,
        )
        output_f = output_dir / f"02_ligand_receptor_connection_{score}.pdf"
        fig.savefig(output_f, dpi=300, bbox_inches="tight")
        plt.close(fig)


    print(f"  ✓ Saved {len(lr_pairs)} ligand-receptor expression plots")


    # Visualization 3: Tumor-Macrophage L-R Pair Expression (PSF Method)
    print("  Generating tumor-macrophage PSF expression plots...")

    for score in lr_pairs:
        print(f"    Processing L-R pair: {score}...")

        fig = plot_tumor_macrophage_lr_expression_2(
            smp_id_fov,
            score,
            df_mac_center,
            rgb_img_lmp1,
            x_min,
            y_min,
            sigma=1000,
        )
        output_f = output_dir / f"03_center_lr_mean_{score}.pdf"
        fig.savefig(output_f, dpi=300, bbox_inches="tight")
        plt.close(fig)
    print(f"  ✓ Saved {len(lr_pairs) * 1} tumor-macrophage PSF expression plots")


print(f"\n{'=' * 60}")
print("All visualizations completed!")
print(f"Total samples processed: {len(smp_id_fov_list)}")
print(f"Output directory: {output_root}")
print(f"{'=' * 60}")

# %%
