# %% Import libraries ==========
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from matplotlib.patches import Circle
from scipy.spatial import cKDTree
from tqdm import tqdm

# %% Configuration ==========

adata_f = (
    Path(__file__).parent.parent
    / "06_figure_6_correlation"
    / "codex_cosmx_anndata.h5ad"
)
print(adata_f)

anno_lmp1_f = Path(__file__).parent / "codex_lmp1_status.csv"

output_dir = Path(__file__).parent

RADIUS = 40  # px, pixel size = 0.5 um / px

# %% Load LMP1 annotation ==========

anno_lmp1 = pd.read_csv(anno_lmp1_f)

anno_lmp1["id"] = (
    anno_lmp1["TMA"].astype(str)
    + "_"
    + anno_lmp1["cropName"].astype(str)
    + "_"
    + anno_lmp1["cellLabel"].astype(str)
)

anno_lmp1.rename(columns={"LMP1_status_sC": "lmp1"}, inplace=True)
anno_lmp1["lmp1"] = (anno_lmp1["lmp1"] == "LMP1+").astype(int)
anno_lmp1.set_index("id", inplace=True)

print(anno_lmp1["lmp1"].value_counts())

# %% Load adata and label tumor ==========

adata = sc.read_h5ad(adata_f)
adata_ebv = adata[adata.obs.ebv == "EBV+"]

adata_ebv.obs = adata_ebv.obs.join(anno_lmp1["lmp1"], how="left")

adata_b = adata_ebv[
    adata_ebv.obs.annotation.isin(["B cell"]) & adata_ebv.obs.lmp1.notna(),
].copy()
adata_m = adata_ebv[adata_ebv.obs.annotation.isin(["Macrophage"])].copy()

adata_b.obs["annotation_2"] = None
adata_b.obs.loc[adata_b.obs.lmp1 == 0, "annotation_2"] = "lmp1n_tumor"
adata_b.obs.loc[adata_b.obs.lmp1 == 1, "annotation_2"] = "lmp1p_tumor"

adata_m.obs["annotation_2"] = "macrophage"

adata_ebv = sc.concat([adata_b, adata_m])

print(adata_ebv.obs.annotation_2.value_counts())

# %% Spatial analysis ==========

# Create sample ID (smp_id = data_tag + core_id) for sample-wise spatial analysis
adata_ebv.obs["smp_id"] = (
    adata_ebv.obs["data_tag"].astype(str) + "_" + adata_ebv.obs["core_id"].astype(str)
)
smp_ids = adata_ebv.obs.smp_id.unique()

# Store all neighbor pairs
neighbor_records = []  # List of dicts: {smp_id, cell_id_mac, cell_id_tumor}

for smp_id in tqdm(smp_ids, desc="Processing samples"):
    # Get current sample data
    adata_smp = adata_ebv[adata_ebv.obs.smp_id == smp_id]

    # Separate macrophages and tumor cells
    mac_mask = adata_smp.obs.annotation_2 == "macrophage"
    tumor_mask = adata_smp.obs.annotation_2.isin(["lmp1p_tumor", "lmp1n_tumor"])

    adata_mac = adata_smp[mac_mask]
    adata_tumor = adata_smp[tumor_mask]

    # Skip if no macrophages or tumor cells in this sample
    if len(adata_mac) == 0 or len(adata_tumor) == 0:
        continue

    # Get spatial coordinates (x_cent, y_cent)
    mac_coords = adata_mac.obs[["x_cent", "y_cent"]].values
    tumor_coords = adata_tumor.obs[["x_cent", "y_cent"]].values

    # Build KDTree for tumor cells for efficient spatial queries
    tree = cKDTree(tumor_coords)

    # Query neighboring tumor cells for each macrophage
    for mac_idx, mac_id in enumerate(adata_mac.obs_names):
        mac_coord = mac_coords[mac_idx]

        # Find tumor cells within RADIUS
        indices = tree.query_ball_point(mac_coord, r=RADIUS)

        # Skip macrophages with no neighboring tumor cells
        # Note: These macrophages will retain their original annotation_2 = "macrophage"
        if len(indices) == 0:
            continue

        # Record all macrophage-tumor pairs
        for tumor_idx in indices:
            tumor_id = adata_tumor.obs_names[tumor_idx]
            neighbor_records.append(
                {
                    "smp_id": smp_id,
                    "cell_id_mac": mac_id,
                    "cell_id_tumor": tumor_id,
                }
            )


# Create neighbor DataFrame from collected records
df_neighbors = pd.DataFrame(neighbor_records)
# Merge tumor cell annotations into neighbor pairs
df_neighbors = df_neighbors.merge(
    adata_ebv.obs[["annotation_2"]], left_on="cell_id_tumor", right_index=True
)

# Count neighboring tumor cells by type for each macrophage
df_neighbors_mac = (
    df_neighbors.groupby(["cell_id_mac", "annotation_2"])
    .size()
    .reset_index(name="n_neighbors")
    .pivot(index="cell_id_mac", columns="annotation_2", values="n_neighbors")
    .fillna(0)
)

# Initialize categorization columns
df_neighbors_mac["tag"] = None
df_neighbors_mac["annotation_3"] = None

# Categorize macrophages based on neighboring tumor types
mask_mixed = (df_neighbors_mac["lmp1p_tumor"] > 0) & (
    df_neighbors_mac["lmp1n_tumor"] > 0
)
mask_lmp1p = (df_neighbors_mac["lmp1p_tumor"] > 0) & (
    df_neighbors_mac["lmp1n_tumor"] == 0
)
mask_lmp1n = (df_neighbors_mac["lmp1p_tumor"] == 0) & (
    df_neighbors_mac["lmp1n_tumor"] > 0
)

# Assign detailed tags
df_neighbors_mac.loc[mask_mixed, "tag"] = "mixed"
df_neighbors_mac.loc[mask_lmp1p, "tag"] = "only_lmp1p"
df_neighbors_mac.loc[mask_lmp1n, "tag"] = "only_lmp1n"

# Assign annotation_3 labels
# Note: Mixed macrophages are labeled as lmp1p_macrophage (LMP1+ dominant)
df_neighbors_mac.loc[mask_mixed, "annotation_3"] = "lmp1p_macrophage"
df_neighbors_mac.loc[mask_lmp1p, "annotation_3"] = "lmp1p_macrophage"
df_neighbors_mac.loc[mask_lmp1n, "annotation_3"] = "lmp1n_macrophage"

print(df_neighbors_mac.tag.value_counts())
print()
print(df_neighbors_mac.annotation_3.value_counts())

# %% Assign annotation_3 labels ==========
# Assign annotation_3 to adata_ebv.obs
# Initialize with annotation_2
adata_ebv.obs["annotation_3"] = adata_ebv.obs["annotation_2"]

# Update macrophages with neighborhood-based labels
mask_mac = adata_ebv.obs.index.isin(df_neighbors_mac.index)
adata_ebv.obs.loc[mask_mac, "annotation_3"] = df_neighbors_mac.loc[
    adata_ebv.obs.index[mask_mac], "annotation_3"
]

# Mark tumor cells that are neighbors of macrophages with "_hopped" suffix
mask_tumor = adata_ebv.obs.index.isin(df_neighbors.cell_id_tumor)
adata_ebv.obs.loc[mask_tumor, "annotation_3"] = (
    adata_ebv.obs.loc[mask_tumor, "annotation_2"] + "_hopped"
)

print(adata_ebv.obs.annotation_3.value_counts())

# %% Visualization: Spatial distribution of cells by annotation_3 ==========
smp_id = "DFCI_c12"

plot_circle = True  # Whether to draw circles around macrophages

df_plot = adata_ebv.obs[adata_ebv.obs.smp_id == smp_id].copy()

# Define fixed color mapping for consistency across plots
color_map = {
    "lmp1n_tumor": "#1f77b4",  # Blue
    "lmp1p_tumor": "#ff7f0e",  # Orange
    "lmp1n_tumor_hopped": "#2ca02c",  # Green
    "lmp1p_tumor_hopped": "#d62728",  # Red
    "lmp1n_macrophage": "#9467bd",  # Purple
    "lmp1p_macrophage": "#8c564b",  # Brown
    "macrophage": "#e377c2",  # Pink
}

fig, axes = plt.subplots(2, 2, figsize=(12, 12))
axes = axes.flatten()

ax = axes[0]
df_plot_sm = df_plot[
    df_plot.annotation_3.isin(["lmp1n_tumor_hopped", "lmp1n_macrophage"])
]
for annotation, group in df_plot_sm.groupby("annotation_3"):
    color = color_map.get(annotation, "#7f7f7f")  # Default gray
    ax.scatter(group["x_cent"], group["y_cent"], label=annotation, alpha=0.6, c=color)
# Draw circles around macrophages
if plot_circle:
    df_mac = df_plot_sm[df_plot_sm.annotation_3 == "lmp1n_macrophage"]
    for _, row in df_mac.iterrows():
        circle = Circle(
            (row["x_cent"], row["y_cent"]),
            RADIUS,
            fill=False,
            edgecolor="#9467bd",
            linewidth=0.5,
            alpha=0.3,
        )
        ax.add_patch(circle)
ax.set_xlabel("X coordinate")
ax.set_ylabel("Y coordinate")
ax.set_title(f"LMP1- cells in sample {smp_id}")
ax.legend()
ax.set_aspect("equal")

ax = axes[1]
df_plot_sm = df_plot[
    df_plot.annotation_3.isin(["lmp1p_tumor_hopped", "lmp1p_macrophage"])
]
for annotation, group in df_plot_sm.groupby("annotation_3"):
    color = color_map.get(annotation, "#7f7f7f")
    ax.scatter(group["x_cent"], group["y_cent"], label=annotation, alpha=0.6, c=color)
# Draw circles around macrophages
if plot_circle:
    df_mac = df_plot_sm[df_plot_sm.annotation_3 == "lmp1p_macrophage"]
    for _, row in df_mac.iterrows():
        circle = Circle(
            (row["x_cent"], row["y_cent"]),
            RADIUS,
            fill=False,
            edgecolor="#8c564b",
            linewidth=0.5,
            alpha=0.3,
        )
        ax.add_patch(circle)
ax.set_xlabel("X coordinate")
ax.set_ylabel("Y coordinate")
ax.set_title(f"LMP1+ cells in sample {smp_id}")
ax.legend()
ax.set_aspect("equal")

ax = axes[2]
df_plot_sm = df_plot[
    ~df_plot.annotation_3.isin(
        [
            "lmp1p_tumor_hopped",
            "lmp1p_macrophage",
            "lmp1n_tumor_hopped",
            "lmp1n_macrophage",
        ]
    )
]
for annotation, group in df_plot_sm.groupby("annotation_3"):
    color = color_map.get(annotation, "#7f7f7f")
    ax.scatter(group["x_cent"], group["y_cent"], label=annotation, alpha=0.6, c=color)
# Draw circles around macrophages
if plot_circle:
    df_mac = df_plot_sm[df_plot_sm.annotation_3 == "macrophage"]
    for _, row in df_mac.iterrows():
        circle = Circle(
            (row["x_cent"], row["y_cent"]),
            RADIUS,
            fill=False,
            edgecolor="#e377c2",
            linewidth=0.5,
            alpha=0.3,
        )
        ax.add_patch(circle)
ax.set_xlabel("X coordinate")
ax.set_ylabel("Y coordinate")
ax.set_title(f"Other cells in sample {smp_id}")
ax.legend()
ax.set_aspect("equal")

ax = axes[3]
df_plot_sm = df_plot
for annotation, group in df_plot_sm.groupby("annotation_3"):
    color = color_map.get(annotation, "#7f7f7f")
    ax.scatter(group["x_cent"], group["y_cent"], label=annotation, alpha=0.6, c=color)
# Draw circles around all macrophages
if plot_circle:
    df_mac_all = df_plot[
        df_plot.annotation_3.isin(
            ["lmp1n_macrophage", "lmp1p_macrophage", "macrophage"]
        )
    ]
    for _, row in df_mac_all.iterrows():
        mac_type = row["annotation_3"]
        circle_color = color_map.get(mac_type, "#7f7f7f")
        circle = Circle(
            (row["x_cent"], row["y_cent"]),
            RADIUS,
            fill=False,
            edgecolor=circle_color,
            linewidth=0.5,
            alpha=0.3,
        )
        ax.add_patch(circle)
ax.set_xlabel("X coordinate")
ax.set_ylabel("Y coordinate")
ax.set_title(f"All cells in sample {smp_id}")
ax.legend()
ax.set_aspect("equal")


# %% Save results ==========
output_dir.mkdir(exist_ok=True)

adata_ebv.write_h5ad(output_dir / "01_adata_with_spatial.h5ad")
df_neighbors.to_csv(output_dir / "01_macrophage_tumor_neighbors.csv", index=False)
df_neighbors_mac.to_csv(output_dir / "01_macrophage_categories.csv")

print(f"\nSaved results to: {output_dir}")
print("  - 01_adata_with_spatial.h5ad: AnnData object with annotation_3")
print("  - 01_macrophage_tumor_neighbors.csv: Macrophage-tumor neighbor pairs")
print("  - 01_macrophage_categories.csv: Macrophage categories and counts")
