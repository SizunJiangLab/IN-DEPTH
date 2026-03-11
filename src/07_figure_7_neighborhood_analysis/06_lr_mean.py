"""
Ligand-Receptor (LR) Mean Expression Analysis

This script calculates mean ligand-receptor expression values for:
1. Tumor-Macrophage cell pairs at single-cell resolution
2. Macrophage-centered aggregations grouped by LMP1 status

Input:
    - 01_adata_with_spatial.h5ad: Spatial transcriptomics data
    - 01_macrophage_tumor_neighbors.csv: Neighbor cell pairs
    - 01_macrophage_categories.csv: Macrophage annotations
    - 05_squidpy_res_df_merged.csv: Squidpy LR interaction results

Output:
    - 06_tumor_macrophage_lr_mean.feather: Single-cell pair LR means (arithmetic mean)
    - 06_tumor_macrophage_lr_mean_pair.feather: Single-cell pair LR means (zero if ligand OR receptor is zero)
    - 06_tumor_macrophage_lr_expr.feather: Single-cell pair ligand and receptor expression values
    - 06_macrophage_center_lr_mean.feather: Macrophage-centered aggregated LR means
    - 06_macrophage_center_lr_mean_pair.feather: Macrophage-centered aggregated LR means (pair version)
"""

# %% Import libraries ==========
from pathlib import Path

import pandas as pd
import scanpy as sc
from tqdm import tqdm

TQDM_FORMAT = "{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"


# %% 01. Data Loading and Preprocessing ==========

print("=" * 80)
print("STEP 1: Loading and Preprocessing Data")
print("=" * 80)

data_dir = Path(__file__).parent

# Load spatial transcriptomics data
print(f"\nLoading AnnData from: {data_dir / '01_adata_with_spatial.h5ad'}")
adata = sc.read_h5ad(data_dir / "01_adata_with_spatial.h5ad")
print(f"  Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Normalize and log-transform expression data
print("\nNormalizing expression data...")
sc.pp.normalize_total(adata, target_sum=1e4)  # CPM normalization (counts per 10k)
sc.pp.log1p(adata)  # Log(x + 1) transformation
print("  ✓ Normalized to 10,000 counts per cell and log-transformed")

# Load neighbor relationships and annotations
print(
    f"\nLoading neighbor relationships from: {data_dir / '01_macrophage_tumor_neighbors.csv'}"
)
df_neighbors = pd.read_csv(data_dir / "01_macrophage_tumor_neighbors.csv")
print(f"  Loaded {len(df_neighbors)} tumor-macrophage neighbor pairs")

print(
    f"\nLoading macrophage categories from: {data_dir / '01_macrophage_categories.csv'}"
)
df_neighbors_mac = pd.read_csv(data_dir / "01_macrophage_categories.csv")
print(f"  Loaded {len(df_neighbors_mac)} macrophage annotations")

# Merge annotations and rename columns for clarity
df_neighbors = df_neighbors.merge(
    df_neighbors_mac[["cell_id_mac", "annotation_3"]], on="cell_id_mac"
)
df_neighbors.rename(
    columns={
        "annotation_2": "annotation_tumor",
        "annotation_3": "annotation_mac",
    },
    inplace=True,
)

print("\nTumor cell type distribution:")
print(df_neighbors.annotation_tumor.value_counts())
print("\nMacrophage cell type distribution:")
print(df_neighbors.annotation_mac.value_counts())

# Add spatial coordinates for tumor and macrophage cells
print("\nAdding spatial coordinates...")
df_neighbors = (
    df_neighbors.merge(
        adata.obs[["x_cent", "y_cent"]], left_on="cell_id_tumor", right_index=True
    )
    .rename(columns={"x_cent": "x_tumor", "y_cent": "y_tumor"})
    .merge(adata.obs[["x_cent", "y_cent"]], left_on="cell_id_mac", right_index=True)
    .rename(columns={"x_cent": "x_mac", "y_cent": "y_mac"})
    .assign(
        x_lr=lambda df: (df["x_tumor"] + df["x_mac"]) / 2,  # Midpoint X coordinate
        y_lr=lambda df: (df["y_tumor"] + df["y_mac"]) / 2,  # Midpoint Y coordinate
        pair_id=lambda df: (
            df["cell_id_tumor"] + "::" + df["cell_id_mac"]
        ),  # Unique pair identifier
    )
    .set_index("pair_id")
)
print(f"  ✓ Added spatial coordinates for {len(df_neighbors)} cell pairs")

# Load Squidpy ligand-receptor interaction results
print(f"\nLoading Squidpy LR results from: {data_dir / '05_squidpy_res_df_merged.csv'}")
df_squidpy = pd.read_csv(data_dir / "05_squidpy_res_df_merged.csv")
print(f"  Loaded {len(df_squidpy)} LR interaction records")

# Filter for tumor-to-macrophage interactions
df_squidpy_t2m = df_squidpy[df_squidpy["interaction"] == "tumor::macrophage"].copy()
print(f"  Filtered to {len(df_squidpy_t2m)} tumor-to-macrophage interactions")


# %% 02. Tumor-Macrophage LR Mean Calculation ==========

print("\n" + "=" * 80)
print("STEP 2: Calculating Ligand-Receptor Mean Expression (Single-Cell Pairs)")
print("=" * 80)

lr_pairs = df_squidpy_t2m["lr_pair"].unique()
print(f"\nProcessing {len(lr_pairs)} unique LR pairs...")

lr_mean_list = []
lr_mean_pair_list = []
lr_expr_list = []

for lr_pair in tqdm(lr_pairs, desc="Processing LR pairs", bar_format=TQDM_FORMAT):
    # Parse ligand and receptor gene names
    source, target = lr_pair.split(
        "::"
    )  # source=ligand (tumor), target=receptor (macrophage)

    # Extract expression values for the ligand (in tumor cells) and receptor (in macrophage cells)
    expr_source = adata[df_neighbors.cell_id_tumor, source].X.toarray()
    expr_target = adata[df_neighbors.cell_id_mac, target].X.toarray()

    # Calculate lr_mean: arithmetic mean of ligand and receptor expression
    lr_mean = (expr_source + expr_target) / 2

    # Calculate lr_mean_pair: same as lr_mean, but set to 0 if EITHER ligand OR receptor is zero
    # This stricter metric requires BOTH genes to be expressed for a non-zero value
    lr_mean_pair = lr_mean.copy()
    lr_mean_pair[(expr_source == 0) | (expr_target == 0)] = 0

    # Store results as DataFrames with pair_id index
    df_lr_mean_pair = pd.DataFrame(
        lr_mean_pair, index=df_neighbors.index, columns=[lr_pair]
    )
    df_lr_mean = pd.DataFrame(lr_mean, index=df_neighbors.index, columns=[lr_pair])
    df_lr_expr = pd.DataFrame(
        {
            f"{source}": expr_source.flatten(),
            f"{target}": expr_target.flatten(),
        },
        index=df_neighbors.index,
    )

    lr_mean_list.append(df_lr_mean)
    lr_mean_pair_list.append(df_lr_mean_pair)
    lr_expr_list.append(df_lr_expr)

# Concatenate all LR pairs into wide-format DataFrames
df_lr_mean = pd.concat(lr_mean_list, axis=1)
df_lr_mean_pair = pd.concat(lr_mean_pair_list, axis=1)
df_lr_expr = pd.concat(lr_expr_list, axis=1)

# Remove duplicate columns if any
df_lr_expr = df_lr_expr.loc[:, ~df_lr_expr.columns.duplicated()]

print(
    f"\n✓ Calculated LR means for {len(lr_pairs)} LR pairs across {len(df_neighbors)} cell pairs"
)
print("  - df_lr_mean: arithmetic mean (ligand + receptor) / 2")
print("  - df_lr_mean_pair: same as lr_mean but zero if ligand OR receptor is zero")
print("  - df_lr_expr: individual ligand and receptor expression values")
print(f"  - Shape (lr_mean/pair): {df_lr_mean.shape}")
print(f"  - Shape (lr_expr): {df_lr_expr.shape}")


# %% 03. Annotate and Save Single-Cell Pair Results ==========

print("\n" + "=" * 80)
print("STEP 3: Annotating LMP1 Status and Saving Single-Cell Results")
print("=" * 80)


def annotate_lmp1(df_neighbors_lr, value_columns):
    """
    Annotate cell pairs with LMP1 status.

    Assigns "lmp1p" if both tumor and macrophage cells are LMP1-positive,
    "lmp1n" if both are LMP1-negative, otherwise None.

    Parameters
    ----------
    df_neighbors_lr : pd.DataFrame
        DataFrame containing neighbor pairs with annotations and LR expression values
    value_columns : list
        List of column names containing expression values (LR pairs or individual genes)

    Returns
    -------
    pd.DataFrame
        Annotated DataFrame with selected columns and lmp1 status
    """
    # Initialize LMP1 status column
    df_neighbors_lr["lmp1"] = None

    # Annotate LMP1-positive pairs (both tumor and macrophage are LMP1+)
    mask_lmp1p = (df_neighbors_lr.annotation_mac == "lmp1p_macrophage") & (
        df_neighbors_lr.annotation_tumor == "lmp1p_tumor"
    )
    df_neighbors_lr.loc[mask_lmp1p, "lmp1"] = "lmp1p"

    # Annotate LMP1-negative pairs (both tumor and macrophage are LMP1-)
    mask_lmp1n = (df_neighbors_lr.annotation_mac == "lmp1n_macrophage") & (
        df_neighbors_lr.annotation_tumor == "lmp1n_tumor"
    )
    df_neighbors_lr.loc[mask_lmp1n, "lmp1"] = "lmp1n"

    # Select and reorder columns: metadata first, then LR expression values
    meta_col = [
        "smp_id",
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
    df_neighbors_lr = df_neighbors_lr[meta_col + value_columns]

    return df_neighbors_lr


# Join neighbor metadata with LR expression values
print("\nAnnotating LMP1 status for single-cell pairs...")
df_neighbors_lr = df_neighbors.join(df_lr_mean)
df_neighbors_lr = annotate_lmp1(df_neighbors_lr, lr_pairs.tolist())

df_neighbors_lr_pair = df_neighbors.join(df_lr_mean_pair)
df_neighbors_lr_pair = annotate_lmp1(df_neighbors_lr_pair, lr_pairs.tolist())

df_neighbors_lr_expr = df_neighbors.join(df_lr_expr)
df_neighbors_lr_expr = annotate_lmp1(df_neighbors_lr_expr, df_lr_expr.columns.tolist())

print("\nLMP1 status distribution (lr_mean):")
print(df_neighbors_lr.lmp1.value_counts(dropna=False))

# Save single-cell pair results
print("\nSaving single-cell pair results...")
output_file_1 = data_dir / "06_tumor_macrophage_lr_mean.feather"
df_neighbors_lr.to_feather(output_file_1)
print(f"  ✓ Saved: {output_file_1}")
print(f"    Shape: {df_neighbors_lr.shape}")

output_file_2 = data_dir / "06_tumor_macrophage_lr_mean_pair.feather"
df_neighbors_lr_pair.to_feather(output_file_2)
print(f"  ✓ Saved: {output_file_2}")
print(f"    Shape: {df_neighbors_lr_pair.shape}")

output_file_3 = data_dir / "06_tumor_macrophage_lr_expr.feather"
df_neighbors_lr_expr.to_feather(output_file_3)
print(f"  ✓ Saved: {output_file_3}")
print(f"    Shape: {df_neighbors_lr_expr.shape}")

# %% 04. Macrophage-Centered LR Mean Aggregation ==========

print("\n" + "=" * 80)
print("STEP 4: Aggregating LR Means per Macrophage Cell")
print("=" * 80)

# Define metadata columns for grouping by macrophage cell
meta_col = ["smp_id", "cell_id_mac", "x_mac", "y_mac", "lmp1"]

# Aggregate lr_mean: average LR expression across all tumor neighbors of each macrophage
print("\nAggregating lr_mean by macrophage cell...")
df_neighbors_lr_lmp1p = df_neighbors_lr[df_neighbors_lr.lmp1 == "lmp1p"]
df_lmp1p_mac_center = (
    df_neighbors_lr_lmp1p.groupby(meta_col)[lr_pairs].mean().reset_index()
)
print(f"  LMP1+ macrophages: {len(df_lmp1p_mac_center)}")

df_neighbors_lr_lmp1n = df_neighbors_lr[df_neighbors_lr.lmp1 == "lmp1n"]
df_lmp1n_mac_center = (
    df_neighbors_lr_lmp1n.groupby(meta_col)[lr_pairs].mean().reset_index()
)
print(f"  LMP1- macrophages: {len(df_lmp1n_mac_center)}")

df_mac_center_lr = pd.concat([df_lmp1p_mac_center, df_lmp1n_mac_center], axis=0)
print(f"  Total aggregated macrophages (lr_mean): {len(df_mac_center_lr)}")

# Aggregate lr_mean_pair: same aggregation for the stricter version (requires both genes expressed)
print("\nAggregating lr_mean_pair by macrophage cell...")
df_neighbors_lr_pair_lmp1p = df_neighbors_lr_pair[df_neighbors_lr_pair.lmp1 == "lmp1p"]
df_lmp1p_mac_center = (
    df_neighbors_lr_pair_lmp1p.groupby(meta_col)[lr_pairs].mean().reset_index()
)

df_neighbors_lr_pair_lmp1n = df_neighbors_lr_pair[df_neighbors_lr_pair.lmp1 == "lmp1n"]
df_lmp1n_mac_center = (
    df_neighbors_lr_pair_lmp1n.groupby(meta_col)[lr_pairs].mean().reset_index()
)
df_mac_center_lr_pair = pd.concat([df_lmp1p_mac_center, df_lmp1n_mac_center], axis=0)
print(f"  Total aggregated macrophages (lr_mean_pair): {len(df_mac_center_lr_pair)}")

# Save macrophage-centered results
print("\nSaving macrophage-centered results...")
output_file_3 = data_dir / "06_macrophage_center_lr_mean.feather"
df_mac_center_lr.to_feather(output_file_3)
print(f"  ✓ Saved: {output_file_3}")
print(f"    Shape: {df_mac_center_lr.shape}")

output_file_4 = data_dir / "06_macrophage_center_lr_mean_pair.feather"
df_mac_center_lr_pair.to_feather(output_file_4)
print(f"  ✓ Saved: {output_file_4}")
print(f"    Shape: {df_mac_center_lr_pair.shape}")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
