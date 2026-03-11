# %% Import libraries ==========
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from tqdm import tqdm

TQDM_FORMAT = "{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"


# %% Load data ==========

output_dir = Path(__file__).parent
adata_ebv = sc.read_h5ad(output_dir / "01_adata_with_spatial.h5ad")

print("Data loaded:")
print(adata_ebv)
print("\nAnnotation distribution:")
print(adata_ebv.obs["annotation_3"].value_counts())

# %% Preprocess data ==========
# Normalize to counts per 10,000 (CPM-like normalization) and log-transform
# Required for ligand-receptor analysis

sc.pp.normalize_total(adata_ebv, target_sum=1e4)
sc.pp.log1p(adata_ebv)

# %% Prepare spatial coordinates with sample offsets ==========
# Create a unified spatial coordinate system for all samples
# This is necessary because each sample has its own coordinate system (0-based)
# We arrange samples in a grid layout to avoid overlapping coordinates

samples = adata_ebv.obs["smp_id"].unique()
print(f"  Number of regions: {len(samples)}")

# Calculate offset for each core (arrange cores in a grid)
max_x = adata_ebv.obs["x_cent"].max()
max_y = adata_ebv.obs["y_cent"].max()
offset_x_step = max_x * 1.2  # 20% padding between cores
offset_y_step = max_y * 1.2

n_cols = int(np.ceil(np.sqrt(len(samples))))  # Square grid layout

# Initialize spatial coordinates array with same shape as obs
# IMPORTANT: This ensures coordinates are in the same order as adata_ebv.obs
adata_ebv.obsm["spatial"] = np.zeros((adata_ebv.n_obs, 2))

for idx, sample in enumerate(samples):
    sample_mask = adata_ebv.obs["smp_id"] == sample

    # Calculate grid position
    row = idx // n_cols
    col = idx % n_cols

    # Calculate offset
    offset_x = col * offset_x_step
    offset_y = row * offset_y_step

    # Assign coordinates directly to correct positions using mask
    # This maintains the original obs order
    adata_ebv.obsm["spatial"][sample_mask, 0] = (
        adata_ebv.obs.loc[sample_mask, "x_cent"].values + offset_x
    )
    adata_ebv.obsm["spatial"][sample_mask, 1] = (
        adata_ebv.obs.loc[sample_mask, "y_cent"].values + offset_y
    )


# %% Define helper functions ==========


def process_ligrec_results(res):
    """
    Process Squidpy ligand-receptor results into tidy format.

    Converts the nested dictionary structure from sq.gr.ligrec() into a
    long-format DataFrame where each row represents one ligand-receptor
    interaction between a source and target cell type.

    Args:
        res: Result dictionary from sq.gr.ligrec() containing 'means' and 'pvalues'

    Returns:
        DataFrame with columns: id, lr_pair, interaction, source, target,
                                ligand, receptor, lr_means, pvalue
    """
    lr_means = res["means"]
    lr_pvals = res["pvalues"]

    # Get interaction names (ligand-receptor pairs)
    interactions = lr_means.columns

    results_list = []

    for source_target in tqdm(
        lr_means.index,
        desc="Processing squidpy results",
        bar_format=TQDM_FORMAT,
    ):
        source, target = source_target

        for interaction in interactions:
            mean_val = lr_means.loc[source_target, interaction]
            pval = lr_pvals.loc[source_target, interaction]

            # Parse ligand and receptor from interaction name (tuple format)
            ligand, receptor = interaction

            results_list.append(
                {
                    "id": f"{source}::{target}_{ligand}::{receptor}",
                    "lr_pair": f"{source}::{target}",
                    "interaction": f"{ligand}::{receptor}",
                    "source": source,
                    "target": target,
                    "ligand": ligand,
                    "receptor": receptor,
                    "lr_means": mean_val,
                    "pvalue": pval,
                }
            )

    results_df = pd.DataFrame(results_list)

    return results_df


# %% LMP1- analysis ==========
# Analyze ligand-receptor interactions between LMP1- tumor cells and macrophages
# Uses permutation testing (n=1000) to assess statistical significance

adata_n = adata_ebv[
    adata_ebv.obs.annotation_3.isin(["lmp1n_tumor_hopped", "lmp1n_macrophage"])
]

res_n = sq.gr.ligrec(
    adata_n,
    n_perms=1000,
    cluster_key="annotation_3",
    copy=True,
    use_raw=False,
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"},
    seed=42,
)


# %% LMP1+ analysis ==========
# Analyze ligand-receptor interactions between LMP1+ tumor cells and macrophages
# Same parameters as LMP1- analysis for direct comparison

adata_p = adata_ebv[
    adata_ebv.obs.annotation_3.isin(["lmp1p_tumor_hopped", "lmp1p_macrophage"])
]

res_p = sq.gr.ligrec(
    adata_p,
    n_perms=1000,
    cluster_key="annotation_3",
    copy=True,
    use_raw=False,
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"},
    seed=42,
)

# %% Process and save results ==========
res_df_n = process_ligrec_results(res_n)
res_df_p = process_ligrec_results(res_p)

# Clean up cell type names in ligand/receptor columns
# Remove prefix (lmp1n_/lmp1p_) and suffix (_hopped) for cleaner comparison
for col in ["ligand", "receptor"]:
    res_df_p[col] = res_df_p[col].str.replace("lmp1p_", "")
    res_df_p[col] = res_df_p[col].str.replace("_hopped", "")

    res_df_n[col] = res_df_n[col].str.replace("lmp1n_", "")
    res_df_n[col] = res_df_n[col].str.replace("_hopped", "")

# Recreate interaction and id columns with cleaned names
# This ensures consistency between LMP1+ and LMP1- results for merging
for df in [res_df_n, res_df_p]:
    df["interaction"] = df["ligand"] + "::" + df["receptor"]
    df["id"] = df["lr_pair"] + "_" + df["interaction"]

# Save individual results
res_df_n.to_csv(output_dir / "05_squidpy_res_df_lmp1n.csv", index=False)
res_df_p.to_csv(output_dir / "05_squidpy_res_df_lmp1p.csv", index=False)


# %% Merge LMP1+ and LMP1- results for comparison ==========
# Create wide-format table for direct comparison and differential analysis

# Select only key columns and rename for clarity
res_df_n_sm = res_df_n.rename(
    columns={"lr_means": "lr_means_lmp1n", "pvalue": "pvalue_lmp1n"}
)[["id", "lr_means_lmp1n", "pvalue_lmp1n"]]
res_df_p_sm = res_df_p.rename(
    columns={"lr_means": "lr_means_lmp1p", "pvalue": "pvalue_lmp1p"}
)[["id", "lr_means_lmp1p", "pvalue_lmp1p"]]

# Merge results using outer join to keep all interactions
res_df_merged = res_df_n_sm.merge(res_df_p_sm, on="id", how="outer")

# Re-extract lr_pair and interaction from id for convenience
res_df_merged[["lr_pair", "interaction"]] = res_df_merged["id"].str.split(
    "_", n=1, expand=True
)
res_df_merged[["ligand", "receptor"]] = res_df_merged["interaction"].str.split(
    "::", expand=True
)
res_df_merged[["source", "target"]] = res_df_merged["lr_pair"].str.split(
    "::", expand=True
)

# Calculate log2 fold change: LMP1+ vs LMP1-
# Positive values indicate higher interaction strength in LMP1+
res_df_merged["log2fc"] = np.log2(
    res_df_merged["lr_means_lmp1p"] / res_df_merged["lr_means_lmp1n"]
)

res_df_merged.to_csv(output_dir / "05_squidpy_res_df_merged.csv", index=False)
