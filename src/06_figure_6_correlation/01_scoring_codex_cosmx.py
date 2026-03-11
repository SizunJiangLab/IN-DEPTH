# env: rsc

# %% Import libraries ==========
from pathlib import Path

import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc

# Constants
TQDM_FORMAT = "{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
SEPARATOR = "=" * 80

# %% Setup Output Directory ==========

print(f"\n{SEPARATOR}")
print("Setup Output Directory")
print(SEPARATOR)

input_dir = Path(__file__).parent
print(f"\nInput directory: {input_dir}")

output_dir = input_dir / "codex_cosmx"
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Output directory: {output_dir}")

# %% Load Data ==========

print(f"\n{SEPARATOR}")
print("Load Data")
print(SEPARATOR)

print("\nLoading data...")
adata = sc.read_h5ad(input_dir / "codex_cosmx_anndata.h5ad")

print("Data loaded successfully")
print(f"   - Total cells: {adata.n_obs:,}")
print(f"   - Total genes: {adata.n_vars:,}")


# %% Compute AUCell Scores ==========

print(f"\n{SEPARATOR}")
print("Compute AUCell Scores")
print(SEPARATOR)

genesets = {
    # T cell dysfunction signature genes
    "Dysfunction-T": [
        "CTLA4",
        "HAVCR2",
        "LAG3",
        "PDCD1",
        "TIGIT",
        "ENTPD1",
        "BTLA",
        "CD244",
        "VSIR",
        "CD160",
        "NT5E",
        "ADORA2A",
        "PVRIG",
        "SIGLEC7",
        "SIGLEC9",
    ],
    # C1QC+ macrophage signature genes
    "C1QC-Mac": [
        "C1QA",
        "C1QB",
        "C1QC",
        "ITM2B",
        "HLA-DMB",
        "MS4A6A",
        "CTSC",
        "TBXAS1",
        "TMEM176B",
        "SYNGR2",
        "ARHGDIB",
        "TMEM176A",
        "UCP2",
        "CAPZB",
        "MAF",
        "TREM2",
        "MSR1",
    ],
}

# Convert gene sets to dataframe format required by AUCell
print("\nPreparing gene sets for AUCell scoring...")
print(f"   - Number of gene sets: {len(genesets)}")
for name, genes in genesets.items():
    print(f"   - {name}: {len(genes)} genes")

genesets_df = []
for term, genes in genesets.items():
    for gene in genes:
        genesets_df.append({"source": term, "target": gene})
genesets_df = pd.DataFrame(genesets_df)

print(f"\nProcessing {adata.n_obs:,} cells for AUCell scoring...")

# Transfer to GPU for accelerated computation (optional)
# Uncomment the next line if GPU is available and you want faster computation
# rsc.get.anndata_to_GPU(adata)

# Compute AUCell scores for each gene set
print("Computing AUCell scores...")
adata_result = rsc.dcg.aucell(adata, net=genesets_df, verbose=True, empty=True)

# [WARNING] Provided AnnData contains empty observations, returning repaired object
if adata_result is not None:
    adata = adata_result

# Transfer back to CPU
rsc.get.anndata_to_CPU(adata)
print("   Scores stored in adata.obsm['score_aucell']")

print("\nAUCell scores computed.")


# %% Compute EBV Burden ==========

print(f"\n{SEPARATOR}")
print("Compute EBV Burden")
print(SEPARATOR)

print("\nComputing EBV burden...")

# EBV viral genes for burden calculation
ebv_genes = [
    "BCRF1",
    "BLLF1",
    "BNLF2A",
    "BZLF1",
    "EBNA2",
    "EBNA3A",
    "EBNA3BC",
    "LMP1",
    "LMP2",
    "RPMS1",
]

# Sum expression of all available EBV genes
ebv_genes_available = [g for g in ebv_genes if g in adata.var_names]
print(
    f"\nComputing EBV burden: {len(ebv_genes_available)}/{len(ebv_genes)} genes available"
)

if len(ebv_genes_available) > 0:
    adata.obs["ebv_burden"] = adata[:, ebv_genes_available].X.sum(axis=1)
else:
    print("Warning: No EBV genes found in the dataset")
    adata.obs["ebv_burden"] = 0


# %% Combine Scores (single-cell-level) ==========

print(f"\n{SEPARATOR}")
print("Combine Scores (single-cell-level)")
print(SEPARATOR)

print("\nCombining all computed scores...")

# Combine AUCell scores with EBV metrics
adata.obsm["analysis"] = pd.concat(
    [
        adata.obsm["score_aucell"],
        pd.DataFrame({"ebv_burden": adata.obs["ebv_burden"]}),
    ],
    axis=1,
)

print("\nSummary of all computed scores:")
print(adata.obsm["analysis"].describe())


# %% Save Processed Data ==========

print(f"\n{SEPARATOR}")
print("Save Processed Data")
print(SEPARATOR)

output_f_with_scores = output_dir / "01_adata_with_scores.h5ad"
adata.write_h5ad(output_f_with_scores)

print(f"\nData with all scores saved to: {output_f_with_scores}")
print(f"   - Total cells: {adata.n_obs}")
print(f"   - Total genes: {adata.n_vars}")
print(f"   - Analysis metrics: {adata.obsm['analysis'].shape[1]}")


# %% Calculate Scores (region-level) ==========

print(f"\n{SEPARATOR}")
print("Calculate Scores (region-level)")
print(SEPARATOR)


def compute_region_level_scores(adata):
    """
    Compute region-level summary statistics for cell-type-specific scores.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with:
        - obs['annotation']: cell type annotations
        - obs['region_id']: region identifiers
        - obsm['analysis']: computed scores (Dysfunction-T, C1QC-Mac, ebv_burden)

    Returns
    -------
    pd.DataFrame
        Region-level summary statistics with columns:
        - Dysfunction_CD4T: mean Dysfunction-T score for CD4 T cells
        - Dysfunction_CD8T: mean Dysfunction-T score for CD8 T cells
        - Dysfunction_T: mean Dysfunction-T score for all T cells
        - C1QC-Mac_Macrophage: mean C1QC-Mac score for macrophages
        - ebv_burden_B: mean EBV burden for B cells
    """
    region_results = []

    # CD4 T cells: Dysfunction-T score
    adata_cd4 = adata[adata.obs["annotation"] == "CD4 T"]
    if adata_cd4.n_obs > 0:
        analysis_df = adata_cd4.obsm["analysis"].copy()
        analysis_df["region_id"] = adata_cd4.obs["region_id"].values
        dysfunction_cd4_stats = (
            analysis_df.groupby("region_id")["Dysfunction-T"]
            .mean()
            .reset_index()
            .rename(columns={"Dysfunction-T": "Dysfunction_CD4T"})
            .set_index("region_id")
        )
        region_results.append(dysfunction_cd4_stats)
        print(f"CD4 T cells analyzed: {adata_cd4.n_obs}")

    # CD8 T cells: Dysfunction-T score
    adata_cd8 = adata[adata.obs["annotation"] == "CD8 T"]
    if adata_cd8.n_obs > 0:
        analysis_df = adata_cd8.obsm["analysis"].copy()
        analysis_df["region_id"] = adata_cd8.obs["region_id"].values
        dysfunction_cd8_stats = (
            analysis_df.groupby("region_id")["Dysfunction-T"]
            .mean()
            .reset_index()
            .rename(columns={"Dysfunction-T": "Dysfunction_CD8T"})
            .set_index("region_id")
        )
        region_results.append(dysfunction_cd8_stats)
        print(f"CD8 T cells analyzed: {adata_cd8.n_obs}")

    # All T cells: Dysfunction-T score
    adata_t = adata[adata.obs["annotation"].isin(["CD4 T", "CD8 T"])]
    if adata_t.n_obs > 0:
        analysis_df = adata_t.obsm["analysis"].copy()
        analysis_df["region_id"] = adata_t.obs["region_id"].values
        dysfunction_t_stats = (
            analysis_df.groupby("region_id")["Dysfunction-T"]
            .mean()
            .reset_index()
            .rename(columns={"Dysfunction-T": "Dysfunction_T"})
            .set_index("region_id")
        )
        region_results.append(dysfunction_t_stats)
        print(f"Total T cells analyzed: {adata_t.n_obs}")

    # Macrophages: C1QC-Mac score
    adata_m = adata[adata.obs["annotation"] == "Macrophage"]
    if adata_m.n_obs > 0:
        analysis_df = adata_m.obsm["analysis"].copy()
        analysis_df["region_id"] = adata_m.obs["region_id"].values
        c1qc_mac_stats = (
            analysis_df.groupby("region_id")["C1QC-Mac"]
            .mean()
            .reset_index()
            .set_index("region_id")
            .rename(columns={"C1QC-Mac": "C1QC-Mac_Macrophage"})
        )
        region_results.append(c1qc_mac_stats)
        print(f"Macrophages analyzed: {adata_m.n_obs}")

    # B cells: EBV burden
    adata_b = adata[adata.obs["annotation"] == "B cell"]
    if adata_b.n_obs > 0:
        analysis_df = adata_b.obsm["analysis"].copy()
        analysis_df["region_id"] = adata_b.obs["region_id"].values
        ebv_burden_b_stats = (
            analysis_df.groupby("region_id")["ebv_burden"]
            .mean()
            .reset_index()
            .set_index("region_id")
            .rename(columns={"ebv_burden": "ebv_burden_B"})
        )
        region_results.append(ebv_burden_b_stats)
        print(f"B cells analyzed: {adata_b.n_obs}")

    # Combine all computed statistics
    df_result = pd.concat(region_results, axis=1)

    print("\nRegion-level summary statistics computed:")
    print(f"  - Number of regions: {len(df_result)}")
    print(f"  - Number of metrics: {len(df_result.columns)}")
    print("\nMetrics computed:")
    for col in df_result.columns:
        print(f"  - {col}")

    return df_result


# Core-level scores
print("\n--- Core-level scores ---")

adata.obs["region_id"] = (
    adata.obs["data_tag"].astype(str) + "_" + adata.obs["core_id"].astype(str)
)

metadata_df = adata.obs[["region_id", "ebv"]].drop_duplicates()

print(f"\nTotal cells: {adata.n_obs:,}")
print(f"Total cores: {adata.obs['region_id'].nunique()}")
print(
    metadata_df.groupby("ebv", observed=True).size().reset_index(name="count"),
    "\n",
)

core_scores = compute_region_level_scores(adata).merge(metadata_df, on="region_id")

print("\nFinal table preview:")
print(core_scores.head())

print("\nCore-level scores computed.")

# FOV-level scores
print("\n--- FOV-level scores ---")

adata.obs["region_id"] = (
    adata.obs["data_tag"].astype(str) + "_" + adata.obs["fov_id"].astype(str)
)

metadata_df = adata.obs[["region_id", "ebv"]].drop_duplicates()

print(f"\nTotal cells: {adata.n_obs:,}")
print(f"Total FOVs: {adata.obs['region_id'].nunique()}")
print(
    metadata_df.groupby("ebv", observed=True).size().reset_index(name="count"),
    "\n",
)

fov_scores = compute_region_level_scores(adata).merge(metadata_df, on="region_id")

print("\nFinal table preview:")
print(fov_scores.head())

print("\nFOV-level scores computed.")

# %% Save Region-Level Scores ==========

print(f"\n{SEPARATOR}")
print("Save Region-Level Scores")
print(SEPARATOR)

core_scores.to_csv(output_dir / "01_scores_core.csv")
print(f"\nCore-level scores saved to: {output_dir / '01_scores_core.csv'}")
print(f"   - Rows (cores): {len(core_scores)}")
print(f"   - Columns (metrics): {len(core_scores.columns)}")

fov_scores.to_csv(output_dir / "01_scores_fov.csv")
print(f"\nFOV-level scores saved to: {output_dir / '01_scores_fov.csv'}")
print(f"   - Rows (FOVs): {len(fov_scores)}")
print(f"   - Columns (metrics): {len(fov_scores.columns)}")

print(f"\n{SEPARATOR}")
print("ANALYSIS COMPLETE")
print(SEPARATOR)
