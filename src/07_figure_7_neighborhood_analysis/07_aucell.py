# %% rsc
from pathlib import Path

import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc

# %% Load data ==========

data_dir = Path(__file__).parent

adata = sc.read_h5ad(data_dir / "01_adata_with_spatial.h5ad")

print("Data loaded:")
print(adata)
print("\nAnnotation distribution:")
print(adata.obs["annotation_2"].value_counts())
print("\nTumor environment distribution:")
print(adata.obs["annotation_3"].value_counts())


# %% Read GMT files and prepare gene sets ==========


def read_gmt(gmt_file):
    """Read GMT file and return a dictionary of gene sets."""
    gene_sets = {}
    with open(gmt_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                term_name = parts[0]
                genes = parts[2:]  # Skip the description field
                gene_sets[term_name] = genes
    return gene_sets


# Read GMT files
gmt_files = list(data_dir.glob("./03_gsea_mac_*/gene_sets.gmt"))

all_gene_sets = {}
for gmt_file in gmt_files:
    gene_sets = read_gmt(gmt_file)
    all_gene_sets.update(gene_sets)

print(f"Loaded {len(all_gene_sets)} gene sets from GMT files")

# Select specific gene sets
terms_mac_up = [
    "Ferroptosis WP4313",
    "Cellular Response To Starvation R-HSA-9711097",
    "Terminal Pathway Of Complement R-HSA-166665",
    "ROS And RNS Production In Phagocytes R-HSA-1222556",
    "Negative Regulation Of Leukocyte Activation (GO:0002695)",
    "Negative Regulation Of T Cell Proliferation (GO:0042130)",
]
terms_mac_down = [
    "Negative Regulation Of Phosphatidylinositol 3-Kinase Signaling (GO:0014067)",
    "Negative Regulation Of Focal Adhesion Assembly (GO:0051895)",
    "Inositol Phosphate Metabolism R-HSA-1483249",
    "CD22 Mediated BCR Regulation R-HSA-5690714",
    "Scavenging Of Heme From Plasma R-HSA-2168880",
    "Electron Transport Chain OXPHOS System In Mitochondria WP111",
]
terms_tumor_up = [
    "STAT5 Activation Downstream Of FLT3 ITD Mutants R-HSA-9702518",
    "Positive Regulation Of Telomerase Activity (GO:0051973)",
    "Tumor Necrosis Factor-Mediated Signaling Pathway (GO:0033209)",
    "Interleukin-27 Signaling R-HSA-9020956",
    "Ebstein Barr Virus LMP1 Signaling WP262",
    "NF-kappa B signaling pathway",
]
terms_tumor_down = [
    "Basement Membrane Organization (GO:0071711)",
    "Negative Regulation Of Calcium Ion Transmembrane Transport (GO:1903170)",
    "Glycosaminoglycan Metabolic Process (GO:0030203)",
    "CD22 Mediated BCR Regulation R-HSA-5690714",
    "Positive Regulation Of Calcium Ion Import Across Plasma Membrane (GO:1905665)",
    "Sulfatase And Aromatase Pathway WP5368",
]

terms_mac = terms_mac_up + terms_mac_down
terms_tumor = terms_tumor_up + terms_tumor_down
terms = terms_mac + terms_tumor

selected_genesets = {}

for term in terms:
    if term in all_gene_sets:
        selected_genesets[term] = all_gene_sets[term]
        print(f"Found: {term} with {len(all_gene_sets[term])} genes")
    else:
        print(f"Warning: {term} not found in GMT files")

print(f"\nTotal gene sets: {len(selected_genesets)}")

# Convert to DataFrame format for rapids_singlecell
print("\nPreparing gene sets for AUCell scoring...")
print(f"   - Number of gene sets: {len(selected_genesets)}")
for name, genes in selected_genesets.items():
    print(f"   - {name}: {len(genes)} genes")

genesets_df = []
for term, genes in selected_genesets.items():
    for gene in genes:
        genesets_df.append({"source": term, "target": gene})
genesets_df = pd.DataFrame(genesets_df)

genesets_df_combined = genesets_df.copy()
mask_mac_up = genesets_df_combined["source"].isin(terms_mac_up)
mask_mac_down = genesets_df_combined["source"].isin(terms_mac_down)
mask_tumor_up = genesets_df_combined["source"].isin(terms_tumor_up)
mask_tumor_down = genesets_df_combined["source"].isin(terms_tumor_down)

genesets_df_combined.loc[mask_mac_up, "source"] = "mac_up"
genesets_df_combined.loc[mask_mac_down, "source"] = "mac_down"
genesets_df_combined.loc[mask_tumor_up, "source"] = "tumor_up"
genesets_df_combined.loc[mask_tumor_down, "source"] = "tumor_down"

print(genesets_df_combined["source"].value_counts())

genesets_df = (
    pd.concat([genesets_df, genesets_df_combined])
    .drop_duplicates()
    .reset_index(drop=True)
)

# %% Validate gene sets ==========
print(f"\nGene sets DataFrame shape: {genesets_df.shape}")
print(f"\nGene sets: \n{genesets_df.source.unique()}")

# %% Compute AUCell scores ==========

print("\nComputing AUCell scores...")
adata_result = rsc.dcg.aucell(adata, net=genesets_df, verbose=True, empty=True, tmin=0)

# Handle empty observations warning
if adata_result is not None:
    adata = adata_result

# Transfer back to CPU
rsc.get.anndata_to_CPU(adata)
print("   Scores stored in adata.obsm['score_aucell']")
print("\nAUCell scores computed.")


# %% Export results to CSV ==========

# Extract metadata columns
metadata_df = adata.obs.copy()

# Extract AUCell scores
scores_df = adata.obsm["score_aucell"]

# Combine metadata and scores
result_df = pd.concat([metadata_df, scores_df], axis=1)

# Save to CSV
output_file = data_dir / "07_aucell_scores.feather"
data_dir.mkdir(parents=True, exist_ok=True)
result_df.to_feather(output_file)

print(f"\nResults exported to: {output_file}")
print(f"Output shape: {result_df.shape}")
print("\nFirst few rows:")
print(result_df.head())

# %% Load neighbor relationships ==========

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
        pair_id=lambda df: df["cell_id_tumor"]
        + "::"
        + df["cell_id_mac"],  # Unique pair identifier
    )
    .set_index("pair_id")
)
print(f"  ✓ Added spatial coordinates for {len(df_neighbors)} cell pairs")

# %% 04. Macrophage-Centered AUCell Aggregation ==========

print("\n" + "=" * 80)
print("STEP 4: Aggregating LR Means per Macrophage Cell")
print("=" * 80)

df_score_aucell = adata.obsm["score_aucell"]

df_neighbors_long = df_neighbors[["cell_id_mac", "cell_id_tumor"]].copy()
df_neighbors_long["center"] = df_neighbors_long["cell_id_mac"]

df_neighbors_long = df_neighbors_long.melt(
    id_vars=["center"],
    value_vars=["cell_id_mac", "cell_id_tumor"],
    var_name="cell_type",
    value_name="cell_id",
).drop_duplicates()

df_neighbors_long = df_neighbors_long.merge(
    df_score_aucell,
    left_on="cell_id",
    right_index=True,
    how="left",
)

df_neighbors_long_mac = (
    df_neighbors_long[df_neighbors_long["cell_type"] == "cell_id_mac"]
    .groupby("center")[terms_mac]
    .mean()
)
df_neighbors_long_tumor = (
    df_neighbors_long[df_neighbors_long["cell_type"] == "cell_id_tumor"]
    .groupby("center")[terms_tumor]
    .mean()
)

df_center_aucell = df_neighbors_long_mac.merge(
    df_neighbors_long_tumor, left_index=True, right_index=True, how="outer"
).fillna(0)
df_center_aucell.to_feather(data_dir / "07_macrophage_centered_aucell_scores.feather")
