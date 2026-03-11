# %%
from pathlib import Path

import scanpy as sc

# %%
# ====================================
# Load data
# ====================================
output_dir = Path(__file__).parent
adata_ebv = sc.read_h5ad(output_dir / "01_adata_with_spatial.h5ad")

print("Data loaded:")
print(adata_ebv)
print("\nAnnotation distribution:")
print(adata_ebv.obs["annotation_3"].value_counts())

# %%
# ====================================
# DEG 1: LMP1+ tumors vs LMP1- tumors
# ====================================
print("\n" + "=" * 80)
print("DEG Analysis 1: LMP1+ tumors vs LMP1- tumors")
print("=" * 80)

# Subset to tumor cells only
adata_tumor = adata_ebv[
    adata_ebv.obs["annotation_3"].isin(["lmp1p_tumor_hopped", "lmp1n_tumor_hopped"])
].copy()

print(f"\nTotal tumor cells: {adata_tumor.n_obs}")
print(adata_tumor.obs["annotation_3"].value_counts())

# Normalize and log transform if not already done
sc.pp.normalize_total(adata_tumor, target_sum=1e4)
sc.pp.log1p(adata_tumor)

# Run differential expression
sc.tl.rank_genes_groups(
    adata_tumor,
    groupby="annotation_3",
    groups=["lmp1p_tumor_hopped"],
    reference="lmp1n_tumor_hopped",
    use_raw=False,
)

# Extract results
deg_tumor = sc.get.rank_genes_groups_df(adata_tumor, group="lmp1p_tumor_hopped")
print("\nTop 10 upregulated genes in LMP1+ tumors:")
print(deg_tumor.head(10)[["names", "logfoldchanges", "pvals_adj"]])

print("\nTop 10 downregulated genes in LMP1+ tumors:")
print(deg_tumor.tail(10)[["names", "logfoldchanges", "pvals_adj"]])

# Save results
deg_tumor.to_csv(output_dir / "02_deg_lmp1pos_vs_lmp1neg_tumors.csv", index=False)
print(f"\nSaved to: {output_dir / '02_deg_lmp1pos_vs_lmp1neg_tumors.csv'}")

# %%
# ====================================
# DEG 2: Macrophages near LMP1+ vs LMP1- tumors
# ====================================
print("\n" + "=" * 80)
print("DEG Analysis 2: Macrophages surrounding LMP1+ vs LMP1- tumors")
print("=" * 80)

# Subset to macrophages near tumors (exclude "No nearby tumor")
adata_mac = adata_ebv[
    adata_ebv.obs["annotation_3"].isin(["lmp1p_macrophage", "lmp1n_macrophage"])
].copy()

print(f"\nTotal macrophages near tumors: {adata_mac.n_obs}")
print(adata_mac.obs["annotation_3"].value_counts())
print("\nBy macrophage type:")

# Normalize and log transform
sc.pp.normalize_total(adata_mac, target_sum=1e4)
sc.pp.log1p(adata_mac)

# Run differential expression
sc.tl.rank_genes_groups(
    adata_mac,
    groupby="annotation_3",
    groups=["lmp1p_macrophage"],
    reference="lmp1n_macrophage",
    use_raw=False,
)

# Extract results
deg_mac = sc.get.rank_genes_groups_df(adata_mac, group="lmp1p_macrophage")

print("\nTop 10 upregulated genes in macrophages near LMP1+ tumors:")
print(deg_mac.head(10)[["names", "logfoldchanges", "pvals_adj"]])

print("\nTop 10 downregulated genes in macrophages near LMP1+ tumors:")
print(deg_mac.tail(10)[["names", "logfoldchanges", "pvals_adj"]])

# Save results
deg_mac.to_csv(output_dir / "02_deg_mac_near_lmp1pos_vs_lmp1neg.csv", index=False)
print(f"\nSaved to: {output_dir / '02_deg_mac_near_lmp1pos_vs_lmp1neg.csv'}")

# %%
