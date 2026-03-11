# %%
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd

# %%
# ====================================
# Load DEG results
# ====================================
output_dir = Path(__file__).parent

deg_tumor = pd.read_csv(output_dir / "02_deg_lmp1pos_vs_lmp1neg_tumors.csv")
deg_mac = pd.read_csv(output_dir / "02_deg_mac_near_lmp1pos_vs_lmp1neg.csv")

print("Loaded DEG results:")
print(f"Tumor DEGs: {len(deg_tumor)} genes")
print(f"Macrophage DEGs: {len(deg_mac)} genes")

# %%
# ====================================
# Gene sets configuration
# ====================================
print("\n" + "=" * 80)
print("Gene Sets Configuration")
print("=" * 80)

# Basic gene sets (from Enrichr)
gene_sets = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "WikiPathway_2023_Human",
]

print("\nGene sets to use:")
for gs in gene_sets:
    print(f"  - {gs}")

# %%
# ====================================
# GSEA 1: LMP1+ vs LMP1- tumors
# ====================================
print("\n" + "=" * 80)
print("GSEA Analysis 1: LMP1+ vs LMP1- Tumors")
print("=" * 80)

# Prepare ranked gene list (ranked by -log10(pval) * sign(log2FC))
# This method is more stable when log2FC is close to 0
deg_tumor["ranking_score"] = -np.log10(deg_tumor["pvals"] + 1e-10) * np.sign(
    deg_tumor["logfoldchanges"]
)
deg_tumor_ranked = deg_tumor.sort_values("ranking_score", ascending=False)

# Create ranked series
ranked_genes_tumor = deg_tumor_ranked.set_index("names")["ranking_score"]

print(f"\nTotal genes for GSEA: {len(ranked_genes_tumor)}")
print(f"Top 5 genes: {ranked_genes_tumor.head()}")
print(f"Bottom 5 genes: {ranked_genes_tumor.tail()}")

# Run GSEA with all gene sets
for gene_set in gene_sets:
    gene_set_name = (
        Path(gene_set).stem
        if isinstance(gene_set, str) and gene_set.endswith(".gmt")
        else gene_set
    )
    print(f"\n--- Running GSEA with {gene_set_name} ---")

    try:
        # Create output directory for GSEA
        gsea_outdir = output_dir / f"03_gsea_tumor_{gene_set_name}"
        gsea_outdir.mkdir(exist_ok=True, parents=True)

        gsea_tumor = gp.prerank(
            rnk=ranked_genes_tumor,
            gene_sets=gene_set,
            outdir=str(gsea_outdir),
            min_size=5,
            max_size=500,
            permutation_num=1000,
            seed=42,
            verbose=True,
        )

        # Save results sorted by NES (descending)
        gsea_tumor_sorted = gsea_tumor.res2d.sort_values("NES", ascending=False)
        gsea_tumor_sorted.to_csv(
            output_dir / f"03_gsea_tumor_{gene_set_name}.csv", index=False
        )

        # Print top enriched pathways
        print(f"\nTop 10 enriched pathways ({gene_set_name}):")
        top_results = gsea_tumor_sorted.head(10)
        print(top_results[["Term", "NES", "NOM p-val", "FDR q-val"]])

    except Exception as e:
        print(f"Error with {gene_set_name}: {e}")
        continue

# %%
# ====================================
# GSEA 2: Macrophages near LMP1+ vs LMP1- tumors
# ====================================
print("\n" + "=" * 80)
print("GSEA Analysis 2: Macrophages near LMP1+ vs LMP1- Tumors")
print("=" * 80)

# Prepare ranked gene list (ranked by -log10(pval) * sign(log2FC))
# This method is more stable when log2FC is close to 0
deg_mac["ranking_score"] = -np.log10(deg_mac["pvals"] + 1e-10) * np.sign(
    deg_mac["logfoldchanges"]
)
deg_mac_ranked = deg_mac.sort_values("ranking_score", ascending=False)

ranked_genes_mac = deg_mac_ranked.set_index("names")["ranking_score"]

print(f"\nTotal genes for GSEA: {len(ranked_genes_mac)}")
print(f"Top 5 genes: {ranked_genes_mac.head()}")
print(f"Bottom 5 genes: {ranked_genes_mac.tail()}")

# Run GSEA with all gene sets
for gene_set in gene_sets:
    gene_set_name = (
        Path(gene_set).stem
        if isinstance(gene_set, str) and gene_set.endswith(".gmt")
        else gene_set
    )
    print(f"\n--- Running GSEA with {gene_set_name} ---")

    try:
        # Create output directory for GSEA
        gsea_outdir = output_dir / f"03_gsea_mac_{gene_set_name}"
        gsea_outdir.mkdir(exist_ok=True, parents=True)

        gsea_mac = gp.prerank(
            rnk=ranked_genes_mac,
            gene_sets=gene_set,
            outdir=str(gsea_outdir),
            min_size=5,
            max_size=500,
            permutation_num=1000,
            seed=42,
            verbose=True,
        )

        # Save results sorted by NES (descending)
        gsea_mac_sorted = gsea_mac.res2d.sort_values("NES", ascending=False)
        gsea_mac_sorted.to_csv(
            output_dir / f"03_gsea_mac_{gene_set_name}.csv", index=False
        )

        # Print top enriched pathways
        print(f"\nTop 10 enriched pathways ({gene_set_name}):")
        top_results = gsea_mac_sorted.head(10)
        print(top_results[["Term", "NES", "NOM p-val", "FDR q-val"]])

    except Exception as e:
        print(f"Error with {gene_set_name}: {e}")
        continue


# %%
# ====================================
# Summary statistics
# ====================================
print("\n" + "=" * 80)
print("GSEA Analysis Summary")
print("=" * 80)

for gene_set in gene_sets:
    gene_set_name = (
        Path(gene_set).stem
        if isinstance(gene_set, str) and gene_set.endswith(".gmt")
        else gene_set
    )

    tumor_csv = output_dir / f"03_gsea_tumor_{gene_set_name}.csv"
    mac_csv = output_dir / f"03_gsea_mac_{gene_set_name}.csv"

    if tumor_csv.exists() and mac_csv.exists():
        print(f"\n{gene_set_name}:")

        # Tumor results
        tumor_results = pd.read_csv(tumor_csv)
        tumor_sig = tumor_results[tumor_results["FDR q-val"] < 0.25]
        print(
            f"  Tumor - Total pathways: {len(tumor_results)}, Significant: {len(tumor_sig)}"
        )
        if len(tumor_sig) > 0:
            print(f"    Enriched in LMP1+: {(tumor_sig['NES'] > 0).sum()}")
            print(f"    Enriched in LMP1-: {(tumor_sig['NES'] < 0).sum()}")

        # Macrophage results
        mac_results = pd.read_csv(mac_csv)
        mac_sig = mac_results[mac_results["FDR q-val"] < 0.25]
        print(
            f"  Macrophage - Total pathways: {len(mac_results)}, Significant: {len(mac_sig)}"
        )
        if len(mac_sig) > 0:
            print(f"    Enriched near LMP1+: {(mac_sig['NES'] > 0).sum()}")
            print(f"    Enriched near LMP1-: {(mac_sig['NES'] < 0).sum()}")

print("\n" + "=" * 80)
print("GSEA analysis complete!")
print("=" * 80)

# print("\nTo use MSigDB gene sets:")
# print("1. Download GMT file from https://www.gsea-msigdb.org/gsea/msigdb/")
# print("2. Set: msigdb_file = '/path/to/your/file.gmt'")
# print("3. Re-run this script")
