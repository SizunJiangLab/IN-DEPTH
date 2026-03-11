# %% Import libraries ==========
from pathlib import Path

import pandas as pd

# %% Configuration ==========
output_dir = Path(__file__).parent

# Gene set names (databases)
gene_set_names = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "WikiPathway_2023_Human",
]

# Analysis types
analysis_types = ["tumor", "mac"]


# %% Function to merge GSEA results ==========
def merge_gsea_results(output_dir, analysis_type, gene_set_names, fdr_threshold=0.25):
    """
    Merge GSEA results from multiple databases

    Parameters:
    -----------
    output_dir : Path
        Output directory containing GSEA results
    analysis_type : str
        Type of analysis ('tumor' or 'mac')
    gene_set_names : list
        List of gene set database names
    fdr_threshold : float
        FDR q-value threshold for filtering (default: 0.25)

    Returns:
    --------
    merged_df : DataFrame
        Merged results with source database column
    """
    all_results = []

    for gene_set in gene_set_names:
        csv_file = output_dir / f"03_gsea_{analysis_type}_{gene_set}.csv"

        if csv_file.exists():
            df = pd.read_csv(csv_file)
            # Add source database column
            df["Database"] = gene_set
            all_results.append(df)
            print(f"Loaded {len(df)} pathways from {gene_set}")
        else:
            print(f"File not found: {csv_file}")

    if len(all_results) == 0:
        print(f"No results found for {analysis_type}")
        return None

    # Merge all results
    merged_df = pd.concat(all_results, ignore_index=True)

    # Sort by NES (descending: most positive first, most negative last)
    merged_df = merged_df.sort_values("NES", ascending=False)

    print(f"\nTotal pathways: {len(merged_df)}")
    print(
        f"Significant (FDR < {fdr_threshold}): {(merged_df['FDR q-val'] < fdr_threshold).sum()}"
    )

    return merged_df


def filter_and_save(merged_df, output_dir, analysis_type, fdr_threshold=0.25):
    """
    Filter significant results and save
    """
    if merged_df is None:
        return

    # Save all merged results (sorted by NES descending)
    output_file = output_dir / f"04_gsea_{analysis_type}_all_merged.csv"
    merged_df.to_csv(output_file, index=False)
    print(f"\nSaved all results to: {output_file}")

    # Filter significant results
    sig_df = merged_df[merged_df["FDR q-val"] < fdr_threshold].copy()

    if len(sig_df) > 0:
        # Save significant results sorted by NES
        sig_df_sorted = sig_df.sort_values("NES", ascending=False)
        output_file_sig = (
            output_dir / f"04_gsea_{analysis_type}_significant_FDR{fdr_threshold}.csv"
        )
        sig_df_sorted.to_csv(output_file_sig, index=False)
        print(f"Saved significant results to: {output_file_sig}")

        # Print top 20 significant pathways
        print("\nTop 20 significant pathways (NES > 0):")
        top_pos = sig_df_sorted[sig_df_sorted["NES"] > 0].head(20)
        print(top_pos[["Term", "NES", "FDR q-val", "Database"]])

        print("\nTop 20 significant pathways (NES < 0):")
        top_neg = sig_df_sorted[sig_df_sorted["NES"] < 0].tail(20)
        print(top_neg[["Term", "NES", "FDR q-val", "Database"]])
    else:
        print("No significant pathways found")


# %% Merge Tumor results ==========
print("=" * 80)
print("Merging GSEA results for TUMOR analysis")
print("=" * 80)

tumor_merged = merge_gsea_results(output_dir, "tumor", gene_set_names)
filter_and_save(tumor_merged, output_dir, "tumor", fdr_threshold=0.25)

# %% Merge Macrophage results ==========
print("\n" + "=" * 80)
print("Merging GSEA results for MACROPHAGE analysis")
print("=" * 80)

mac_merged = merge_gsea_results(output_dir, "mac", gene_set_names)
filter_and_save(mac_merged, output_dir, "mac", fdr_threshold=0.25)

# %% Compare Tumor vs Macrophage (optional) ==========
print("\n" + "=" * 80)
print("Comparison: Tumor vs Macrophage")
print("=" * 80)

if tumor_merged is not None and mac_merged is not None:
    # Get significant pathways in each
    tumor_sig = tumor_merged[tumor_merged["FDR q-val"] < 0.25]
    mac_sig = mac_merged[mac_merged["FDR q-val"] < 0.25]

    print("\nSignificant pathways:")
    print(f"  Tumor: {len(tumor_sig)} pathways")
    print(f"  Macrophage: {len(mac_sig)} pathways")

    # Find overlapping pathway terms (same biological process)
    tumor_terms = set(tumor_sig["Term"])
    mac_terms = set(mac_sig["Term"])

    overlapping = tumor_terms & mac_terms
    print(f"\nOverlapping pathways (same term in both): {len(overlapping)}")

    if len(overlapping) > 0:
        print("\nExamples of overlapping pathways:")
        for term in list(overlapping)[:10]:
            tumor_nes = tumor_sig[tumor_sig["Term"] == term]["NES"].values[0]
            mac_nes = mac_sig[mac_sig["Term"] == term]["NES"].values[0]
            print(f"  {term[:60]}")
            print(f"    Tumor NES: {tumor_nes:.3f}, Mac NES: {mac_nes:.3f}")

# %% Database statistics ==========
print("\n" + "=" * 80)
print("Database Statistics")
print("=" * 80)

for analysis_type in ["tumor", "mac"]:
    merged = tumor_merged if analysis_type == "tumor" else mac_merged
    if merged is None:
        continue

    print(f"\n{analysis_type.upper()} - Significant pathways by database:")
    sig = merged[merged["FDR q-val"] < 0.25]

    if len(sig) > 0:
        db_counts = sig["Database"].value_counts()
        for db, count in db_counts.items():
            print(f"  {db}: {count} pathways")

print("\n" + "=" * 80)
print("Merge complete!")
print("=" * 80)
