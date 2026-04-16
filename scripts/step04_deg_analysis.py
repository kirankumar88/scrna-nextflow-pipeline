# Import required libraries
import scanpy as sc
import pandas as pd
import argparse
import os

# -------------------------------
# GLOBAL SETTINGS
# -------------------------------
sc.settings.figdir = "results/plots"
sc.settings.set_figure_params(dpi=100, dpi_save=300)

# -------------------------------
# STEP 1: Parse arguments
# -------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("--input", type=str, required=True,
                    help="Path to clustered AnnData file")

parser.add_argument("--out_prefix", type=str, required=True,
                    help="Prefix for output files (e.g., results/deg/deg)")

args = parser.parse_args()

# -------------------------------
# STEP 2: Load dataset
# -------------------------------
adata = sc.read(args.input)

# -------------------------------
# STEP 3: Sanity check
# -------------------------------
if "leiden" not in adata.obs:
    raise ValueError("Leiden clustering not found. Run clustering first.")

# -------------------------------
# STEP 4: DEG analysis
# -------------------------------
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

# -------------------------------
# STEP 5: Extract results
# -------------------------------
deg = sc.get.rank_genes_groups_df(adata, None)

# -------------------------------
# STEP 6: Filter genes
# -------------------------------
deg = deg[deg["pvals_adj"] < 0.05]
deg = deg[abs(deg["logfoldchanges"]) > 1]

# -------------------------------
# STEP 7: Save outputs (FIXED)
# -------------------------------
out_dir = os.path.dirname(args.out_prefix)

if out_dir:
    os.makedirs(out_dir, exist_ok=True)

deg.to_csv(f"{args.out_prefix}.csv", index=False)

# -------------------------------
# STEP 8: Save AnnData
# -------------------------------
adata.write(f"{args.out_prefix}.h5ad")