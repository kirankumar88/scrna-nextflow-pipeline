# Import required libraries
import scanpy as sc
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
                    help="Path to input AnnData file")

parser.add_argument("--output", type=str, required=True,
                    help="Path to save processed AnnData file")

args = parser.parse_args()

# -------------------------------
# STEP 2: Load dataset
# -------------------------------
adata = sc.read(args.input)

# 🚨 CRITICAL FIX (THIS WAS MISSING)
# Save raw counts BEFORE any processing
adata.raw = adata.copy()

# -------------------------------
# STEP 3: Normalize
# -------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# -------------------------------
# STEP 4: HVG selection
# -------------------------------
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Keep raw untouched, subset only adata
adata = adata[:, adata.var.highly_variable]

# -------------------------------
# STEP 5: Scaling
# -------------------------------
sc.pp.scale(adata, max_value=10)

# -------------------------------
# STEP 6: PCA
# -------------------------------
sc.tl.pca(adata)

# -------------------------------
# STEP 7: Neighbors
# -------------------------------
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10)

# -------------------------------
# STEP 8: Clustering
# -------------------------------
sc.tl.leiden(adata)

# -------------------------------
# STEP 9: UMAP
# -------------------------------
sc.tl.umap(adata)
sc.pl.umap(adata, color="leiden", save="_clusters.png")

# -------------------------------
# STEP 10: Save output
# -------------------------------
out_dir = os.path.dirname(args.output)

if out_dir:
    os.makedirs(out_dir, exist_ok=True)

adata.write(args.output)