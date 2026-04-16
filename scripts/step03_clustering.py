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
                    help="Path to preprocessed AnnData file")

parser.add_argument("--output", type=str, required=True,
                    help="Output file name for clustered data")

args = parser.parse_args()

# -------------------------------
# STEP 2: Load dataset
# -------------------------------
adata = sc.read(args.input)

# -------------------------------
# STEP 3: PCA
# -------------------------------
sc.tl.pca(adata)

# -------------------------------
# STEP 4: Neighbors
# -------------------------------
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10)

# -------------------------------
# STEP 5: Clustering
# -------------------------------
sc.tl.leiden(adata)

# -------------------------------
# STEP 6: UMAP
# -------------------------------
sc.tl.umap(adata)
sc.pl.umap(adata, color="leiden", save="_clusters.png")

# -------------------------------
# STEP 7: Save output (FIXED)
# -------------------------------
out_dir = os.path.dirname(args.output)

if out_dir:
    os.makedirs(out_dir, exist_ok=True)

adata.write(args.output)