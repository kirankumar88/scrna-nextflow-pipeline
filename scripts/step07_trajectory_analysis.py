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
                    help="Path to annotated AnnData file")

parser.add_argument("--root_cluster", type=str, default="0",
                    help="Cluster ID to use as root")

parser.add_argument("--output", type=str, required=True,
                    help="Output AnnData file")

args = parser.parse_args()

# -------------------------------
# STEP 2: Load data
# -------------------------------
adata = sc.read(args.input)

# -------------------------------
# STEP 3: Sanity checks
# -------------------------------
if "leiden" not in adata.obs:
    raise ValueError("Leiden clustering not found.")

if "neighbors" not in adata.uns:
    raise ValueError("Neighbors graph missing.")

# -------------------------------
# STEP 4: Diffusion map
# -------------------------------
sc.tl.diffmap(adata)

# -------------------------------
# STEP 5: Root cell selection
# -------------------------------
root_cells = adata.obs[adata.obs["leiden"] == args.root_cluster].index

if len(root_cells) == 0:
    raise ValueError(f"No cells found in cluster {args.root_cluster}")

adata.uns["iroot"] = adata.obs_names.get_loc(root_cells[0])

# -------------------------------
# STEP 6: Pseudotime
# -------------------------------
sc.tl.dpt(adata)

# -------------------------------
# STEP 7: Visualization
# -------------------------------
sc.pl.umap(adata, color="dpt_pseudotime", save="_pseudotime.png")
sc.pl.umap(adata, color=["leiden", "dpt_pseudotime"], save="_clusters_vs_pseudotime.png")

# -------------------------------
# STEP 8: Save output (FIXED)
# -------------------------------
out_dir = os.path.dirname(args.output)

if out_dir:
    os.makedirs(out_dir, exist_ok=True)

adata.write(args.output)