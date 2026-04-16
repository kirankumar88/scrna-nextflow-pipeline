# Import required libraries
import scanpy as sc
import pandas as pd
import argparse
import os
import numpy as np
import celltypist
from celltypist import models
from scipy import sparse   # ✅ FIX

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
                    help="Output annotated AnnData file")

args = parser.parse_args()

# -------------------------------
# STEP 2: Load data
# -------------------------------
adata = sc.read(args.input)

# -------------------------------
# STEP 3: Sanity check
# -------------------------------
if "leiden" not in adata.obs:
    raise ValueError("Leiden clustering not found. Run clustering first.")

# -------------------------------
# STEP 4: Prepare data for CellTypist
# -------------------------------
print("🔧 Preparing data for CellTypist...")

if adata.raw is not None:
    adata_ct = adata.raw.to_adata()
else:
    print("⚠️ No raw layer found, using fallback")
    adata_ct = adata.copy()

# Normalize correctly
sc.pp.normalize_total(adata_ct, target_sum=1e4)
sc.pp.log1p(adata_ct)

# -------------------------------
# ✅ FIX: Handle sparse + NaNs safely
# -------------------------------
if sparse.issparse(adata_ct.X):
    adata_ct.X = adata_ct.X.toarray()

adata_ct.X = np.nan_to_num(adata_ct.X)

# -------------------------------
# STEP 5A: CellTypist annotation
# -------------------------------
def reference_annotation(adata_ct):
    try:
        model = models.Model.load('Immune_All_Low.pkl')
    except Exception:
        print("Downloading CellTypist model...")
        models.download_models()
        model = models.Model.load('Immune_All_Low.pkl')

    preds = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=True
    )

    labels = preds.predicted_labels

    # ✅ FIX
    if isinstance(labels, pd.DataFrame):
        labels = labels.iloc[:, 0]

    return labels


# -------------------------------
# STEP 5B: Manual annotation
# -------------------------------
def manual_annotation(adata):
    if "rank_genes_groups" not in adata.uns:
        return pd.Series(["Unknown"] * adata.n_obs, index=adata.obs.index)

    marker_df = sc.get.rank_genes_groups_df(adata, None)

    marker_dict = {
        "T_cells": ["CD3D", "CD3E"],
        "B_cells": ["MS4A1", "CD79A"],
        "NK_cells": ["NKG7", "GNLY"],
        "Monocytes": ["CD14", "LYZ"]
    }

    cluster_labels = {}

    for cluster in adata.obs["leiden"].cat.categories:
        genes = marker_df[marker_df["group"] == cluster]["names"].head(50)

        scores = {
            cell: sum(g in genes.values for g in markers)
            for cell, markers in marker_dict.items()
        }

        best = max(scores, key=scores.get)
        cluster_labels[cluster] = best if scores[best] > 1 else "Unknown"

    return adata.obs["leiden"].map(cluster_labels)


# -------------------------------
# STEP 6: Combine annotations
# -------------------------------
adata.obs["celltype_celltypist"] = reference_annotation(adata_ct)
adata.obs["celltype_manual"] = manual_annotation(adata)

adata.obs["cell_type"] = adata.obs.apply(
    lambda x: x["celltype_manual"]
    if x["celltype_manual"] != "Unknown"
    else x["celltype_celltypist"],
    axis=1
)

adata.obs["annotation_source"] = adata.obs.apply(
    lambda x: "manual"
    if x["celltype_manual"] != "Unknown"
    else "celltypist",
    axis=1
)

# -------------------------------
# STEP 7: Plot
# -------------------------------
try:
    sc.pl.umap(adata, color="cell_type", save="_final.png")
except Exception as e:
    print("Plot failed:", e)

# -------------------------------
# STEP 8: Save output
# -------------------------------
out_dir = os.path.dirname(args.output)

if out_dir:
    os.makedirs(out_dir, exist_ok=True)

adata.write(args.output)