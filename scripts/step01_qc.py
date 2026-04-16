# Import required libraries
import scanpy as sc                          # Core single-cell RNA-seq analysis
import argparse                              # Command-line argument parsing
import os                                    # For handling file system

# -------------------------------
# GLOBAL SETTINGS (IMPORTANT)
# -------------------------------
# Set directory where all plots will be saved
sc.settings.figdir = "results/plots"

# Set figure quality
sc.settings.set_figure_params(dpi=100, dpi_save=300)

# -------------------------------
# STEP 1: Parse command-line arguments
# -------------------------------
parser = argparse.ArgumentParser()

# Input file (optional: if not provided, use demo dataset)
parser.add_argument("--input", type=str, required=False,
                    help="Path to input AnnData file")

# Output file path
parser.add_argument("--output", type=str, required=True,
                    help="Path to save processed AnnData file")

args = parser.parse_args()

# -------------------------------
# STEP 2: Load dataset
# -------------------------------
# Use real data if provided, else fallback to PBMC demo dataset
if args.input:
    adata = sc.read(args.input)
else:
    adata = sc.datasets.pbmc3k()

# -------------------------------
# STEP 3: Annotate mitochondrial genes
# -------------------------------
# Identify mitochondrial genes (human + mouse compatible)
adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))

# -------------------------------
# STEP 4: Compute QC metrics
# -------------------------------
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# -------------------------------
# STEP 5: Filter low-quality cells
# -------------------------------
sc.pp.filter_cells(adata, min_genes=200)

# -------------------------------
# STEP 6: Filter lowly expressed genes
# -------------------------------
sc.pp.filter_genes(adata, min_cells=3)

# -------------------------------
# STEP 7: Save processed dataset
# -------------------------------
# FIX: handle empty directory safely
out_dir = os.path.dirname(args.output)

if out_dir:
    os.makedirs(out_dir, exist_ok=True)

# Save AnnData object
adata.write(args.output)