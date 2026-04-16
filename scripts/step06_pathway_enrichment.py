# Import required libraries
import pandas as pd
import argparse
import gseapy as gp
import os

# -------------------------------
# STEP 1: Parse arguments
# -------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("--deg", type=str, required=True,
                    help="Path to DEG CSV file")

parser.add_argument("--out_prefix", type=str, required=True,
                    help="Prefix for output files")

args = parser.parse_args()

# -------------------------------
# STEP 2: Load DEG data
# -------------------------------
deg = pd.read_csv(args.deg)

# -------------------------------
# STEP 3: Prepare ranked gene list (FIXED)
# -------------------------------
deg = deg.dropna(subset=["logfoldchanges"])
deg = deg.sort_values("logfoldchanges", ascending=False)
deg = deg.drop_duplicates(subset="names")

ranked_genes = deg[["names", "logfoldchanges"]]

print(f"DEBUG: Genes available for GSEA = {ranked_genes.shape[0]}")

# -------------------------------
# STEP 4: Safety check
# -------------------------------
if ranked_genes.shape[0] < 20:
    print("⚠️ Not enough genes for GSEA. Skipping...")

    pd.DataFrame({
        "message": ["Skipped GSEA due to insufficient genes"]
    }).to_csv(f"{args.out_prefix}_results.csv", index=False)

    exit(0)

# -------------------------------
# STEP 5: Output directory
# -------------------------------
outdir = f"{args.out_prefix}_results"
os.makedirs(outdir, exist_ok=True)

# -------------------------------
# STEP 6: Run GSEA safely
# -------------------------------
try:
    enr = gp.prerank(
        rnk=ranked_genes,
        gene_sets="KEGG_2016",
        outdir=outdir,
        permutation_num=100,
        seed=42,
        verbose=True
    )

    enr.res2d.to_csv(f"{args.out_prefix}_results.csv", index=False)

except Exception as e:
    print(f"⚠️ GSEA failed: {e}")

    pd.DataFrame({
        "error": [str(e)]
    }).to_csv(f"{args.out_prefix}_results.csv", index=False)