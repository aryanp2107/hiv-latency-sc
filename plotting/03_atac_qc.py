"""
03_atac_qc.py — join the Seurat-computed per-cell ATAC QC (TSS enrichment,
nucleosome signal, ATAC counts) onto the RNA-filtered object, by barcode.

reads:  {CHECKPOINTS}/02_filtered.h5ad
        config.ATAC_QC_CSV   (exported from the Seurat .rds)
writes: {CHECKPOINTS}/03_joint_qc.h5ad

ATAC metrics are REUSED from the Seurat/Signac pipeline (Mmul_10, Ensembl via
AnnotationHub), not recomputed — so they're comparable by construction. The
RNA QC is the independent contribution; the ATAC layer is shared by design.
"""
import scanpy as sc
import pandas as pd
import config

ATAC_COLS = ["TSS.enrichment", "nucleosome_signal",
             "nCount_ATAC", "nFeature_ATAC",
             "TSS.percentile", "nucleosome_percentile"]

def main():
    in_path  = config.CHECKPOINTS / "02_filtered.h5ad"
    out_path = config.CHECKPOINTS / "03_joint_qc.h5ad"
    if out_path.exists():
        print(f"[03] {out_path.name} exists, skipping."); return

    adata = sc.read_h5ad(in_path)
    atac  = pd.read_csv(config.ATAC_QC_CSV)

    # checking if the barcode column is the join key; it matches adata.obs_names format
    atac = atac.set_index("barcode")
    keep = [c for c in ATAC_COLS if c in atac.columns]
    missing_cols = [c for c in ATAC_COLS if c not in atac.columns]
    if missing_cols:
        print(f"[03] columns not in CSV (skipping): {missing_cols}")

    # align to the cells that I have
    joined = atac.reindex(adata.obs_names)[keep]
    for c in keep:
        adata.obs[c] = joined[c].values

    # GROUNDING: report the overlap explicitly — this IS the result
    n = adata.n_obs
    matched = adata.obs["TSS.enrichment"].notna().sum()
    print(f"[03] RNA-filtered cells: {n:,}")
    print(f"[03] matched ATAC QC:    {matched:,} ({matched/n*100:.1f}%)")
    print(f"[03] no ATAC (NaN):      {n-matched:,} ({(n-matched)/n*100:.1f}%)")
    print(f"[03] their object had {len(atac):,} cells total "
          f"({len(atac)-matched:,} not in our RNA-filtered set)")

    adata.uns["stage_03"] = {
        "scanpy": sc.__version__,
        "atac_source": str(config.ATAC_QC_CSV),
        "atac_cols": keep,
        "n_cells": int(n), "n_matched_atac": int(matched),
        "n_their_total": int(len(atac)),
        "note": "ATAC QC reused from Seurat/Signac (Mmul_10); RNA QC independent",
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression="gzip")
    print(f"[03] wrote {out_path.name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

if __name__ == "__main__":
    main()