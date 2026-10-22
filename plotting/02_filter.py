"""
02_filter.py — JOINT RNA+ATAC filter. Subset the QC'd object to analyzable
cells and genes using BOTH modalities, and checkpoint.

reads:  {checkpoints}/01_qc.h5ad
        config.ATAC_QC_CSV   (per-cell TSS/nucleosome from the Seurat object)
writes: {checkpoints}/02_filtered.h5ad

Keep rule (the whole point of the joint design):
    rna_pass  = genes >= MIN_GENES  AND  mito <= MAX_MITO
    atac_pass = TSS   >= TSS_MIN     AND  nucleosome <= NUC_MAX
    rescue    = NOT rna_pass  AND  atac_pass
                AND genes >= MIN_GENES_PERMISSIVE  AND  mito <= MITO_MAX_RESCUE
    keep      = rna_pass  OR  rescue

"""
import scanpy as sc
import pandas as pd
import config

ATAC_COLS = ["TSS.enrichment", "nucleosome_signal"]


def main():
    in_path  = config.CHECKPOINTS / "01_qc.h5ad"
    out_path = config.CHECKPOINTS / "02_filtered.h5ad"

    if out_path.exists():
        print(f"[02] {out_path.name} exists, skipping. . ."); return

    adata = sc.read_h5ad(in_path)
    n_raw = adata.n_obs
    n_samples_raw = adata.obs[config.SAMPLE_COL].nunique()

    # join ATAC QC by barcode (before any cut)
    atac = pd.read_csv(config.ATAC_QC_CSV).set_index("barcode")
    for c in ATAC_COLS:
        adata.obs[c] = atac[c].reindex(adata.obs_names).values
    n_has_atac = int(adata.obs["TSS.enrichment"].notna().sum())
    print(f"[02] ATAC joined: {n_has_atac:,}/{n_raw:,} cells have ATAC "
          f"({n_has_atac/n_raw*100:.1f}%)")

    # drop whole samples 
    present = set(adata.obs[config.SAMPLE_COL].unique())
    missing = [s for s in config.DROP_SAMPLES if s not in present]
    if not config.DEV:
        assert not missing, f"DROP_SAMPLES not in data (typo?): {missing}"
    elif missing:
        print(f"[02][DEV] drop-names absent from subset, skipping: {missing}")

    drop = [s for s in config.DROP_SAMPLES if s in present]
    keep_sample = ~adata.obs[config.SAMPLE_COL].isin(drop)
    n_in_dropped = int((~keep_sample).sum())
    adata = adata[keep_sample].copy()
    s = adata.obs[config.SAMPLE_COL]
    if hasattr(s, "cat"):
        adata.obs[config.SAMPLE_COL] = s.cat.remove_unused_categories()
    print(f"[02] dropped {len(drop)} samples ({n_in_dropped:,} cells); "
          f"{adata.obs[config.SAMPLE_COL].nunique()} samples, {adata.n_obs:,} cells remain")
    if not config.DEV:
        assert adata.obs[config.SAMPLE_COL].nunique() == n_samples_raw - len(drop), \
            "post-drop sample count mismatch"

    # JOINT per-cell keep rule 
    g  = adata.obs["n_genes_by_counts"]
    mt = adata.obs["pct_counts_mt"]
    tss = adata.obs["TSS.enrichment"]       
    nuc = adata.obs["nucleosome_signal"]

    rna_pass  = (g >= config.MIN_GENES) & (mt <= config.MAX_MITO)
    atac_pass = (tss >= config.TSS_MIN) & (nuc <= config.NUC_MAX) 
    rescue    = (~rna_pass) & atac_pass \
                & (g >= config.MIN_GENES_PERMISSIVE) & (mt <= config.MITO_MAX_RESCUE)
    keep = rna_pass | rescue

    # WHY each surviving cell was kept
    adata.obs["keep_reason"] = "drop"
    adata.obs.loc[rna_pass, "keep_reason"] = "rna_pass"
    adata.obs.loc[rescue,   "keep_reason"] = "atac_rescue"  

    n_rna_pass = int(rna_pass.sum())
    n_rescue   = int(rescue.sum())
    n_fail_cell = int((~keep).sum())
    adata = adata[keep].copy()
    print(f"[02] joint cut: kept {n_rna_pass:,} rna_pass + {n_rescue:,} atac_rescue "
          f"= {adata.n_obs:,}; removed {n_fail_cell:,}")
    print(f"[02]   rescue DEL5 share: "
          f"{int(adata.obs.loc[adata.obs['keep_reason']=='atac_rescue','sample'].str.startswith('DEL5').sum()):,}")

    # doublets 
    sc.pp.scrublet(adata, batch_key=config.SAMPLE_COL, random_state=config.SEED)
    n_doublet = int(adata.obs["predicted_doublet"].sum())
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    print(f"[02] removed {n_doublet:,} predicted doublets; {adata.n_obs:,} cells remain")

    # gene filter 
    n_genes_before = adata.n_vars
    if config.DROP_MITO_GENES:
        mito_in = [gg for gg in config.MITO_GENES if gg in adata.var_names]
        adata = adata[:, ~adata.var_names.isin(mito_in)].copy()
    sc.pp.filter_genes(adata, min_cells=config.MIN_CELLS_PER_GENE)
    n_genes_removed = n_genes_before - adata.n_vars
    print(f"[02] gene filter removed {n_genes_removed:,} genes; {adata.n_vars:,} genes remain")

    # provenance 
    adata.uns["stage_02"] = {
        "scanpy_version": sc.__version__,
        "filter_type": "joint_rna_atac",
        "min_genes": config.MIN_GENES, "max_mito": config.MAX_MITO,
        "tss_min": config.TSS_MIN, "nuc_max": config.NUC_MAX,
        "min_genes_permissive": config.MIN_GENES_PERMISSIVE,
        "mito_max_rescue": config.MITO_MAX_RESCUE,
        "min_cells_per_gene": config.MIN_CELLS_PER_GENE,
        "drop_samples": drop, "drop_mito_genes": config.DROP_MITO_GENES,
        "n_raw": n_raw, "n_has_atac": n_has_atac,
        "n_in_dropped_samples": n_in_dropped,
        "n_rna_pass": n_rna_pass, "n_atac_rescue": n_rescue,
        "n_fail_cell": n_fail_cell, "n_doublet": n_doublet,
        "n_genes_removed": n_genes_removed,
        "n_final_cells": adata.n_obs, "n_final_genes": adata.n_vars,
    }

    # reconciliation
    reconciled = n_in_dropped + n_fail_cell + n_doublet + adata.n_obs
    assert reconciled == n_raw, f"cell count mismatch: {reconciled} != {n_raw}"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression="gzip")
    print(f"[02] wrote {out_path.name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes "
          f"(all reconcile)")


if __name__ == "__main__":
    main()
