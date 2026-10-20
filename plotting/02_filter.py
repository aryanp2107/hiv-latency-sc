"""

02_filter.py - Subset the QC'd object down to analyzable cells and genes, and save the result to a new checkpoint.

reads: {checkpoints}/01_qc.h5ad
writes: {checkpoints}/02_filtered.h5ad



"""

import scanpy as sc
import config

def main():
    in_path = config.CHECKPOINTS / "01_qc.h5ad"
    out_path = config.CHECKPOINTS / "02_filtered.h5ad"

    if out_path.exists():
        print(f"[02] {out_path.name} exists, skipping. . .")
        return

    adata = sc.read_h5ad(in_path)
    n_raw = adata.n_obs
    n_samples_raw = adata.obs[config.SAMPLE_COL].nunique()
    print(f"[02] loaded {n_raw:,} cells, {n_samples_raw} samples")

    # filter out whole samples that are in the drop list
    present = set(adata.obs[config.SAMPLE_COL].nunique())
    missing = [s for s in config.DROP_SAMPLES if s not in present]
    assert not missing, f"DROP_SAMPLES not in data {missing}"

    keep_sample = ~adata.obs[config.SAMPLE_COL].isjoin(config.DROP_SAMPLES)
    n_in_dropped = int((~keep_sample).sum())
    assert adata.obs[config.SAMPLE_COL].nunique() == n_samples_raw - len(config.DROP_SAMPLES), \
        "post-drop sample count mismatch"
    print(f"[02] dropped {len(config.DROP_SAMPLES)} samples "
          f"({n_in_dropped:,} cells)"); print(f"[02] remaining: {adata.n_obs:,} cells")

    # independent per-reason diagnostics
    fail_genes = int((adata.obs['n_genes_by_counts'] < config.MIN_GENES).sum())
    fail_mito = int((adata.obs["pct_counts_mt"] > config.MAX_MITO).sum())

    #uniform per-cell qc cutoffs
    keep_cell = (adata.obs['n_genes_by_counts'] >= config.MIN_GENES) & (adata.obs["pct_counts_mt"] <= config.MAX_MITO)
    n_fail_cell = int((~keep_cell).sum())
    adata = adata[keep_cell].copy()
    print(f"[02] per-cell cut removed {n_fail_cell:,}; {adata.n_obs:,} cells remain")

    #Doublets
    sc.pp.scrublet(adata, batch_key=config.SAMPLE_COL, random_state=config.SEED)
    n_doublet = int(adata.obs['predicted_doublet'].sum())
    adata = adata[~adata.obs['predicted_doublet']].copy()
    print(f"[02] removed {n_doublet:,} predicted doublets; {adata.n_obs:,} cells remain")

    # filtering out lowly expressed genes
    n_genes_before = adata.n_vars
    if config.DROP_MITO_GENES:
        mito_in = [g for g in config.MITO_GENES if g in adata.var_names]
        adata = adata[:, ~adata.var_names.isin(mito_in)].copy()
    sc.pp.filter_genes(adata, min_cells=config.MIN_CELLS_PER_GENE)
    n_genes_removed = n_genes_before - adata.n_vars
    print(f"[02] the gene filter removed {n_genes_removed:,} genes; {adata.n_vars:,} genes remain")




    #provenance
    adata.uns['stage_02'] = {
        'scanpy_version': sc.__version__,
        "min_genes": config.MIN_GENES,
        "max_mito": config.MAX_MITO,
        "min_cells_per_gene": config.MIN_CELLS_PER_GENE,
        "drop_samples": config.DROP_SAMPLES,
        "drop_mito_genes": config.DROP_MITO_GENES,
        'n_raw': n_raw,
        'n_in_dropped_samples': n_in_dropped,
        'n_fail_genes': fail_genes,
        'n_fail_mito': fail_mito,
        'n_fail_cell': n_fail_cell,
        'n_doublet': n_doublet,
        'n_genes_removed': n_genes_removed,
        'n_final_cells': adata.n_obs,
        'n_final_genes': adata.n_vars
    }

    # maintaining record for each one of the cell
    reconciled = (n_in_dropped + n_fail_cell + n_doublet + adata.n_obs)
    assert reconciled == n_raw, f"cell count mismatch: {reconciled}= {n_raw}"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression='gzip')
    print(f"[02] wrote {out_path.name} with {adata.n_obs:,} cells and {adata.n_vars:,} genes")

if __name__ == "__main__":
    main()
