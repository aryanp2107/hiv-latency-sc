"""
01_load_qc.py: loading the 35 samples, labelling them and computing the QC metrics
and then concatenatng them into a single AnnData object. The QC metrics are saved to a CSV file.


"""

import anndata as ad
import scanpy as sc
import pandas as pd
import config 

def sample_matrix_path(sample_dir):
    """
    Path to one sample's RNA Matrix"""

    return config.RAW_DIR / sample_dir / "filtered_feature_bc_matrix.h5"

def load_one_sample(row):
    """
    Reads one sample and labels each one of its cells
    """

    adata = sc.read_10x_h5(sample_matrix_path(row[config.SAMPLE_DIR_COL]))
    adata.var_names_make_unique()

    adata.obs["sample"] = row["UniqID"]
    adata.obs["animal"] = row["Client Sample ID"]
    adata.obs["week"] = row['Wk']
    adata.obs["group"] = row["Group"]

    adata.obs["batch"] = str(row['Batch'])

    #making barcodes gloablly unique before merging
    adata.obs_names = row["UniqID"] + "_" + adata.obs_names

    if config.DEV and adata.n_obs > config.DEV_MAX_CELLS:
        sc.pp.subsample(adata, n_obs=config.DEV_MAX_CELLS, random_state=config.SEED)
    return adata

def main():
    out_path = config.CHECKPOINTS / "01_qc.h5ad"
    if out_path.exists():
        print(f"[01] {out_path.name} exists, skipping. . . ")
        return

    meta = pd.read_csv(config.META_PATH)
    assert len(meta) == 35, f"expected 35 samples, got {len(meta)}"

    # deriving Batch from Admera submission plate (06 -> 1, 07 -> 2) 
    plate = meta["Admera ID"].str.split("-").str[1]     
    meta["Batch"] = plate.map({"06": "1", "07": "2"})
    assert meta["Batch"].notna().all(), \
        f"unmapped Admera plates: {sorted(plate[meta['Batch'].isna()].unique())}"
    bc = meta["Batch"].value_counts().to_dict()
    assert bc == {"1": 24, "2": 11}, f"batch split {bc} != expected 24/11"
    print(f"[01] batch from Admera plate: {bc}")

    # DEV mode: only load a few samples

    if config.DEV:
        meta = meta.head(config.DEV_N_SAMPLES)
        print(f"[01][DEV] subset to {len(meta)} samples, "
              f"<= {config.DEV_MAX_CELLS} cells each (count checks relaxed)")
 
    print(f"[01] loading {len(meta)} samples")

        
    # loading each sample and checking its count
    per_sample = []
    for _, row in meta.iterrows():
        a = load_one_sample(row)
        expected = int(row["#cells"])
        flag = "" if a.n_obs == expected else f"expected {expected}, got {a.n_obs}"
        print(f"[01] {row['UniqID']}: {a.n_obs} cells {flag}")
        per_sample.append(a)

    # concatenating all samples into one AnnData object
    adata = ad.concat(per_sample, join="outer", merge="same")
    adata.obs_names_make_unique()
    print(f"[01] concatenated: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # grounding: grand total == sum of sheet #cells (skip in DEV — deliberately subsampled)
    if not config.DEV:
        expected_total = int(meta["#cells"].sum())
        assert adata.n_obs == expected_total, \
            f"total {adata.n_obs} cells, expected {expected_total}"

    # flagging macaque mitochondrial genes
    mito_in = [g for g in config.MITO_GENES if g in adata.var_names]
    missing = [g for g in config.MITO_GENES if g not in adata.var_names]
    assert mito_in, "No mitochondrial genes foound - check the gene names in the matrix"
    if missing:
        print(f"[01] mito genes not found: {missing}")
    adata.var["mt"] = adata.var_names.isin(mito_in)
    print(f"[01] flagged {len(mito_in)} mitochondrial genes")

    # computing QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    print(f"[01] median pct_counts_mt = {adata.obs['pct_counts_mt'].median():.2f}%")

    adata.uns["stage_01"] ={
        "scanpy": sc.__version__, "anndata": ad.__version__,
        "n_samples": int(meta["UniqID"].nunique()), "n_cells": int(adata.n_obs),
        "mito_genes": len(mito_in), "missing_mito_genes": missing,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path, compression="gzip")
    print(f"[01] saved {out_path.name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

if __name__ == "__main__":
    main()