"""
05_atac_quant.py — build the cell x consensus-peak ATAC matrix for all 26 kept samples.

reads:  per-sample atac_fragments.tsv.gz  +  consensus_peaks.bed
        02_filtered.h5ad (for the kept-cell barcodes, the alignment anchor)
writes: {checkpoints}/05_atac_peaks.h5ad   (cells x 231,794 peaks, raw counts)

Heavy step (reads 26 gigabyte fragment files) -> run in tmux.
Barcodes are prefixed UniqID_ to match the RNA object so the modalities align for WNN.
"""
import snapatac2 as snap
import scanpy as sc, anndata as ad
import pandas as pd
from pathlib import Path
import config

RAW       = Path("/scratch/maid/DP2_PBMC/PBMC_rawData")
CONSENSUS = str(config.CHECKPOINTS / "consensus_peaks.bed")
FAI       = "/scratch/references/refdata-cellranger-arc-Mmul_10-2023-A/fasta/genome.fa.fai"


def main():
    out_path = config.CHECKPOINTS / "05_atac_peaks.h5ad"
    if out_path.exists():
        print(f"[05] {out_path.name} exists, skipping."); return

    # macaque chrom sizes from the reference .fai (same ref the data was aligned to)
    fai = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1],
                      names=["chrom", "length"], dtype={"chrom": str})
    chrom_sizes = dict(zip(fai["chrom"], fai["length"]))

    # kept barcodes per sample, from the RNA object (alignment anchor)
    rna = sc.read_h5ad(config.CHECKPOINTS / "02_filtered.h5ad", backed="r")
    samples = sorted(rna.obs["sample"].unique())
    kept_by_sample = {s: set(rna.obs_names[rna.obs["sample"] == s]) for s in samples}
    print(f"[05] {len(samples)} kept samples; consensus = {CONSENSUS}")

    per_sample, match_log = [], {}
    for i, s in enumerate(samples, 1):
        a = snap.pp.import_fragments(str(RAW / s / "atac_fragments.tsv.gz"),
                                     chrom_sizes=chrom_sizes, sorted_by_barcode=False)
        a.obs_names = [s + "_" + bc for bc in a.obs_names]
        keep = [bc for bc in a.obs_names if bc in kept_by_sample[s]]
        match_log[s] = (len(keep), len(kept_by_sample[s]))
        print(f"[05] ({i}/{len(samples)}) {s}: imported {a.n_obs:,} | "
              f"matched {len(keep):,}/{len(kept_by_sample[s]):,}", flush=True)
        a = a[keep].copy()
        per_sample.append(snap.pp.make_peak_matrix(a, peak_file=CONSENSUS))

    atac = ad.concat(per_sample)
    print(f"[05] ATAC matrix: {atac.n_obs:,} cells x {atac.n_vars:,} peaks")

    # grounding: total matched cells should be close to RNA cell count
    total_matched = sum(m for m, _ in match_log.values())
    total_kept    = sum(t for _, t in match_log.values())
    print(f"[05] matched {total_matched:,} / {total_kept:,} kept RNA cells "
          f"({total_matched/total_kept*100:.1f}%)")

    atac.uns["stage_05"] = {
        "snapatac2": snap.__version__,
        "n_consensus_peaks": int(atac.n_vars),
        "n_cells": int(atac.n_obs),
        "reference": "refdata-cellranger-arc-Mmul_10-2023-A",
        "peak_method": "pyranges overlap-merge, main chromosomes only",
        "total_matched_cells": int(total_matched),
        "total_kept_rna_cells": int(total_kept),
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    atac.write_h5ad(out_path, compression="gzip")
    print(f"[05] wrote {out_path.name}: {atac.n_obs:,} x {atac.n_vars:,}")


if __name__ == "__main__":
    main()