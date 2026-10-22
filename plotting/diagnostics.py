"""
diagnostics.py 

Joins the Seurat ATAC QC (TSS enrichment, nucleosome signal) onto the FULL
01_qc object BEFORE any cell cut, then produces the plots/tables that set the
joint keep-rule thresholds — especially the n_genes-vs-TSS view that reveals
the quiescent (low-RNA, intact-chromatin) population the RNA-only cut discards.

 Outputs: PNGs + printed tables.

"""
from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import config

# ---- candidate thresholds — TUNE from these plots, not final ---------------
MIN_GENES_STRICT     = 200    # current RNA-only floor
MIN_GENES_PERMISSIVE = 100    # proposed rescue floor for ATAC-sound cells
TSS_MIN              = 3.5     # ATAC "good nucleus" floor
NUC_MAX              = 2.0     # ATAC over-digestion ceiling
MITO_MAX_RESCUE      = 15 

PLOTS = Path("plots"); PLOTS.mkdir(exist_ok=True)


def main():
    adata = sc.read_h5ad(config.CHECKPOINTS / "01_qc.h5ad")

    # ---- join ATAC metrics by barcode (NO cutting) ------------------------
    atac = pd.read_csv(config.ATAC_QC_CSV).set_index("barcode")
    for c in ["TSS.enrichment", "nucleosome_signal"]:
        adata.obs[c] = atac[c].reindex(adata.obs_names).values
    obs = adata.obs

    has_atac = obs["TSS.enrichment"].notna()
    print(f"cells total: {len(obs):,} | with ATAC: {has_atac.sum():,} "
          f"({has_atac.mean()*100:.1f}%) | NaN ATAC: {(~has_atac).sum():,}")

    # ---- ATAC coverage per sample (quantifies the NaN fallback) ----------
    cov = obs.groupby("sample").agg(
        n_cells=("TSS.enrichment", "size"),
        n_atac=("TSS.enrichment", "count"))          # count() ignores NaN
    cov["pct_atac"] = (cov["n_atac"] / cov["n_cells"] * 100).round(1)
    print("\nATAC coverage per sample:")
    print(cov.sort_values("pct_atac").to_string())

    # ---- quadrant table: RNA-pass x ATAC-pass (cells WITH ATAC) ----------
    a = obs[has_atac].copy()
    rna_pass  = a["n_genes_by_counts"] >= MIN_GENES_STRICT
    atac_pass = (a["TSS.enrichment"] >= TSS_MIN) & (a["nucleosome_signal"] <= NUC_MAX)
    rescue = (~rna_pass) & atac_pass \
         & (a["n_genes_by_counts"] >= MIN_GENES_PERMISSIVE) \
         & (a["pct_counts_mt"] <= MITO_MAX_RESCUE)

    print("\nquadrants (cells with ATAC):")
    print(f"  RNA-pass & ATAC-pass : {(rna_pass & atac_pass).sum():>9,}  keep (good both)")
    print(f"  RNA-pass & ATAC-fail : {(rna_pass & ~atac_pass).sum():>9,}  good RNA, poor chromatin")
    print(f"  RNA-fail & ATAC-pass : {((~rna_pass) & atac_pass).sum():>9,}  QUIESCENT candidates")
    print(f"  RNA-fail & ATAC-fail : {((~rna_pass) & (~atac_pass)).sum():>9,}  drop (bad both)")
    print(f"  -> rescued by OR-rule (genes>={MIN_GENES_PERMISSIVE} & ATAC-pass): {rescue.sum():,}")
    print(f"     of which DEL5: {a.loc[rescue, 'sample'].str.startswith('DEL5').sum():,}")

    # ---- PLOT 1: THE view — n_genes vs TSS enrichment --------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    hb = ax.hexbin(a["n_genes_by_counts"], a["TSS.enrichment"],
                   xscale="log", gridsize=60, bins="log", cmap="viridis", mincnt=1)
    ax.axvline(MIN_GENES_STRICT, color="red", ls="--", lw=1, label=f"RNA floor {MIN_GENES_STRICT}")
    ax.axvline(MIN_GENES_PERMISSIVE, color="orange", ls="--", lw=1, label=f"rescue floor {MIN_GENES_PERMISSIVE}")
    ax.axhline(TSS_MIN, color="white", ls=":", lw=1.2, label=f"TSS min {TSS_MIN}")
    ax.set(xlabel="n_genes_by_counts (log)", ylabel="TSS enrichment",
           title="Quiescent cloud = left of red (low genes) + above white (intact chromatin)")
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(hb, label="log10(cells)")
    fig.tight_layout(); fig.savefig(PLOTS / "genes_vs_tss.png", dpi=130); plt.close(fig)

    # ---- PLOT 2: high-TSS cells — quiescent (low mito) or dying (high mito)? ----
    fig, ax = plt.subplots(figsize=(8, 6))
    hb = ax.hexbin(a["TSS.enrichment"], a["pct_counts_mt"],
                   gridsize=60, bins="log", cmap="magma", mincnt=1)
    ax.axvline(TSS_MIN, color="cyan", ls=":", lw=1.2)
    ax.set(xlabel="TSS enrichment", ylabel="pct_counts_mt",
           title="High-TSS cells should be LOW mito (quiescent, not dying)")
    fig.colorbar(hb, label="log10(cells)")
    fig.tight_layout(); fig.savefig(PLOTS / "tss_vs_mito.png", dpi=130); plt.close(fig)

    # ---- PLOT 3 & 4: per-sample ATAC distributions ----------------------
    order = a.groupby("sample", observed=True)["TSS.enrichment"].median().sort_values().index
    for metric, thr in [("TSS.enrichment", TSS_MIN), ("nucleosome_signal", NUC_MAX)]:
        fig, ax = plt.subplots(figsize=(12, 5))
        sns.violinplot(data=a, x="sample", y=metric, order=order, inner="quartile", ax=ax)
        ax.axhline(thr, color="red", ls="--", lw=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
        ax.set_title(f"{metric} per sample (red = threshold {thr})")
        fig.tight_layout(); fig.savefig(PLOTS / f"{metric}_per_sample.png", dpi=130); plt.close(fig)

    print(f"\nsaved 4 plots to {PLOTS.resolve()}/")


if __name__ == "__main__":
    main()