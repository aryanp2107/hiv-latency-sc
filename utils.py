"""
utils.py — shared helpers imported by stages and notebooks.
"""
from pathlib import Path
import anndata as ad


def audit(checkpoints=None):
    """Print one row per checkpoint: size on disk, shape, and provenance.

    Reads each .h5ad in BACKED mode (obs/var/uns load, the big X matrix stays
    on disk), so this is instant even at 300k cells.
    """
    import config
    cp = Path(checkpoints or config.CHECKPOINTS)
    files = sorted(cp.glob("*.h5ad"))
    if not files:
        print(f"(no checkpoints in {cp})"); return

    print(f"tracking: {cp}")
    print(f"{'file':<22}{'size':>7}  {'cells':>9} x {'genes':<7}  stage / scanpy")
    print("-" * 70)
    for f in files:
        size_mb = f.stat().st_size / 1e6
        a = ad.read_h5ad(f, backed="r")                 
        stages = [k for k in a.uns_keys() if k.startswith("stage_")]
        prov = a.uns[stages[-1]] if stages else {}
        tag = f"{stages[-1] if stages else '?':<9} scanpy {prov.get('scanpy', '?')}"
        warn = "  <-- suspiciously small" if size_mb < 1 else ""
        print(f"{f.name:<22}{size_mb:>6.0f}M  {a.n_obs:>9,} x {a.n_vars:<7,}  {tag}{warn}")
        a.file.close()