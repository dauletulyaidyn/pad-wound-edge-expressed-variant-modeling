#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run a minimal Scanpy pipeline (QC + HVG + PCA/UMAP + Leiden clustering + marker tables)."
    )
    ap.add_argument("--starsolo-dir", default="starsolo", help="Directory containing per-sample STARsolo outputs.")
    ap.add_argument("--metadata", default="data/metadata_samples.csv", help="CSV with columns sample,condition.")
    ap.add_argument("--outdir", default="results/expression", help="Output directory.")
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-genes", type=int, default=9000)
    ap.add_argument("--hvg", type=int, default=2000)
    ap.add_argument("--neighbors", type=int, default=15)
    ap.add_argument("--pcs", type=int, default=30)
    ap.add_argument("--leiden-resolution", type=float, default=1.0)
    args = ap.parse_args()

    try:
        import scanpy as sc
        import numpy as np
    except Exception as e:
        raise SystemExit(
            "Missing Scanpy.\n"
            "Recommended (WSL/Ubuntu):\n"
            "  bash scripts/install_micromamba_env.sh\n"
            "  micromamba activate scrna-variant-ml\n"
            "Then run:\n"
            "  python scripts/run_expression_scanpy.py\n"
            f"\nOriginal import error: {e}"
        )

    root = Path(__file__).resolve().parents[1]
    starsolo = (root / args.starsolo_dir).resolve()
    md = pd.read_csv(root / args.metadata)
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    samples = md["sample"].astype(str).tolist()
    sample_to_cond = dict(zip(md["sample"].astype(str), md["condition"].astype(str)))

    adatas = []
    for s in samples:
        mtx_dir = starsolo / s / "Solo.out" / "Gene" / "filtered"
        if not mtx_dir.exists():
            print(f"[skip] missing STARsolo matrix: {mtx_dir}")
            continue
        a = sc.read_10x_mtx(str(mtx_dir), var_names="gene_symbols", make_unique=True)
        a.obs["sample_id"] = s
        a.obs["condition"] = sample_to_cond.get(s, "unknown")
        adatas.append(a)

    if not adatas:
        raise SystemExit("No STARsolo matrices found. Run alignment first (pipeline.py star).")

    adata = adatas[0].concatenate(*adatas[1:], batch_key="batch", batch_categories=None, index_unique=None)

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    adata = adata[adata.obs["n_genes_by_counts"] <= args.max_genes].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.hvg,
        flavor="seurat_v3",
        batch_key="sample_id" if "sample_id" in adata.obs.columns else None,
    )
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=args.neighbors, n_pcs=args.pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.leiden_resolution, key_added="leiden_clusters")
    adata.obs["final_clusters"] = adata.obs["leiden_clusters"].astype(str)

    sc.tl.rank_genes_groups(adata, groupby="final_clusters", method="wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(outdir / "markers_by_cluster.csv", index=False)

    if adata.obs["condition"].nunique() > 1:
        sc.tl.rank_genes_groups(adata, groupby="condition", method="wilcoxon")
        deg = sc.get.rank_genes_groups_df(adata, group=None)
        deg.to_csv(outdir / "markers_by_condition.csv", index=False)

    sc.pl.umap(
        adata,
        color=["final_clusters", "sample_id", "condition", "n_genes_by_counts"],
        wspace=0.4,
        show=False,
    )
    import matplotlib.pyplot as plt

    plt.savefig(outdir / "umap_overview.png", dpi=200)
    plt.close()

    metrics = {
        "n_cells": int(adata.n_obs),
        "n_features_hvg": int(adata.n_vars),
        "n_samples": int(adata.obs["sample_id"].nunique()) if "sample_id" in adata.obs.columns else 0,
        "n_clusters": int(adata.obs["final_clusters"].nunique()),
        "min_genes": int(args.min_genes),
        "max_genes": int(args.max_genes),
        "hvg": int(args.hvg),
        "neighbors": int(args.neighbors),
        "pcs": int(args.pcs),
        "leiden_resolution": float(args.leiden_resolution),
    }
    pd.DataFrame([metrics]).to_csv(outdir / "metrics.csv", index=False)

    adata.write(outdir / "clustered_adata.h5ad")
    print(f"Wrote: {outdir / 'clustered_adata.h5ad'}")
    print(f"Wrote: {outdir / 'umap_overview.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

