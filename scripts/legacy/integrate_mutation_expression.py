#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_variant_gene_long(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip")


def build_gene_burden(df: pd.DataFrame) -> pd.DataFrame:
    df = df.dropna(subset=["gene_name", "sample", "chrom", "pos", "ref", "alt"]).copy()
    df["variant_key"] = (
        df["chrom"].astype(str)
        + ":"
        + df["pos"].astype(str)
        + ":"
        + df["ref"].astype(str)
        + ">"
        + df["alt"].astype(str)
    )
    dedup = df.drop_duplicates(subset=["sample", "gene_name", "variant_key"])
    mat = (
        dedup.groupby(["sample", "gene_name"])["variant_key"]
        .count()
        .unstack(fill_value=0)
        .sort_index()
    )
    return mat


def write_gene_set(path: Path, genes: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(genes) + ("\n" if genes else ""), encoding="utf-8")


def top_markers(markers: pd.DataFrame, group_col: str, gene_col: str, top_n: int = 200) -> dict[str, list[str]]:
    score_col = "scores" if "scores" in markers.columns else None
    lfc_col = "logfoldchanges" if "logfoldchanges" in markers.columns else None
    sort_col = score_col or lfc_col
    out: dict[str, list[str]] = {}
    for g, sub in markers.groupby(group_col):
        s = sub.copy()
        if sort_col:
            s[sort_col] = pd.to_numeric(s[sort_col], errors="coerce")
            s = s.sort_values(sort_col, ascending=False)
        genes = s[gene_col].dropna().astype(str).tolist()
        out[str(g)] = list(dict.fromkeys(genes))[:top_n]
    return out


def find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Integrate mutation-derived gene sets with expression marker genes (clusters/conditions). "
            "Outputs are hypothesis-generating (no medical causal claims)."
        )
    )
    ap.add_argument("--variant-gene-long", required=True, help="Path to variant_gene_long*.tsv.gz (cohort preferred).")
    ap.add_argument("--metadata", required=True, help="CSV with columns sample,condition.")
    ap.add_argument("--expression-dir", default="results/expression", help="Directory containing markers CSVs.")
    ap.add_argument("--outdir", default="results/integration", help="Output directory.")
    ap.add_argument("--top-genes", type=int, default=200, help="Gene set size (default: 200).")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    vgl = (root / args.variant_gene_long).resolve()
    md_path = (root / args.metadata).resolve()
    expr_dir = (root / args.expression_dir).resolve()
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    gene_sets_dir = outdir / "gene_sets"
    gene_sets_dir.mkdir(parents=True, exist_ok=True)

    if not vgl.exists():
        raise SystemExit(f"Missing: {vgl}")
    if not md_path.exists():
        raise SystemExit(f"Missing: {md_path}")

    md = pd.read_csv(md_path)
    if not {"sample", "condition"}.issubset(set(md.columns)):
        raise SystemExit("metadata must contain columns: sample,condition")

    df = load_variant_gene_long(vgl)
    burden = build_gene_burden(df)
    burden.to_csv(outdir / "gene_burden_per_sample.csv.gz", index=True, compression="gzip")

    # Align metadata with burden samples
    md = md[["sample", "condition"]].copy()
    md["sample"] = md["sample"].astype(str)
    md["condition"] = md["condition"].astype(str)
    md = md[md["sample"].isin(burden.index)].copy()

    # Gene sets
    total = burden.sum(axis=0).sort_values(ascending=False)
    top_all = total.head(int(args.top_genes)).index.astype(str).tolist()
    write_gene_set(gene_sets_dir / "top_mutated_genes.txt", top_all)

    gene_sets_summary = [
        {
            "set": "top_mutated_genes",
            "n_genes": len(top_all),
            "description": f"Top {len(top_all)} genes by total unique expressed-variant count across samples.",
            "file": str((gene_sets_dir / "top_mutated_genes.txt").relative_to(outdir)),
        }
    ]

    # Group-differential burden gene sets (WE vs UWE or any two labels)
    conds = sorted(md["condition"].unique().tolist())
    if len(conds) == 2:
        c0, c1 = conds
        s0 = md.loc[md["condition"] == c0, "sample"].tolist()
        s1 = md.loc[md["condition"] == c1, "sample"].tolist()
        m0 = burden.loc[s0].mean(axis=0)
        m1 = burden.loc[s1].mean(axis=0)
        delta = (m1 - m0).sort_values(ascending=False)

        up1 = delta.head(int(args.top_genes)).index.astype(str).tolist()
        up0 = delta.tail(int(args.top_genes)).index.astype(str).tolist()[::-1]

        write_gene_set(gene_sets_dir / f"{c1}_higher_mutation_burden.txt", up1)
        write_gene_set(gene_sets_dir / f"{c0}_higher_mutation_burden.txt", up0)

        gene_sets_summary.extend(
            [
                {
                    "set": f"{c1}_higher_mutation_burden",
                    "n_genes": len(up1),
                    "description": f"Top {len(up1)} genes by mean burden difference ({c1} - {c0}).",
                    "file": str((gene_sets_dir / f"{c1}_higher_mutation_burden.txt").relative_to(outdir)),
                },
                {
                    "set": f"{c0}_higher_mutation_burden",
                    "n_genes": len(up0),
                    "description": f"Top {len(up0)} genes by mean burden difference ({c0} - {c1}).",
                    "file": str((gene_sets_dir / f"{c0}_higher_mutation_burden.txt").relative_to(outdir)),
                },
            ]
        )

        delta_df = (
            pd.DataFrame(
                {
                    "gene": delta.index.astype(str),
                    "mean_" + c0: m0.reindex(delta.index).to_numpy(),
                    "mean_" + c1: m1.reindex(delta.index).to_numpy(),
                    "delta_" + c1 + "_minus_" + c0: delta.to_numpy(),
                }
            )
            .reset_index(drop=True)
            .sort_values("delta_" + c1 + "_minus_" + c0, ascending=False)
        )
        delta_df.to_csv(outdir / "mutation_burden_delta_by_condition.csv", index=False)

    pd.DataFrame(gene_sets_summary).to_csv(outdir / "gene_sets_summary.csv", index=False)

    # Overlap with expression markers (best-effort)
    overlaps = []
    gene_sets = {p.stem: [ln.strip() for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()] for p in gene_sets_dir.glob("*.txt")}

    markers_cluster = expr_dir / "markers_by_cluster.csv"
    if markers_cluster.exists():
        m = pd.read_csv(markers_cluster)
        group_col = find_col(m, ["group", "clusters", "cluster", "final_clusters"])
        gene_col = find_col(m, ["names", "gene", "gene_name"])
        if group_col and gene_col:
            top_by_cluster = top_markers(m, group_col=group_col, gene_col=gene_col, top_n=500)
            for cluster, genes in top_by_cluster.items():
                gset = set(genes)
                for set_name, set_genes in gene_sets.items():
                    overlaps.append(
                        {
                            "source": "markers_by_cluster",
                            "group": cluster,
                            "set": set_name,
                            "n_markers_top500": len(gset),
                            "n_set": len(set_genes),
                            "n_overlap": len(gset.intersection(set(set_genes))),
                        }
                    )

    markers_cond = expr_dir / "markers_by_condition.csv"
    if markers_cond.exists():
        m = pd.read_csv(markers_cond)
        group_col = find_col(m, ["group", "condition"])
        gene_col = find_col(m, ["names", "gene", "gene_name"])
        if group_col and gene_col:
            top_by_cond = top_markers(m, group_col=group_col, gene_col=gene_col, top_n=500)
            for cond, genes in top_by_cond.items():
                gset = set(genes)
                for set_name, set_genes in gene_sets.items():
                    overlaps.append(
                        {
                            "source": "markers_by_condition",
                            "group": cond,
                            "set": set_name,
                            "n_markers_top500": len(gset),
                            "n_set": len(set_genes),
                            "n_overlap": len(gset.intersection(set(set_genes))),
                        }
                    )

    if overlaps:
        pd.DataFrame(overlaps).sort_values(["source", "group", "set"]).to_csv(outdir / "marker_overlap.csv", index=False)

    # Minimal report
    lines = [
        "# Integration outputs",
        "",
        f"- Variant-gene input: `{vgl}`",
        f"- Expression dir: `{expr_dir}`",
        f"- Gene sets: `{gene_sets_dir}`",
        "",
        "Notes:",
        "- These results are hypothesis-generating links between mutation-derived gene sets and expression markers.",
        "- Without matched normal + population filtering + consequence annotation, 'somatic-like' is not confirmed somatic.",
    ]
    (outdir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote: {outdir / 'gene_sets_summary.csv'}")
    if (outdir / "marker_overlap.csv").exists():
        print(f"Wrote: {outdir / 'marker_overlap.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
