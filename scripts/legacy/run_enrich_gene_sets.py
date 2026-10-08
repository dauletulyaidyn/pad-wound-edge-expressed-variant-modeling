#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


def gprofiler_enrich(genes: list[str], organism: str = "hsapiens") -> list[dict[str, Any]]:
    import urllib.request

    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {
        "organism": organism,
        "query": genes,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"],
        "user_threshold": 0.05,
    }
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as resp:
        res = json.loads(resp.read().decode("utf-8"))
    return res.get("result", [])


def main() -> int:
    ap = argparse.ArgumentParser(description="Run g:Profiler enrichment for every *.txt gene set in a directory.")
    ap.add_argument("--gene-sets-dir", required=True, help="Directory containing *.txt (one gene per line).")
    ap.add_argument("--outdir", required=True, help="Output directory (subdir per gene set).")
    ap.add_argument("--organism", default="hsapiens", help="g:Profiler organism (default: hsapiens).")
    ap.add_argument("--top", type=int, default=20, help="Top terms to plot per gene set (default: 20).")
    args = ap.parse_args()

    gene_sets_dir = Path(args.gene_sets_dir)
    if not gene_sets_dir.exists():
        raise SystemExit(f"Missing: {gene_sets_dir}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows_index = []
    for genes_path in sorted(gene_sets_dir.glob("*.txt")):
        set_name = genes_path.stem
        genes = [ln.strip() for ln in genes_path.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
        genes = list(dict.fromkeys(genes))  # stable unique
        if len(genes) < 5:
            print(f"[skip] {set_name}: need >=5 genes (got {len(genes)})")
            continue

        out_set = outdir / set_name
        out_set.mkdir(parents=True, exist_ok=True)

        res = gprofiler_enrich(genes, organism=args.organism)
        if not res:
            print(f"[none] {set_name}: no significant results")
            continue

        df = pd.DataFrame(res)
        df.to_csv(out_set / "enrichment_results.csv", index=False)

        # plot best-effort
        try:
            import matplotlib.pyplot as plt
            import numpy as np

            d = df.sort_values("p_value", ascending=True).head(args.top).copy()
            d["score"] = -np.log10(d["p_value"].astype(float))
            d = d.iloc[::-1]
            plt.figure(figsize=(10, max(4, 0.25 * len(d) + 1)))
            plt.barh(d["name"], d["score"])
            plt.xlabel("-log10(p-value)")
            plt.title(f"Enrichment: {set_name} (top {min(args.top, len(d))})")
            plt.tight_layout()
            plt.savefig(out_set / "enrichment_top_terms.png", dpi=200)
            plt.close()
        except Exception as e:
            print(f"[warn] {set_name}: plotting skipped: {e}")

        top_row = df.sort_values("p_value", ascending=True).head(1)
        if len(top_row) == 1:
            r = top_row.iloc[0].to_dict()
            rows_index.append(
                {
                    "set": set_name,
                    "n_genes": len(genes),
                    "top_term": r.get("name"),
                    "top_source": r.get("source"),
                    "top_p_value": r.get("p_value"),
                }
            )
        print(f"[ok] {set_name}: wrote {out_set / 'enrichment_results.csv'}")

    if rows_index:
        pd.DataFrame(rows_index).sort_values("top_p_value").to_csv(outdir / "enrichment_index.csv", index=False)
        print(f"Wrote: {outdir / 'enrichment_index.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
