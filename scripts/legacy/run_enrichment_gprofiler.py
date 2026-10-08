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
    ap = argparse.ArgumentParser(description="Run pathway/GO enrichment using g:Profiler API.")
    ap.add_argument("--genes", required=True, help="Path to a text file with one gene symbol per line.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--organism", default="hsapiens", help="g:Profiler organism (default: hsapiens).")
    ap.add_argument("--top", type=int, default=20, help="Top terms to plot (default: 20).")
    args = ap.parse_args()

    genes_path = Path(args.genes)
    genes = [ln.strip() for ln in genes_path.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
    genes = list(dict.fromkeys(genes))  # stable unique
    if len(genes) < 5:
        raise SystemExit(f"Need >=5 genes for enrichment. Got: {len(genes)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = gprofiler_enrich(genes, organism=args.organism)
    if not rows:
        print("No enrichment results.")
        return 0

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "enrichment_results.csv", index=False)

    # Plot top terms by -log10(p)
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        d = df.sort_values("p_value", ascending=True).head(args.top).copy()
        d["score"] = -np.log10(d["p_value"].astype(float))
        d = d.iloc[::-1]
        plt.figure(figsize=(10, max(4, 0.25 * len(d) + 1)))
        plt.barh(d["name"], d["score"])
        plt.xlabel("-log10(p-value)")
        plt.title(f"Enrichment (top {min(args.top,len(d))})")
        plt.tight_layout()
        plt.savefig(outdir / "enrichment_top_terms.png", dpi=200)
        plt.close()
    except Exception as e:
        print(f"[warn] plotting skipped: {e}")

    print(f"Wrote: {outdir / 'enrichment_results.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
