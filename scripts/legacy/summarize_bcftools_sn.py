#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import os
from pathlib import Path


def parse_sn_file(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("SN"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            key = parts[2].strip()
            val = parts[3].strip()
            out[key] = val
    return out


def infer_sample_and_kind(filename: str) -> tuple[str, str]:
    base = os.path.basename(filename)
    sample = base.split(".", 1)[0]
    kind = "unknown"
    if ".filtered.novel.cohort.vcf.stats.summary.txt" in base:
        kind = "filtered_novel_cohort"
    elif ".filtered.novel.vcf.stats.summary.txt" in base:
        kind = "filtered_novel"
    elif ".filtered.stats.summary.txt" in base:
        kind = "filtered"
    return sample, kind


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize bcftools stats SN lines into a CSV table.")
    ap.add_argument("--glob", required=True, help="Glob for *.stats.summary.txt files.")
    ap.add_argument("--out", required=True, help="Output CSV path.")
    args = ap.parse_args()

    paths = sorted(glob.glob(args.glob))
    if not paths:
        raise SystemExit(f"No files matched: {args.glob}")

    rows: list[dict[str, str]] = []
    keys: set[str] = set()
    for p in paths:
        sample, kind = infer_sample_and_kind(p)
        sn = parse_sn_file(Path(p))
        row = {"sample": sample, "kind": kind, **sn}
        rows.append(row)
        keys.update(sn.keys())

    header = ["sample", "kind"] + sorted(keys)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

