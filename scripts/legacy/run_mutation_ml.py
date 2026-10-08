#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import gzip
import csv
import warnings

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (
    StratifiedKFold,
    RepeatedStratifiedKFold,
    cross_validate,
    permutation_test_score,
)
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif


BASE = Path(__file__).resolve().parents[1]


def read_metadata(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path)


def load_variant_gene_long(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip")


def build_gene_matrix(df: pd.DataFrame) -> pd.DataFrame:
    # Count unique variants per gene per sample.
    df = df.dropna(subset=["gene_name"])
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


def run_unsupervised(X: np.ndarray, samples: list[str], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xs = scaler.fit_transform(X)

    pca = PCA(n_components=2, random_state=0)
    Xp = pca.fit_transform(Xs)

    # Save PCA coords
    pd.DataFrame({"sample": samples, "PC1": Xp[:, 0], "PC2": Xp[:, 1]}).to_csv(
        outdir / "pca_coords.csv", index=False
    )
    with (outdir / "pca_explained_variance.txt").open("w", encoding="utf-8") as f:
        f.write(f"explained_variance_ratio: {pca.explained_variance_ratio_.tolist()}\n")

    # KMeans sweep
    rows = []
    for k in range(2, min(7, len(samples))):
        km = KMeans(n_clusters=k, n_init=20, random_state=0)
        labels = km.fit_predict(Xs)
        sil = silhouette_score(Xs, labels) if len(set(labels)) > 1 else float("nan")
        rows.append({"k": k, "silhouette": sil})
        pd.DataFrame({"sample": samples, "cluster": labels}).to_csv(
            outdir / f"kmeans_k{k}_clusters.csv", index=False
        )
    pd.DataFrame(rows).to_csv(outdir / "kmeans_silhouette.csv", index=False)


def _supervised_models(n_features: int) -> list[tuple[str, object]]:
    k = min(50, max(1, n_features))
    return [
        (
            "logreg_l2",
            Pipeline(
                steps=[
                    ("var", VarianceThreshold()),
                    ("scaler", StandardScaler(with_mean=True, with_std=True)),
                    ("clf", LogisticRegression(max_iter=5000)),
                ]
            ),
        ),
        (
            "logreg_l1",
            Pipeline(
                steps=[
                    ("var", VarianceThreshold()),
                    ("scaler", StandardScaler(with_mean=True, with_std=True)),
                    ("clf", LogisticRegression(max_iter=5000, penalty="l1", solver="liblinear")),
                ]
            ),
        ),
        (
            "linear_svc",
            Pipeline(
                steps=[
                    ("var", VarianceThreshold()),
                    ("scaler", StandardScaler(with_mean=True, with_std=True)),
                    ("clf", LinearSVC(dual="auto")),
                ]
            ),
        ),
        ("rf_500", RandomForestClassifier(n_estimators=500, random_state=0)),
        # More conservative variants (dimensionality reduction inside CV)
        (
            "logreg_pca10",
            Pipeline(
                steps=[
                    ("var", VarianceThreshold()),
                    ("scaler", StandardScaler(with_mean=True, with_std=True)),
                    ("pca", PCA(n_components=min(10, max(2, n_features)), random_state=0)),
                    ("clf", LogisticRegression(max_iter=5000)),
                ]
            ),
        ),
        (
            "logreg_kbest",
            Pipeline(
                steps=[
                    ("var", VarianceThreshold()),
                    ("scaler", StandardScaler(with_mean=True, with_std=True)),
                    ("kbest", SelectKBest(score_func=f_classif, k=k)),
                    ("clf", LogisticRegression(max_iter=5000)),
                ]
            ),
        ),
    ]


def _summarize_cv(res: dict[str, np.ndarray]) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, vals in res.items():
        if not key.startswith("test_"):
            continue
        metric = key.replace("test_", "")
        out[f"{metric}_mean"] = float(np.mean(vals))
        out[f"{metric}_std"] = float(np.std(vals))
        out[f"{metric}_q025"] = float(np.quantile(vals, 0.025))
        out[f"{metric}_q975"] = float(np.quantile(vals, 0.975))
    return out


def load_baseline_features_from_bcftools_sn(path: Path, samples: list[str]) -> tuple[np.ndarray, list[str]] | None:
    if not path.exists():
        return None
    df = pd.read_csv(path)
    if "sample" not in df.columns:
        return None

    cols = [c for c in ["number of records:", "number of SNPs:", "number of indels:"] if c in df.columns]
    if not cols:
        return None

    df = df.set_index("sample").reindex(samples)
    X = df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(dtype=float)
    return X, cols


def run_supervised(
    X: np.ndarray,
    samples: list[str],
    metadata: pd.DataFrame,
    outdir: Path,
    repeats: int,
    permutations: int,
    baseline_bcftools_sn: Path,
) -> None:
    md = metadata.set_index("sample").reindex(samples)
    if "condition" not in md.columns:
        return
    y = md["condition"].astype("category")
    if y.isna().any():
        return
    if y.nunique() < 2:
        return

    outdir.mkdir(parents=True, exist_ok=True)

    n_splits = min(5, int(y.value_counts().min()))
    if n_splits < 2:
        return
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
    rcv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=max(1, repeats), random_state=0)

    scoring = {"accuracy": "accuracy", "balanced_accuracy": "balanced_accuracy"}
    if y.nunique() == 2:
        scoring["roc_auc"] = "roc_auc"

    models = _supervised_models(n_features=X.shape[1])

    # Single CV (kept for backwards compatibility)
    rows = []
    for name, model in models:
        res = cross_validate(model, X, y, cv=cv, scoring=scoring, n_jobs=None)
        row = {"model": name, "cv": f"StratifiedKFold(n_splits={n_splits})"}
        row.update(_summarize_cv(res))
        rows.append(row)

    pd.DataFrame(rows).sort_values("balanced_accuracy_mean", ascending=False).to_csv(
        outdir / "supervised_metrics.csv", index=False
    )

    # Repeated CV (more stable estimates on small n)
    rows_rep = []
    for name, model in models:
        res = cross_validate(model, X, y, cv=rcv, scoring=scoring, n_jobs=None)
        row = {
            "model": name,
            "cv": f"RepeatedStratifiedKFold(n_splits={n_splits}, n_repeats={max(1,repeats)})",
        }
        row.update(_summarize_cv(res))
        rows_rep.append(row)

    pd.DataFrame(rows_rep).sort_values("balanced_accuracy_mean", ascending=False).to_csv(
        outdir / "supervised_metrics_repeated.csv", index=False
    )

    # Permutation test (leakage sanity check). Use single CV for speed.
    if permutations > 0:
        model = dict(models)["logreg_l2"]
        score_name = "balanced_accuracy"
        score, perm_scores, pvalue = permutation_test_score(
            model,
            X,
            y,
            scoring=score_name,
            cv=cv,
            n_permutations=permutations,
            n_jobs=None,
            random_state=0,
        )
        pd.DataFrame({"perm_score": perm_scores}).to_csv(
            outdir / "permutation_scores.csv", index=False
        )
        with (outdir / "permutation_test.txt").open("w", encoding="utf-8") as f:
            f.write(f"scoring: {score_name}\n")
            f.write(f"cv: StratifiedKFold(n_splits={n_splits})\n")
            f.write(f"n_permutations: {permutations}\n")
            f.write(f"observed_score: {float(score):.6f}\n")
            f.write(f"perm_mean: {float(np.mean(perm_scores)):.6f}\n")
            f.write(f"perm_std: {float(np.std(perm_scores)):.6f}\n")
            f.write(f"pvalue: {float(pvalue):.6f}\n")

    # Baseline sanity check using only simple per-sample variant counts (bcftools SN)
    baseline = load_baseline_features_from_bcftools_sn(baseline_bcftools_sn, samples)
    if baseline is not None:
        Xb, feat_names = baseline
        model = Pipeline(
            steps=[
                ("var", VarianceThreshold()),
                ("scaler", StandardScaler(with_mean=True, with_std=True)),
                ("clf", LogisticRegression(max_iter=5000)),
            ]
        )
        res = cross_validate(model, Xb, y, cv=rcv, scoring=scoring, n_jobs=None)
        row = {"model": "logreg_l2", "features": "+".join(feat_names)}
        row.update(_summarize_cv(res))
        pd.DataFrame([row]).to_csv(outdir / "baseline_qc_metrics.csv", index=False)
    return


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--variant-gene-long",
        default=str(BASE / "results" / "tables" / "variant_gene_long.tsv.gz"),
    )
    ap.add_argument(
        "--metadata",
        default=str(BASE / "data" / "metadata_samples.csv"),
    )
    ap.add_argument(
        "--outdir",
        default=str(BASE / "results" / "ml"),
    )
    ap.add_argument(
        "--repeats",
        type=int,
        default=30,
        help="Number of repeats for RepeatedStratifiedKFold (default: 30).",
    )
    ap.add_argument(
        "--permutations",
        type=int,
        default=200,
        help="Permutation test count (default: 200, 0 disables).",
    )
    ap.add_argument(
        "--baseline-bcftools-sn",
        default=str(BASE / "results" / "tables" / "bcftools_sn.filtered_novel_cohort.csv"),
        help="Baseline feature source (default: bcftools SN summary for cohort-filtered VCFs).",
    )
    args = ap.parse_args()

    # Keep output clean and avoid drowning runs in expected warnings on tiny datasets.
    warnings.filterwarnings(
        "ignore",
        category=UserWarning,
        module=r"sklearn\.feature_selection\._univariate_selection",
    )
    warnings.filterwarnings(
        "ignore",
        message=r"invalid value encountered in divide",
        category=RuntimeWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message=r"divide by zero encountered in divide",
        category=RuntimeWarning,
    )

    variant_gene_long = Path(args.variant_gene_long)
    outdir = Path(args.outdir)

    df = load_variant_gene_long(variant_gene_long)
    mat = build_gene_matrix(df)
    outdir.mkdir(parents=True, exist_ok=True)
    mat.to_csv(outdir / "gene_matrix.csv")

    samples = mat.index.tolist()
    X = mat.to_numpy(dtype=float)
    run_unsupervised(X, samples, outdir)

    metadata = read_metadata(Path(args.metadata))
    if metadata is not None:
        # Supervised metrics
        run_supervised(
            X,
            samples,
            metadata,
            outdir,
            repeats=args.repeats,
            permutations=args.permutations,
            baseline_bcftools_sn=Path(args.baseline_bcftools_sn),
        )

        # Feature ranking for interpretability (logreg_l2 on standardized X)
        md = metadata.set_index("sample").reindex(samples)
        if "condition" in md.columns:
            y = md["condition"].astype("category")
            if not y.isna().any() and y.nunique() == 2:
                scaler = StandardScaler(with_mean=True, with_std=True)
                Xs = scaler.fit_transform(X)
                clf = LogisticRegression(max_iter=5000)
                clf.fit(Xs, y)
                coefs = clf.coef_[0]
                genes = mat.columns.to_list()
                df_coef = pd.DataFrame({"gene": genes, "coef": coefs})
                df_coef["abs_coef"] = df_coef["coef"].abs()
                df_coef.sort_values("abs_coef", ascending=False).to_csv(
                    outdir / "logreg_top_genes.csv", index=False
                )

        # Write a short report for the paper draft.
        rep_path = outdir / "ml_report.txt"
        try:
            rep = []
            rep.append(f"samples: {len(samples)}")
            rep.append(f"features: {int(X.shape[1])}")
            rep.append(f"repeats: {int(args.repeats)}")
            rep.append(f"permutations: {int(args.permutations)}")
            sm = outdir / "supervised_metrics_repeated.csv"
            if sm.exists():
                dfm = pd.read_csv(sm).sort_values("balanced_accuracy_mean", ascending=False)
                best = dfm.iloc[0].to_dict()
                rep.append("best_model_repeated: " + str(best.get("model")))
                rep.append("best_balanced_accuracy_mean: " + f"{float(best.get('balanced_accuracy_mean', float('nan'))):.6f}")
                rep.append("best_balanced_accuracy_q025: " + f"{float(best.get('balanced_accuracy_q025', float('nan'))):.6f}")
                rep.append("best_balanced_accuracy_q975: " + f"{float(best.get('balanced_accuracy_q975', float('nan'))):.6f}")
            pt = outdir / "permutation_test.txt"
            if pt.exists():
                rep.append("--- permutation_test ---")
                rep.extend([line.rstrip('\n') for line in pt.read_text(encoding='utf-8').splitlines() if line.strip()])
            bq = outdir / "baseline_qc_metrics.csv"
            if bq.exists():
                dfb = pd.read_csv(bq)
                rep.append("--- baseline_qc_metrics ---")
                rep.append(dfb.to_string(index=False))
            rep.append(
                "note: With n=14, perfect scores can happen even without leakage; use repeated CV + permutation test and avoid strong generalization claims."
            )
            rep_path.write_text("\n".join(rep) + "\n", encoding="utf-8")
        except Exception:
            # Report generation should never fail the run.
            pass

    print(f"Wrote ML outputs to {outdir}")


if __name__ == "__main__":
    main()
