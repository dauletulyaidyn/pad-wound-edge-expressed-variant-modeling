# Revision 2026-11-21 (manuscript rev2)

The revised manuscript is titled *Sample identity as a source of data leakage in machine learning on scRNA-seq expressed-variant features: a genotype audit and exact permutation test of grouped validation*. This folder holds its supplementary materials, figures and analysis scripts.

| Path | What it is |
|---|---|
| `supplementary/KBTU_rev2_Supplementary_Materials.docx` | Notes S1–S2, Tables S1–S5 |
| `supplementary/KBTU_S*.csv` | the same tables as CSV (S1–S6) |
| `figures/` | Figures 1–5 (600 dpi) |
| `scripts/` | light_00–05 and the shared model definitions |

The models (`light_03`, `light_05`) must run with scikit-learn 1.9.0 and XGBoost 3.2.0. The full pipeline is in https://github.com/dauletulyaidyn/somatic-like-scrna-variants-ml-pipeline.

## Key results

- Runs cluster by GEO sample, not by condition: 1-NN identifies the sample of 14/14 runs under all four similarity measures.
- Sc_Skin1 (GSM5364333) and Sc_W2 (GSM5364338) are genotype-identical (concordance 0.997 at depth ≥ 20), so the 7 samples come from 6 people.
- No classifier in any grouped design is significant on the exact test over all 35 labelings (minimum p = 0.200).
- Chance balanced accuracy under leave-one-group-out validation is below 0.5 (null mean 0.28–0.48).
