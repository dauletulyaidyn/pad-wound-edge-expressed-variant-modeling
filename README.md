# Machine Learning-Based Gene-Burden Modeling from scRNA-seq Expressed Variants for Wound-Edge Analysis in Peripheral Artery Disease

Wound-edge expressed-variant modeling, pathway interpretation and cellular context.

Article-specific research support repository, version 1.0.2.

**Status: historical source-artifact package; manuscript numerical results are not independently reproduced.**

The manuscript-specific numerical results have not been reproduced by this package. Archived source files and manuscript-reported claims are separated. Public GSE176415 metadata contain four unwounded-skin and three wound biological samples; the fourteen runs comprise eight UWE and six WE observations. Legacy labels disagree with this mapping. The archived CV uses run-level splits and fixed hyperparameters, not donor/GSM-grouped nested tuning. Cohort-frequency filtering uses the entire cohort. Archived PCA values, coefficients and some classifier summaries differ from the manuscript. No matched DNA or verified PAD clinical status is established by this dataset. All expressed variants and biological interpretations are exploratory. Reanalysis with reconciled labels and group-aware training-only preprocessing is required before publication claims can be considered reproducible.

## Contents

- `supplementary/`: Word supplement, supplementary methods, numbered Tables S1–S3 and source-result tables.
- `archive/`: unchanged historical matrices, metrics, coefficients, metadata and dependency definitions.
- `data/metadata_source_reconciled.csv`: source-derived run labels and GSM groups for future reanalysis.
- `scripts/legacy/`: original pipeline code retained for provenance.
- `provenance.json`: SHA-256 manifest and relative PAD source paths.
- `REPRODUCIBILITY.md`: unresolved claim-to-artifact differences.

## Verify the package

```bash
python scripts/validate_package.py
```

This uses only Python standard-library modules. Success indicates file integrity, not scientific validation.

## Reanalysis prerequisites

Use Linux/WSL with STAR, samtools, bcftools, BEDTools and Python scientific packages described in `archive/requirements-wsl.txt` and `archive/environment.yml`. The historical scripts expect a conventional repository root one level above `scripts`; here they are preserved under `scripts/legacy`, so inspect their `BASE`/`ROOT` variables and pass explicit input/output paths or reconstruct their original layout. Download source FASTQ and compatible GRCh38/GENCODE references separately. Do not run legacy labels as a confirmatory WE/UWE analysis. Preserve GSM groups during evaluation and fit recurrent-locus filtering inside training folds.

## Relationship to the original project

Derived from local PAD source materials. Original general project: https://github.com/dauletulyaidyn/somatic-like-scrna-variants-ml-pipeline. This repository has a separate article-specific name, description and support scope. It does not represent an independent cohort or a new validation study.

## Sources

- [GSE176415](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176415)
- [PRJNA736095](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA736095)

No manuscript is published here. Code/data licensing has not been newly assigned; consult original owners and source terms.

## Journal and supplied template

Target: [Herald of KBTU](https://vestnik.kbtu.edu.kz/jour/about/submissions#authorGuidelines). Article/template number: 2. Supplied template: `2. шаблон_Econ and bus paper_template.docx + 2. шаблон_Author details_in 3 languages.docx`. The earlier ETASR/MDPI appearance came from preliminary manuscript files and is not the journal destination. The corrected local manuscripts retain the supplied filled journal versions and now cite this support repository and its supplementary package.

See `JOURNAL_ALIGNMENT.md` for boundaries and submission-file handling.

## Journal attachment contents

`supplementary/` contains only the journal supplement DOCX and Tables S1-S3 (CSV). Historical output copies are kept in `archive/submission_excluded`; the editable supplement source is in `docs/`. Submission ZIP files exclude source code, repository documentation, provenance manifests and historical archives. The manuscript explicitly cites all three supplementary tables and methods.
