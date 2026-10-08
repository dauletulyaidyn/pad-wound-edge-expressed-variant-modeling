# Supplementary Methods S1. Data provenance and analysis units

This supplement supports the manuscript: Machine Learning-Based Gene-Burden Modeling from scRNA-seq Expressed Variants for Wound-Edge Analysis in Peripheral Artery Disease. The source data are PRJNA736095 / GSE176415. Table S1 provides run identifiers, GEO biological sample identifiers, source-derived conditions and archived labels. Biological samples and sequencing runs are different units. The original PAD project remains the provenance source; this repository packages the article-specific support materials.

# Supplementary Methods S2. Expressed-variant and gene-burden workflow

The archived code covers splice-aware STARsolo alignment, bcftools expressed-variant discovery, cohort-common-position exclusion, GENCODE exon overlap and unique sample-by-gene variant counts. Manuscript acceptance thresholds are DP >= 10, QUAL >= 30, alternate depth >= 3 and VAF >= 0.10. Exact command expressions and annotation defaults must be checked in the archived scripts before reuse. RNA calls are not confirmed somatic or germline variants. FASTQ/BAM/VCF/reference inputs are not bundled.

# Supplementary Methods S3. Classification and evaluation

The archived Python code compares linear SVC, L1/L2 logistic regression, PCA10 and KBest logistic regression, and a 500-tree random forest. Scaling and variance filtering are in sklearn pipelines. Saved repeated CV summaries are provided for baseline and cohort feature sets. The archived permutation test concerns logistic regression under single stratified CV; it does not establish model-specific permutation results for all classifiers. Repeated run-level resampling does not estimate independent-patient performance.

# Supplementary Note S1. Reproducibility limitations

The manuscript-specific numerical results have not been reproduced by this package. Archived source files and manuscript-reported claims are separated. Public GSE176415 metadata contain four unwounded-skin and three wound biological samples; the fourteen runs comprise eight UWE and six WE observations. Legacy labels disagree with this mapping. The archived CV uses run-level splits and fixed hyperparameters, not donor/GSM-grouped nested tuning. Cohort-frequency filtering uses the entire cohort. Archived PCA values, coefficients and some classifier summaries differ from the manuscript. No matched DNA or verified PAD clinical status is established by this dataset. All expressed variants and biological interpretations are exploratory. Reanalysis with reconciled labels and group-aware training-only preprocessing is required before publication claims can be considered reproducible.

# Table S1. Run metadata reconciliation

Machine-readable table: supplementary/tables/Table_S1_run_metadata_reconciliation.csv. Each row records the run, GSM, BioSample, source tissue label, legacy analysis label and agreement status. Source-derived labels are also provided in data/metadata_source_reconciled.csv; these are not claimed to be the labels used for archived results.

# Table S2. Manuscript-reported numerical claims

Machine-readable table: supplementary/tables/Table_S2_manuscript_reported_claims.csv. Values and text are transcribed from the supplied manuscript and are explicitly not new computational evidence.

# Table S3. Claim-to-artifact reconciliation

Machine-readable table: supplementary/tables/Table_S3_claim_artifact_reconciliation.csv. This table records differences in replication, classifier metrics, PCA, coefficient provenance and validation design. Archive files retain their original values; no missing predictions or significant tests are synthesized.

# Supplementary Data S1. Archived computational outputs

The archive/baseline and archive/cohort folders contain saved gene-burden matrices, classifier summaries, PCA coordinates and variance, logistic-regression coefficients, and permutation outputs. supplementary/tables contains local and cohort-filtered variant summaries and historical publication-package summaries. These outputs are historical artifacts with unresolved label provenance.

# Supplementary Data S2. Pathways and cellular context

supplementary/enrichment contains complete saved enrichment CSV files. supplementary/tables includes gene-set sizes, cluster-marker lists, broad cell annotations and historical logistic-regression rankings. Gene expression and cluster overlap provide contextual evidence; they do not assign run-level variants to individual cells. Cellular counts and effect directions reported in the manuscript remain subject to the claim-to-artifact reconciliation in Table S3.

# Supplementary Code S1. Source code and integrity verification

Run python scripts/validate_package.py from the repository root to verify the SHA-256 manifest, code syntax, and metadata structure. The validator checks packaging integrity; it does not validate classifier performance. Archived scripts in scripts/legacy require Linux/WSL and their original layout or explicit input/output paths. They are supplied for inspection and reconstruction, not as a tested one-command manuscript reproduction.

# Data and code availability

Article-specific repository: https://github.com/dauletulyaidyn/pad-wound-edge-expressed-variant-modeling. Source sequencing accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176415. The previous general project is https://github.com/dauletulyaidyn/somatic-like-scrna-variants-ml-pipeline. This repository has a different article-specific name and scope. SHA-256 provenance is recorded in provenance.json. No journal publication, DOI or scientific replication is asserted.
