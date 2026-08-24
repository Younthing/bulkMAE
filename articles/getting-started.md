# A stateless bulk transcriptomics workflow with bulkMAE

Most backend-specific sections remain an unevaluated cookbook because
the engines are optional and several resources are online or
cache-backed. The construction and local QC chunks are deliberately
executable during vignette builds; the package test suite executes the
complete offline public-API chains.

## 1. Data contract

`bulkMAE` uses a `MultiAssayExperiment` as the public input because it
already defines the three relationships needed by bulk studies:

1.  `experiments()` stores `SummarizedExperiment` leaves;
2.  `colData()` stores primary subject/sample metadata;
3.  `sampleMap()` relates assay columns to primary samples, including
    replicated measurements or assay-specific names.

The wrappers call `getWithColData()` before analysis. Consequently, a
backend sees an ordinary assay and metadata table whose rows and columns
are aligned. No result is written back to the MAE.

Public functions use `<family>_<method>()` or `<family>_<operation>()`
names. Thus `de_`, `enrich_`, `score_`, `activity_`, `coexpr_`, and
`surv_` each form one discoverable autocomplete family. The complete
0.3-to-0.4 rename table is installed as
`inst/guides/naming-migration-zh.md`.

## 2. Construct an MAE

``` r

library(bulkMAE)

mae <- mae_simulate(n_features = 100, n_samples = 12, seed = 1)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'

# The same public constructor accepts a user's matrix and aligned metadata.
from_matrix <- mae_from_matrix(
  expression = mae_pull_assay(mae, "rna", "counts"),
  samples = mae_samples(mae, "rna"),
  row_data = mae_feature_data(mae, "rna"),
  experiment = "rna",
  assay = "counts"
)

mae_validate(mae, "rna")
mae_assays(mae, "rna")
#> [1] "counts"         "log_expression" "tpm"
```

For Salmon, kallisto, or RSEM estimates, use
[`import_tximport()`](https://younthing.github.io/bulkMAE/reference/import_tximport.md)
to create the leaf. Do not round abundance-scaled values and pass them
to DESeq2 unless the DESeq2/tximport documentation explicitly supports
that construction.

## 3. Quality control

Use a transformed, approximately homoscedastic assay for distances and
sample correlations. Use raw counts for library depth and expression
filtering.

``` r

metrics <- qc_library(mae, "rna", assay = "counts")
pca_local <- reduce_pca(mae, "rna", assay = "log_expression", top_n = 50)

head(metrics)
#>     sample library_size detected_features zero_fraction
#> 1 sample01         9648               100             0
#> 2 sample02        22341               100             0
#> 3 sample03        39286               100             0
#> 4 sample04        13595               100             0
#> 5 sample05        15097               100             0
#> 6 sample06        53445               100             0
head(pca_local$x[, seq_len(2), drop = FALSE])
#>                 PC1        PC2
#> sample01 -5.7251950  0.3683962
#> sample02 -1.2488415 -3.4133992
#> sample03  7.5177991  6.0149035
#> sample04 -6.9218010  0.6164581
#> sample05  0.1968213  2.1009616
#> sample06  6.5691377 -5.5654282
```

``` r

keep <- filter_expr(
  mae,
  "rna",
  group = "condition"
)

qc_library(mae, "rna")

vst <- transform_vst(mae, "rna", blind = TRUE)
vst_mae <- mae_add_assay(mae, "rna", value = vst, name = "vst")

pca <- reduce_pca(vst_mae, "rna", "vst", top_n = 50)
qc_outliers(pca)
```

[`qc_outliers()`](https://younthing.github.io/bulkMAE/reference/qc_outliers.md)
is only a review aid. It must not be used as an automatic sample-removal
rule. Investigate library metrics, metadata, and experimental records
before excluding a sample or relabelling it.

## 4. Differential expression

### DESeq2

``` r

dds <- de_deseq2(
  mae,
  "rna",
  design = ~ batch + condition
)

de <- de_deseq2_results(
  dds,
  contrast = c("condition", "treated", "control")
)
```

### edgeR quasi-likelihood

``` r

edge_test <- de_edger(
  mae,
  "rna",
  formula = ~ batch + condition,
  coef = "conditiontreated"
)
```

### limma-voom

``` r

design <- de_design(mae, "rna", ~ batch + condition)
contrast <- de_contrast(
  mae,
  "rna",
  ~ batch + condition,
  contrasts = "conditiontreated"
)

limma_fit <- de_voom(
  mae,
  "rna",
  formula = ~ batch + condition,
  contrasts = contrast
)
```

### Repeated measurements

``` r

dream_fit <- de_dream(
  mae,
  "rna",
  formula = ~ condition + (1 | subject_id)
)
```

Time, dose, interaction, and paired designs are model-formula problems.
The wrapper does not convert strings such as “paired” into hidden design
matrices; write the intended formula and inspect its coefficients.

## 5. Batch and unwanted variation

ComBat and ComBat-seq produce corrected matrices for visualization or
methods that explicitly need corrected input. For differential
expression, modelling a known batch in the design is usually more
transparent than replacing counts.

``` r

combat <- adjust_combat(
  vst_mae,
  "rna",
  assay = "vst",
  batch = "batch",
  preserve = ~ condition
)

sv <- adjust_sva(
  vst_mae,
  "rna",
  assay = "vst",
  full = ~ batch + condition,
  null = ~ batch
)

sv_columns <- adjust_covariates(sv)
adjusted_mae <- mae_add_sample_data(vst_mae, "rna", sv_columns)
```

Surrogate variables and RUV factors are added explicitly to a returned
MAE copy before they are named in a downstream formula. `bulkMAE` does
not silently refit a model or invent control genes.

## 6. Functional interpretation

``` r

ranks <- de_ranks(de)
selected <- de_selected(de, fdr = 0.05, min_abs_effect = 1)

# A fully offline synthetic set for this simulated MAE.
ranked_ids <- names(ranks)
synthetic_sets <- gene_sets_prepare(list(
  leading = head(ranked_ids, 25),
  trailing = tail(ranked_ids, 25)
))
fgsea_result <- enrich_fgsea(ranks, synthetic_sets)

# For a real human dataset, this online/cache-backed helper retrieves MSigDB.
hallmark_sets <- gene_sets_msigdb(
  species = "human",
  collection = "H",
  id_type = "gene_symbol"
)

# The following mapping is likewise for real Ensembl-keyed results.
org_db <- annotation_orgdb("human")
id_mapping <- annotate_ids(
  names(selected),
  database = org_db,
  from = "ENSEMBL",
  to = "ENTREZID",
  multi_values = "CharacterList"
)
selected_entrez <- annotate_rekey(
  names(selected)[selected],
  id_mapping,
  source = "source_id",
  target = "target_id",
  duplicates = "first"
)
tested_entrez <- annotate_rekey(
  names(selected),
  id_mapping,
  source = "source_id",
  target = "target_id",
  duplicates = "first"
)
go_result <- enrich_go(
  genes = selected_entrez,
  method = "ora",
  org_db = org_db,
  universe = tested_entrez
)

feature_data <- mae_feature_data(mae, "rna")
tested_gene_lengths <- setNames(
  feature_data$gene_length,
  rownames(feature_data)
)
goseq_result <- enrich_goseq(
  selected = selected,
  genome = "hg38",
  id = "ensGene",
  bias = tested_gene_lengths,
  plot_fit = FALSE
)

gsva_scores <- score_gsva(
  vst_mae,
  "rna",
  gene_sets = synthetic_sets,
  assay = "vst"
)
pathway_mae <- mae_add_experiment(
  vst_mae,
  value = gsva_scores,
  name = "pathway",
  assay_name = "score",
  source_experiment = "rna"
)

pathway_activity <- activity_progeny(vst_mae, "rna", assay = "vst")
tf_activity <- activity_tf(vst_mae, "rna", assay = "vst")
```

The gene universe for ORA should be the set of genes that could have
been selected after filtering, not every gene in an annotation database.

### Choosing among gene-set tests

Use ORA for a selected gene list, preranked GSEA for a genome-wide
statistic, CAMERA/FRY/mroast when the expression matrix and design
should remain in the test, and GSVA/ssGSEA/singscore for per-sample
scores.

``` r

camera_result <- enrich_camera(
  vst_mae,
  "rna",
  gene_sets = synthetic_sets,
  formula = ~ batch + condition,
  contrast = contrast[, 1],
  assay = "vst"
)

reactome_result <- enrich_reactome(
  ranks = ranks,
  method = "gsea",
  organism = "human"
)
sample_scores <- score_singscore(
  vst_mae,
  "rna",
  up_set = interferon_up,
  down_set = interferon_down,
  assay = "vst"
)
```

### decoupleR is a multi-method engine

[`activity_decouple()`](https://younthing.github.io/bulkMAE/reference/activity_decouple.md)
is not limited to transcription factors. The official engine also
exposes enrichment-style algorithms and network-aware activity methods.
The network remains an explicit input so its direction and weights are
visible.

``` r

activity_methods()

resource <- activity_resource("PROGENy", organism = "human")
activities <- activity_decouple(
  vst_mae,
  "rna",
  network = resource,
  assay = "vst",
  statistics = c("mlm", "ulm", "wsum"),
  method_args = list(NULL, NULL, NULL),
  source = "pathway",
  target = "genesymbol",
  mor = "weight",
  consensus = TRUE
)
```

Use `statistics = c("aucell", "fgsea", "gsva", "ora")` only when those
statistics match the question and input scale. Use `mlm`, `ulm`,
`viper`, `wmean`, or `wsum` when regulator-target signs or weights are
part of the scientific model. Confirm current names with
[`activity_methods()`](https://younthing.github.io/bulkMAE/reference/activity_methods.md)
because the backend, rather than bulkMAE, owns the method registry.

## 7. Networks and subtypes

``` r

nmf_fit <- cluster_nmf(mae, "rna", assay = "log_expression", rank = 2:6)
nmf_classes <- cluster_nmf_classes(nmf_fit, what = "samples")

consensus_fit <- cluster_consensus(
  vst_mae,
  "rna",
  assay = "vst",
  max_k = 6,
  repetitions = 1000,
  plot = NULL
)
consensus_diagnostics <- cluster_consensus_diagnostics(
  consensus_fit,
  plot = NULL
)
classes_k3 <- cluster_consensus_classes(consensus_fit, k = 3)

wgcna_fit <- coexpr_wgcna(
  vst_mae,
  "rna",
  assay = "vst",
  power = 6,
  top_n = 5000
)
module_colors <- coexpr_modules(wgcna_fit)

genie_weights <- network_genie3(
  vst_mae,
  "rna",
  assay = "vst",
  regulators = transcription_factors
)
genie_edges <- network_genie3_links(genie_weights, top = 1000)
```

Choose WGCNA’s soft threshold from the study data; the example value is
not a universal default. Likewise, decide the number of molecular
subtypes using stability and biological/clinical validation rather than
a single index.
[`coexpr_wgcna()`](https://younthing.github.io/bulkMAE/reference/coexpr_wgcna.md)
follows `goodSamplesGenes()`: unusable samples and features are removed
automatically, a warning lists the change, and the returned native WGCNA
result retains the quality object and removed names for audit.

### Differential co-expression (DCA)

In this package DCA means differential co-expression analysis. It asks
whether gene-gene correlations differ between two groups; it is distinct
from differential expression.

``` r

dca <- coexpr_differential(
  vst_mae,
  "rna",
  group = "condition",
  contrast = c("control", "treated"),
  assay = "vst",
  features = highly_variable_genes
)
```

Correlation testing scales roughly with the square of the feature count.
Apply a biologically defensible, outcome-independent feature filter
before DCA and report all correlation and FDR thresholds.

### Time courses and differential usage

[`de_masigpro()`](https://younthing.github.io/bulkMAE/reference/de_masigpro.md)
is the specialized path for multi-group polynomial time-course or
dose-response models. Its `edesign` must follow maSigPro’s documented
`Time`, `Replicates`, and group-indicator layout.
[`dtu_diffsplice()`](https://younthing.github.io/bulkMAE/reference/dtu_diffsplice.md)
operates after a native edgeR or limma exon/transcript-level fit.

``` r

time_result <- de_masigpro(
  time_mae,
  "rna",
  edesign = de_masigpro_design(
    time_mae,
    "rna",
    time = "timepoint",
    replicate = "replicate",
    group = "condition"
  ),
  assay = "counts",
  degree = 2,
  counts = TRUE
)

splice_result <- dtu_diffsplice(
  exon_fit,
  gene_id = mae_feature_data(exon_mae, "rna")$gene_id,
  feature_id = rownames(mae_feature_data(exon_mae, "rna"))
)
```

## 8. Tissue deconvolution

``` r

reference_sce <- deconv_reference(
  counts = single_cell_counts,
  cell_data = single_cell_metadata,
  cell_type = "cell_type",
  sample = "donor"
)

fractions <- deconv(
  tpm_mae,
  "rna",
  assay = "tpm",
  method = "quantiseq"
)

music_fit <- deconv_music(
  mae,
  "rna",
  assay = "counts",
  sc_reference = reference_sce,
  clusters = "cell_type",
  samples = "sample_id"
)

cibersortx_table <- deconv_cibersortx_input(tpm_mae, "rna")
```

The CIBERSORTx helper stops at input preparation because execution
occurs in an external service or container. Check data-governance
requirements before any upload. BayesPrism and MuSiC require an
appropriate single-cell reference; the MAE holds the bulk data, not the
reference-specific modelling choices.

## 9. Clinical models and drug reversal

Feature selection must occur inside each resampling fold. Selecting
genes on the full cohort before cross-validation leaks outcome
information.

``` r

cox_fit <- surv_cox(
  vst_mae,
  "rna",
  formula = surv_formula(
    time = "survival_time",
    event = "event",
    predictors = c("age", "risk_score")
  )
)

km_fit <- surv_km(
  vst_mae,
  "rna",
  formula = surv_formula(
    time = "survival_time",
    event = "event",
    predictors = "condition"
  )
)

roc_fit <- surv_roc(
  vst_mae,
  "rna",
  time = "survival_time",
  event = "event",
  marker = "risk_score",
  times = c(300, 600, 900)
)

penalized <- surv_penalized(
  vst_mae,
  "rna",
  time = "survival_time",
  event = "event",
  assay = "vst",
  folds = 5
)

entrez_ranks <- annotate_rekey(
  ranks,
  id_mapping,
  source = "source_id",
  target = "target_id",
  duplicates = "max_abs"
)
query <- drug_query(entrez_ranks, n = 150)
drug_lincs_databases()
lincs_result <- drug_lincs(query, reference_database = "lincs2")
```

For disease reversal, prioritize negative connectivity scores. Interpret
them with cell line, dose, time, mechanism, toxicity, and external
validation; a connectivity score is not clinical evidence.

## 10. Reproducibility

Record the MAE creation code, package versions, design formulas,
contrasts, feature filters, gene identifier release, pathway database
release, random seeds, and the exact external reference used for
deconvolution or LINCS. The package stays stateless so these choices
remain visible in the analysis script.
