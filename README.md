# bulkMAE

`bulkMAE` is a deliberately small adapter layer for bulk transcriptomics in R.
Every assay-level workflow starts from a `MultiAssayExperiment` (MAE), selects
one `SummarizedExperiment` leaf, aligns it with primary sample metadata through
`sampleMap`, and calls an established backend package.

The public API is organized as `<family>_<method>()` or
`<family>_<operation>()`. For example, typing `de_`, `enrich_`, `score_`,
`coexpr_`, or `surv_` exposes one coherent analysis family in autocomplete.

The package does **not** keep analysis history, mutate the input MAE, invent a
result registry, or wrap backend result objects in a new class. DESeq2 returns a
`DESeqDataSet`, edgeR returns its native fit/test object, limma returns an
`MArrayLM`, and so on.

## Contract

```r
result <- analysis_function(
  x = mae,
  experiment = "rna",
  assay = "counts",
  ...
)
```

- `x` is always the MAE containing samples and assays.
- `experiment` is explicit; hidden “active assays” are not used.
- `assay` is explicit whenever more than one scale is plausible.
- samples are aligned with `MultiAssayExperiment::getWithColData()`.
- transformed data and model fits are returned; they are not written back.

## Input completeness

Starting from a feature-by-sample matrix and aligned sample metadata, the
public API can construct an MAE (`mae_from_matrix()`), add assays, experiments,
sample covariates, and feature annotations, and convert native DE results into
ranks, selections, or effect tables. `mae_simulate()` supplies an installed,
fully synthetic MAE for examples and offline smoke tests.

The package also provides the bridges that previously forced users to assemble
backend objects by hand: `de_design()`/`de_contrast()`,
`de_masigpro_design()`, `adjust_covariates()`, `gene_sets_prepare()`,
`annotation_orgdb()`, `annotate_rekey()`, `deconv_reference()`,
`surv_formula()`, and native-result accessors such as `de_table()` and
`coexpr_modules()`.

Expression and sample metadata alone cannot determine study facts or reference
biology. Real survival outcomes, pairing/time semantics, RUV negative controls,
locked signature parameters, transcript/effective lengths, a tissue-matched
single-cell reference, validation cohorts, and cross-study effect estimates
must come from the study or a declared reference. BioMart, KEGG, STRING,
OmniPath/MSigDB, and LINCS are online or cache-backed; CIBERSORTx execution is
external. These are explicit input/resource boundaries, not values that
`bulkMAE` silently guesses or replaces with simulation.

## Coverage

| Layer | Main functions | Native backends |
|---|---|---|
| Access/import | `mae_from_matrix()`, `mae_add_experiment()`, `mae_add_sample_data()`, `import_tximport()` | MAE, SE, tximport |
| Annotation/QC | `annotation_orgdb()`, `annotate_rekey()`, `annotate_gene_lengths()`, `filter_expr()`, `reduce_pca()` | AnnotationDbi, biomaRt, edgeR, limma |
| Preprocessing | `normalize_tmm()`, `transform_vst()`, `adjust_combat()`, `adjust_batch()` | edgeR, DESeq2, sva, limma |
| Differential | `de_deseq2()`, `de_edger()`, `de_limma()`, `de_dream()`, `de_masigpro()` | DESeq2, edgeR, limma, dream, maSigPro |
| Splicing/co-expression | `dtu_diffsplice()`, `coexpr_differential()` | limma, edgeR, diffcoexp |
| Gene sets/pathways | `gene_sets_prepare()`, `gene_sets_read_gmt()`, `gene_sets_msigdb()`, `enrich_fgsea()`, `score_gsva()` | msigdbr, clusterProfiler, goseq, fgsea, GSVA |
| Regulatory activity | `activity_decouple()`, `activity_progeny()`, `activity_tf()` | decoupleR |
| Systems/subtypes | `cluster_consensus()`, `cluster_nmf()`, `coexpr_wgcna()`, `coexpr_preservation()`, `network_genie3()` | ConsensusClusterPlus, NMF, WGCNA, GENIE3 |
| Deconvolution | `deconv_reference()`, `deconv()`, `deconv_music()`, `deconv_bayesprism()` | SingleCellExperiment, immunedeconv, MuSiC, BayesPrism |
| Clinical | `score_signature()`, `surv_formula()`, `surv_cox()`, `meta_collect()`, `meta_effect()`, `drug_lincs()` | survival, glmnet, timeROC, metafor, signatureSearch |

`activity_decouple()` is intentionally cross-cutting. Its `statistics` argument can
select enrichment-style algorithms (`aucell`, `fgsea`, `gsva`, `ora`) or
network-aware activity estimators (`mlm`, `ulm`, `viper`, `wmean`, `wsum`), and
can ask decoupleR for a consensus. Use `activity_methods()` to inspect the
methods supported by the installed decoupleR release. For signed or weighted
methods, map the resource's weight column explicitly with `mor = "weight"` (or
the corresponding column name).

`coexpr_wgcna()` automatically applies `goodSamplesGenes()` filtering. It warns
when samples or features are removed and appends their names plus the native
quality result to the returned WGCNA list; it does not pause for confirmation.

## Minimal example

```r
library(bulkMAE)

set.seed(1)
sample_data <- data.frame(
  condition = factor(rep(c("control", "treated"), each = 4)),
  batch = factor(rep(c("A", "B"), times = 4)),
  row.names = sprintf("sample%02d", seq_len(8))
)
counts <- matrix(
  rnbinom(120 * 8, mu = 100, size = 5),
  nrow = 120,
  dimnames = list(
    sprintf("gene%03d", seq_len(120)),
    rownames(sample_data)
  )
)
feature_data <- data.frame(
  gene_length = seq(500, 2500, length.out = nrow(counts)),
  row.names = rownames(counts)
)

mae <- mae_from_matrix(
  expression = counts,
  samples = sample_data,
  row_data = feature_data,
  experiment = "rna",
  assay = "counts"
)

tpm <- normalize_tpm(mae, "rna", lengths = "gene_length")
mae_with_tpm <- mae_add_assay(mae, "rna", tpm, name = "tpm")
mae_small <- mae_subset_features(
  mae_with_tpm,
  "rna",
  features = rownames(counts)[seq_len(50)]
)

mae_assays(mae_with_tpm, "rna")
dim(mae_pull_assay(mae_small, "rna", "tpm"))
colSums(mae_pull_assay(mae_with_tpm, "rna", "tpm"))
```

These helpers return modified copies: `mae` remains the raw-count input, while
`mae_with_tpm` and `mae_small` make each added or filtered data state explicit.

Read the [complete walkthrough](vignettes/getting-started.Rmd), the
[method-selection guide](inst/guides/expanded-methods-zh.md), and the
[input-completeness and resource-boundary audit](inst/guides/input-completeness-zh.md),
plus the [implementation-to-documentation audit map](inst/guides/official-sources.md).
Users upgrading from 0.3 or earlier should also read the
[0.4 naming migration map](inst/guides/naming-migration-zh.md).
After installation, the audit files are also available under
`system.file("guides", package = "bulkMAE")`.

## Dependency policy

Only MAE/SE infrastructure is imported. Analysis engines are optional and
checked when their wrapper is called. This lets users install only the methods
needed for a project instead of forcing one large, conflict-prone environment.

MuSiC, immunedeconv, and BayesPrism are declared optional dependencies but are
not available from every standard Bioconductor/CRAN repository. Their wrappers
resolve the installed backend only when called and restore any temporary search
path compatibility changes; install them from their official repositories and
record the commit or release in the project lockfile.

See the [Chinese pak dependency installation guide](inst/guides/dependency-installation-zh.md)
for the complete CRAN/Bioconductor backend list, verified GitHub package
specifications, and notes about resources that package installation cannot
provide.
