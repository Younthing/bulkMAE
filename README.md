# bulkMAE

<!-- badges: start -->
[![R-CMD-check](https://github.com/Younthing/bulkMAE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Younthing/bulkMAE/actions/workflows/R-CMD-check.yaml)
[![Full backend check](https://github.com/Younthing/bulkMAE/actions/workflows/full-backend-check.yaml/badge.svg)](https://github.com/Younthing/bulkMAE/actions/workflows/full-backend-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Younthing/bulkMAE/graph/badge.svg)](https://app.codecov.io/gh/Younthing/bulkMAE)
[![pkgdown](https://github.com/Younthing/bulkMAE/actions/workflows/pkgdown.yaml/badge.svg)](https://younthing.github.io/bulkMAE/)
[![GitHub release](https://img.shields.io/github/v/release/Younthing/bulkMAE)](https://github.com/Younthing/bulkMAE/releases/latest)
<!-- badges: end -->

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

## Installation

Install the current GitHub Release (`v0.4.0`):

```r
# install.packages("pak")
pak::pak("Younthing/bulkMAE@v0.4.0")
```

Or install the source tarball attached to that
[release](https://github.com/Younthing/bulkMAE/releases/tag/v0.4.0):

```r
install.packages(
  "https://github.com/Younthing/bulkMAE/releases/download/v0.4.0/bulkMAE_0.4.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

The development tree on `main` is:

```r
pak::pak("Younthing/bulkMAE")
```

Analysis backends remain optional and are installed per project. See
[Dependency policy](#dependency-policy) and the
[Chinese pak dependency installation guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/dependency-installation-zh.md).

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

## Publication-size plots

Every `plot_*()` helper returns a standard ggplot object and carries a
recommended physical size tuned to that plot family. Continue styling with
ordinary ggplot2 syntax, then let `plot_save()` export at that final size:

```r
qc_plot <- plot_qc_library(qc_library(mae, "rna", "counts")) +
  labs(title = "Library QC") +
  theme(legend.position = "bottom")

plot_save("library-qc.pdf", qc_plot)
```

The default is `scale = 1`; PDF output uses Cairo when available. Override
either dimension explicitly when a journal requires an exact layout:

```r
plot_save("library-qc.pdf", qc_plot, width = 8.5, height = 11, units = "cm")
```

`theme_bulkmae()` uses publication-oriented 6 pt base text, and data labels
drawn by bulkMAE use point units explicitly. A theme controls styling but not
the physical graphics device: RStudio's plot pane and knitr chunks therefore
do not read the recommended dimensions. Use `plot_save()` for exact files, and
set chunk `fig.width` / `fig.height` explicitly when document previews also
need a fixed aspect ratio.

Specialized enrichment constructors keep method-specific inputs and statistics
explicit:

```r
plot_gsea_classic(gsea_result, term = "HALLMARK_INFLAMMATORY_RESPONSE")
plot_gsea_ridge(gsea_result, terms = pathway_ids)

plot_ora_bubble(ora_result, terms = pathway_ids)
plot_ora_network(ora_result, terms = pathway_ids, feature_values = effects)
plot_ora_radial(
  ora_result,
  terms = pathway_ids,
  feature_values = effects,
  label_features = c("IL6", "CXCL8")
)
```

Every selection is explicit. A native clusterProfiler GSEA object carries its
ranked vector, gene sets, and weighting exponent; a tabular fgsea result does
not, so classic plots require those analysis inputs to be supplied rather than
reconstructing them heuristically. ORA bubbles distinguish rich factor, gene
ratio, and fold enrichment. Community-network term edges use Jaccard overlap
from complete enriched-feature membership, while radial term edges encode the
number of shared enriched features.

Read the [Chinese executable getting-started guide](https://younthing.github.io/bulkMAE/articles/getting-started.html), the
[airway QC and paired differential-expression tutorial](https://younthing.github.io/bulkMAE/articles/airway-qc-de.html), the
[airway GO ORA、GSEA 与绘图教程](https://younthing.github.io/bulkMAE/articles/enrichment-analysis.html), the
[method-selection guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/expanded-methods-zh.md), and the
[input-completeness and resource-boundary audit](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/input-completeness-zh.md),
plus the [implementation-to-documentation audit map](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/official-sources.md).
Users upgrading from 0.3 or earlier should also read the
[0.4 naming migration map](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/naming-migration-zh.md).
After installation, the audit files are also available under
`system.file("guides", package = "bulkMAE")`.

## Continuous integration

Pull requests run one Ubuntu R-release CMD check against hard dependencies.
After merge, `main` runs macOS/Windows/Ubuntu-devel portability checks, the
full-backend job, coverage, and pkgdown. Those heavier workflows also run
nightly or from **Actions → Run workflow**. Network-backed tests remain
opt-in so transient failures from BioMart, KEGG, STRING, OmniPath, or LINCS
do not block ordinary changes.

The [pkgdown website](https://younthing.github.io/bulkMAE/) is rebuilt from
`main` and from a published Release, then deployed to `gh-pages`.

## Dependency policy

Only MAE/SE infrastructure and ggplot2 are imported. Analysis engines are optional and
checked when their wrapper is called. This lets users install only the methods
needed for a project instead of forcing one large, conflict-prone environment.

MuSiC, immunedeconv, and BayesPrism are declared optional dependencies but are
not available from every standard Bioconductor/CRAN repository. Their wrappers
resolve the installed backend only when called and restore any temporary search
path compatibility changes; install them from their official repositories and
record the commit or release in the project lockfile.

See the [Chinese pak dependency installation guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/dependency-installation-zh.md)
for the complete CRAN/Bioconductor backend list, verified GitHub package
specifications, and notes about resources that package installation cannot
provide.
