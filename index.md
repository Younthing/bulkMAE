# bulkMAE

[![R-CMD-check](https://github.com/Younthing/bulkMAE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Younthing/bulkMAE/actions/workflows/R-CMD-check.yaml)
[![Full backend
check](https://github.com/Younthing/bulkMAE/actions/workflows/full-backend-check.yaml/badge.svg)](https://github.com/Younthing/bulkMAE/actions/workflows/full-backend-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Younthing/bulkMAE/graph/badge.svg)](https://app.codecov.io/gh/Younthing/bulkMAE)
[![pkgdown](https://github.com/Younthing/bulkMAE/actions/workflows/pkgdown.yaml/badge.svg)](https://younthing.github.io/bulkMAE/)

`bulkMAE` is a deliberately small adapter layer for bulk transcriptomics
in R. Every assay-level workflow starts from a `MultiAssayExperiment`
(MAE), selects one `SummarizedExperiment` leaf, aligns it with primary
sample metadata through `sampleMap`, and calls an established backend
package.

The public API is organized as `<family>_<method>()` or
`<family>_<operation>()`. For example, typing `de_`, `enrich_`,
`score_`, `coexpr_`, or `surv_` exposes one coherent analysis family in
autocomplete.

The package does **not** keep analysis history, mutate the input MAE,
invent a result registry, or wrap backend result objects in a new class.
DESeq2 returns a `DESeqDataSet`, edgeR returns its native fit/test
object, limma returns an `MArrayLM`, and so on.

## Contract

``` r

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
- samples are aligned with
  [`MultiAssayExperiment::getWithColData()`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-helpers.html).
- transformed data and model fits are returned; they are not written
  back.

## Input completeness

Starting from a feature-by-sample matrix and aligned sample metadata,
the public API can construct an MAE
([`mae_from_matrix()`](https://younthing.github.io/bulkMAE/reference/mae_from_matrix.md)),
add assays, experiments, sample covariates, and feature annotations, and
convert native DE results into ranks, selections, or effect tables.
[`mae_simulate()`](https://younthing.github.io/bulkMAE/reference/mae_simulate.md)
supplies an installed, fully synthetic MAE for examples and offline
smoke tests.

The package also provides the bridges that previously forced users to
assemble backend objects by hand:
[`de_design()`](https://younthing.github.io/bulkMAE/reference/de_design.md)/[`de_contrast()`](https://younthing.github.io/bulkMAE/reference/de_contrast.md),
[`de_masigpro_design()`](https://younthing.github.io/bulkMAE/reference/de_masigpro_design.md),
[`adjust_covariates()`](https://younthing.github.io/bulkMAE/reference/adjust_covariates.md),
[`gene_sets_prepare()`](https://younthing.github.io/bulkMAE/reference/gene_sets_prepare.md),
[`annotation_orgdb()`](https://younthing.github.io/bulkMAE/reference/annotation_orgdb.md),
[`annotate_rekey()`](https://younthing.github.io/bulkMAE/reference/annotate_rekey.md),
[`deconv_reference()`](https://younthing.github.io/bulkMAE/reference/deconv_reference.md),
[`surv_formula()`](https://younthing.github.io/bulkMAE/reference/surv_formula.md),
and native-result accessors such as
[`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md)
and
[`coexpr_modules()`](https://younthing.github.io/bulkMAE/reference/coexpr_modules.md).

Expression and sample metadata alone cannot determine study facts or
reference biology. Real survival outcomes, pairing/time semantics, RUV
negative controls, locked signature parameters, transcript/effective
lengths, a tissue-matched single-cell reference, validation cohorts, and
cross-study effect estimates must come from the study or a declared
reference. BioMart, KEGG, STRING, OmniPath/MSigDB, and LINCS are online
or cache-backed; CIBERSORTx execution is external. These are explicit
input/resource boundaries, not values that `bulkMAE` silently guesses or
replaces with simulation.

## Coverage

| Layer | Main functions | Native backends |
|----|----|----|
| Access/import | [`mae_from_matrix()`](https://younthing.github.io/bulkMAE/reference/mae_from_matrix.md), [`mae_add_experiment()`](https://younthing.github.io/bulkMAE/reference/mae_add_experiment.md), [`mae_add_sample_data()`](https://younthing.github.io/bulkMAE/reference/mae_add_sample_data.md), [`import_tximport()`](https://younthing.github.io/bulkMAE/reference/import_tximport.md) | MAE, SE, tximport |
| Annotation/QC | [`annotation_orgdb()`](https://younthing.github.io/bulkMAE/reference/annotation_orgdb.md), [`annotate_rekey()`](https://younthing.github.io/bulkMAE/reference/annotate_rekey.md), [`annotate_gene_lengths()`](https://younthing.github.io/bulkMAE/reference/annotate_gene_lengths.md), [`filter_expr()`](https://younthing.github.io/bulkMAE/reference/filter_expr.md), [`reduce_pca()`](https://younthing.github.io/bulkMAE/reference/reduce_pca.md) | AnnotationDbi, biomaRt, edgeR, limma |
| Preprocessing | [`normalize_tmm()`](https://younthing.github.io/bulkMAE/reference/normalize_tmm.md), [`transform_vst()`](https://younthing.github.io/bulkMAE/reference/transform_vst.md), [`adjust_combat()`](https://younthing.github.io/bulkMAE/reference/adjust_combat.md), [`adjust_batch()`](https://younthing.github.io/bulkMAE/reference/adjust_batch.md) | edgeR, DESeq2, sva, limma |
| Differential | [`de_deseq2()`](https://younthing.github.io/bulkMAE/reference/de_deseq2.md), [`de_edger()`](https://younthing.github.io/bulkMAE/reference/de_edger.md), [`de_limma()`](https://younthing.github.io/bulkMAE/reference/de_limma.md), [`de_dream()`](https://younthing.github.io/bulkMAE/reference/de_dream.md), [`de_masigpro()`](https://younthing.github.io/bulkMAE/reference/de_masigpro.md) | DESeq2, edgeR, limma, dream, maSigPro |
| Splicing/co-expression | [`dtu_diffsplice()`](https://younthing.github.io/bulkMAE/reference/dtu_diffsplice.md), [`coexpr_differential()`](https://younthing.github.io/bulkMAE/reference/coexpr_differential.md) | limma, edgeR, diffcoexp |
| Gene sets/pathways | [`gene_sets_prepare()`](https://younthing.github.io/bulkMAE/reference/gene_sets_prepare.md), [`gene_sets_read_gmt()`](https://younthing.github.io/bulkMAE/reference/gene_sets_read_gmt.md), [`gene_sets_msigdb()`](https://younthing.github.io/bulkMAE/reference/gene_sets_msigdb.md), [`enrich_fgsea()`](https://younthing.github.io/bulkMAE/reference/enrich_fgsea.md), [`score_gsva()`](https://younthing.github.io/bulkMAE/reference/score_gsva.md) | msigdbr, clusterProfiler, goseq, fgsea, GSVA |
| Regulatory activity | [`activity_decouple()`](https://younthing.github.io/bulkMAE/reference/activity_decouple.md), [`activity_progeny()`](https://younthing.github.io/bulkMAE/reference/activity_progeny.md), [`activity_tf()`](https://younthing.github.io/bulkMAE/reference/activity_tf.md) | decoupleR |
| Systems/subtypes | [`cluster_consensus()`](https://younthing.github.io/bulkMAE/reference/cluster_consensus.md), [`cluster_nmf()`](https://younthing.github.io/bulkMAE/reference/cluster_nmf.md), [`coexpr_wgcna()`](https://younthing.github.io/bulkMAE/reference/coexpr_wgcna.md), [`coexpr_preservation()`](https://younthing.github.io/bulkMAE/reference/coexpr_preservation.md), [`network_genie3()`](https://younthing.github.io/bulkMAE/reference/network_genie3.md) | ConsensusClusterPlus, NMF, WGCNA, GENIE3 |
| Deconvolution | [`deconv_reference()`](https://younthing.github.io/bulkMAE/reference/deconv_reference.md), [`deconv()`](https://younthing.github.io/bulkMAE/reference/deconv.md), [`deconv_music()`](https://younthing.github.io/bulkMAE/reference/deconv_music.md), [`deconv_bayesprism()`](https://younthing.github.io/bulkMAE/reference/deconv_bayesprism.md) | SingleCellExperiment, immunedeconv, MuSiC, BayesPrism |
| Clinical | [`score_signature()`](https://younthing.github.io/bulkMAE/reference/score_signature.md), [`surv_formula()`](https://younthing.github.io/bulkMAE/reference/surv_formula.md), [`surv_cox()`](https://younthing.github.io/bulkMAE/reference/surv_cox.md), [`meta_collect()`](https://younthing.github.io/bulkMAE/reference/meta_collect.md), [`meta_effect()`](https://younthing.github.io/bulkMAE/reference/meta_effect.md), [`drug_lincs()`](https://younthing.github.io/bulkMAE/reference/drug_lincs.md) | survival, glmnet, timeROC, metafor, signatureSearch |

[`activity_decouple()`](https://younthing.github.io/bulkMAE/reference/activity_decouple.md)
is intentionally cross-cutting. Its `statistics` argument can select
enrichment-style algorithms (`aucell`, `fgsea`, `gsva`, `ora`) or
network-aware activity estimators (`mlm`, `ulm`, `viper`, `wmean`,
`wsum`), and can ask decoupleR for a consensus. Use
[`activity_methods()`](https://younthing.github.io/bulkMAE/reference/activity_methods.md)
to inspect the methods supported by the installed decoupleR release. For
signed or weighted methods, map the resource’s weight column explicitly
with `mor = "weight"` (or the corresponding column name).

[`coexpr_wgcna()`](https://younthing.github.io/bulkMAE/reference/coexpr_wgcna.md)
automatically applies `goodSamplesGenes()` filtering. It warns when
samples or features are removed and appends their names plus the native
quality result to the returned WGCNA list; it does not pause for
confirmation.

## Minimal example

``` r

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

These helpers return modified copies: `mae` remains the raw-count input,
while `mae_with_tpm` and `mae_small` make each added or filtered data
state explicit.

Read the [complete
walkthrough](https://younthing.github.io/bulkMAE/articles/getting-started.html),
the [method-selection
guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/expanded-methods-zh.md),
and the [input-completeness and resource-boundary
audit](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/input-completeness-zh.md),
plus the [implementation-to-documentation audit
map](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/official-sources.md).
Users upgrading from 0.3 or earlier should also read the [0.4 naming
migration
map](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/naming-migration-zh.md).
After installation, the audit files are also available under
`system.file("guides", package = "bulkMAE")`.

## Continuous integration

Pull requests and the default branch run a cross-platform R CMD check
against the package’s hard dependencies. A separate Ubuntu job installs
every optional offline analysis backend with `pak`, while the coverage
workflow exercises the same full backend set. Network-backed tests
remain opt-in so transient failures from BioMart, KEGG, STRING,
OmniPath, or LINCS do not block ordinary changes.

The [pkgdown website](https://younthing.github.io/bulkMAE/) is rebuilt
from `main` and published from the generated `gh-pages` branch.

## Dependency policy

Only MAE/SE infrastructure is imported. Analysis engines are optional
and checked when their wrapper is called. This lets users install only
the methods needed for a project instead of forcing one large,
conflict-prone environment.

MuSiC, immunedeconv, and BayesPrism are declared optional dependencies
but are not available from every standard Bioconductor/CRAN repository.
Their wrappers resolve the installed backend only when called and
restore any temporary search path compatibility changes; install them
from their official repositories and record the commit or release in the
project lockfile.

See the [Chinese pak dependency installation
guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/dependency-installation-zh.md)
for the complete CRAN/Bioconductor backend list, verified GitHub package
specifications, and notes about resources that package installation
cannot provide.
