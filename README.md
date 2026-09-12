# bulkMAE

<!-- badges: start -->
[![R-CMD-check](https://github.com/Younthing/bulkMAE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Younthing/bulkMAE/actions/workflows/R-CMD-check.yaml)
[![Full backend check](https://github.com/Younthing/bulkMAE/actions/workflows/full-backend-check.yaml/badge.svg)](https://github.com/Younthing/bulkMAE/actions/workflows/full-backend-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Younthing/bulkMAE/graph/badge.svg)](https://app.codecov.io/gh/Younthing/bulkMAE)
[![pkgdown](https://github.com/Younthing/bulkMAE/actions/workflows/pkgdown.yaml/badge.svg)](https://younthing.github.io/bulkMAE/)
[![GitHub release](https://img.shields.io/github/v/release/Younthing/bulkMAE)](https://github.com/Younthing/bulkMAE/releases/latest)
[![License: Artistic-2.0](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://opensource.org/licenses/Artistic-2.0)
<!-- badges: end -->

Stateless adapters from
[`MultiAssayExperiment`](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html)
to established bulk transcriptomics methods.

Each analysis call names an MAE, one `SummarizedExperiment` leaf, and an
assay. `bulkMAE` aligns samples through `sampleMap` and returns the backend's
native object. It does not keep analysis history, mutate the input MAE, or wrap
results in a new class. DESeq2 still returns a `DESeqDataSet`. edgeR still
returns its fit and test objects. limma still returns an `MArrayLM`.

The public API is `<family>_<method>()` or `<family>_<operation>()`. Typing
`de_`, `enrich_`, `score_`, `coexpr_`, or `surv_` shows one family in
autocomplete.

Read the [pkgdown site](https://younthing.github.io/bulkMAE/) for reference
pages and executable tutorials.

## Install

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

Analysis backends stay optional. Install only the engines a project uses. See
[Dependencies](#dependencies) and the
[Chinese pak installation guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/dependency-installation-zh.md).

## Quick start

`mae_simulate()` builds a small synthetic MAE for examples. The annotations
and outcomes have no biological meaning.

```r
library(bulkMAE)

mae <- mae_simulate(n_features = 120, n_samples = 8, seed = 1)
mae_validate(mae, "rna")
mae_assays(mae, "rna")
```

Import a real feature-by-sample matrix with aligned sample metadata:

```r
mae <- mae_from_matrix(
  expression = counts,
  samples = sample_data,
  row_data = feature_data,
  experiment = "rna",
  assay = "counts"
)
```

Helpers such as `mae_add_assay()` and `mae_subset_features()` return modified
copies. The input MAE stays unchanged.

## Analysis families

| Layer | Main functions | Native backends |
|---|---|---|
| Access and import | `mae_from_matrix()`, `mae_add_experiment()`, `mae_add_sample_data()`, `import_tximport()` | [MAE](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html), [SE](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html), [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) |
| Annotation and QC | `annotation_orgdb()`, `annotate_rekey()`, `annotate_gene_lengths()`, `filter_expr()`, `reduce_pca()` | [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html), [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) |
| Preprocessing | `normalize_tmm()`, `transform_vst()`, `adjust_combat()`, `adjust_batch()` | [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [sva](https://bioconductor.org/packages/release/bioc/html/sva.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) |
| Differential | `de_deseq2()`, `de_edger()`, `de_limma()`, `de_dream()`, `de_masigpro()` | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [dream](https://bioconductor.org/packages/release/bioc/html/variancePartition.html), [maSigPro](https://bioconductor.org/packages/release/bioc/html/maSigPro.html) |
| Splicing and co-expression | `dtu_diffsplice()`, `coexpr_differential()` | [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [diffcoexp](https://bioconductor.org/packages/release/bioc/html/diffcoexp.html) |
| Gene sets and pathways | `gene_sets_prepare()`, `gene_sets_read_gmt()`, `gene_sets_msigdb()`, `enrich_fgsea()`, `score_gsva()` | [msigdbr](https://cran.r-project.org/package=msigdbr), [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [goseq](https://bioconductor.org/packages/release/bioc/html/goseq.html), [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html), [GSVA](https://bioconductor.org/packages/release/bioc/html/GSVA.html) |
| Regulatory activity | `activity_decouple()`, `activity_progeny()`, `activity_tf()` | [decoupleR](https://bioconductor.org/packages/release/bioc/html/decoupleR.html) |
| Systems/subtypes | `cluster_consensus()`, `cluster_nmf()`, `coexpr_wgcna()`, `coexpr_preservation()`, `network_genie3()` | [ConsensusClusterPlus](https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html), [NMF](https://cran.r-project.org/package=NMF), [WGCNA](https://cran.r-project.org/package=WGCNA), [GENIE3](https://bioconductor.org/packages/release/bioc/html/GENIE3.html) |
| Deconvolution | `deconv_reference()`, `deconv()`, `deconv_music()`, `deconv_bayesprism()` | [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html), [immunedeconv](https://omnideconv.org/immunedeconv/), [MuSiC](https://xuranw.github.io/MuSiC/), [BayesPrism](https://github.com/Danko-Lab/BayesPrism) |
| Clinical | `score_signature()`, `surv_formula()`, `surv_cox()`, `meta_collect()`, `meta_effect()`, `drug_lincs()` | [survival](https://cran.r-project.org/package=survival), [glmnet](https://cran.r-project.org/package=glmnet), [timeROC](https://cran.r-project.org/package=timeROC), [metafor](https://wviechtb.github.io/metafor/), [signatureSearch](https://bioconductor.org/packages/release/bioc/html/signatureSearch.html) |

Each backend name links to its official package page. Official pages for
wrappers outside this table are
[RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html),
[ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html),
[singscore](https://bioconductor.org/packages/release/bioc/html/singscore.html),
[STRINGdb](https://bioconductor.org/packages/release/bioc/html/STRINGdb.html),
[CIBERSORTx](https://cibersortx.stanford.edu/),
[uwot](https://cran.r-project.org/package=uwot),
[Rtsne](https://cran.r-project.org/package=Rtsne), and
[OmnipathR](https://bioconductor.org/packages/release/bioc/html/OmnipathR.html).

`activity_decouple()` can run enrichment-style statistics (`aucell`, `fgsea`,
`gsva`, `ora`) or network-aware estimators (`mlm`, `ulm`, `viper`, `wmean`,
`wsum`), and can ask decoupleR for a consensus. Use `activity_methods()` to
list methods in the installed decoupleR release. For signed or weighted
methods, map the resource weight column with `mor = "weight"` or the matching
column name.

`coexpr_wgcna()` applies `goodSamplesGenes()` before the network fit. It warns
when samples or features are removed and appends their names plus the native
quality result to the returned list.

## Call shape

```r
result <- analysis_function(
  x = mae,
  experiment = "rna",
  assay = "counts",
  ...
)
```

- `x` is the MAE that holds samples and assays.
- `experiment` is always named. There is no hidden active assay.
- `assay` is named whenever more than one scale is plausible.
- Samples are aligned with `MultiAssayExperiment::getWithColData()`.
- Fits and transformed matrices are returned. They are not written back.

Bridges that used to be assembled by hand include `de_design()`,
`de_contrast()`, `de_masigpro_design()`, `adjust_covariates()`,
`gene_sets_prepare()`, `annotation_orgdb()`, `annotate_rekey()`,
`deconv_reference()`, `surv_formula()`, `de_table()`, and `coexpr_modules()`.

## What you supply

Expression and sample metadata do not determine study facts or reference
biology. Survival outcomes, pairing and time semantics, RUV negative controls,
locked signature weights, transcript or effective lengths, a tissue-matched
single-cell reference, validation cohorts, and cross-study effect estimates
must come from the study or a declared resource.

BioMart, KEGG, STRING, OmniPath, MSigDB, and LINCS are online or cache-backed.
CIBERSORTx execution is external. `bulkMAE` does not fill those gaps with
simulation.

## Plots

<img src="man/figures/readme-gallery.png" alt="Volcano, heatmap, ORA community network, and GSEA ridge plots from bulkMAE">

## Documentation

- [Getting started](https://younthing.github.io/bulkMAE/articles/getting-started.html)
  (Chinese executable guide)
- [airway QC and paired differential expression](https://younthing.github.io/bulkMAE/articles/airway-qc-de.html)
- [airway GO ORA, GSEA, and plots](https://younthing.github.io/bulkMAE/articles/enrichment-analysis.html)
- P0 stubs: [co-expression modules](https://younthing.github.io/bulkMAE/articles/coexpression-modules.html),
  [molecular subtyping](https://younthing.github.io/bulkMAE/articles/molecular-subtyping.html),
  [immune deconvolution](https://younthing.github.io/bulkMAE/articles/immune-deconvolution.html)
- [Method-selection guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/expanded-methods-zh.md)
- [Input-completeness audit](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/input-completeness-zh.md)
- [Adapter-to-source map](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/official-sources.md)
- [0.4 naming migration](https://younthing.github.io/bulkMAE/articles/naming-migration.html)
  for upgrades from 0.3 or earlier. Call `bulkmae_rename()` to map an old name.

Installed copies also live under `system.file("guides", package = "bulkMAE")`.

## Tests and CI

Pull requests run one Ubuntu R-release CMD check against hard dependencies.
After merge, `main` runs macOS, Windows, and Ubuntu-devel portability checks,
the full-backend job, coverage, pkgdown, and BiocCheck. Those heavier workflows also run
nightly or from **Actions → Run workflow**. Network-backed tests stay opt-in
so BioMart, KEGG, STRING, OmniPath, or LINCS outages do not block ordinary
changes.

The [pkgdown site](https://younthing.github.io/bulkMAE/) rebuilds from `main`
and from a published Release, then deploys to `gh-pages`.

## Dependencies

Hard imports are MAE and SE infrastructure plus ggplot2. Analysis engines are
checked when their wrapper is called.

[MuSiC](https://xuranw.github.io/MuSiC/),
[immunedeconv](https://omnideconv.org/immunedeconv/), and
[BayesPrism](https://github.com/Danko-Lab/BayesPrism) are optional and are
not on every standard Bioconductor or CRAN repository. Install them from
their official sources and record the commit or release in the project
lockfile.

The [Chinese pak installation guide](https://github.com/Younthing/bulkMAE/blob/main/inst/guides/dependency-installation-zh.md)
lists CRAN and Bioconductor backends, verified GitHub specifications, and
resources that package installation cannot provide.

## License

`bulkMAE` is released under the [Artistic License 2.0](https://opensource.org/licenses/Artistic-2.0).
