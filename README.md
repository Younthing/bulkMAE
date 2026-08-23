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

## Coverage

| Layer | Main functions | Native backends |
|---|---|---|
| Access/import | `mae_pull_experiment()`, `mae_create()`, `import_tximport()` | MAE, SE, tximport |
| Annotation/QC | `annotate_ids()`, `filter_expr()`, `reduce_pca()`, `reduce_umap()` | AnnotationDbi, edgeR, limma, uwot |
| Preprocessing | `normalize_tmm()`, `transform_vst()`, `adjust_combat()`, `adjust_batch()` | edgeR, DESeq2, sva, limma |
| Differential | `de_deseq2()`, `de_edger()`, `de_limma()`, `de_dream()`, `de_masigpro()` | DESeq2, edgeR, limma, dream, maSigPro |
| Splicing/co-expression | `dtu_diffsplice()`, `coexpr_differential()` | limma, edgeR, diffcoexp |
| Gene sets/pathways | `enrich_ora()`, `enrich_goseq()`, `enrich_fgsea()`, `enrich_camera()`, `score_gsva()`, `score_singscore()` | clusterProfiler, goseq, fgsea, limma, GSVA, singscore |
| Regulatory activity | `activity_decouple()`, `activity_progeny()`, `activity_tf()` | decoupleR |
| Systems/subtypes | `cluster_consensus()`, `cluster_nmf()`, `coexpr_wgcna()`, `coexpr_preservation()`, `network_genie3()` | ConsensusClusterPlus, NMF, WGCNA, GENIE3 |
| Deconvolution | `deconv()`, `deconv_music()`, `deconv_bayesprism()` | immunedeconv, MuSiC, BayesPrism |
| Clinical | `score_signature()`, `surv_cox()`, `surv_roc()`, `ml_glmnet()`, `meta_effect()`, `drug_lincs()` | survival, glmnet, timeROC, metafor, signatureSearch |

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

rna <- mae_create_experiment(
  assays = counts,
  col_data = sample_table,
  assay_name = "counts"
)

mae <- mae_create(
  experiments = list(rna = rna),
  col_data = sample_table
)

dds <- de_deseq2(mae, "rna", design = ~ batch + condition)
res <- de_deseq2_results(dds, contrast = c("condition", "treated", "control"))
```

Read the [complete walkthrough](vignettes/getting-started.Rmd), the
[method-selection guide](inst/guides/expanded-methods-zh.md), and the
[implementation-to-documentation audit map](inst/guides/official-sources.md).
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
load the installed namespace only when called; install them from their official
repositories and record the commit or release in the project lockfile.
