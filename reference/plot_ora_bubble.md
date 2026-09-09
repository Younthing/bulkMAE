# Plot explicitly selected ORA terms as bubbles

Displays one of three distinct over-representation metrics without
treating them as interchangeable: rich factor is overlap divided by the
background term size, gene ratio is overlap divided by the input list
size, and fold enrichment is gene ratio divided by background ratio.
Point area represents overlap count and colour represents finite
`-log10()` evidence; adjusted evidence is labelled compactly as `adj p`.

## Usage

``` r
plot_ora_bubble(
  result,
  terms,
  x = c("rich_factor", "gene_ratio", "fold_enrichment"),
  p_value = c("auto", "adjusted", "raw"),
  term_labels = NULL
)
```

## Arguments

- result:

  A clusterProfiler `enrichResult` or an ORA result data frame
  containing `ID`, `Description`, `GeneRatio`, `BgRatio`, `pvalue`,
  `p.adjust`, and `Count`.

- terms:

  Unique term identifiers to display, in the desired order. Terms are
  never selected automatically.

- x:

  ORA metric displayed on the horizontal axis.

- p_value:

  Evidence column to display. `"auto"` uses adjusted p-values only when
  every selected term has one; otherwise the whole plot uses raw
  p-values.

- term_labels:

  Optional complete term-named character vector describing the selected
  `terms`.

## Value

An unprinted standard ggplot object carrying recommended physical
dimensions for
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md).

## See also

Other plotting:
[`plot_assay_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_assay_heatmap.md),
[`plot_de_ma()`](https://younthing.github.io/bulkMAE/reference/plot_de_ma.md),
[`plot_de_volcano()`](https://younthing.github.io/bulkMAE/reference/plot_de_volcano.md),
[`plot_embedding()`](https://younthing.github.io/bulkMAE/reference/plot_embedding.md),
[`plot_gsea_classic()`](https://younthing.github.io/bulkMAE/reference/plot_gsea_classic.md),
[`plot_gsea_ridge()`](https://younthing.github.io/bulkMAE/reference/plot_gsea_ridge.md),
[`plot_ora_network()`](https://younthing.github.io/bulkMAE/reference/plot_ora_network.md),
[`plot_ora_radial()`](https://younthing.github.io/bulkMAE/reference/plot_ora_radial.md),
[`plot_qc_correlation()`](https://younthing.github.io/bulkMAE/reference/plot_qc_correlation.md),
[`plot_qc_library()`](https://younthing.github.io/bulkMAE/reference/plot_qc_library.md),
[`plot_qc_outliers()`](https://younthing.github.io/bulkMAE/reference/plot_qc_outliers.md),
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md),
[`theme_bulkmae()`](https://younthing.github.io/bulkMAE/reference/theme_bulkmae.md)
