# Plot GSEA rank-metric distributions as ridges

Draws density ridges and gene-position barcodes for explicitly selected
terms. This visualizes where member genes lie in the ranked metric; it
does not refit GSEA or select terms.

## Usage

``` r
plot_gsea_ridge(
  result,
  terms,
  ranks = NULL,
  gene_sets = NULL,
  membership = c("leading_edge", "gene_set"),
  term_labels = NULL,
  show_statistics = FALSE
)
```

## Arguments

- result:

  A native clusterProfiler gseaResult or a GSEA result table.

- terms:

  Unique explicit term identifiers, in top-to-bottom order.

- ranks:

  Optional decreasing, gene-named, finite numeric rank vector.

- gene_sets:

  Optional uniquely named list mapping terms to genes.

- membership:

  Genes used for each density: leading_edge or gene_set.

- term_labels:

  Optional complete term-named labels for selected terms.

- show_statistics:

  Show available NES and p-value text.

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
[`plot_lincs_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_heatmap.md),
[`plot_lincs_overlap()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_overlap.md),
[`plot_lincs_rank()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_rank.md),
[`plot_lincs_waterfall()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_waterfall.md),
[`plot_ora_bubble()`](https://younthing.github.io/bulkMAE/reference/plot_ora_bubble.md),
[`plot_ora_network()`](https://younthing.github.io/bulkMAE/reference/plot_ora_network.md),
[`plot_ora_radial()`](https://younthing.github.io/bulkMAE/reference/plot_ora_radial.md),
[`plot_qc_correlation()`](https://younthing.github.io/bulkMAE/reference/plot_qc_correlation.md),
[`plot_qc_library()`](https://younthing.github.io/bulkMAE/reference/plot_qc_library.md),
[`plot_qc_outliers()`](https://younthing.github.io/bulkMAE/reference/plot_qc_outliers.md),
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md),
[`theme_bulkmae()`](https://younthing.github.io/bulkMAE/reference/theme_bulkmae.md)
