# Plot selected assay features as a heatmap

Plot selected assay features as a heatmap

## Usage

``` r
plot_assay_heatmap(
  x,
  experiment,
  assay,
  features,
  scale = c("row", "none"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = NULL,
  feature_label = NULL
)
```

## Arguments

- x:

  A `MultiAssayExperiment`.

- experiment:

  Experiment name.

- assay:

  Assay name or index.

- features:

  Explicit feature identifiers to plot.

- scale:

  Standardize each feature row or retain assay values.

- cluster_rows, cluster_columns:

  Reorder features or samples by hierarchical clustering. Dendrograms
  are not drawn.

- column_split:

  Optional sample metadata column or complete sample-named vector used
  to split the x axis.

- feature_label:

  Optional feature metadata column or complete feature-named labels.
  Duplicate display labels are made unique.

## Value

An unprinted standard ggplot object carrying recommended physical
dimensions for
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md).

## See also

Other plotting:
[`plot_de_ma()`](https://younthing.github.io/bulkMAE/reference/plot_de_ma.md),
[`plot_de_volcano()`](https://younthing.github.io/bulkMAE/reference/plot_de_volcano.md),
[`plot_embedding()`](https://younthing.github.io/bulkMAE/reference/plot_embedding.md),
[`plot_gsea_classic()`](https://younthing.github.io/bulkMAE/reference/plot_gsea_classic.md),
[`plot_gsea_ridge()`](https://younthing.github.io/bulkMAE/reference/plot_gsea_ridge.md),
[`plot_ora_bubble()`](https://younthing.github.io/bulkMAE/reference/plot_ora_bubble.md),
[`plot_ora_network()`](https://younthing.github.io/bulkMAE/reference/plot_ora_network.md),
[`plot_ora_radial()`](https://younthing.github.io/bulkMAE/reference/plot_ora_radial.md),
[`plot_qc_correlation()`](https://younthing.github.io/bulkMAE/reference/plot_qc_correlation.md),
[`plot_qc_library()`](https://younthing.github.io/bulkMAE/reference/plot_qc_library.md),
[`plot_qc_outliers()`](https://younthing.github.io/bulkMAE/reference/plot_qc_outliers.md),
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md),
[`theme_bulkmae()`](https://younthing.github.io/bulkMAE/reference/theme_bulkmae.md)
