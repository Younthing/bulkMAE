# Plot a differential-expression volcano plot

Plot a differential-expression volcano plot

## Usage

``` r
plot_de_volcano(
  result,
  coef = NULL,
  fdr = 0.05,
  min_abs_effect = 1,
  label_features = NULL,
  feature_labels = NULL,
  symmetric_x = TRUE
)
```

## Arguments

- result:

  Input accepted by
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md).

- coef:

  Optional coefficient passed to
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md).

- fdr:

  Maximum adjusted p-value.

- min_abs_effect:

  Minimum absolute effect estimate.

- label_features:

  Feature identifiers to label; no features are selected for labels
  automatically.

- feature_labels:

  Optional complete feature-named labels or a named feature-label column
  in a data frame result.

- symmetric_x:

  Centre the x axis on zero with equal negative and positive limits. Set
  to `FALSE` to use ggplot2's automatic x range.

## Value

An unprinted standard ggplot object carrying recommended physical
dimensions for
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md).

## See also

Other plotting:
[`plot_assay_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_assay_heatmap.md),
[`plot_de_ma()`](https://younthing.github.io/bulkMAE/reference/plot_de_ma.md),
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
