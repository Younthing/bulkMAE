# Plot a differential-expression MA plot

Plot a differential-expression MA plot

## Usage

``` r
plot_de_ma(
  result,
  coef = NULL,
  fdr = 0.05,
  min_abs_effect = 1,
  label_features = NULL,
  feature_labels = NULL,
  abundance = NULL,
  abundance_column = NULL,
  abundance_transform = c("auto", "identity", "log10")
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

- abundance:

  Optional complete feature-named abundance vector.

- abundance_column:

  Optional abundance column in a data-frame result.

- abundance_transform:

  Display abundance unchanged, on a base-10 log axis, or choose the
  backend convention automatically. Zero abundances are omitted from a
  log10 plot; negative abundances are rejected.

## Value

An unprinted standard ggplot object carrying recommended physical
dimensions for
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md).

## See also

Other plotting:
[`plot_assay_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_assay_heatmap.md),
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
