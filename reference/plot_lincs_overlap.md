# Plot query up/down overlap counts for LINCS hits

Uses the `N_upset` and `N_downset` columns returned by
`signatureSearch::gess_lincs()`. Those counts are the query genes that
overlap the extreme ranks of each reference signature. They are not an
independent enrichment test.

## Usage

``` r
plot_lincs_overlap(result, n = 12L, score = NULL, compounds = NULL)
```

## Arguments

- result:

  A `gessResult` or LINCS table, usually from
  [`drug_lincs()`](https://younthing.github.io/bulkMAE/reference/drug_lincs.md)
  or
  [`drug_lincs_table()`](https://younthing.github.io/bulkMAE/reference/drug_lincs_table.md).

- n:

  Number of already-ranked rows to draw, taking the most extreme
  absolute scores when the table is longer.

- score:

  Optional score column forwarded to
  [`drug_lincs_table()`](https://younthing.github.io/bulkMAE/reference/drug_lincs_table.md).

- compounds:

  Optional compound identifiers. When supplied, rows are restricted to
  those names before ranking.

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
[`plot_lincs_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_heatmap.md),
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
