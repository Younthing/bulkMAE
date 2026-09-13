# Plot one classic GSEA running-score profile

Displays the running enrichment score, gene-set hit positions with a
rank heat band, and the ranked metric as three intrinsic parts of one
standard ggplot. The parts share a rank-position axis and use fixed
relative panel heights of 0.50, 0.20, and 0.30. This function does not
call private enrichplot helpers or rerun GSEA.

## Usage

``` r
plot_gsea_classic(
  result,
  term,
  ranks = NULL,
  gene_sets = NULL,
  exponent = NULL,
  term_label = NULL
)
```

## Arguments

- result:

  A native clusterProfiler gseaResult or a GSEA result table.

- term:

  One explicit term identifier to display.

- ranks:

  Optional decreasing, gene-named, finite numeric rank vector.

- gene_sets:

  Optional uniquely named list mapping terms to genes.

- exponent:

  Optional finite non-negative GSEA weighting exponent.

- term_label:

  Optional single non-empty display label for term.

## Value

An unprinted standard ggplot object carrying recommended physical
dimensions for
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md).

## Details

A native clusterProfiler gseaResult supplies its ranked vector, gene
sets, exponent, and leading-edge membership automatically. For an fgsea
or generic GSEA table, ranks, gene_sets, and exponent are all required.
Those inputs must reproduce the reported enrichment score.

## See also

Other plotting:
[`plot_assay_heatmap()`](https://younthing.github.io/bulkMAE/reference/plot_assay_heatmap.md),
[`plot_de_ma()`](https://younthing.github.io/bulkMAE/reference/plot_de_ma.md),
[`plot_de_volcano()`](https://younthing.github.io/bulkMAE/reference/plot_de_volcano.md),
[`plot_embedding()`](https://younthing.github.io/bulkMAE/reference/plot_embedding.md),
[`plot_gsea_ridge()`](https://younthing.github.io/bulkMAE/reference/plot_gsea_ridge.md),
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
