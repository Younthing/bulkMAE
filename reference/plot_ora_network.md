# Plot an ORA term-feature community network

Builds a deterministic community layout from explicitly selected ORA
terms and their enriched-feature membership. Term-to-term edges
represent the Jaccard coefficient of the complete enriched-feature sets,
not ontology or full-pathway similarity. A shared feature is drawn once
inside every selected term community that contains it, matching the
reference community grammar. Optional feature values control every
corresponding visual node's size and alpha by absolute magnitude. Values
do not affect term selection.

## Usage

``` r
plot_ora_network(
  result,
  terms,
  membership = NULL,
  feature_values = NULL,
  features = NULL,
  term_labels = NULL,
  feature_labels = NULL,
  label_features = NULL,
  p_value = c("auto", "adjusted", "raw")
)
```

## Arguments

- result:

  A clusterProfiler `enrichResult` or compatible ORA result data frame.
  See
  [`plot_ora_bubble()`](https://younthing.github.io/bulkMAE/reference/plot_ora_bubble.md).

- terms:

  Unique term identifiers to display, in the desired order.

- membership:

  Optional term-named list of enriched feature identifiers. When
  omitted, membership is read from the result's `geneID` column. An
  explicit list must describe the selected terms exactly.

- feature_values:

  Optional feature-named finite numeric vector. It must cover every
  displayed feature; additional values are ignored.

- features:

  Optional explicit feature identifiers to display. Supply a character
  vector for one global feature set, or a term-named list to select
  displayed membership independently within each selected term. `NULL`
  displays all enriched features. Selection never changes overlap
  statistics computed from the complete membership.

- term_labels:

  Optional complete term-named labels.

- feature_labels:

  Optional feature-named labels covering every displayed feature.
  Additional labels are ignored.

- label_features:

  Explicit displayed feature identifiers to label. `NULL` labels no
  feature nodes.

- p_value:

  Evidence column used for term-node size.

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
[`plot_lincs_overlap()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_overlap.md),
[`plot_lincs_rank()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_rank.md),
[`plot_lincs_waterfall()`](https://younthing.github.io/bulkMAE/reference/plot_lincs_waterfall.md),
[`plot_ora_bubble()`](https://younthing.github.io/bulkMAE/reference/plot_ora_bubble.md),
[`plot_ora_radial()`](https://younthing.github.io/bulkMAE/reference/plot_ora_radial.md),
[`plot_qc_correlation()`](https://younthing.github.io/bulkMAE/reference/plot_qc_correlation.md),
[`plot_qc_library()`](https://younthing.github.io/bulkMAE/reference/plot_qc_library.md),
[`plot_qc_outliers()`](https://younthing.github.io/bulkMAE/reference/plot_qc_outliers.md),
[`plot_save()`](https://younthing.github.io/bulkMAE/reference/plot_save.md),
[`theme_bulkmae()`](https://younthing.github.io/bulkMAE/reference/theme_bulkmae.md)
