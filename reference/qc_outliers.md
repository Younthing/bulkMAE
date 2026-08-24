# Flag samples with unusually large robust PCA distances

This is a transparent screening heuristic, not a test that automatically
excludes samples. Distances are calculated from median/MAD-standardized
PC scores and compared with a chi-squared quantile.

## Usage

``` r
qc_outliers(
  pca,
  components = seq_len(min(5L, ncol(pca$x))),
  probability = 0.99
)
```

## Arguments

- pca:

  A `prcomp` object, usually from
  [`reduce_pca()`](https://younthing.github.io/bulkMAE/reference/reduce_pca.md).

- components:

  Principal components to include.

- probability:

  Chi-squared cutoff probability.

## Value

A data frame containing robust distance, cutoff, and flag.
