# Run sample principal component analysis

Run sample principal component analysis

## Usage

``` r
reduce_pca(
  x,
  experiment,
  assay,
  features = NULL,
  top_n = NULL,
  center = TRUE,
  scale. = FALSE
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- features:

  Optional feature subset.

- top_n:

  Optional number of most variable features to retain.

- center, scale.:

  Passed to [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html).

## Value

A native `prcomp` object.
