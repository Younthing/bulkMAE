# Calculate sample-to-sample correlations

Calculate sample-to-sample correlations

## Usage

``` r
qc_correlation(
  x,
  experiment,
  assay,
  method = "pearson",
  use = "pairwise.complete.obs",
  features = NULL
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- method:

  Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

- use:

  Missing-value policy passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

- features:

  Optional feature subset.

## Value

A sample-by-sample correlation matrix.
