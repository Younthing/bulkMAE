# Fit a repeated-measures model with dream

Fit a repeated-measures model with dream

## Usage

``` r
de_dream(
  x,
  experiment,
  formula,
  assay = "counts",
  filter = TRUE,
  normalize_method = "TMM",
  contrast = NULL,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Formula containing fixed and optional random effects.

- assay:

  Assay name or one-based assay index.

- filter:

  Remove low-expression features with
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).

- normalize_method:

  Library normalization method.

- contrast:

  Optional contrast matrix passed as `L` to
  `variancePartition::dream()`.

- ...:

  Additional arguments passed to `variancePartition::dream()`.

## Value

A native, empirical-Bayes moderated `MArrayLM` object.
