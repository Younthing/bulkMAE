# Adjust a continuous assay with ComBat

Adjust a continuous assay with ComBat

## Usage

``` r
adjust_combat(
  x,
  experiment,
  batch,
  assay,
  preserve = NULL,
  parametric = TRUE,
  mean_only = FALSE,
  reference_batch = NULL,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- batch:

  Metadata column name or batch vector.

- assay:

  Assay name or one-based assay index.

- preserve:

  Optional formula for biological covariates to preserve.

- parametric:

  Use parametric empirical Bayes priors.

- mean_only:

  Adjust batch means but not scales.

- reference_batch:

  Optional reference batch.

- ...:

  Additional arguments passed to `sva::ComBat()`.

## Value

A corrected feature-by-sample matrix.
