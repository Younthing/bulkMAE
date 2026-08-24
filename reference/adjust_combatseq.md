# Adjust integer counts with ComBat-seq

Adjust integer counts with ComBat-seq

## Usage

``` r
adjust_combatseq(
  x,
  experiment,
  batch,
  assay = "counts",
  group = NULL,
  covariates = NULL,
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

- group:

  Optional biological group column or vector to preserve.

- covariates:

  Optional formula for additional covariates to preserve.

- ...:

  Additional arguments passed to `sva::ComBat_seq()`.

## Value

A corrected count matrix.
