# Remove batch effects from continuous expression data with limma

This convenience method is intended for visualization or exploratory
analysis. Statistical models should usually include batch in their
design.

## Usage

``` r
adjust_batch(
  x,
  experiment,
  batch,
  assay,
  preserve = ~1,
  batch2 = NULL,
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

- preserve:

  Formula describing biological effects to preserve.

- batch2:

  Optional second batch column or vector.

- covariates:

  Optional numeric matrix, or metadata column names, for continuous
  nuisance covariates.

- ...:

  Additional arguments passed to `limma::removeBatchEffect()`.

## Value

A corrected feature-by-sample matrix.
