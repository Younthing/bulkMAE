# Fit a cross-validated penalized Cox model

Fit a cross-validated penalized Cox model

## Usage

``` r
surv_penalized(
  x,
  experiment,
  time,
  event,
  assay,
  features = NULL,
  alpha = 1,
  folds = 10,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- time, event:

  Metadata columns containing follow-up time and event.

- assay:

  Assay name or one-based assay index.

- features:

  Optional expression feature subset.

- alpha:

  Elastic-net mixing parameter.

- folds:

  Number of cross-validation folds.

- seed:

  Optional local random seed used while constructing folds.

- ...:

  Additional arguments passed to `glmnet::cv.glmnet()`.

## Value

A native `cv.glmnet` object.
