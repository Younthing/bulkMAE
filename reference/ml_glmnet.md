# Fit a cross-validated glmnet prediction model

Fit a cross-validated glmnet prediction model

## Usage

``` r
ml_glmnet(
  x,
  experiment,
  outcome,
  assay,
  features = NULL,
  family = "binomial",
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

- outcome:

  Metadata column containing the outcome.

- assay:

  Assay name or one-based assay index.

- features:

  Optional expression feature subset.

- family:

  glmnet model family.

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
