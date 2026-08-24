# Construct a model design matrix from MAE sample metadata

Construct a model design matrix from MAE sample metadata

## Usage

``` r
de_design(x, experiment, formula)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Model formula evaluated against aligned sample metadata.

## Value

A numeric design matrix with rows aligned to assay samples.
