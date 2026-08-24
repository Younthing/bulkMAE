# Validate the package's MultiAssayExperiment contract

Checks only the structural assumptions shared by the stateless wrappers:
`x` is a
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html),
the selected experiment exists, and its columns can be aligned to
primary samples.

## Usage

``` r
mae_validate(x, experiment = NULL)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

## Value

`x`, invisibly.
