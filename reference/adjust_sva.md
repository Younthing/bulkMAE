# Estimate surrogate variables with sva

Estimate surrogate variables with sva

## Usage

``` r
adjust_sva(
  x,
  experiment,
  assay,
  full,
  null = ~1,
  n_surrogates = NULL,
  method = "irw",
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- full:

  Full biological model formula.

- null:

  Null model formula.

- n_surrogates:

  Optional number of surrogate variables. When omitted, `sva::num.sv()`
  is used.

- method:

  Method passed to `sva::sva()`.

- ...:

  Additional arguments passed to `sva::sva()`.

## Value

A native `sva` result list.
