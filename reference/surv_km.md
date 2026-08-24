# Fit Kaplan-Meier survival curves

Fit Kaplan-Meier survival curves

## Usage

``` r
surv_km(x, experiment, formula, ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Survival formula evaluated against aligned sample metadata.

- ...:

  Additional arguments passed to
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html).

## Value

A native `survfit` object.
