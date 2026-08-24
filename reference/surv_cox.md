# Fit a Cox proportional-hazards model

Fit a Cox proportional-hazards model

## Usage

``` r
surv_cox(x, experiment, formula, assay = NULL, features = NULL, ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Cox model formula. Clinical variables are read from MAE sample
  metadata.

- assay:

  Assay name or one-based assay index.

- features:

  Optional expression features appended to the model data.

- ...:

  Additional arguments passed to
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).

## Value

A native `coxph` fit.
