# Construct a namespace-safe survival formula

This convenience helper lets users name MAE metadata columns without
first attaching or calling the `survival` package. More complex
time-dependent or stratified formulas can still be supplied directly to
[`surv_km()`](https://younthing.github.io/bulkMAE/reference/surv_km.md)
and
[`surv_cox()`](https://younthing.github.io/bulkMAE/reference/surv_cox.md).

## Usage

``` r
surv_formula(time, event, predictors = NULL)
```

## Arguments

- time, event:

  Metadata column names containing follow-up time and event.

- predictors:

  Optional right-hand-side variable names. These must be sample-metadata
  columns for
  [`surv_km()`](https://younthing.github.io/bulkMAE/reference/surv_km.md).
  For
  [`surv_cox()`](https://younthing.github.io/bulkMAE/reference/surv_cox.md),
  expression features are also available when the same names are
  supplied through its `features` argument.

## Value

A formula with a
[`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
response.
