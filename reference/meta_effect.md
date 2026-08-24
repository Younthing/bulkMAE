# Fit a univariate or meta-regression model

Fit a univariate or meta-regression model

## Usage

``` r
meta_effect(
  effects,
  standard_errors = NULL,
  moderators = NULL,
  method = "REML",
  ...
)
```

## Arguments

- effects:

  Numeric effect estimates, or a data frame containing `effect` and
  `standard_error` columns, such as the output of
  [`meta_collect()`](https://younthing.github.io/bulkMAE/reference/meta_collect.md).

- standard_errors:

  Numeric standard errors. Leave `NULL` when `effects` is a data frame.

- moderators:

  Optional moderator matrix or formula accepted by `metafor::rma.uni()`.

- method:

  Meta-analysis estimator.

- ...:

  Additional arguments passed to `metafor::rma.uni()`.

## Value

A native `rma.uni` fit.
