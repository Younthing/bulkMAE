# Collect one feature's effect estimates across differential analyses

Collect one feature's effect estimates across differential analyses

## Usage

``` r
meta_collect(results, feature, coef = NULL)
```

## Arguments

- results:

  A named list of inputs accepted by
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md).

- feature:

  One feature identifier present in every result.

- coef:

  Optional coefficient shared by all results, or a list with one
  coefficient per result.

## Value

A data frame with `study`, `feature_id`, `effect`, and `standard_error`,
accepted directly by
[`meta_effect()`](https://younthing.github.io/bulkMAE/reference/meta_effect.md).
