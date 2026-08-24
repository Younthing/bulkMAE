# Convert native differential-expression results to a stable table

The helper does not replace native result objects. It extracts the small
set of fields needed by downstream enrichment, drug-query, and
meta-analysis functions while retaining feature identifiers explicitly.
For edgeR likelihood-ratio and quasi-likelihood F tests, `statistic` is
a directional signed square-root score for ranking only. Such tests do
not expose a coefficient standard error, so `standard_error` remains
`NA` and
[`meta_collect()`](https://younthing.github.io/bulkMAE/reference/meta_collect.md)
rejects them instead of manufacturing an invalid value. A DESeq2
likelihood-ratio statistic is retained in the table but marked as
non-directional;
[`de_ranks()`](https://younthing.github.io/bulkMAE/reference/de_ranks.md)
requires `column = "effect"` for that result, or a Wald result when a
directional test statistic is desired. Generic data-frame `F` and `LR`
columns are likewise retained as non-directional because their degrees
of freedom cannot be inferred safely.

## Usage

``` r
de_table(
  result,
  coef = NULL,
  feature_column = NULL,
  effect_column = NULL,
  standard_error_column = NULL,
  statistic_column = NULL,
  p_value_column = NULL,
  adjusted_p_value_column = NULL
)
```

## Arguments

- result:

  A `DESeqResults`, edgeR `DGELRT`/`DGEExact`, limma/dream `MArrayLM`,
  or data-frame-like result.

- coef:

  Coefficient name or index for an `MArrayLM` with more than one
  coefficient.

- feature_column, effect_column, standard_error_column,
  statistic_column:

  Optional feature, effect, standard-error, and statistic column names
  for a generic data frame. Common backend names are detected
  automatically.

- p_value_column, adjusted_p_value_column:

  Optional raw and adjusted p-value column names for a generic data
  frame. Common backend names are detected automatically.

## Value

A data frame with `feature_id`, `effect`, `standard_error`, `statistic`,
`p_value`, and `adjusted_p_value` columns.
