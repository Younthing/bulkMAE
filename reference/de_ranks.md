# Extract finite named ranks from differential-expression results

Extract finite named ranks from differential-expression results

## Usage

``` r
de_ranks(result, column = c("statistic", "effect"), decreasing = TRUE, ...)
```

## Arguments

- result:

  Input accepted by
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md).

- column:

  Normalized
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md)
  column used for ranking.

- decreasing:

  Sort ranks from high to low.

- ...:

  Additional arguments passed to
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md).

## Value

A finite named numeric vector suitable for `enrich_*()` and
[`drug_query()`](https://younthing.github.io/bulkMAE/reference/drug_query.md).
