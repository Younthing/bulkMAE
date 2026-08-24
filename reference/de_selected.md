# Select differential features as a named logical vector

The returned vector covers every tested feature and can be passed
directly to
[`enrich_goseq()`](https://younthing.github.io/bulkMAE/reference/enrich_goseq.md).
Its names also define the tested universe for ORA.

## Usage

``` r
de_selected(
  result,
  fdr = 0.05,
  min_abs_effect = 0,
  direction = c("both", "up", "down"),
  ...
)
```

## Arguments

- result:

  Input accepted by
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md).

- fdr:

  Maximum adjusted p-value.

- min_abs_effect:

  Minimum absolute effect size.

- direction:

  Select both directions, positive effects, or negative effects.

- ...:

  Additional arguments passed to
  [`de_table()`](https://younthing.github.io/bulkMAE/reference/de_table.md).

## Value

A named logical vector over all result features.
