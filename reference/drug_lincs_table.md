# Extract a ranked LINCS connectivity table

Converts a native signatureSearch `gessResult` from
[`drug_lincs()`](https://younthing.github.io/bulkMAE/reference/drug_lincs.md),
or a data frame with the same columns, into a ranked table. The helper
does not recompute weighted connectivity scores. The score column used
for ranking is the one
[`drug_lincs()`](https://younthing.github.io/bulkMAE/reference/drug_lincs.md)
actually stored (`NCS` by default; `Tau` only when it was requested and
is finite).

## Usage

``` r
drug_lincs_table(result, score = NULL)
```

## Arguments

- result:

  A `gessResult` or a data frame with LINCS/CMap-style columns.

- score:

  Optional score column to rank by. When `NULL`, prefer a finite `NCS`
  column, then `Tau`, `WTCS`, or `NCSct`.

## Value

A data frame ranked by the resolved score, with `score_column`,
`score_label`, and `direction` attached. `direction` is the sign of the
score (`Reverse`, `Mimic`, or `Unrelated`), not a wet-lab result.
