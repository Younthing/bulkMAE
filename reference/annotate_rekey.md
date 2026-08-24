# Apply an identifier mapping to common gene-level values

A source identifier that maps to multiple targets is expanded to every
target. When multiple input rows or elements then share a target, the
explicitly selected `duplicates` rule is applied. Exact duplicate
mapping rows are removed before expansion so they cannot multiply a
value.

## Usage

``` r
annotate_rekey(
  values,
  mapping,
  source = "source_id",
  target = "target_id",
  duplicates = c("sum", "mean", "max_abs", "first"),
  drop_unmapped = TRUE
)
```

## Arguments

- values:

  A numeric gene-by-sample matrix, a named numeric vector, or a
  character vector of identifiers.

- mapping:

  A data frame containing source and target identifier columns.

- source, target:

  Column names in `mapping`.

- duplicates:

  Aggregation rule for numeric values sharing a target.

- drop_unmapped:

  Drop source identifiers without a non-missing target. When `FALSE`,
  their original source identifiers are retained.

## Value

A target-keyed numeric matrix, named numeric vector, or character
vector, matching the input value type.
