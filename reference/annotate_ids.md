# Map feature identifiers with an AnnotationDbi database

This is a thin wrapper around `AnnotationDbi::mapIds()`. It deliberately
returns a mapping table and does not alter row names in the MAE.

## Usage

``` r
annotate_ids(ids, database, from, to, multi_values = "first")
```

## Arguments

- ids:

  Character vector of source identifiers.

- database:

  An installed `AnnotationDb` object, such as `org.Hs.eg.db`.

- from:

  Source key type.

- to:

  Destination column.

- multi_values:

  Policy passed to `multiVals`.

## Value

A long data frame with atomic `source_id` and `target_id` columns. When
`multi_values` returns a list or `CharacterList`, one-to-many mappings
are expanded to repeated source rows so the result can be passed
directly to
[`annotate_rekey()`](https://younthing.github.io/bulkMAE/reference/annotate_rekey.md).
