# Look up the 0.4 name for a pre-0.4 function

Version 0.4 removed the old exported aliases so autocomplete shows one
API. Use this helper when an old script still calls a 0.3 name.

## Usage

``` r
bulkmae_rename(old_name = NULL)
```

## Arguments

- old_name:

  A single pre-0.4 function name, with or without `()`. Use `NULL` to
  return the full map.

## Value

A replacement string such as
[`score_gsva()`](https://younthing.github.io/bulkMAE/reference/score_gsva.md)
or `enrich_go(method = "ora")`. If `old_name` is `NULL`, a data frame
with `old_name` and `new_name` columns.

## Examples

``` r
bulkmae_rename("run_gsva")
#> [1] "score_gsva()"
bulkmae_rename("run_go_ora")
#> [1] "enrich_go(method = \"ora\")"
```
