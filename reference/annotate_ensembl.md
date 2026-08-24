# Remove Ensembl version suffixes

Remove Ensembl version suffixes

## Usage

``` r
annotate_ensembl(ids)
```

## Arguments

- ids:

  Character vector of identifiers.

## Value

A character vector with a terminal dot and digits removed.

## Examples

``` r
annotate_ensembl(c("ENSG00000141510.18", "TP53"))
#> [1] "ENSG00000141510" "TP53"           
```
