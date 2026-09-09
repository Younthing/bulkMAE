# Run preranked gene-set enrichment with fgsea

Run preranked gene-set enrichment with fgsea

## Usage

``` r
enrich_fgsea(ranks, pathways, min_size = 15, max_size = 500, ...)
```

## Arguments

- ranks:

  Named numeric gene-level statistics.

- pathways:

  Named list of gene sets.

- min_size, max_size:

  Gene-set size limits.

- ...:

  Additional arguments passed to
  [`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html).

## Value

A native fgsea result table.
