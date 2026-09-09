# Run over-representation analysis against arbitrary gene sets

Run over-representation analysis against arbitrary gene sets

## Usage

``` r
enrich_ora(genes, gene_sets, universe = NULL, ...)
```

## Arguments

- genes:

  Gene identifiers of interest.

- gene_sets:

  Named list of gene vectors, or a two-column term-to-gene data frame.

- universe:

  Optional background gene universe.

- ...:

  Additional arguments passed to
  [`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).

## Value

A native `enrichResult` object.
