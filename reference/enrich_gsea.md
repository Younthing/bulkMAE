# Run preranked enrichment against arbitrary gene sets with clusterProfiler

Run preranked enrichment against arbitrary gene sets with
clusterProfiler

## Usage

``` r
enrich_gsea(ranks, gene_sets, ...)
```

## Arguments

- ranks:

  Named numeric gene-level statistics.

- gene_sets:

  Named list of gene sets, or a term-to-gene data frame.

- ...:

  Additional arguments passed to
  [`clusterProfiler::GSEA()`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html).

## Value

A native `gseaResult` object.
