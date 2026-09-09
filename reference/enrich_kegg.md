# Run KEGG enrichment

Run KEGG enrichment

## Usage

``` r
enrich_kegg(
  genes = NULL,
  ranks = NULL,
  method = c("ora", "gsea"),
  organism = "hsa",
  key_type = "kegg",
  universe = NULL,
  ...
)
```

## Arguments

- genes:

  Gene identifiers of interest. Required only for `method = "ora"`.

- ranks:

  Named numeric gene-level statistics. Required only for
  `method = "gsea"`.

- method:

  Enrichment method: over-representation analysis (`ora`) or preranked
  GSEA (`gsea`).

- organism:

  KEGG organism code, for example `hsa` or `mmu`.

- key_type:

  Identifier type accepted by the clusterProfiler backend.

- universe:

  Optional background gene universe for ORA.

- ...:

  Additional arguments passed to
  [`clusterProfiler::enrichKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/enrichKEGG.html)
  or
  [`clusterProfiler::gseKEGG()`](https://rdrr.io/pkg/clusterProfiler/man/gseKEGG.html).

## Value

A native `enrichResult` or `gseaResult` object.
