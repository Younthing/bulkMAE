# Run Reactome enrichment

Run Reactome enrichment

## Usage

``` r
enrich_reactome(
  genes = NULL,
  ranks = NULL,
  method = c("ora", "gsea"),
  organism = "human",
  universe = NULL,
  ...
)
```

## Arguments

- genes:

  Entrez gene identifiers of interest. Required only for
  `method = "ora"`.

- ranks:

  Named numeric Entrez-level statistics. Required only for
  `method = "gsea"`.

- method:

  Enrichment method: over-representation analysis (`ora`) or preranked
  GSEA (`gsea`).

- organism:

  Organism accepted by `ReactomePA::enrichPathway()`.

- universe:

  Optional background gene universe for ORA.

- ...:

  Additional arguments passed to `ReactomePA::enrichPathway()` or
  `ReactomePA::gsePathway()`.

## Value

A native `enrichResult` or `gseaResult` object.
