# Run Gene Ontology enrichment

Run Gene Ontology enrichment

## Usage

``` r
enrich_go(
  genes = NULL,
  ranks = NULL,
  method = c("ora", "gsea"),
  org_db,
  key_type = "ENTREZID",
  ontology = "BP",
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

- org_db:

  An installed organism annotation database.

- key_type:

  Identifier type used in `genes` or `ranks`.

- ontology:

  GO ontology: `BP`, `MF`, `CC`, or `ALL`.

- universe:

  Optional background gene universe for ORA.

- ...:

  Additional arguments passed to `clusterProfiler::enrichGO()` or
  `clusterProfiler::gseGO()`.

## Value

A native `enrichResult` or `gseaResult` object.
