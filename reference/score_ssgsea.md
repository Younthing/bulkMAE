# Calculate single-sample GSEA scores

Uses the current `GSVA::ssgseaParam()` parameter-object API.

## Usage

``` r
score_ssgsea(
  x,
  experiment,
  gene_sets,
  assay,
  min_size = 1,
  max_size = Inf,
  normalize = TRUE,
  alpha = 0.25,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- gene_sets:

  Named list of gene sets.

- assay:

  Assay name or one-based assay index.

- min_size, max_size:

  Gene-set size limits.

- normalize:

  Normalize ssGSEA scores by their range.

- alpha:

  Tail-weight exponent used by ssGSEA.

- verbose:

  Show backend progress.

- ...:

  Additional arguments passed to `GSVA::gsvaParam()`.

## Value

The native matrix-like object returned by `GSVA::gsva()`.
