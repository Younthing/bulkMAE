# Calculate GSVA scores

Uses the current parameter-object API: `GSVA::gsvaParam()` followed by
`GSVA::gsva()`.

## Usage

``` r
score_gsva(
  x,
  experiment,
  gene_sets,
  assay,
  kcdf = "auto",
  min_size = 1,
  max_size = Inf,
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

- kcdf:

  Kernel used by GSVA.

- min_size, max_size:

  Gene-set size limits.

- verbose:

  Show backend progress.

- ...:

  Additional arguments passed to `GSVA::gsvaParam()`.

## Value

The native matrix-like object returned by `GSVA::gsva()`.
