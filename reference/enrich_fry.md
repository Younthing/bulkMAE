# Run FRY rotation gene-set testing

Run FRY rotation gene-set testing

## Usage

``` r
enrich_fry(
  x,
  experiment,
  gene_sets,
  formula,
  contrast,
  assay,
  weights = NULL,
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

- formula:

  Model formula used to construct the design matrix.

- contrast:

  Contrast accepted by
  [`limma::camera()`](https://rdrr.io/pkg/limma/man/camera.html).

- assay:

  Assay name or one-based assay index.

- weights:

  Optional observation-level precision-weight matrix, or the name of an
  assay containing one. It must have the same dimensions and dimnames as
  `assay`. When omitted, no voom precision weights are used.

- ...:

  Additional arguments passed to
  [`limma::camera()`](https://rdrr.io/pkg/limma/man/camera.html).

## Value

The native FRY result data frame.
