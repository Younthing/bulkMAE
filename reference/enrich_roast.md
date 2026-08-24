# Run mroast rotation gene-set testing

Run mroast rotation gene-set testing

## Usage

``` r
enrich_roast(
  x,
  experiment,
  gene_sets,
  formula,
  contrast,
  assay,
  weights = NULL,
  rotations = 1999,
  seed = NULL,
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

  Contrast accepted by `limma::camera()`.

- assay:

  Assay name or one-based assay index.

- weights:

  Optional observation-level precision-weight matrix, or the name of an
  assay containing one. It must have the same dimensions and dimnames as
  `assay`. When omitted, no voom precision weights are used.

- rotations:

  Number of rotations passed as `nrot`.

- seed:

  Optional random seed.

- ...:

  Additional arguments passed to `limma::mroast()`.

## Value

The native mroast result data frame.
