# Apply limma-voom transformation and precision weighting

Apply limma-voom transformation and precision weighting

## Usage

``` r
transform_voom(
  x,
  experiment,
  formula,
  assay = "counts",
  normalize_method = "TMM",
  plot = FALSE,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Fixed-effect model formula.

- assay:

  Assay name or one-based assay index.

- normalize_method:

  edgeR library normalization method.

- plot:

  Whether `voom` should draw its diagnostic plot.

- ...:

  Additional arguments passed to
  [`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html).

## Value

A native limma `EList`.
