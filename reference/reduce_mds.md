# Run sample multidimensional scaling

Uses [`limma::plotMDS()`](https://rdrr.io/pkg/limma/man/plotMDS.html)
with `plot = FALSE` and returns the native MDS object without opening a
graphics device.

## Usage

``` r
reduce_mds(
  x,
  experiment,
  assay,
  features = NULL,
  top = 500,
  gene_selection = "pairwise",
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- features:

  Optional feature subset.

- top:

  Number of leading features used by `plotMDS`.

- gene_selection:

  Feature-selection rule accepted by `plotMDS`.

- ...:

  Additional arguments passed to
  [`limma::plotMDS()`](https://rdrr.io/pkg/limma/man/plotMDS.html).

## Value

A native limma MDS object.
