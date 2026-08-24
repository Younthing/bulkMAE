# Construct a one-experiment MAE from a matrix

Construct a one-experiment MAE from a matrix

## Usage

``` r
mae_from_matrix(
  expression,
  samples,
  row_data = NULL,
  experiment = "rna",
  assay = "counts"
)
```

## Arguments

- expression:

  A feature-by-sample matrix-like expression object with unique feature
  and sample names.

- samples:

  Sample metadata whose row names exactly match the expression columns.

- row_data:

  Optional feature metadata whose row names exactly match the expression
  rows.

- experiment:

  Name of the experiment to create.

- assay:

  Name of the expression assay.

## Value

A
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html).
