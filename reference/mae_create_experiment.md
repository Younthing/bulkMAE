# Construct a SummarizedExperiment leaf

Construct a SummarizedExperiment leaf

## Usage

``` r
mae_create_experiment(assays, col_data, row_data = NULL, assay_name = "counts")
```

## Arguments

- assays:

  A matrix, matrix-like object, or named list of assays.

- col_data:

  Sample metadata coercible to
  [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html).

- row_data:

  Optional feature metadata coercible to
  [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html).

- assay_name:

  Name used when `assays` is a single matrix.

## Value

A
[SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

## Examples

``` r
counts <- matrix(
  1:6,
  nrow = 3,
  dimnames = list(paste0("gene", 1:3), c("sample1", "sample2"))
)
samples <- data.frame(
  condition = c("control", "treated"),
  row.names = colnames(counts)
)
mae_create_experiment(counts, samples)
#> class: SummarizedExperiment 
#> dim: 3 2 
#> metadata(0):
#> assays(1): counts
#> rownames(3): gene1 gene2 gene3
#> rowData names(0):
#> colnames(2): sample1 sample2
#> colData names(1): condition
```
