# Construct a MultiAssayExperiment

Construct a MultiAssayExperiment

## Usage

``` r
mae_create(experiments, col_data, sample_map = NULL)
```

## Arguments

- experiments:

  A named list or `ExperimentList` of experiment leaves.

- col_data:

  Primary sample metadata with unique row names.

- sample_map:

  Optional explicit MAE sample map. When omitted, experiment column
  names must map directly to primary sample names.

## Value

A
[MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html).

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
rna <- mae_create_experiment(counts, samples)
mae_create(list(rna = rna), samples)
#> A MultiAssayExperiment object of 1 listed
#>  experiment with a user-defined name and respective class.
#>  Containing an ExperimentList class object of length 1:
#>  [1] rna: SummarizedExperiment with 3 rows and 2 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
