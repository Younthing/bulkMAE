# Extract an experiment aligned to primary sample data

Uses
[`MultiAssayExperiment::getWithColData()`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment-helpers.html)
so replicated or differently named assay columns are mapped through
`sampleMap` before analysis.

## Usage

``` r
mae_pull_experiment(x, experiment)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

## Value

A `SummarizedExperiment` with aligned `colData`.
