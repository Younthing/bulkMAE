# Calculate sample-level library metrics

Calculate sample-level library metrics

## Usage

``` r
qc_library(x, experiment, assay = "counts", detected_above = 0)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- detected_above:

  A feature is detected when its value exceeds this threshold.

## Value

A data frame with one row per sample.
