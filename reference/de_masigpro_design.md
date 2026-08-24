# Build a maSigPro experimental-design table from sample metadata

Build a maSigPro experimental-design table from sample metadata

## Usage

``` r
de_masigpro_design(x, experiment, time, replicate, group)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- time:

  Metadata column containing numeric time or dose.

- replicate:

  Metadata column identifying biological replicates.

- group:

  One categorical metadata column, or multiple existing binary
  indicator-column names.

## Value

A sample-aligned data frame containing `Time`, `Replicates`, and one
binary column per group.
