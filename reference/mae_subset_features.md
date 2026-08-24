# Subset features in one experiment

This pure helper returns a modified copy of `x`, leaving all other
experiments and the `sampleMap` unchanged.

## Usage

``` r
mae_subset_features(x, experiment, features)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- features:

  Unique feature names to retain, in the requested order.

## Value

A `MultiAssayExperiment` copy with one experiment subset by row.
