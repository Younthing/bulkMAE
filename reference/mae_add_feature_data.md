# Add aligned feature metadata to one experiment

This pure helper updates the raw experiment leaf and returns a modified
copy of `x`. It never changes feature identifiers or assay values.

## Usage

``` r
mae_add_feature_data(x, experiment, value, name = NULL, overwrite = FALSE)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- value:

  A feature-named atomic vector, or a feature-by-covariate data frame or
  matrix with feature row names.

- name:

  Optional name for a single metadata column. It is required for vector
  `value` and may rename a one-column data frame or matrix.

- overwrite:

  Replace metadata columns already present in the experiment.

## Value

A `MultiAssayExperiment` copy containing the feature metadata.
