# Add aligned sample metadata to one experiment

This pure helper updates `colData` on the raw experiment leaf and
returns a modified copy of `x`. New columns may not duplicate primary
MAE metadata; `overwrite` applies only to columns already stored on the
raw leaf.

## Usage

``` r
mae_add_sample_data(x, experiment, value, name = NULL, overwrite = FALSE)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- value:

  A sample-named atomic vector, or a sample-by-covariate data frame or
  matrix with assay-sample row names.

- name:

  Optional name for a single metadata column. It is required for vector
  `value` and may rename a one-column data frame or matrix.

- overwrite:

  Replace metadata columns already present on the raw experiment leaf.

## Value

A `MultiAssayExperiment` copy containing the sample metadata.
