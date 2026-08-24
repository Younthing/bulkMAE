# Add an experiment while preserving the MAE sample map

A matrix is wrapped in a `SummarizedExperiment`. With neither
`source_experiment` nor `sample_map`, its column names must be primary
MAE sample identifiers. `source_experiment` copies an existing
experiment's map when the new value has the same assay-column
identifiers.

## Usage

``` r
mae_add_experiment(
  x,
  value,
  name,
  assay_name = "value",
  source_experiment = NULL,
  sample_map = NULL,
  overwrite = FALSE
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- value:

  A feature-by-sample matrix-like object or a `SummarizedExperiment`.

- name:

  Name for the new experiment.

- assay_name:

  Assay name used when `value` is a matrix-like object.

- source_experiment:

  Optional existing experiment whose sample map is copied for the new
  experiment.

- sample_map:

  Optional map for the new experiment, with `primary` and `colname`
  columns and an optional `assay` column equal to `name`.

- overwrite:

  Replace an existing experiment of the same name.

## Value

A `MultiAssayExperiment` copy containing the new experiment.
