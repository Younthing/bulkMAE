# Infer a gene regulatory network with GENIE3

Infer a gene regulatory network with GENIE3

## Usage

``` r
network_genie3(
  x,
  experiment,
  assay,
  regulators = NULL,
  targets = NULL,
  trees = 1000,
  cores = 1,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- regulators:

  Optional regulator names or indices.

- targets:

  Optional target names or indices.

- trees:

  Number of trees.

- cores:

  Number of worker cores.

- seed:

  Optional random seed.

- ...:

  Additional arguments passed to `GENIE3::GENIE3()`.

## Value

A native GENIE3 weight matrix.
