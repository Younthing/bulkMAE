# Infer transcription-factor activity with decoupleR

Infer transcription-factor activity with decoupleR

## Usage

``` r
activity_tf(
  x,
  experiment,
  assay,
  organism = "human",
  resource = c("collectri", "dorothea"),
  min_size = 5,
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

- organism:

  Organism accepted by the selected network resource.

- resource:

  TF network: current CollecTRI or DoRothEA.

- min_size:

  Minimum regulon size passed to `decoupleR::run_ulm()`.

- ...:

  Additional arguments passed to `decoupleR::run_ulm()`.

## Value

The native decoupleR result table.
