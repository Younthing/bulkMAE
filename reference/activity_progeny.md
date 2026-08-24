# Infer pathway activity with PROGENy and decoupleR

Infer pathway activity with PROGENy and decoupleR

## Usage

``` r
activity_progeny(
  x,
  experiment,
  assay,
  organism = "human",
  top = 500,
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

  Organism accepted by `decoupleR::get_progeny()`.

- top:

  Number of responsive genes retained per pathway.

- min_size:

  Minimum regulon size passed to `decoupleR::run_mlm()`.

- ...:

  Additional arguments passed to `decoupleR::run_mlm()`.

## Value

The native decoupleR result table.
