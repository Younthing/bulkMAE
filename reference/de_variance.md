# Partition expression variance among model terms

Partition expression variance among model terms

## Usage

``` r
de_variance(x, experiment, formula, assay, ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Variance-partition formula. Random effects use `(1|group)`.

- assay:

  Assay name or one-based assay index.

- ...:

  Additional arguments passed to
  `variancePartition::fitExtractVarPartModel()`.

## Value

The native variance-partition result.
