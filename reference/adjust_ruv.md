# Estimate unwanted factors with RUVg

Estimate unwanted factors with RUVg

## Usage

``` r
adjust_ruv(x, experiment, controls, k = 1, assay = "counts", ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- controls:

  Control-feature names or indices passed as `cIdx`.

- k:

  Number of unwanted factors.

- assay:

  Assay name or one-based assay index.

- ...:

  Additional arguments passed to `RUVSeq::RUVg()`.

## Value

The native object returned by `RUVSeq::RUVg()`.
