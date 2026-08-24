# Extract sample or feature classes from an NMF fit

Extract sample or feature classes from an NMF fit

## Usage

``` r
cluster_nmf_classes(fit, what = c("samples", "features"))
```

## Arguments

- fit:

  A native result returned by
  [`cluster_nmf()`](https://younthing.github.io/bulkMAE/reference/cluster_nmf.md).

- what:

  Return sample or feature classes.

## Value

A named class vector from the NMF backend.
