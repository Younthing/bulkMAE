# Run non-negative matrix factorization

Run non-negative matrix factorization

## Usage

``` r
cluster_nmf(
  x,
  experiment,
  assay,
  rank,
  method = "brunet",
  runs = 30,
  seed = 1,
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

- rank:

  Factorization rank or vector of candidate ranks.

- method:

  NMF algorithm.

- runs:

  Number of runs.

- seed:

  Random seed or NMF seed method.

- ...:

  Additional arguments passed to `NMF::nmf()`.

## Value

A native NMF fit object.
