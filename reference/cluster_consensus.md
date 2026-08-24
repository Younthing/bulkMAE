# Run consensus clustering

Run consensus clustering

## Usage

``` r
cluster_consensus(
  x,
  experiment,
  assay,
  max_k = 6,
  repetitions = 1000,
  item_fraction = 0.8,
  feature_fraction = 1,
  algorithm = "hc",
  distance = "pearson",
  seed = 1,
  title = tempdir(),
  plot = NULL,
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

- max_k:

  Maximum number of clusters.

- repetitions:

  Number of resampling repetitions.

- item_fraction, feature_fraction:

  Sampling fractions.

- algorithm:

  Clustering algorithm.

- distance:

  Distance metric.

- seed:

  Random seed.

- title:

  Output title or directory used by the backend.

- plot:

  Plot format accepted by ConsensusClusterPlus, or `NULL`.

- ...:

  Additional arguments passed to
  `ConsensusClusterPlus::ConsensusClusterPlus()`.

## Value

The native ConsensusClusterPlus result list.
