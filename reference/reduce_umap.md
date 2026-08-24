# Run UMAP on samples

Run UMAP on samples

## Usage

``` r
reduce_umap(
  x,
  experiment,
  assay,
  features = NULL,
  top_n = NULL,
  neighbors = NULL,
  components = 2,
  metric = "euclidean",
  min_distance = 0.01,
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

- features:

  Optional feature subset.

- top_n:

  Optional number of most variable features to retain.

- neighbors:

  Number of nearest neighbors. When `NULL`, uses the smaller of 15 and
  the number of samples minus one.

- components:

  Embedding dimensions.

- metric:

  Distance metric.

- min_distance:

  Minimum embedding distance.

- seed:

  Optional random seed.

- ...:

  Additional arguments passed to `uwot::umap()`.

## Value

The native sample-by-component matrix returned by uwot.
