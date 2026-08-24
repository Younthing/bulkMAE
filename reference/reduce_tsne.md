# Run t-SNE on samples

Run t-SNE on samples

## Usage

``` r
reduce_tsne(
  x,
  experiment,
  assay,
  features = NULL,
  top_n = NULL,
  dimensions = 2,
  perplexity = NULL,
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

- dimensions:

  Embedding dimensions.

- perplexity:

  Perplexity passed to `Rtsne::Rtsne()`. When `NULL`, uses the largest
  safe integer no greater than 30.

- seed:

  Optional random seed.

- ...:

  Additional arguments passed to `Rtsne::Rtsne()`.

## Value

A native `Rtsne` result object.
