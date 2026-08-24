# Calculate consensus-clustering stability diagnostics

Calculate consensus-clustering stability diagnostics

## Usage

``` r
cluster_consensus_diagnostics(results, title = tempdir(), plot = NULL, ...)
```

## Arguments

- results:

  Result list returned by
  [`cluster_consensus()`](https://younthing.github.io/bulkMAE/reference/cluster_consensus.md).

- title:

  Output title or directory used by ConsensusClusterPlus.

- plot:

  Plot format accepted by `ConsensusClusterPlus::calcICL()`, or `NULL`
  to suppress plots.

- ...:

  Additional arguments passed to `ConsensusClusterPlus::calcICL()`.

## Value

The native list containing cluster and item consensus tables.
