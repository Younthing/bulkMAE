# Convert GENIE3 weights to a ranked edge list

Convert GENIE3 weights to a ranked edge list

## Usage

``` r
network_genie3_links(weights, threshold = 0, top = NULL)
```

## Arguments

- weights:

  Weight matrix returned by
  [`network_genie3()`](https://younthing.github.io/bulkMAE/reference/network_genie3.md).

- threshold:

  Optional minimum edge weight.

- top:

  Optional maximum number of edges.

## Value

The native data frame returned by `GENIE3::getLinkList()`.
