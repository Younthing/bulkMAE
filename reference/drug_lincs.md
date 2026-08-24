# Search a LINCS reference database for connected signatures

Search a LINCS reference database for connected signatures

## Usage

``` r
drug_lincs(
  query,
  reference_database,
  workers = 1,
  sort_by = "NCS",
  tau = FALSE,
  annotations = TRUE,
  ...
)
```

## Arguments

- query:

  A list with `upset` and `downset`, usually from
  [`drug_query()`](https://younthing.github.io/bulkMAE/reference/drug_query.md).

- reference_database:

  Path or identifier accepted by `signatureSearch::qSig()`.

- workers:

  Number of workers.

- sort_by:

  LINCS score used to rank results.

- tau:

  Calculate the standardized Tau score.

- annotations:

  Add compound annotations when available.

- ...:

  Additional arguments passed to `signatureSearch::gess_lincs()`.

## Value

A native `gessResult` object.
