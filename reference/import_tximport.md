# Import transcript abundance estimates with tximport

The returned object is an experiment leaf, not a stateful analysis
object. Add it to a `MultiAssayExperiment` with
[`mae_create()`](https://younthing.github.io/bulkMAE/reference/mae_create.md).

## Usage

``` r
import_tximport(
  files,
  col_data,
  type = "salmon",
  tx2gene = NULL,
  counts_from_abundance = "no",
  ...
)
```

## Arguments

- files:

  Named character vector of quantification files.

- col_data:

  Sample metadata; row names must equal `names(files)`.

- type:

  Quantifier type accepted by `tximport::tximport()`.

- tx2gene:

  Optional transcript-to-gene mapping.

- counts_from_abundance:

  Value passed to `countsFromAbundance`.

- ...:

  Additional arguments passed to `tximport::tximport()`.

## Value

A `SummarizedExperiment` with available `counts`, `abundance`, and
`length` assays.
