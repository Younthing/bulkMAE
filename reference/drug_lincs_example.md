# Diagnostic toy LINCS connectivity table

Returns a small table with the same columns as
`signatureSearch::gess_lincs()` (`pert`, `cell`, `type`, `trend`,
`WTCS`, `NCS`, `Tau`, `N_upset`, `N_downset`). Scores are synthetic.
They are not a download of CMap, LINCS L1000, or any commercial
reference, and they are not a wet-lab result.

## Usage

``` r
drug_lincs_example()
```

## Value

A data frame labelled `diagnostic_toy`.

## Details

Use this helper when no local HDF5 reference is available. The
production path remains
[`drug_lincs()`](https://younthing.github.io/bulkMAE/reference/drug_lincs.md)
against an explicit local or cached database.
