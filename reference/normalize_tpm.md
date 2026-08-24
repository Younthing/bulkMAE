# Convert counts to transcripts per million

Gene lengths may be supplied directly, read from one `rowData` column,
or read from an aligned assay such as the sample-specific
effective-length assay produced by
[`import_tximport()`](https://younthing.github.io/bulkMAE/reference/import_tximport.md).
No external preprocessing step is required.

## Usage

``` r
normalize_tpm(x, experiment, lengths = NULL, assay = "counts")
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- lengths:

  A positive numeric vector, a feature-named numeric vector, a
  feature-by-sample matrix, or the name of a `rowData` column or assay.
  When `NULL`, an assay named `length` (as produced by
  [`import_tximport()`](https://younthing.github.io/bulkMAE/reference/import_tximport.md))
  is used.

- assay:

  Assay name or one-based assay index.

## Value

A finite feature-by-sample TPM matrix whose columns sum to one million.
