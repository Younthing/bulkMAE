# Add an aligned assay to one experiment

This pure helper returns a modified copy of `x`. It accepts a matrix, a
single-assay `SummarizedExperiment` such as a DESeq2 transformation, or
a limma `EList` containing an `E` matrix. Feature and sample identifiers
are aligned by name before insertion.

## Usage

``` r
mae_add_assay(x, experiment, value, name, overwrite = FALSE)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- value:

  Matrix-like values, a single-assay `SummarizedExperiment`, or an
  `EList` with an `E` matrix.

- name:

  Name for the new assay.

- overwrite:

  Replace an existing assay of the same name.

## Value

A `MultiAssayExperiment` copy containing the added assay.
