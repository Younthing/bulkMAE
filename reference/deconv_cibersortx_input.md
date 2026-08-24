# Prepare a CIBERSORTx mixture table

CIBERSORTx is an external service or container rather than an R backend.
This helper prepares its required gene-by-sample table without uploading
or writing data.

## Usage

``` r
deconv_cibersortx_input(
  x,
  experiment,
  assay = "tpm",
  gene_column = "GeneSymbol"
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- gene_column:

  Name of the first identifier column.

## Value

A data frame ready to write as a tab-delimited mixture file.
