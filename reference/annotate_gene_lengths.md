# Retrieve representative transcript lengths from Ensembl

This online convenience function queries BioMart for transcript lengths
and summarizes multiple annotated transcripts per input gene. These
annotation lengths are not sample-specific effective lengths; quantified
RNA-seq data should preferentially use the `length` assay from
[`import_tximport()`](https://younthing.github.io/bulkMAE/reference/import_tximport.md).

## Usage

``` r
annotate_gene_lengths(
  ids,
  dataset = "hsapiens_gene_ensembl",
  id_type = "ensembl_gene_id",
  summary = c("median", "mean", "max"),
  strip_version = TRUE,
  biomart = "genes",
  mirror = NULL,
  version = NULL
)
```

## Arguments

- ids:

  Unique feature identifiers.

- dataset:

  Ensembl dataset, for example `hsapiens_gene_ensembl`.

- id_type:

  BioMart gene identifier attribute and filter.

- summary:

  How to summarize multiple positive transcript lengths per gene.

- strip_version:

  Remove terminal numeric Ensembl version suffixes before querying.

- biomart, mirror, version:

  Passed to
  [`annotate_biomart()`](https://younthing.github.io/bulkMAE/reference/annotate_biomart.md).

## Value

A feature-aligned data frame with `source_id`, `query_id`,
`gene_length`, and `n_transcripts`. Missing genes retain `NA` lengths.
