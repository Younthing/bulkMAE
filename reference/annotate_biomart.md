# Query Ensembl through biomaRt

Query Ensembl through biomaRt

## Usage

``` r
annotate_biomart(
  values,
  attributes,
  filter,
  dataset,
  biomart = "genes",
  mirror = NULL,
  version = NULL
)
```

## Arguments

- values:

  Values supplied to the selected filter.

- attributes:

  Ensembl attributes to return.

- filter:

  Ensembl filter name.

- dataset:

  Ensembl dataset, for example `hsapiens_gene_ensembl`.

- biomart:

  Ensembl BioMart name.

- mirror:

  Optional Ensembl mirror accepted by `biomaRt::useEnsembl()`.

- version:

  Optional archived Ensembl version.

## Value

A data frame returned by `biomaRt::getBM()`.
