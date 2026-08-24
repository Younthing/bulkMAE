# List supported signatureSearch reference databases

The databases are distributed through Bioconductor ExperimentHub and are
downloaded into its local cache when first requested. This function
reports metadata only; it does not download a database.

## Usage

``` r
drug_lincs_databases()
```

## Value

A data frame describing supported database identifiers, human Entrez
gene IDs, stored value types, ExperimentHub records, and access mode.
