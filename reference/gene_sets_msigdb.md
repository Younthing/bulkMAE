# Retrieve MSigDB gene sets through msigdbr

Retrieve MSigDB gene sets through msigdbr

## Usage

``` r
gene_sets_msigdb(
  species = "human",
  collection = "H",
  subcollection = NULL,
  id_type = c("gene_symbol", "ncbi_gene", "ensembl_gene"),
  db_species = "HS",
  ...
)
```

## Arguments

- species:

  Target species accepted by `msigdbr::msigdbr()`.

- collection:

  MSigDB collection code.

- subcollection:

  Optional MSigDB subcollection code.

- id_type:

  Identifier column returned in each set.

- db_species:

  MSigDB source database species.

- ...:

  Additional arguments passed to `msigdbr::msigdbr()`.

## Value

A named gene-set list with MSigDB source, version, species, collection,
subcollection, and identifier metadata stored as attributes.
