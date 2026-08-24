# Resolve a local organism annotation database

Resolve a local organism annotation database

## Usage

``` r
annotation_orgdb(organism = c("human", "mouse"), package = NULL)
```

## Arguments

- organism:

  Supported organism name.

- package:

  Optional installed `OrgDb` package name. When omitted, `org.Hs.eg.db`
  or `org.Mm.eg.db` is selected from `organism`.

## Value

The installed `OrgDb` object exported by `package`.
