# Standardize gene sets as a named list

Standardize gene sets as a named list

## Usage

``` r
gene_sets_prepare(x, term = NULL, gene = NULL, min_size = 1, max_size = Inf)
```

## Arguments

- x:

  A named list of gene vectors, or a long data frame. For a data frame,
  the first two columns are used when `term` and `gene` are omitted.

- term, gene:

  Column name or one-based column index for terms and genes.

- min_size, max_size:

  Inclusive unique-gene size limits.

## Value

A uniquely named list of unique character gene identifiers.
