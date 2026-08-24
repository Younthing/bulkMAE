# Read gene sets from a GMT file

Read gene sets from a GMT file

## Usage

``` r
gene_sets_read_gmt(file, ...)
```

## Arguments

- file:

  Path to a Gene Matrix Transposed (`.gmt`) file.

- ...:

  Size filters passed to
  [`gene_sets_prepare()`](https://younthing.github.io/bulkMAE/reference/gene_sets_prepare.md).

## Value

A named gene-set list. Set descriptions and the normalized source path
are stored in `description` and `source_file` attributes.
