# Filter lowly expressed features with edgeR

Filter lowly expressed features with edgeR

## Usage

``` r
filter_expr(
  x,
  experiment,
  assay = "counts",
  group = NULL,
  design = NULL,
  min_count = 10,
  min_total_count = 15,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- group:

  Optional metadata column or grouping vector.

- design:

  Optional design matrix. Supply either `group` or `design`.

- min_count:

  Minimum count passed to `edgeR::filterByExpr()`.

- min_total_count:

  Minimum total count passed to `edgeR::filterByExpr()`.

- ...:

  Additional arguments passed to `edgeR::filterByExpr()`.

## Value

A named logical vector, one value per feature.
