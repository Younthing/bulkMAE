# Extract aligned adjustment covariates

Supports an SVA result containing `sv`, an RUV matrix result containing
`W`, an RUVSeq `SeqExpressionSet` with `W_` phenotype columns, or an
already sample-by-covariate matrix/data frame. Sample names are never
inferred from row positions.

## Usage

``` r
adjust_covariates(result, columns = NULL)
```

## Arguments

- result:

  An adjustment result or sample-by-covariate object.

- columns:

  Optional unique column names or indices to retain.

## Value

A finite numeric data frame with unique sample row names.
