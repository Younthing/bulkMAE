# Fit a limma-voom differential-expression model

Fit a limma-voom differential-expression model

## Usage

``` r
de_voom(
  x,
  experiment,
  formula,
  assay = "counts",
  filter = TRUE,
  normalize_method = "TMM",
  contrasts = NULL,
  trend = FALSE,
  voom_plot = FALSE,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- formula:

  Fixed-effect model formula.

- assay:

  Assay name or one-based assay index.

- filter:

  Remove low-expression features with
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).

- normalize_method:

  Library normalization method.

- contrasts:

  Optional contrast matrix passed to
  [`limma::contrasts.fit()`](https://rdrr.io/pkg/limma/man/contrasts.fit.html).

- trend:

  Passed to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html).

- voom_plot:

  Draw the voom mean-variance diagnostic plot.

- ...:

  Additional arguments passed to
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html).

## Value

A native `MArrayLM` object.
