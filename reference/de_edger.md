# Fit an edgeR quasi-likelihood model

Fit an edgeR quasi-likelihood model

## Usage

``` r
de_edger(
  x,
  experiment,
  formula,
  assay = "counts",
  filter = TRUE,
  normalize_method = "TMM",
  robust = TRUE,
  coef = NULL,
  contrast = NULL,
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

  Remove low-expression features with `edgeR::filterByExpr()`.

- normalize_method:

  Library normalization method.

- robust:

  Use robust empirical Bayes dispersion estimation.

- coef:

  Coefficient index or name for `edgeR::glmQLFTest()`.

- contrast:

  Optional numeric contrast for `edgeR::glmQLFTest()`.

- ...:

  Additional arguments passed to `edgeR::glmQLFit()`.

## Value

A native `DGEGLM` fit when neither `coef` nor `contrast` is given;
otherwise a native `DGELRT` test object.
