# Test differential exon or transcript usage from a fitted limma/edgeR model

Test differential exon or transcript usage from a fitted limma/edgeR
model

## Usage

``` r
dtu_diffsplice(
  fit,
  gene_id,
  feature_id = NULL,
  robust = FALSE,
  coef = NULL,
  contrast = NULL,
  ...
)
```

## Arguments

- fit:

  A native `MArrayLM` or `DGEGLM` model fit.

- gene_id:

  Gene identifier for every fitted feature.

- feature_id:

  Optional exon or transcript identifier.

- robust:

  Use robust empirical Bayes estimation.

- coef, contrast:

  Coefficient or contrast for the edgeR method. Supply at most one. For
  limma, apply contrasts to `fit` before calling this function; limma's
  `diffSplice()` calculates statistics for the fitted coefficients.

- ...:

  Additional arguments passed to the selected backend's `diffSplice()`
  method.

## Value

A native limma or edgeR differential-splicing result.

## Details

`gene_id` and `feature_id` are aligned to the fitted rows when they are
named. This is particularly useful after expression filtering. Supply an
untested `DGEGLM` to the edgeR method; choose its test with `coef` or
`contrast`. For limma, apply the desired contrast with
[`limma::contrasts.fit()`](https://rdrr.io/pkg/limma/man/contrasts.fit.html)
before calling this function when a contrast, rather than a coefficient,
is required.
