# Fit a limma model to continuous expression or score data

Use this adapter for microarray intensities, log-expression assays,
pathway scores, or other approximately continuous outcomes. For counts,
use
[`de_voom()`](https://younthing.github.io/bulkMAE/reference/de_voom.md)
instead.

## Usage

``` r
de_limma(x, experiment, formula, assay, contrasts = NULL, trend = FALSE, ...)
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

- contrasts:

  Optional contrast matrix passed to `limma::contrasts.fit()`.

- trend:

  Passed to `limma::eBayes()`.

- ...:

  Additional arguments passed to `limma::eBayes()`.

## Value

A native `MArrayLM` object.
