# Fit a DESeq2 differential-expression model

Fit a DESeq2 differential-expression model

## Usage

``` r
de_deseq2(
  x,
  experiment,
  design,
  assay = "counts",
  test = c("Wald", "LRT"),
  reduced = NULL,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- design:

  DESeq2 design formula.

- assay:

  Assay name or one-based assay index.

- test:

  Wald or likelihood-ratio test.

- reduced:

  Reduced formula required for a likelihood-ratio test.

- ...:

  Additional arguments passed to
  [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html).

## Value

A fitted native `DESeqDataSet`.
