# Estimate DESeq2 size factors

Estimate DESeq2 size factors

## Usage

``` r
normalize_deseq(x, experiment, assay = "counts", type = "ratio", ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- type:

  Size-factor estimator passed to `DESeq2::estimateSizeFactors()`.

- ...:

  Additional arguments passed to `DESeq2::estimateSizeFactors()`.

## Value

A native `DESeqDataSet` with estimated size factors.
