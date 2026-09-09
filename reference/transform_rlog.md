# Apply the DESeq2 regularized-log transformation

Apply the DESeq2 regularized-log transformation

## Usage

``` r
transform_rlog(
  x,
  experiment,
  assay = "counts",
  blind = TRUE,
  design = ~1,
  fit_type = "parametric",
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

- blind:

  Whether the transformation should ignore the design.

- design:

  Design used when `blind = FALSE`.

- fit_type:

  Dispersion fit type passed to
  [`DESeq2::vst()`](https://rdrr.io/pkg/DESeq2/man/vst.html).

- ...:

  Additional arguments passed to
  [`DESeq2::vst()`](https://rdrr.io/pkg/DESeq2/man/vst.html).

## Value

A native `DESeqTransform` object.
