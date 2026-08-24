# Apply the DESeq2 variance-stabilizing transformation

Apply the DESeq2 variance-stabilizing transformation

## Usage

``` r
transform_vst(
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

  Dispersion fit type passed to `DESeq2::vst()`.

- ...:

  Additional arguments passed to `DESeq2::vst()`.

## Value

A native `DESeqTransform` object.
