# Calculate TMM normalization factors

Calculate TMM normalization factors

## Usage

``` r
normalize_tmm(x, experiment, assay = "counts", method = "TMM", ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- assay:

  Assay name or one-based assay index.

- method:

  Normalization method passed to `edgeR::normLibSizes()`.

- ...:

  Additional arguments passed to `edgeR::normLibSizes()`.

## Value

A native edgeR::DGEList containing normalization factors.
