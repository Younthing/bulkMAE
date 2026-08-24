# Deconvolve bulk expression with immunedeconv

The selected assay must contain non-log TPM-like expression with HGNC
gene symbols as row names, as required by immunedeconv. xCell returns
relative enrichment scores rather than cell fractions. EPIC and
quanTIseq produce proportion-like estimates, but their output scales and
supported cell types remain method-specific and should not be pooled as
if they were interchangeable measurements.

## Usage

``` r
deconv(x, experiment, method, assay = "tpm", ...)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- method:

  Method supported by `immunedeconv::deconvolute()`, such as `epic`,
  `quantiseq`, or `xcell`.

- assay:

  Assay name or one-based assay index.

- ...:

  Additional arguments passed to `immunedeconv::deconvolute()`.

## Value

The native immunedeconv result table.
