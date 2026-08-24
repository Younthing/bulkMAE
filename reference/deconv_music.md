# Deconvolve bulk counts with MuSiC

Deconvolve bulk counts with MuSiC

## Usage

``` r
deconv_music(
  x,
  experiment,
  sc_reference,
  clusters = "cell_type",
  samples = "sample_id",
  assay = "counts",
  ...,
  reference_assay = NULL
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- sc_reference:

  Single-cell reference accepted by `MuSiC::music_prop()`.

- clusters:

  Cell-type annotation column in the reference. The default matches the
  canonical column created by
  [`deconv_reference()`](https://younthing.github.io/bulkMAE/reference/deconv_reference.md).

- samples:

  Sample identifier column in the reference. The default matches the
  canonical column created by
  [`deconv_reference()`](https://younthing.github.io/bulkMAE/reference/deconv_reference.md).

- assay:

  Assay name or one-based assay index.

- ...:

  Additional arguments passed to `MuSiC::music_prop()`.

- reference_assay:

  Raw-count assay in `sc_reference`. `NULL` selects `counts`, or the
  sole assay when no `counts` assay is present.

## Value

The native MuSiC result list. Its weighted estimates are cell-type
proportions; they are not on the same scale as xCell enrichment scores.
