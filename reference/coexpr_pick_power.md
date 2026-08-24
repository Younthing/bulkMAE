# Evaluate candidate WGCNA soft-thresholding powers

Runs WGCNA's scale-free topology diagnostics directly from an MAE assay
and applies the same automatic sample/feature quality filtering used by
[`coexpr_wgcna()`](https://younthing.github.io/bulkMAE/reference/coexpr_wgcna.md).
The returned fit indices support, but do not replace, the scientific
choice of a soft-thresholding power.

## Usage

``` r
coexpr_pick_power(
  x,
  experiment,
  assay,
  powers = c(seq_len(10L), seq(12L, 20L, by = 2L)),
  top_n = NULL,
  network_type = "signed",
  r_squared = 0.85,
  verbose = 0,
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

- powers:

  Positive candidate powers.

- top_n:

  Optional number of most variable genes.

- network_type:

  WGCNA network type. Use the same value for module preservation.

- r_squared:

  Target scale-free topology fit passed as `RsquaredCut`.

- verbose:

  Backend verbosity.

- ...:

  Additional arguments passed to `WGCNA::pickSoftThreshold()`.

## Value

The native `pickSoftThreshold()` list, augmented with `bulkMAEQuality`,
`bulkMAERemovedSamples`, and `bulkMAERemovedFeatures` elements.
