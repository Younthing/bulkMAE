# Detect co-expression modules with WGCNA

Detect co-expression modules with WGCNA

## Usage

``` r
coexpr_wgcna(
  x,
  experiment,
  assay,
  power,
  top_n = NULL,
  network_type = "signed",
  tom_type = NULL,
  min_module_size = 30,
  merge_cut_height = 0.25,
  numeric_labels = TRUE,
  seed = 1,
  verbose = 2,
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

- power:

  Soft-thresholding power chosen for the dataset.

- top_n:

  Optional number of most variable genes.

- network_type:

  WGCNA network type. Use the same value for module preservation.

- tom_type:

  Optional topological-overlap type. By default it follows
  `network_type` (`unsigned` stays unsigned; signed variants use signed
  TOM).

- min_module_size:

  Minimum module size.

- merge_cut_height:

  Module merging threshold.

- numeric_labels:

  Return numeric module labels.

- seed:

  Random seed used by WGCNA without changing the caller's random number
  generator state.

- verbose:

  Backend verbosity.

- ...:

  Additional arguments passed to `WGCNA::blockwiseModules()`.

## Value

The native WGCNA module result list, augmented with `bulkMAEQuality`,
`bulkMAERemovedSamples`, and `bulkMAERemovedFeatures` elements
documenting automatic quality filtering.
