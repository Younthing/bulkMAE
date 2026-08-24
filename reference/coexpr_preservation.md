# Test preservation of reference WGCNA modules in another cohort

The two cohorts remain explicit MAE inputs. Features are intersected and
ordered identically before constructing the `multiData` and `multiColor`
lists required by `WGCNA::modulePreservation()`.

## Usage

``` r
coexpr_preservation(
  reference,
  reference_experiment,
  test,
  test_experiment,
  module_colors,
  reference_assay,
  test_assay,
  permutations = 200,
  network_type = "signed",
  seed = 1,
  verbose = 2,
  save_permuted_statistics = FALSE,
  ...
)
```

## Arguments

- reference, test:

  Reference and test `MultiAssayExperiment` objects.

- reference_experiment, test_experiment:

  Experiment names.

- module_colors:

  Named vector assigning every reference feature to a WGCNA module color
  or label.

- reference_assay, test_assay:

  Assays on comparable transformed scales.

- permutations:

  Number of permutations.

- network_type:

  WGCNA network type.

- seed:

  Random seed passed to WGCNA.

- verbose:

  Backend verbosity.

- save_permuted_statistics:

  Save permutation statistics to a file. The default avoids hidden
  filesystem output.

- ...:

  Additional arguments passed to `WGCNA::modulePreservation()`.

## Value

The native WGCNA module-preservation result list, augmented with a
`bulkMAEInputQuality` element documenting feature intersections and
automatic sample/feature filtering.
