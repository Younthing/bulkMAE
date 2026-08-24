# Run one or more decoupleR methods

This general adapter deliberately exposes both enrichment-style methods
(`aucell`, `fgsea`, `gsva`, and `ora`) and network-aware activity
methods (`mlm`, `ulm`, `viper`, `wmean`, and `wsum`). With
`statistics = NULL`, decoupleR runs its documented default methods and
can calculate a consensus.

## Usage

``` r
activity_decouple(
  x,
  experiment,
  network,
  assay,
  statistics = NULL,
  method_args = list(NULL),
  consensus = TRUE,
  consensus_statistics = NULL,
  source = "source",
  target = "target",
  mor = NULL,
  min_size = 5,
  include_time = FALSE,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- network:

  Long-format network or gene-set table.

- assay:

  Assay name or one-based assay index.

- statistics:

  Method names accepted by `decoupleR::decouple()`, or `NULL` for its
  defaults.

- method_args:

  List of method-specific argument lists.

- consensus:

  Calculate a consensus score.

- consensus_statistics:

  Optional decoupleR result-statistic names used for consensus scoring,
  for example `norm_mlm` and `norm_ulm`. `NULL` preserves the backend
  default.

- source, target:

  Column names identifying regulators/sets and targets.

- mor:

  Optional column containing signed interaction weights. It is mapped to
  decoupleR's conventional `mor` column without discarding the original
  column.

- min_size:

  Minimum number of targets per source.

- include_time:

  Include execution time in the result.

- ...:

  Additional arguments passed to `decoupleR::decouple()`.

## Value

The native long-format table returned by `decoupleR::decouple()`.
