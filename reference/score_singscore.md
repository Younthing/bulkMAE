# Score samples with rank-based singscore

Score samples with rank-based singscore

## Usage

``` r
score_singscore(
  x,
  experiment,
  up_set,
  down_set = NULL,
  assay,
  ties_method = "min",
  known_direction = TRUE,
  ...
)
```

## Arguments

- x:

  A `MultiAssayExperiment` object.

- experiment:

  Optional experiment name.

- up_set:

  Features expected to be up-regulated.

- down_set:

  Optional features expected to be down-regulated.

- assay:

  Assay name or one-based assay index.

- ties_method:

  Ranking method passed to `singscore::rankGenes()`.

- known_direction:

  Whether signature direction is known.

- ...:

  Additional arguments passed to `singscore::simpleScore()`.

## Value

The native data frame returned by `singscore::simpleScore()`.
