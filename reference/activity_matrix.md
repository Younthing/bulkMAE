# Convert long-format activity results to a source-by-sample matrix

Convert long-format activity results to a source-by-sample matrix

## Usage

``` r
activity_matrix(
  result,
  value = "score",
  source = "source",
  sample = "condition",
  statistic = NULL,
  statistic_column = "statistic"
)
```

## Arguments

- result:

  A data-frame-like result returned by
  [`activity_decouple()`](https://younthing.github.io/bulkMAE/reference/activity_decouple.md),
  [`activity_progeny()`](https://younthing.github.io/bulkMAE/reference/activity_progeny.md),
  or
  [`activity_tf()`](https://younthing.github.io/bulkMAE/reference/activity_tf.md).

- value, source, sample:

  Column names containing activity values, sources, and sample
  identifiers.

- statistic:

  Optional statistic/method value to retain when a result contains
  multiple methods.

- statistic_column:

  Column containing statistic/method names.

## Value

A numeric source-by-sample matrix suitable for
[`mae_add_experiment()`](https://younthing.github.io/bulkMAE/reference/mae_add_experiment.md).
